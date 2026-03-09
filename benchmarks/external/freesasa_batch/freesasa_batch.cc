// freesasa_batch.cc
// Batch process PDB files with FreeSASA C API and file-level parallelism
//
// Usage: freesasa_batch <input_dir> <output_dir> [--n-threads=N] [--n-points=N] [--algorithm=sr|lr]
//
// Build (from benchmarks/external/):
//   c++ -O3 -std=c++17 -I freesasa/src \
//     -o freesasa_batch/freesasa_batch freesasa_batch/freesasa_batch.cc \
//     freesasa/src/libfreesasa.a -lpthread

#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <pthread.h>
#include <string>
#include <vector>

extern "C" {
#include "freesasa.h"
}

namespace fs = std::filesystem;

struct WorkItem {
    fs::path input;
    fs::path output;
};

struct WorkerContext {
    const std::vector<WorkItem>* items;
    std::atomic<size_t>* next_index;
    int n_points;
    freesasa_algorithm algorithm;
    std::atomic<int>* success_count;
    std::atomic<int>* fail_count;
};

static bool has_pdb_extension(const fs::path& path) {
    std::string ext = path.extension().string();
    for (char& c : ext) c = std::tolower(c);
    return ext == ".pdb";
}

static std::vector<fs::path> find_pdb_files(const fs::path& dir) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && has_pdb_extension(entry.path())) {
            files.push_back(entry.path());
        }
    }
    return files;
}

static bool process_pdb_file(const fs::path& input, const fs::path& output,
                             int n_points, freesasa_algorithm algorithm) {
    FILE* in = fopen(input.c_str(), "r");
    if (!in) {
        fprintf(stderr, "Warning: cannot open %s: %s\n", input.c_str(), strerror(errno));
        return false;
    }

    freesasa_structure* structure =
        freesasa_structure_from_pdb(in, &freesasa_default_classifier, 0);
    fclose(in);

    if (!structure) {
        fprintf(stderr, "Warning: failed to parse structure from %s\n", input.c_str());
        return false;
    }

    freesasa_parameters params = freesasa_default_parameters;
    params.n_threads = 1;
    params.alg = algorithm;
    if (algorithm == FREESASA_LEE_RICHARDS) {
        params.lee_richards_n_slices = n_points;
    } else {
        params.shrake_rupley_n_points = n_points;
    }

    freesasa_result* result = freesasa_calc_structure(structure, &params);
    freesasa_structure_free(structure);
    if (!result) {
        fprintf(stderr, "Warning: SASA calculation failed for %s\n", input.c_str());
        return false;
    }

    FILE* out = fopen(output.c_str(), "w");
    if (!out) {
        fprintf(stderr, "Warning: cannot write %s: %s\n", output.c_str(), strerror(errno));
        freesasa_result_free(result);
        return false;
    }

    fprintf(out, "%.2f\n", result->total);
    fclose(out);
    freesasa_result_free(result);
    return true;
}

static void* worker_thread(void* arg) {
    auto* ctx = static_cast<WorkerContext*>(arg);
    size_t total = ctx->items->size();

    while (true) {
        size_t idx = ctx->next_index->fetch_add(1);
        if (idx >= total) break;

        const auto& item = (*ctx->items)[idx];
        if (process_pdb_file(item.input, item.output, ctx->n_points, ctx->algorithm)) {
            ctx->success_count->fetch_add(1);
        } else {
            ctx->fail_count->fetch_add(1);
        }
    }

    return nullptr;
}

static bool parse_int_arg(const char* arg, const char* prefix, int* out) {
    size_t prefix_len = strlen(prefix);
    if (strncmp(arg, prefix, prefix_len) != 0) return false;
    char* end;
    errno = 0;
    long val = strtol(arg + prefix_len, &end, 10);
    if (*end != '\0' || errno != 0 || val <= 0 || val > 100000) return false;
    *out = static_cast<int>(val);
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dir> <output_dir> [--n-threads=N] [--n-points=N]"
                     " [--algorithm=sr|lr]\n";
        return 1;
    }

    fs::path input_dir(argv[1]);
    fs::path output_dir(argv[2]);
    int n_threads = 1;
    int n_points = 100;
    freesasa_algorithm algorithm = FREESASA_SHRAKE_RUPLEY;

    for (int i = 3; i < argc; i++) {
        if (parse_int_arg(argv[i], "--n-threads=", &n_threads)) continue;
        if (parse_int_arg(argv[i], "--n-points=", &n_points)) continue;
        if (strncmp(argv[i], "--algorithm=", 12) == 0) {
            const char* alg = argv[i] + 12;
            if (strcmp(alg, "lr") == 0 || strcmp(alg, "lee-richards") == 0) {
                algorithm = FREESASA_LEE_RICHARDS;
            } else if (strcmp(alg, "sr") == 0 || strcmp(alg, "shrake-rupley") == 0) {
                algorithm = FREESASA_SHRAKE_RUPLEY;
            } else {
                std::cerr << "Error: Unknown algorithm: " << alg
                          << " (use sr or lr)\n";
                return 1;
            }
            continue;
        }
        std::cerr << "Error: Unknown argument: " << argv[i] << "\n";
        return 1;
    }

    if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
        std::cerr << "Error: Invalid input directory: " << input_dir << "\n";
        return 1;
    }

    try {
        fs::create_directories(output_dir);
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error: Cannot create output directory: " << e.what() << "\n";
        return 1;
    }

    auto files = find_pdb_files(input_dir);
    if (files.empty()) {
        std::cerr << "No PDB files found in " << input_dir << "\n";
        return 1;
    }

    // Build work items
    std::vector<WorkItem> items;
    items.reserve(files.size());
    for (const auto& file : files) {
        items.push_back({file, output_dir / (file.stem().string() + ".txt")});
    }

    freesasa_set_verbosity(FREESASA_V_NOWARNINGS);

    // Clamp thread count to file count
    int actual_threads = n_threads;
    if (actual_threads > static_cast<int>(items.size())) {
        actual_threads = static_cast<int>(items.size());
    }

    std::atomic<size_t> next_index{0};
    std::atomic<int> success_count{0};
    std::atomic<int> fail_count{0};

    WorkerContext ctx{&items, &next_index, n_points, algorithm, &success_count, &fail_count};

    if (actual_threads <= 1) {
        worker_thread(&ctx);
    } else {
        std::vector<pthread_t> threads(actual_threads);
        int created = 0;
        for (int i = 0; i < actual_threads; i++) {
            int rc = pthread_create(&threads[i], nullptr, worker_thread, &ctx);
            if (rc != 0) {
                std::cerr << "Error: Failed to create thread " << i << "\n";
                next_index.store(items.size());
                for (int j = 0; j < created; j++) {
                    pthread_join(threads[j], nullptr);
                }
                return 1;
            }
            created++;
        }
        for (int i = 0; i < actual_threads; i++) {
            pthread_join(threads[i], nullptr);
        }
    }

    const char* alg_name = (algorithm == FREESASA_LEE_RICHARDS) ? "lr" : "sr";
    std::cout << "Done: " << success_count.load() << " succeeded, "
              << fail_count.load() << " failed ("
              << actual_threads << " threads, " << n_points << " points, "
              << alg_name << ").\n";

    return fail_count.load() > 0 ? 1 : 0;
}
