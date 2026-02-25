// sasa_batch.cpp
// Batch process CIF/PDB files with FreeSASA C API
//
// Usage: sasa_batch <input_folder> <output_folder>
// Compile: g++ -std=c++17 -O2 -o sasa_batch sasa_batch.cpp -lfreesasa -ljson-c -lxml2

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

extern "C" {
#include <freesasa.h>
}

namespace fs = std::filesystem;

bool has_structure_extension(const fs::path& path) {
    std::string ext = path.extension().string();
    for (char& c : ext) c = std::tolower(c);
    return ext == ".pdb" || ext == ".cif";
}

std::vector<fs::path> find_structure_files(const fs::path& dir) {
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file() && has_structure_extension(entry.path())) {
            files.push_back(entry.path());
        }
    }
    return files;
}

bool process_file(const fs::path& input, const fs::path& output, int n_points) {
    FILE* in = fopen(input.c_str(), "r");
    if (!in) return false;

    freesasa_structure* structure = freesasa_structure_from_pdb(
        in, &freesasa_default_classifier, 0);
    fclose(in);
    if (!structure) return false;

    freesasa_parameters params = freesasa_default_parameters;
    params.n_threads = 1;
    params.alg = FREESASA_SHRAKE_RUPLEY;
    params.shrake_rupley_n_points = n_points;

    freesasa_result* result = freesasa_calc_structure(structure, &params);
    freesasa_structure_free(structure);
    if (!result) return false;

    FILE* out = fopen(output.c_str(), "w");
    if (!out) { freesasa_result_free(result); return false; }

    fprintf(out, "%.2f\n", result->total);
    fclose(out);
    freesasa_result_free(result);

    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <input_folder> <output_folder> [n_points]\n";
        return 1;
    }

    fs::path input_dir(argv[1]), output_dir(argv[2]);
    int n_points = (argc >= 4) ? std::stoi(argv[3]) : 100;

    if (!fs::exists(input_dir) || !fs::is_directory(input_dir)) {
        std::cerr << "Error: Invalid input folder: " << input_dir << "\n";
        return 1;
    }

    fs::create_directories(output_dir);
    auto files = find_structure_files(input_dir);

    if (files.empty()) {
        std::cerr << "No PDB or CIF files found in " << input_dir << "\n";
        return 1;
    }

    freesasa_set_verbosity(FREESASA_V_NOWARNINGS);

    int success = 0, failed = 0;
    for (const auto& file : files) {
        fs::path out = output_dir / (file.stem().string() + ".json");
        if (process_file(file, out, n_points)) ++success; else ++failed;
    }

    std::cout << "Done: " << success << " succeeded, " << failed << " failed.\n";
    return failed > 0 ? 1 : 0;
}
