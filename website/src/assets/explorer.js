// SASA explorer: search the AFDB dimer and SwissProt tables with DuckDB-WASM.
// The Parquet files live next to the site (data/explorer/); a table is copied
// into DuckDB's memory the first time it is searched.

const DUCKDB = "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.32.0/+esm";
const DATA = new URL("../data/explorer/", import.meta.url);
const PAGE_SIZE = 50;
// The full tables are the Parquet downloads; a CSV holds the matching rows up to this
const CSV_LIMIT = 100000;

const area = (v) => (v == null ? "–" : Math.round(v).toLocaleString("en-US"));
const score = (v, digits = 2) => (v == null ? "–" : Number(v).toFixed(digits));
const pct = (part, whole) => (part == null || !whole ? "–" : `${Math.round((100 * part) / whole)}%`);

const esc = (s) =>
  String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

const uniprot = (acc) =>
  acc ? `<a class="ext" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(acc)}" rel="noopener">${esc(acc)}</a>` : "–";
const afdb = (id) => `<a class="ext" href="https://alphafold.ebi.ac.uk/entry/${encodeURIComponent(id)}" rel="noopener">${esc(id)}</a>`;
const gene = (g) => `<b>${esc(g || "–")}</b>`;

// A column: key, header, cell formatter, SQL sort expression (null: not
// sortable), text alignment, and whether it is shown before the reader picks.
// `fixed` columns are always shown.
const col = (key, label, fmt, sort, opts = {}) => ({ key, label, fmt, sort, num: true, show: true, ...opts });

const DIMER_AREAS = [
  col("sasa_a", "SASA A", (r) => area(r.sasa_a), "sasa_a"),
  col("sasa_b", "SASA B", (r) => area(r.sasa_b), "sasa_b"),
  col("sasa_ab", "SASA AB", (r) => area(r.sasa_ab), "sasa_ab"),
  col("delta_sasa", "ΔSASA", (r) => `<b>${area(r.delta_sasa)}</b>`, "delta_sasa"),
  col("bsa", "BSA", (r) => area(r.bsa), "bsa", { show: false }),
  col("delta_sasa_polar", "ΔSASA polar", (r) => area(r.delta_sasa_polar), "delta_sasa_polar", { show: false }),
  col("delta_sasa_apolar", "ΔSASA apolar", (r) => area(r.delta_sasa_apolar), "delta_sasa_apolar", { show: false }),
  col("apolar_pct", "Apolar %", (r) => pct(r.delta_sasa_apolar, r.delta_sasa), "delta_sasa_apolar / nullif(delta_sasa, 0)"),
  col("ipsae_min", "ipSAEmin", (r) => score(r.ipsae_min), "ipsae_min"),
  col("iptm", "ipTM", (r) => score(r.iptm), "iptm"),
  col("pdockq", "pDockQ", (r) => score(r.pdockq), "pdockq", { show: false }),
  col("plddt", "pLDDT", (r) => score(r.plddt, 1), "plddt", { show: false }),
  col("high_quality", "HQ", (r) => (r.high_quality ? '<span class="explorer__hq">HQ</span>' : ""), "high_quality"),
];

// Most dimer models are only in the AFDB/NVIDIA FTP release, so the ID is not a link
const MODEL = col("model_id", "Model", (r) => `<span class="explorer__id">${esc(r.model_id)}</span>`, "model_id", { num: false, fixed: true });

// A range filter: key, label, SQL expression, slider step and value formatter.
const range = (key, label, expr, step, fmt) => ({ key, label, expr, step, fmt });
const int = (v) => Math.round(v).toLocaleString("en-US");
const DIMER_FILTERS = [
  range("delta_sasa", "ΔSASA (Å²)", "delta_sasa", 50, int),
  range("apolar_pct", "Apolar %", "100 * delta_sasa_apolar / nullif(delta_sasa, 0)", 1, (v) => `${Math.round(v)}%`),
  range("ipsae_min", "ipSAEmin", "ipsae_min", 0.01, (v) => v.toFixed(2)),
  range("iptm", "ipTM", "iptm", 0.01, (v) => v.toFixed(2)),
  range("plddt", "pLDDT", "plddt", 1, (v) => v.toFixed(0)),
];

// Each dataset: its file, model count, the columns of the result table, the
// columns a search term is matched in, its range filters and its first order.
const DATASETS = {
  heterodimers: {
    label: "Heterodimers",
    file: "afdb_heterodimers.parquet",
    count: 531605,
    filters: DIMER_FILTERS,
    search: "UniProt accession, gene, AF model ID or protein name",
    match: ["accession_a", "accession_b", "gene_a", "gene_b", "model_id", "description_a", "description_b"],
    columns: [
      MODEL,
      col("gene_a", "Gene A", (r) => gene(r.gene_a), "gene_a", { num: false, fixed: true }),
      col("accession_a", "UniProt A", (r) => uniprot(r.accession_a), "accession_a", { num: false }),
      col("gene_b", "Gene B", (r) => gene(r.gene_b), "gene_b", { num: false, fixed: true }),
      col("accession_b", "UniProt B", (r) => uniprot(r.accession_b), "accession_b", { num: false }),
      ...DIMER_AREAS,
    ],
    sort: ["delta_sasa", "desc"],
    dimer: true,
  },
  homodimers: {
    label: "Homo-dimers",
    file: "afdb_homodimers.parquet",
    count: 91637,
    filters: DIMER_FILTERS,
    search: "UniProt accession, gene, AF model ID or protein name",
    match: ["accession_a", "gene_a", "model_id", "description_a"],
    columns: [
      MODEL,
      col("gene_a", "Gene", (r) => gene(r.gene_a), "gene_a", { num: false, fixed: true }),
      col("accession_a", "UniProt", (r) => uniprot(r.accession_a), "accession_a", { num: false }),
      ...DIMER_AREAS,
    ],
    sort: ["delta_sasa", "desc"],
    dimer: true,
  },
  swissprot: {
    label: "SwissProt",
    file: "swissprot.parquet",
    count: 550122,
    filters: [
      range("n_residues", "Residues", "n_residues", 1, int),
      range("sasa", "SASA (Å²)", "sasa", 50, int),
      range("apolar_pct", "Apolar %", "100 * sasa_apolar / nullif(sasa, 0)", 1, (v) => `${Math.round(v)}%`),
      range("plddt", "pLDDT", "plddt", 1, (v) => v.toFixed(0)),
    ],
    search: "UniProt accession, protein name or organism",
    match: ["accession", "model_id", "description", "organism"],
    columns: [
      col("accession", "UniProt", (r) => uniprot(r.accession), "accession", { num: false, fixed: true }),
      col("model_id", "AFDB", (r) => afdb(r.model_id), "model_id", { num: false }),
      col("protein", "Protein", (r) => `${esc(r.description)}<br><span class="explorer__sub">${esc(r.organism)}</span>`, "description", { num: false, fixed: true }),
      col("n_residues", "Residues", (r) => Number(r.n_residues).toLocaleString("en-US"), "n_residues"),
      col("sasa", "SASA", (r) => `<b>${area(r.sasa)}</b>`, "sasa"),
      col("sasa_polar", "Polar", (r) => area(r.sasa_polar), "sasa_polar"),
      col("sasa_apolar", "Apolar", (r) => area(r.sasa_apolar), "sasa_apolar"),
      col("apolar_pct", "Apolar %", (r) => pct(r.sasa_apolar, r.sasa), "sasa_apolar / nullif(sasa, 0)"),
      col("per_residue", "SASA / residue", (r) => (r.n_residues ? score(r.sasa / Number(r.n_residues), 1) : "–"), "sasa / nullif(n_residues, 0)"),
      col("plddt", "pLDDT", (r) => score(r.plddt, 1), "plddt", { show: false }),
    ],
    sort: ["sasa", "desc"],
    dimer: false,
  },
};

const state = {
  dataset: "heterodimers",
  hq: false,
  // Narrowed range filters of the current dataset: key -> [low, high]
  ranges: {},
  query: "",
  sort: DATASETS.heterodimers.sort,
  offset: 0,
};

// The columns a reader chose to see, kept per dataset in this browser only. The
// choice is stored as the columns turned on and off, so a column added later
// starts as its default.
function shownColumns(ds) {
  let saved = null;
  try {
    saved = JSON.parse(localStorage.getItem(`zsasa-explorer-cols-${ds.name}`));
  } catch {
    saved = null;
  }
  const choice = saved && typeof saved === "object" ? saved : {};
  return ds.columns.filter((c) => c.fixed || (choice[c.key] ?? c.show));
}

function saveColumns(ds, choice) {
  try {
    localStorage.setItem(`zsasa-explorer-cols-${ds.name}`, JSON.stringify(choice));
  } catch {
    // Storage may be unavailable (private window); the choice lasts for this page
  }
}

function readHash() {
  const params = new URLSearchParams(location.hash.slice(1));
  if (params.get("set") in DATASETS) state.dataset = params.get("set");
  const ds = DATASETS[state.dataset];
  state.query = params.get("q") ?? "";
  const [key, dir] = (params.get("sort") ?? "").split("-");
  state.sort = ds.columns.some((c) => c.key === key && c.sort) && ["asc", "desc"].includes(dir) ? [key, dir] : ds.sort;
  state.hq = params.get("hq") === "1";
  // r=delta_sasa:1000~5000,ipsae_min:0.6~1
  state.ranges = {};
  for (const item of (params.get("r") ?? "").split(",")) {
    const [k, span] = item.split(":");
    const [lo, hi] = (span ?? "").split("~").map(Number);
    if (ds.filters.some((f) => f.key === k) && Number.isFinite(lo) && Number.isFinite(hi) && lo <= hi) {
      state.ranges[k] = [lo, hi];
    }
  }
}

function writeHash() {
  const params = new URLSearchParams({ set: state.dataset });
  if (state.query) params.set("q", state.query);
  params.set("sort", state.sort.join("-"));
  if (state.hq && DATASETS[state.dataset].dimer) params.set("hq", "1");
  const ranges = Object.entries(state.ranges).map(([k, [lo, hi]]) => `${k}:${lo}~${hi}`);
  if (ranges.length) params.set("r", ranges.join(","));
  history.replaceState(null, "", `#${params}`);
}

function orderBy(ds) {
  const [key, dir] = state.sort;
  const column = ds.columns.find((c) => c.key === key) ?? ds.columns.find((c) => c.key === ds.sort[0]);
  return `${column.sort} ${dir === "asc" ? "ASC" : "DESC"} NULLS LAST, model_id`;
}

// The WHERE clause of the current state and its parameters, in order.
function where() {
  const ds = DATASETS[state.dataset];
  const clauses = [];
  const params = [];
  const q = state.query.trim();
  if (q) {
    clauses.push("contains(search_text, ?)");
    params.push(q.toLowerCase());
  }
  if (ds.dimer && state.hq) clauses.push("high_quality");
  for (const f of ds.filters) {
    const span = state.ranges[f.key];
    if (!span) continue;
    clauses.push(`(${f.expr}) BETWEEN ? AND ?`);
    params.push(span[0], span[1]);
  }
  return [clauses.length ? `WHERE ${clauses.join(" AND ")}` : "", params];
}

let conn = null;

async function connect(status) {
  status.textContent = "Loading DuckDB…";
  const duckdb = await import(DUCKDB);
  const bundle = await duckdb.selectBundle(duckdb.getJsDelivrBundles());
  const workerUrl = URL.createObjectURL(
    new Blob([`importScripts("${bundle.mainWorker}");`], { type: "text/javascript" }),
  );
  const db = new duckdb.AsyncDuckDB(new duckdb.VoidLogger(), new Worker(workerUrl));
  await db.instantiate(bundle.mainModule, bundle.pthreadWorker);
  URL.revokeObjectURL(workerUrl);
  for (const ds of Object.values(DATASETS)) {
    await db.registerFileURL(ds.file, new URL(ds.file, DATA).href, duckdb.DuckDBDataProtocol.HTTP, false);
  }
  return db.connect();
}

// Copy a dataset into DuckDB's memory the first time it is searched, with the
// searched columns joined into one lower-case column. Over HTTP every search
// would read and scan those columns again.
async function load(ds, status) {
  if (ds.loading) return ds.loading;
  status.textContent = `Loading the ${ds.label} table…`;
  const text = `lower(concat_ws(' ', ${ds.match.join(", ")})) AS search_text`;
  ds.loading = (async () => {
    await conn.query(`CREATE TABLE ${ds.name} AS SELECT *, ${text} FROM read_parquet('${ds.file}')`);
    ds.table = ds.name;
    // The slider ends: the range of each filter, on its step
    const bounds = ds.filters.map((f) => `min(${f.expr}) AS "${f.key}_lo", max(${f.expr}) AS "${f.key}_hi"`).join(", ");
    const [row] = (await conn.query(`SELECT ${bounds} FROM ${ds.name}`)).toArray().map((r) => r.toJSON());
    // Rounded to the decimals of the step, so 0.94 does not come out as 0.9400000000000001
    const onStep = (v, f) => Number(v.toFixed((String(f.step).split(".")[1] ?? "").length));
    ds.bounds = Object.fromEntries(
      ds.filters.map((f) => [
        f.key,
        [
          onStep(Math.floor(Number(row[`${f.key}_lo`]) / f.step) * f.step, f),
          onStep(Math.ceil(Number(row[`${f.key}_hi`]) / f.step) * f.step, f),
        ],
      ]),
    );
  })();
  return ds.loading;
}

async function run(sql, params) {
  const stmt = await conn.prepare(sql);
  try {
    const result = await stmt.query(...params);
    return result.toArray().map((row) => row.toJSON());
  } finally {
    await stmt.close();
  }
}

function el(tag, attrs = {}, html = "") {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
  node.innerHTML = html;
  return node;
}

function field(label, attrs, value, onInput) {
  const wrap = el("label", { class: "explorer__ctl" });
  wrap.append(el("span", { class: "label" }, esc(label)));
  const input = el("input", attrs);
  input.value = value;
  input.addEventListener("input", () => onInput(input.value));
  wrap.append(input);
  return wrap;
}

function detailTable(r) {
  const rows = [
    ["Chain A", r.sasa_a, r.sasa_a_polar, r.sasa_a_apolar],
    ["Chain B", r.sasa_b, r.sasa_b_polar, r.sasa_b_apolar],
    ["Complex AB", r.sasa_ab, r.sasa_ab_polar, r.sasa_ab_apolar],
    ["ΔSASA (A + B − AB)", r.delta_sasa, r.delta_sasa_polar, r.delta_sasa_apolar],
    ["BSA (ΔSASA / 2)", r.bsa, r.delta_sasa_polar / 2, r.delta_sasa_apolar / 2],
  ];
  const body = rows
    .map(
      ([name, total, polar, apolar]) =>
        `<tr><th scope="row">${name}</th><td class="num">${area(total)}</td><td class="num">${area(polar)}</td>` +
        `<td class="num">${area(apolar)}</td><td class="num">${pct(apolar, total)}</td></tr>`,
    )
    .join("");
  const chain = (label, desc, code, n) => `<p><b>${label}</b> ${esc(desc)} (${esc(code)}, ${Number(n)} residues)</p>`;
  return (
    `<div class="explorer__detail">${chain("A", r.description_a, r.entry_name_a, r.n_residues_a)}` +
    `${chain("B", r.description_b, r.entry_name_b, r.n_residues_b)}` +
    `<table><thead><tr><th>Area (Å²)</th><th class="num">Total</th><th class="num">Polar</th>` +
    `<th class="num">Apolar</th><th class="num">Apolar %</th></tr></thead><tbody>${body}</tbody></table>` +
    `<p class="explorer__sub">pLDDT ${score(r.plddt, 1)} · pDockQ ${score(r.pdockq)} · backbone clashes ${r.n_clash_backbone ?? "–"}</p></div>`
  );
}

function toCsv(rows) {
  if (!rows.length) return "";
  const keys = Object.keys(rows[0]);
  const cell = (v) => {
    const s = typeof v === "bigint" ? v.toString() : String(v ?? "");
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  };
  return [keys.join(","), ...rows.map((r) => keys.map((k) => cell(r[k])).join(","))].join("\n");
}

// A two-thumb slider over `bounds` for one range filter. `value` is the
// narrowed range or undefined; `onChange` gets the new range, or undefined when
// both thumbs are back at the ends.
function rangeSlider(f, bounds, value, onChange) {
  const [min, max] = bounds;
  const wrap = el("div", { class: "explorer__range-field" });
  const head = el("div", { class: "explorer__range-head" });
  const out = el("output", { class: "explorer__range-value" });
  head.append(el("span", { class: "label" }, esc(f.label)), out);
  const track = el("div", { class: "explorer__slider" });
  const fill = el("div", { class: "explorer__slider-fill" });
  const attrs = { type: "range", min: String(min), max: String(max), step: String(f.step) };
  const low = el("input", { ...attrs, "aria-label": `${f.label}, minimum` });
  const high = el("input", { ...attrs, "aria-label": `${f.label}, maximum` });
  low.value = String(value?.[0] ?? min);
  high.value = String(value?.[1] ?? max);
  track.append(fill, low, high);
  wrap.append(head, track);
  const show = () => {
    const [a, b] = [Number(low.value), Number(high.value)];
    const span = max - min || 1;
    fill.style.left = `${(100 * (a - min)) / span}%`;
    fill.style.right = `${(100 * (max - b)) / span}%`;
    out.textContent = `${f.fmt(a)} – ${f.fmt(b)}`;
    wrap.classList.toggle("is-set", narrowed(a, b));
  };
  const narrowed = (a, b) => a > min + f.step / 2 || b < max - f.step / 2;
  const changed = (moved, other, keepBelow) => {
    // The thumbs do not cross
    if (keepBelow ? Number(moved.value) > Number(other.value) : Number(moved.value) < Number(other.value)) {
      moved.value = other.value;
    }
    show();
    const [a, b] = [Number(low.value), Number(high.value)];
    onChange(narrowed(a, b) ? [a, b] : undefined);
  };
  low.addEventListener("input", () => changed(low, high, true));
  high.addEventListener("input", () => changed(high, low, false));
  show();
  return wrap;
}

async function main() {
  const root = document.querySelector("[data-explorer]");
  if (!root) return;
  for (const [name, ds] of Object.entries(DATASETS)) ds.name = name;
  readHash();
  root.classList.add("explorer");
  root.innerHTML = "";

  // The datasets get a row of their own; more may be added later
  const tabs = el("div", { class: "explorer__tabs", role: "tablist", "aria-label": "Dataset" });
  const top = el("div", { class: "explorer__top" });
  const sliders = el("div", { class: "explorer__ranges" });
  const options = el("div", { class: "explorer__options" });
  const picker = el("fieldset", { class: "explorer__columns" });
  const status = el("p", { class: "explorer__status", role: "status" });
  const tableWrap = el("div", { class: "table explorer__table" });
  const pager = el("div", { class: "explorer__pager" });
  const search = el("div", { class: "explorer__panel explorer__panel--search", role: "search" });
  search.append(top, sliders, options, picker);
  // The count of the results, and the CSV of them at its right
  const resultsHead = el("div", { class: "explorer__results-head" });
  const results = el("div", { class: "explorer__panel explorer__panel--results" });
  results.append(resultsHead, tableWrap, pager);
  root.append(tabs, search, results);

  let timer = null;
  const refresh = () => {
    state.offset = 0;
    clearTimeout(timer);
    timer = setTimeout(update, 250);
  };
  let lastRows = [];

  const csv = el("button", { type: "button", class: "explorer__btn explorer__btn--primary" }, "Download CSV");
  csv.disabled = true;
  csv.addEventListener("click", () => downloadCsv());
  resultsHead.append(status, csv);

  function buildTabs() {
    tabs.replaceChildren();
    for (const [key, ds] of Object.entries(DATASETS)) {
      const tab = el(
        "button",
        { type: "button", role: "tab", class: "explorer__tab", "aria-selected": String(key === state.dataset) },
        `${esc(ds.label)}<span>${ds.count.toLocaleString("en-US")}</span>`,
      );
      tab.addEventListener("click", () => {
        if (key === state.dataset) return;
        state.dataset = key;
        state.sort = ds.sort;
        state.ranges = {};
        buildControls();
        refresh();
      });
      tabs.append(tab);
    }
  }

  function buildControls() {
    const ds = DATASETS[state.dataset];
    buildTabs();
    top.replaceChildren(
      field("Search", { type: "search", placeholder: ds.search, class: "explorer__search" }, state.query, (v) => {
        state.query = v;
        refresh();
      }),
    );
    buildSliders();
    options.replaceChildren();
    if (ds.dimer) {
      const hq = el("label", { class: "explorer__check" });
      const box = el("input", { type: "checkbox" });
      box.checked = state.hq;
      box.addEventListener("change", () => {
        state.hq = box.checked;
        refresh();
      });
      hq.append(box, document.createTextNode(" AFDB high quality only"));
      options.append(hq);
    }
    const reset = el("button", { type: "button", class: "explorer__btn" }, "Reset filters");
    reset.addEventListener("click", () => {
      state.ranges = {};
      state.hq = false;
      state.query = "";
      buildControls();
      refresh();
    });
    options.append(reset);
    buildPicker();
  }

  // The sliders need the range of each filter, known once the table is loaded
  function buildSliders() {
    const ds = DATASETS[state.dataset];
    sliders.replaceChildren();
    if (!ds.bounds) return;
    for (const f of ds.filters) {
      const bounds = ds.bounds[f.key];
      const value = state.ranges[f.key];
      if (value) state.ranges[f.key] = [Math.max(value[0], bounds[0]), Math.min(value[1], bounds[1])];
      sliders.append(
        rangeSlider(f, bounds, state.ranges[f.key], (span) => {
          if (span) state.ranges[f.key] = span;
          else delete state.ranges[f.key];
          refresh();
        }),
      );
    }
  }

  async function downloadCsv() {
    const ds = DATASETS[state.dataset];
    if (!ds.table) return;
    const [clause, params] = where();
    csv.disabled = true;
    const label = csv.textContent;
    csv.textContent = "Preparing CSV…";
    try {
      const rows = await run(
        `SELECT * EXCLUDE (search_text) FROM ${ds.table} ${clause} ORDER BY ${orderBy(ds)} LIMIT ${CSV_LIMIT}`,
        params,
      );
      const link = el("a", {
        href: URL.createObjectURL(new Blob([toCsv(rows)], { type: "text/csv" })),
        download: `zsasa-${ds.name}.csv`,
      });
      link.click();
      URL.revokeObjectURL(link.href);
    } finally {
      csv.disabled = false;
      csv.textContent = label;
    }
  }

  function csvLabel(n) {
    if (!n) return "Download CSV";
    const rows = Math.min(n, CSV_LIMIT).toLocaleString("en-US");
    return n > CSV_LIMIT ? `Download CSV (first ${rows} rows)` : `Download CSV (${rows} rows)`;
  }

  function buildPicker() {
    const ds = DATASETS[state.dataset];
    const shown = new Set(shownColumns(ds).map((c) => c.key));
    picker.replaceChildren(el("legend", { class: "label" }, "Columns"));
    for (const c of ds.columns.filter((c) => !c.fixed)) {
      const chip = el("label", { class: "explorer__chip" });
      const box = el("input", { type: "checkbox" });
      box.checked = shown.has(c.key);
      box.addEventListener("change", () => {
        const choice = Object.fromEntries(
          [...picker.querySelectorAll("input")].map((b, i) => [ds.columns.filter((c) => !c.fixed)[i].key, b.checked]),
        );
        saveColumns(ds, choice);
        render(ds, lastRows);
      });
      chip.append(box, document.createTextNode(c.label));
      picker.append(chip);
    }
  }

  let generation = 0;
  async function update() {
    const mine = ++generation;
    const ds = DATASETS[state.dataset];
    writeHash();
    const [clause, params] = where();
    try {
      const fresh = !ds.table;
      await load(ds, status);
      if (mine !== generation) return;
      if (fresh || !sliders.childElementCount) buildSliders();
      status.textContent = "Searching…";
      const summary = ds.dimer ? "median(delta_sasa)" : "median(sasa)";
      // One connection runs one statement at a time
      const [stats] = await run(`SELECT count(*) AS n, ${summary} AS median_value FROM ${ds.table} ${clause}`, params);
      const rows = await run(
        `SELECT * EXCLUDE (search_text) FROM ${ds.table} ${clause} ORDER BY ${orderBy(ds)} LIMIT ${PAGE_SIZE} OFFSET ${state.offset}`,
        params,
      );
      if (mine !== generation) return;
      const n = Number(stats.n);
      const median = ds.dimer ? "median ΔSASA" : "median SASA";
      status.textContent = n
        ? `${n.toLocaleString("en-US")} ${n === 1 ? "model matches" : "models match"} · ${median} ${area(stats.median_value)} Å²`
        : "No model matches.";
      lastRows = rows;
      csv.textContent = csvLabel(n);
      csv.disabled = !n;
      render(ds, rows);
      renderPager(n);
    } catch (err) {
      if (mine !== generation) return;
      status.textContent = `The query failed: ${err.message}`;
    }
  }

  function render(ds, rows) {
    const columns = shownColumns(ds);
    const [sortKey, sortDir] = state.sort;
    const head = columns
      .map((c) => {
        const cls = c.num ? ' class="num"' : "";
        if (!c.sort) return `<th${cls}>${esc(c.label)}</th>`;
        const active = c.key === sortKey;
        const ariaSort = active ? (sortDir === "asc" ? "ascending" : "descending") : "none";
        const mark = active ? (sortDir === "asc" ? " ▲" : " ▼") : "";
        return `<th${cls} aria-sort="${ariaSort}"><button type="button" class="explorer__sort" data-key="${c.key}">${esc(c.label)}${mark}</button></th>`;
      })
      .join("");
    const body = rows
      .map(
        (r, i) =>
          `<tr data-row="${i}"${ds.dimer ? ' class="explorer__row" tabindex="0" aria-expanded="false"' : ""}>` +
          columns.map((c) => `<td${c.num ? ' class="num"' : ""}>${c.fmt(r)}</td>`).join("") +
          "</tr>",
      )
      .join("");
    tableWrap.innerHTML = `<table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;

    // A first click sorts numbers from the largest and text from A
    for (const button of tableWrap.querySelectorAll(".explorer__sort")) {
      button.addEventListener("click", () => {
        const key = button.dataset.key;
        const column = ds.columns.find((c) => c.key === key);
        const dir = state.sort[0] === key ? (state.sort[1] === "desc" ? "asc" : "desc") : column.num ? "desc" : "asc";
        state.sort = [key, dir];
        state.offset = 0;
        update();
      });
    }

    if (!ds.dimer) return;
    const toggle = (tr) => {
      const next = tr.nextElementSibling;
      if (next?.classList.contains("explorer__open")) {
        next.remove();
        tr.setAttribute("aria-expanded", "false");
        return;
      }
      const detail = el("tr", { class: "explorer__open" });
      detail.append(el("td", { colspan: String(columns.length) }, detailTable(rows[Number(tr.dataset.row)])));
      tr.after(detail);
      tr.setAttribute("aria-expanded", "true");
    };
    for (const tr of tableWrap.querySelectorAll(".explorer__row")) {
      tr.addEventListener("click", (e) => {
        if (!e.target.closest("a")) toggle(tr);
      });
      tr.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          toggle(tr);
        }
      });
    }
  }

  function renderPager(n) {
    pager.replaceChildren();
    if (!n) return;
    const first = state.offset + 1;
    const last = Math.min(state.offset + PAGE_SIZE, n);
    const prev = el("button", { type: "button", class: "explorer__btn" }, "← Previous");
    prev.disabled = state.offset === 0;
    prev.addEventListener("click", () => {
      state.offset = Math.max(0, state.offset - PAGE_SIZE);
      update();
    });
    const next = el("button", { type: "button", class: "explorer__btn" }, "Next →");
    next.disabled = last >= n;
    next.addEventListener("click", () => {
      state.offset += PAGE_SIZE;
      update();
    });
    pager.append(
      prev,
      el("span", { class: "explorer__range" }, `${first.toLocaleString("en-US")}–${last.toLocaleString("en-US")} of ${n.toLocaleString("en-US")}`),
      next,
    );
  }


  buildControls();
  try {
    conn = await connect(status);
    const probe = await fetch(new URL(DATASETS[state.dataset].file, DATA), { method: "HEAD" });
    if (!probe.ok) {
      status.textContent = "The data files are not part of this build of the site.";
      return;
    }
  } catch (err) {
    status.textContent = `DuckDB could not start in this browser: ${err.message}`;
    return;
  }
  // A link to another search on this page changes only the hash; a link to a
  // section (#download) is not a search and leaves the explorer as it is
  window.addEventListener("hashchange", () => {
    if (!new URLSearchParams(location.hash.slice(1)).has("set")) return;
    readHash();
    state.offset = 0;
    buildControls();
    update();
  });
  update();
}

main();
