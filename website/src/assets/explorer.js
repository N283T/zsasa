// SASA explorer: search the AFDB dimer and SwissProt tables with DuckDB-WASM.
// The Parquet files live next to the site (data/explorer/); DuckDB reads them
// with HTTP range requests, so a query fetches only the columns it needs.

const DUCKDB = "https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@1.32.0/+esm";
const DATA = new URL("../data/explorer/", import.meta.url);
const PAGE_SIZE = 50;
const CSV_LIMIT = 10000;

const area = (v) => (v == null ? "–" : Math.round(v).toLocaleString("en-US"));
const score = (v, digits = 2) => (v == null ? "–" : Number(v).toFixed(digits));
const pct = (part, whole) => (part == null || !whole ? "–" : `${Math.round((100 * part) / whole)}%`);

const esc = (s) =>
  String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);

const uniprot = (acc) =>
  acc ? `<a class="ext" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(acc)}" rel="noopener">${esc(acc)}</a>` : "–";
const afdb = (id) => `<a class="ext" href="https://alphafold.ebi.ac.uk/entry/${encodeURIComponent(id)}" rel="noopener">${esc(id)}</a>`;

const partner = (gene, acc) => `<b>${esc(gene || "–")}</b><br><span class="explorer__sub">${uniprot(acc)}</span>`;

// Each dataset: its file, the columns of the result table, how a search term
// matches, the extra filters it offers and the orders a reader can pick.
const DATASETS = {
  dimers: {
    label: "AFDB human dimers",
    file: "afdb_dimers.parquet",
    search: "UniProt accession, gene, AF model ID or protein name",
    match: ["accession_a", "accession_b", "gene_a", "gene_b", "model_id", "description_a", "description_b"],
    columns: [
      ["Model", (r) => `${afdb(r.model_id)}<br><span class="explorer__sub">${r.kind}-dimer</span>`],
      ["Chain A", (r) => partner(r.gene_a, r.accession_a)],
      ["Chain B", (r) => partner(r.gene_b, r.accession_b)],
      ["SASA A", (r) => area(r.sasa_a), "num"],
      ["SASA B", (r) => area(r.sasa_b), "num"],
      ["SASA AB", (r) => area(r.sasa_ab), "num"],
      ["ΔSASA", (r) => `<b>${area(r.delta_sasa)}</b>`, "num"],
      ["BSA", (r) => area(r.bsa), "num"],
      ["Apolar %", (r) => pct(r.delta_sasa_apolar, r.delta_sasa), "num"],
      ["ipSAEmin", (r) => score(r.ipsae_min), "num"],
      ["ipTM", (r) => score(r.iptm), "num"],
      ["pLDDT", (r) => score(r.plddt, 1), "num"],
      ["HQ", (r) => (r.high_quality ? '<span class="explorer__hq">HQ</span>' : ""), "num"],
    ],
    orders: {
      delta: ["ΔSASA, largest first", "delta_sasa DESC"],
      ipsae: ["ipSAEmin, highest first", "ipsae_min DESC NULLS LAST, delta_sasa DESC"],
      apolar: ["Apolar %, highest first", "delta_sasa_apolar / nullif(delta_sasa, 0) DESC NULLS LAST"],
      polar: ["Polar %, highest first", "delta_sasa_polar / nullif(delta_sasa, 0) DESC NULLS LAST"],
      plddt: ["pLDDT, highest first", "plddt DESC NULLS LAST"],
      gene: ["Gene A, A–Z", "gene_a, gene_b, model_id"],
    },
    detail: true,
  },
  swissprot: {
    label: "SwissProt",
    file: "swissprot.parquet",
    search: "UniProt accession, protein name or organism",
    match: ["accession", "model_id", "description", "organism"],
    columns: [
      ["Accession", (r) => `${uniprot(r.accession)}<br><span class="explorer__sub">${afdb(r.model_id)}</span>`],
      ["Protein", (r) => `${esc(r.description)}<br><span class="explorer__sub">${esc(r.organism)}</span>`],
      ["Residues", (r) => Number(r.n_residues).toLocaleString("en-US"), "num"],
      ["SASA", (r) => `<b>${area(r.sasa)}</b>`, "num"],
      ["Polar", (r) => area(r.sasa_polar), "num"],
      ["Apolar", (r) => area(r.sasa_apolar), "num"],
      ["Apolar %", (r) => pct(r.sasa_apolar, r.sasa), "num"],
      ["SASA / residue", (r) => (r.n_residues ? score(r.sasa / Number(r.n_residues), 1) : "–"), "num"],
      ["pLDDT", (r) => score(r.plddt, 1), "num"],
    ],
    orders: {
      sasa: ["SASA, largest first", "sasa DESC"],
      apolar: ["Apolar %, highest first", "sasa_apolar / nullif(sasa, 0) DESC"],
      polar: ["Polar %, highest first", "sasa_polar / nullif(sasa, 0) DESC"],
      per_residue: ["SASA / residue, highest first", "sasa / nullif(n_residues, 0) DESC"],
      plddt: ["pLDDT, highest first", "plddt DESC NULLS LAST"],
      accession: ["Accession, A–Z", "accession, fragment"],
    },
    detail: false,
  },
};

const state = {
  dataset: "dimers",
  kind: "all",
  hq: false,
  minDelta: "",
  minIpsae: "",
  query: "",
  order: "delta",
  offset: 0,
};

function readHash() {
  const params = new URLSearchParams(location.hash.slice(1));
  if (params.get("set") in DATASETS) state.dataset = params.get("set");
  if (["all", "homo", "hetero"].includes(params.get("kind"))) state.kind = params.get("kind");
  state.query = params.get("q") ?? "";
  if (params.get("order") in DATASETS[state.dataset].orders) state.order = params.get("order");
  else state.order = Object.keys(DATASETS[state.dataset].orders)[0];
}

function writeHash() {
  const params = new URLSearchParams({ set: state.dataset });
  if (state.dataset === "dimers" && state.kind !== "all") params.set("kind", state.kind);
  if (state.query) params.set("q", state.query);
  params.set("order", state.order);
  history.replaceState(null, "", `#${params}`);
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
  if (state.dataset === "dimers") {
    if (state.kind !== "all") {
      clauses.push("kind = ?");
      params.push(state.kind);
    }
    if (state.hq) clauses.push("high_quality");
    const minDelta = Number.parseFloat(state.minDelta);
    if (Number.isFinite(minDelta)) {
      clauses.push("delta_sasa >= ?");
      params.push(minDelta);
    }
    const minIpsae = Number.parseFloat(state.minIpsae);
    if (Number.isFinite(minIpsae)) {
      clauses.push("ipsae_min >= ?");
      params.push(minIpsae);
    }
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
  for (const [name, ds] of Object.entries(DATASETS)) {
    await db.registerFileURL(ds.file, new URL(ds.file, DATA).href, duckdb.DuckDBDataProtocol.HTTP, false);
    ds.name = name;
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
  ds.loading = conn
    .query(`CREATE TABLE ${ds.name} AS SELECT *, ${text} FROM read_parquet('${ds.file}')`)
    .then(() => {
      ds.table = ds.name;
    });
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

function seg(label, options, current, onPick) {
  const wrap = el("div", { class: "explorer__ctl" });
  wrap.append(el("span", { class: "label" }, esc(label)));
  const group = el("div", { class: "seg", role: "group", "aria-label": label });
  for (const [value, text] of options) {
    const b = el("button", { type: "button", "aria-pressed": String(value === current) }, esc(text));
    b.addEventListener("click", () => {
      for (const other of group.children) other.setAttribute("aria-pressed", "false");
      b.setAttribute("aria-pressed", "true");
      onPick(value);
    });
    group.append(b);
  }
  wrap.append(group);
  return wrap;
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
  const names = [
    `<p><b>A</b> ${esc(r.description_a)} (${esc(r.entry_name_a)}, ${Number(r.n_residues_a)} residues)</p>`,
    `<p><b>B</b> ${esc(r.description_b)} (${esc(r.entry_name_b)}, ${Number(r.n_residues_b)} residues)</p>`,
  ].join("");
  return (
    `<div class="explorer__detail">${names}` +
    `<table><thead><tr><th>Area (Å²)</th><th class="num">Total</th><th class="num">Polar</th>` +
    `<th class="num">Apolar</th><th class="num">Apolar %</th></tr></thead><tbody>${body}</tbody></table>` +
    `<p class="explorer__sub">pDockQ ${score(r.pdockq)} · backbone clashes ${r.n_clash_backbone ?? "–"}</p></div>`
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

async function main() {
  const root = document.querySelector("[data-explorer]");
  if (!root) return;
  readHash();
  root.classList.add("explorer");
  root.innerHTML = "";

  const controls = el("div", { class: "explorer__controls" });
  const filters = el("div", { class: "explorer__controls" });
  const status = el("p", { class: "explorer__status", role: "status" });
  const tableWrap = el("div", { class: "table explorer__table" });
  const pager = el("div", { class: "explorer__pager" });
  // The search panel stays in view while the results scroll under it
  const search = el("div", { class: "explorer__panel explorer__panel--search", role: "search" });
  search.append(controls, filters);
  const results = el("div", { class: "explorer__panel explorer__panel--results" });
  results.append(status, tableWrap, pager);
  root.append(search, results);

  let timer = null;
  const refresh = (resetOffset = true) => {
    if (resetOffset) state.offset = 0;
    clearTimeout(timer);
    timer = setTimeout(update, 250);
  };

  function buildControls() {
    const ds = DATASETS[state.dataset];
    controls.replaceChildren(
      seg(
        "Data",
        Object.entries(DATASETS).map(([k, v]) => [k, v.label]),
        state.dataset,
        (v) => {
          state.dataset = v;
          state.order = Object.keys(DATASETS[v].orders)[0];
          buildControls();
          refresh();
        },
      ),
      field("Search", { type: "search", placeholder: ds.search, class: "explorer__search" }, state.query, (v) => {
        state.query = v;
        refresh();
      }),
    );
    const order = el("label", { class: "explorer__ctl" });
    order.append(el("span", { class: "label" }, "Sort"));
    const select = el("select");
    for (const [k, [text]] of Object.entries(ds.orders)) {
      const option = el("option", { value: k }, esc(text));
      option.selected = k === state.order;
      select.append(option);
    }
    select.addEventListener("change", () => {
      state.order = select.value;
      refresh();
    });
    order.append(select);
    filters.replaceChildren(order);
    if (state.dataset === "dimers") {
      filters.append(
        seg(
          "Set",
          [
            ["all", "All"],
            ["homo", "Homo"],
            ["hetero", "Hetero"],
          ],
          state.kind,
          (v) => {
            state.kind = v;
            refresh();
          },
        ),
        field("Min ΔSASA (Å²)", { type: "number", min: "0", step: "100", inputmode: "decimal" }, state.minDelta, (v) => {
          state.minDelta = v;
          refresh();
        }),
        field("Min ipSAEmin", { type: "number", min: "0", max: "1", step: "0.1", inputmode: "decimal" }, state.minIpsae, (v) => {
          state.minIpsae = v;
          refresh();
        }),
      );
      const hq = el("label", { class: "explorer__check" });
      const box = el("input", { type: "checkbox" });
      box.checked = state.hq;
      box.addEventListener("change", () => {
        state.hq = box.checked;
        refresh();
      });
      hq.append(box, document.createTextNode(" AFDB high quality only"));
      filters.append(hq);
    }
  }

  let generation = 0;
  async function update() {
    const mine = ++generation;
    const ds = DATASETS[state.dataset];
    writeHash();
    const [clause, params] = where();
    try {
      await load(ds, status);
      if (mine !== generation) return;
      status.textContent = "Searching…";
      const summary = state.dataset === "dimers"
        ? "median(delta_sasa) AS median_value"
        : "median(sasa) AS median_value";
      // One connection runs one statement at a time
      const [stats] = await run(`SELECT count(*) AS n, ${summary} FROM ${ds.table} ${clause}`, params);
      const rows = await run(
        `SELECT * EXCLUDE (search_text) FROM ${ds.table} ${clause} ORDER BY ${ds.orders[state.order][1]} LIMIT ${PAGE_SIZE} OFFSET ${state.offset}`,
        params,
      );
      if (mine !== generation) return;
      const n = Number(stats.n);
      const median = state.dataset === "dimers" ? "median ΔSASA" : "median SASA";
      status.textContent = n
        ? `${n.toLocaleString("en-US")} ${n === 1 ? "model matches" : "models match"} · ${median} ${area(stats.median_value)} Å²`
        : "No model matches.";
      render(ds, rows);
      renderPager(ds, n, clause, params);
    } catch (err) {
      if (mine !== generation) return;
      status.textContent = `The query failed: ${err.message}`;
    }
  }

  function render(ds, rows) {
    const head = ds.columns.map(([name, , cls]) => `<th${cls ? ` class="${cls}"` : ""}>${esc(name)}</th>`).join("");
    const body = rows
      .map(
        (r, i) =>
          `<tr data-row="${i}"${ds.detail ? ' class="explorer__row" tabindex="0" aria-expanded="false"' : ""}>` +
          ds.columns.map(([, fmt, cls]) => `<td${cls ? ` class="${cls}"` : ""}>${fmt(r)}</td>`).join("") +
          "</tr>",
      )
      .join("");
    tableWrap.innerHTML = `<table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
    if (!ds.detail) return;
    const toggle = (tr) => {
      const next = tr.nextElementSibling;
      if (next?.classList.contains("explorer__open")) {
        next.remove();
        tr.setAttribute("aria-expanded", "false");
        return;
      }
      const detail = el("tr", { class: "explorer__open" });
      detail.append(el("td", { colspan: String(ds.columns.length) }, detailTable(rows[Number(tr.dataset.row)])));
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

  function renderPager(ds, n, clause, params) {
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
    const csv = el("button", { type: "button", class: "explorer__btn" }, `Download CSV (first ${Math.min(n, CSV_LIMIT).toLocaleString("en-US")})`);
    csv.addEventListener("click", async () => {
      csv.disabled = true;
      try {
        const rows = await run(
          `SELECT * EXCLUDE (search_text) FROM ${ds.table} ${clause} ORDER BY ${ds.orders[state.order][1]} LIMIT ${CSV_LIMIT}`,
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
      }
    });
    pager.append(
      prev,
      el("span", { class: "explorer__range" }, `${first.toLocaleString("en-US")}–${last.toLocaleString("en-US")} of ${n.toLocaleString("en-US")}`),
      next,
      csv,
    );
  }

  buildControls();
  try {
    conn = await connect(status);
    const probe = await fetch(new URL(DATASETS.dimers.file, DATA), { method: "HEAD" });
    if (!probe.ok) {
      status.textContent = "The data files are not part of this build of the site.";
      return;
    }
  } catch (err) {
    status.textContent = `DuckDB could not start in this browser: ${err.message}`;
    return;
  }
  // A link to another search on this page changes only the hash
  window.addEventListener("hashchange", () => {
    readHash();
    state.offset = 0;
    buildControls();
    update();
  });
  update();
}

main();
