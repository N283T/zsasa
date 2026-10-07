# /// script
# requires-python = ">=3.11"
# dependencies = ["markdown-it-py>=3.0", "pygments>=2.17"]
# ///
"""Static site builder for the zsasa website.

Renders the hand-written landing page and the Markdown docs into ``website/dist``.

    uv run website/build.py

The docs are written in the Docusaurus dialect of Markdown; ``preprocess``
translates the handful of constructs they use (admonitions, tabs, heading ids,
JSX attributes) so the sources do not need to change.
"""

from __future__ import annotations

import html
import json
import posixpath
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

from markdown_it import MarkdownIt
from markdown_it.token import Token
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import TextLexer, get_lexer_by_name
from pygments.util import ClassNotFound

WEB = Path(__file__).resolve().parent
REPO = WEB.parent
DOCS = WEB / "docs"
STATIC = WEB / "static"
SRC = WEB / "src"
TEMPLATES = WEB / "templates"
DIST = WEB / "dist"
DATA = WEB / "data" / "benchmarks"

SITE_URL = "https://n283t.github.io/zsasa/"
BASE_PATH = (
    "/zsasa/"  # prefix the deployed site lives under; hard-coded in some doc links
)
EDIT = "https://github.com/N283T/zsasa/tree/main/website/docs"

# Sidebar: (group label, [doc id | (label, [doc id, ...])]).
NAV: list[tuple[str, list]] = [
    ("Start", ["index", "getting-started"]),
    (
        "Guides",
        [
            "guide/choosing-tool",
            "guide/batch",
            "guide/workflows",
            "guide/classifiers",
            "guide/trajectory",
            "guide/algorithms",
        ],
    ),
    (
        "Reference",
        [
            "cli/commands",
            "cli/input",
            "cli/output",
            (
                "Python API",
                [
                    "python-api/index",
                    "python-api/core",
                    "python-api/classifier",
                    "python-api/analysis",
                    "python-api/xtc",
                ],
            ),
            (
                "Integrations",
                [
                    "integrations/index",
                    "integrations/gemmi",
                    "integrations/biopython",
                    "integrations/biotite",
                    "integrations/mdtraj",
                    "integrations/mdanalysis",
                ],
            ),
            "zig-api/autodoc",
        ],
    ),
    (
        "Benchmarks",
        [
            "benchmarks/index",
            "benchmarks/validation",
            "benchmarks/batch",
            "benchmarks/md",
            "benchmarks/single-file",
        ],
    ),
    ("Project", ["citation", "comparison", "changelog"]),
]
LABELS = {"index": "Introduction"}

FRONTMATTER = re.compile(r"\A---\n(.*?)\n---\n", re.S)
MDX_IMPORT = re.compile(r"^import .+ from ['\"]@theme/.+['\"];?\n", re.M)
ADMONITION = re.compile(
    r"^:::(note|tip|info|warning|caution|danger)(?:\[(.*?)\]|[ \t]+(.+?))?[ \t]*\n(.*?)\n:::[ \t]*$",
    re.M | re.S,
)
TABS = re.compile(r"^<Tabs[^>]*>\n(.*?)^</Tabs>[ \t]*$", re.M | re.S)
TAB_ITEM = re.compile(
    r"^[ \t]*<TabItem\s+([^>]*)>\n(.*?)^[ \t]*</TabItem>[ \t]*$", re.M | re.S
)
HEADING_ID = re.compile(r"\s*\{#([\w-]+)\}\s*$")
CODE_INCLUDE = re.compile(r"<!--code:([\w-]+)-->\n(.*?)<!--/code-->", re.S)
PLACEHOLDER = re.compile(r"\{\{([\w:/.-]+)\}\}")
CHART = re.compile(r'<div data-chart="([\w-]+)/([\w-]+)"></div>')
TABLE = re.compile(r'<div data-table="([\w-]+)/([\w-]+)"></div>')
CHARTS_SCRIPT = '<script src="{root}assets/charts.js" defer></script>'


@dataclass
class Page:
    doc_id: str
    source: Path
    group: str
    label: str = ""
    title: str = ""
    description: str = ""
    body: str = ""
    toc: list[tuple[str, str, str]] = field(default_factory=list)
    anchors: set[str] = field(default_factory=set)
    # (target doc id, fragment, original href) for every link to a heading.
    refs: list[tuple[str, str, str]] = field(default_factory=list)

    @property
    def path(self) -> str:
        """Site-relative URL path, matching the URLs the Docusaurus site used."""
        if self.doc_id == "index":
            return "docs/"
        return f"docs/{self.doc_id.removesuffix('/index')}/"

    @property
    def root(self) -> str:
        return "../" * self.path.count("/")


def flat_nav() -> list[tuple[str, str]]:
    """(doc id, group label) in reading order."""
    order = []
    for group, items in NAV:
        for item in items:
            for doc_id in [item] if isinstance(item, str) else item[1]:
                order.append((doc_id, group))
    return order


def find_source(doc_id: str) -> Path:
    for ext in (".md", ".mdx"):
        path = DOCS / f"{doc_id}{ext}"
        if path.exists():
            return path
    raise SystemExit(f"NAV lists '{doc_id}' but website/docs has no such file")


def slugify(text: str) -> str:
    """GitHub-style heading slug, as Docusaurus generates."""
    return re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")


def render_code(code: str, lang: str) -> str:
    try:
        lexer = get_lexer_by_name(lang)
    except ClassNotFound:
        lexer = TextLexer()
    body = highlight(code, lexer, HtmlFormatter(nowrap=True)).rstrip("\n")
    return (
        '<div class="code">'
        f'<div class="code__lang"><span>{html.escape(lang)}</span>'
        '<button class="copy" type="button" data-copy>Copy</button></div>'
        f"<pre><code>{body}</code></pre></div>\n"
    )


def make_markdown() -> MarkdownIt:
    md = MarkdownIt("commonmark", {"html": True}).enable(["table", "strikethrough"])

    def fence(self, tokens, idx, options, env):
        token = tokens[idx]
        return render_code(token.content, (token.info.strip().split() or ["text"])[0])

    md.add_render_rule("fence", fence)
    md.add_render_rule("table_open", lambda *_: '<div class="table"><table>\n')
    md.add_render_rule("table_close", lambda *_: "</table></div>\n")
    return md


def preprocess(text: str) -> str:
    """Translate the Docusaurus-flavoured syntax the docs use into plain HTML + Markdown."""
    text = FRONTMATTER.sub("", text)
    text = MDX_IMPORT.sub("", text)
    text = text.replace('className="', 'class="')

    def admonition(m: re.Match) -> str:
        kind = m.group(1)
        title = m.group(2) or m.group(3) or kind
        return (
            f'<div class="note note--{kind}">\n'
            f'<p class="note__t">{html.escape(title)}</p>\n\n{m.group(4)}\n\n</div>'
        )

    def tabs(m: re.Match) -> str:
        items = TAB_ITEM.findall(m.group(1))
        labels = [re.search(r'label="([^"]*)"', attrs).group(1) for attrs, _ in items]
        buttons = "".join(
            f'<button type="button" role="tab" aria-selected="{str(i == 0).lower()}"'
            f"{'' if i == 0 else ' tabindex=-1'}>{html.escape(label)}</button>"
            for i, label in enumerate(labels)
        )
        out = [
            f'<div class="tabs" data-tabs>\n<div class="tabs__list" role="tablist">{buttons}</div>'
        ]
        for i, (_, body) in enumerate(items):
            hidden = "" if i == 0 else " hidden"
            out.append(
                f'<div class="tabs__panel" role="tabpanel"{hidden}>\n\n{body.strip()}\n\n</div>'
            )
        return "\n".join(out) + "\n</div>"

    return TABS.sub(tabs, ADMONITION.sub(admonition, text))


def inline_text(inline: Token) -> str:
    return "".join(
        c.content for c in inline.children or [] if c.type in ("text", "code_inline")
    ).strip()


class Site:
    def __init__(self) -> None:
        self.md = make_markdown()
        self.pages: dict[str, Page] = {}
        for doc_id, group in flat_nav():
            self.pages[doc_id] = Page(doc_id, find_source(doc_id), group)
        listed = {p.source for p in self.pages.values()}
        for path in sorted(DOCS.rglob("*.md*")):
            if path not in listed:
                print(
                    f"warning: {path.relative_to(REPO)} is not in NAV and will not be built",
                    file=sys.stderr,
                )
        self.version = re.search(
            r'\.version = "([^"]+)"', (REPO / "build.zig.zon").read_text()
        ).group(1)
        self.errors: list[str] = []
        # Benchmark charts and tables, exported by scripts/export_benchmarks.py.
        self.bench = {
            path.stem: json.loads(path.read_text())
            for path in sorted(DATA.glob("*.json"))
        }

    # ----- links -----

    def resolve(self, href: str, page: Page) -> str:
        """Rewrite a link or image URL found in a doc to a path relative to that page."""
        if href.startswith("pathname://"):
            href = href.removeprefix("pathname://")
        elif re.match(r"^[a-z][a-z0-9+.-]*:|^//", href):
            return href
        target, _, frag = href.partition("#")
        anchor = f"#{frag}" if frag else ""
        if not target:
            page.refs.append((page.doc_id, frag, href))
            return href
        if target.startswith(BASE_PATH):
            return page.root + target.removeprefix(BASE_PATH) + anchor
        if target.startswith("/docs/"):
            doc_id = target.removeprefix("/docs/").rstrip("/")
        elif target.startswith("/"):
            if not (STATIC / target.lstrip("/")).exists():
                self.errors.append(f"{page.doc_id}: missing static file {href}")
            return page.root + target.lstrip("/") + anchor
        else:
            base = posixpath.dirname(page.doc_id)
            doc_id = posixpath.normpath(
                posixpath.join(base, re.sub(r"\.mdx?$", "", target))
            ).rstrip("/")
        for candidate in (doc_id, f"{doc_id}/index"):
            if candidate in self.pages:
                if frag:
                    page.refs.append((candidate, frag, href))
                return page.root + self.pages[candidate].path + anchor
        self.errors.append(f"{page.doc_id}: broken link {href}")
        return href

    # ----- docs -----

    def parse(self, page: Page) -> None:
        text = page.source.read_text()
        meta = m.group(1) if (m := FRONTMATTER.match(text)) else ""
        tokens = self.md.parse(preprocess(text), {})

        for i, token in enumerate(tokens):
            if token.type == "heading_open":
                inline = tokens[i + 1]
                custom = None
                if inline.children and inline.children[-1].type == "text":
                    last = inline.children[-1]
                    if m := HEADING_ID.search(last.content):
                        custom, last.content = m.group(1), last.content[: m.start()]
                text_ = inline_text(inline)
                if token.tag == "h1":
                    page.title = page.title or text_
                    continue
                slug = base = custom or slugify(text_)
                n = 0
                while slug in page.anchors:
                    n += 1
                    slug = f"{base}-{n}"
                page.anchors.add(slug)
                token.attrSet("id", slug)
                if token.tag in ("h2", "h3"):
                    page.toc.append((token.tag, slug, text_))
                link = Token("html_inline", "", 0)
                link.content = f'<a class="header-anchor" href="#{slug}" aria-label="Link to this section">#</a>'
                inline.children.append(link)
            elif token.type == "paragraph_open" and page.title and not page.description:
                page.description = inline_text(tokens[i + 1])[:200]
            for child in token.children or []:
                if child.type == "link_open":
                    child.attrSet("href", self.resolve(child.attrGet("href"), page))
                elif child.type == "image":
                    child.attrSet("src", self.resolve(child.attrGet("src"), page))

        body = self.md.renderer.render(tokens, self.md.options, {})
        # Raw HTML in the docs hard-codes the deployed base path.
        page.body = self.expand(
            re.sub(rf'(src|href)="{BASE_PATH}', rf'\1="{page.root}', body)
        )
        page.anchors |= set(re.findall(r'\bid="([^"]+)"', page.body))

        if m := re.search(r"^title:\s*[\"']?(.+?)[\"']?\s*$", meta, re.M):
            page.title = m.group(1)
        label = re.search(r"^sidebar_label:\s*[\"']?(.+?)[\"']?\s*$", meta, re.M)
        page.label = LABELS.get(page.doc_id) or (
            label.group(1) if label else page.title
        )

    def figure(self, m: re.Match) -> str:
        """Expand ``<div data-chart="page/id"></div>`` into a chart figure with its inline spec."""
        try:
            spec = self.bench[m.group(1)]["charts"][m.group(2)]
        except KeyError:
            self.errors.append(f"unknown chart {m.group(1)}/{m.group(2)}")
            return ""
        payload = json.dumps(spec, ensure_ascii=False, separators=(",", ":")).replace(
            "</", "<\\/"
        )
        return (
            '<figure class="chart">'
            f'<div class="chart__head"><p class="chart__title">{html.escape(spec["title"])}</p>'
            '<div class="chart__controls"></div></div>'
            '<div class="chart__legend"></div>'
            '<div class="chart__body"><noscript>This chart needs JavaScript. '
            "The tables on this page carry the headline values.</noscript></div>"
            '<p class="chart__stats"></p>'
            f'<figcaption class="chart__note">{html.escape(spec.get("note", ""))}</figcaption>'
            '<details class="chart__data"><summary>Data table</summary><div class="table"></div></details>'
            f'<script type="application/json">{payload}</script>'
            "</figure>"
        )

    def table(self, m: re.Match) -> str:
        """Expand ``<div data-table="page/id"></div>`` into a static table."""
        try:
            spec = self.bench[m.group(1)]["tables"][m.group(2)]
        except KeyError:
            self.errors.append(f"unknown table {m.group(1)}/{m.group(2)}")
            return ""
        head = "".join(f"<th>{html.escape(c)}</th>" for c in spec["columns"])
        body = "".join(
            "<tr>"
            + "".join(
                f"<td>{'<strong>' + html.escape(c) + '</strong>' if i == 0 and c.startswith('zsasa') else html.escape(c)}</td>"
                for i, c in enumerate(row)
            )
            + "</tr>"
            for row in spec["rows"]
        )
        return (
            f'<div class="table"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
            f'<p class="table-cap">{html.escape(spec["caption"])}</p>'
        )

    def expand(self, text: str) -> str:
        return TABLE.sub(self.table, CHART.sub(self.figure, text))

    def sidebar(self, current: Page) -> str:
        def link(doc_id: str, label: str | None = None) -> str:
            page = self.pages[doc_id]
            cur = ' aria-current="page"' if page is current else ""
            return f'<a href="{current.root}{page.path}"{cur}>{html.escape(label or page.label)}</a>'

        out = []
        for group, items in NAV:
            out.append(
                f'<div class="side__group"><span class="label">{group}</span><ul>'
            )
            for item in items:
                if isinstance(item, str):
                    out.append(f"<li>{link(item)}</li>")
                else:
                    label, (head, *rest) = item
                    out.append(f"<li>{link(head, label)}<ul>")
                    out.extend(f"<li>{link(child)}</li>" for child in rest)
                    out.append("</ul></li>")
            out.append("</ul></div>")
        return "\n".join(out)

    def fill(self, template: str, values: dict[str, str], root: str) -> str:
        def sub(m: re.Match) -> str:
            key = m.group(1)
            if key.startswith("doc:"):
                return root + self.pages[key[4:]].path
            if key.startswith("fact:"):
                return self.bench["batch"]["facts"][key[5:]]
            return values[key]

        return PLACEHOLDER.sub(sub, template)

    def chrome(self, values: dict[str, str], root: str) -> dict[str, str]:
        """Shared head / nav / footer partials, filled for one page."""
        return {
            name: self.fill((TEMPLATES / f"{name}.html").read_text(), values, root)
            for name in ("head", "nav", "footer")
        }

    def write_doc(self, page: Page) -> None:
        order = [doc_id for doc_id, _ in flat_nav()]
        index = order.index(page.doc_id)
        pager = []
        for cls, word, j in (
            ("prev", "Previous", index - 1),
            ("next", "Next", index + 1),
        ):
            if 0 <= j < len(order):
                other = self.pages[order[j]]
                pager.append(
                    f'<a class="{cls}" href="{page.root}{other.path}"><span class="label">{word}</span>'
                    f"<b>{html.escape(other.title)}</b></a>"
                )

        values = {
            "root": page.root,
            "version": self.version,
            "title": f"{page.label if page.doc_id == 'index' else page.title} · zsasa",
            "desc": html.escape(
                page.description or f"{page.title}: zsasa documentation."
            ),
            "docs_current": ' aria-current="page"',
            "menu_hidden": "",
        }
        values |= self.chrome(values, page.root)
        values |= {
            "crumbs": f"<span>{page.group}</span><span>{html.escape(page.title)}</span>",
            "sidebar": self.sidebar(page),
            "content": page.body,
            "toc": "\n".join(
                f'<li><a class="d{tag[1]}" href="#{anchor}">{html.escape(text)}</a></li>'
                for tag, anchor, text in page.toc
            ),
            "pager": "".join(pager),
            "edit": f"{EDIT}/{page.source.relative_to(DOCS).as_posix()}",
            "scripts": CHARTS_SCRIPT.format(root=page.root)
            if 'class="chart"' in page.body
            else "",
        }
        out_dir = DIST / page.path
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "index.html").write_text(
            self.fill((TEMPLATES / "doc.html").read_text(), values, page.root)
        )

    def write_landing(self) -> None:
        root = "./"
        values = {
            "root": root,
            "version": self.version,
            "title": "zsasa: solvent accessible surface area at proteome scale",
            "desc": "High-performance SASA calculation in Zig, with a CLI and Python bindings.",
            "docs_current": "",
            "menu_hidden": " hidden",
        }
        values |= self.chrome(values, root)
        page = (SRC / "index.html").read_text()
        page = self.expand(
            CODE_INCLUDE.sub(lambda m: render_code(m.group(2), m.group(1)), page)
        )
        (DIST / "index.html").write_text(self.fill(page, values, root))

    def write_extras(self) -> None:
        """404 page (absolute paths: it is served from any depth) and sitemap."""
        values = {
            "root": BASE_PATH,
            "version": self.version,
            "title": "Page not found · zsasa",
            "desc": "Page not found.",
            "docs_current": "",
            "menu_hidden": " hidden",
        }
        values |= self.chrome(values, BASE_PATH)
        page = self.fill((SRC / "404.html").read_text(), values, BASE_PATH)
        (DIST / "404.html").write_text(page)

        urls = [SITE_URL] + [SITE_URL + page.path for page in self.pages.values()]
        entries = "".join(f"<url><loc>{html.escape(url)}</loc></url>" for url in urls)
        (DIST / "sitemap.xml").write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            f'<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">{entries}</urlset>\n'
        )

    def build(self) -> None:
        if DIST.exists():
            shutil.rmtree(DIST)
        shutil.copytree(STATIC, DIST)
        shutil.copytree(SRC / "assets", DIST / "assets", dirs_exist_ok=True)

        for page in self.pages.values():
            self.parse(page)
        for page in self.pages.values():
            for target, frag, href in page.refs:
                if frag not in self.pages[target].anchors:
                    self.errors.append(
                        f"{page.doc_id}: no heading '#{frag}' for link {href}"
                    )
        if self.errors:
            raise SystemExit("broken links:\n  " + "\n  ".join(self.errors))

        self.write_landing()
        if self.errors:
            raise SystemExit("landing page:\n  " + "\n  ".join(self.errors))
        for page in self.pages.values():
            self.write_doc(page)
        self.write_extras()
        print(f"built {1 + len(self.pages)} pages -> {DIST.relative_to(REPO)}")


if __name__ == "__main__":
    Site().build()
