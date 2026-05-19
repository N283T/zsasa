import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
  docs: [
    "index",
    {
      type: "category",
      label: "Getting Started",
      items: ["getting-started"],
    },
    {
      type: "category",
      label: "Guides",
      items: [
        "guide/choosing-tool",
        "guide/batch",
        "guide/workflows",
        "guide/classifiers",
        "guide/trajectory",
        "guide/algorithms",
      ],
    },
    {
      type: "category",
      label: "Reference",
      items: [
        "cli/commands",
        "cli/input",
        "cli/output",
        {
          type: "category",
          label: "Python API",
          link: { type: "doc", id: "python-api/index" },
          items: [
            "python-api/core",
            "python-api/classifier",
            "python-api/analysis",
            "python-api/xtc",
          ],
        },
        {
          type: "category",
          label: "Integrations",
          items: [
            "integrations/index",
            "integrations/gemmi",
            "integrations/biopython",
            "integrations/biotite",
            "integrations/mdtraj",
            "integrations/mdanalysis",
          ],
        },
        "zig-api/autodoc",
      ],
    },
    {
      type: "category",
      label: "Benchmarks",
      items: [
        "benchmarks/index",
        "benchmarks/single-file",
        "benchmarks/batch",
        "benchmarks/md",
        "benchmarks/validation",
      ],
    },
    {
      type: "category",
      label: "Project",
      items: ["comparison", "changelog"],
    },
  ],
};

export default sidebars;
