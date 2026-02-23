import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
  docs: [
    "index",
    "getting-started",
    {
      type: "category",
      label: "User Guide",
      items: [
        "guide/algorithms",
        "guide/classifiers",
        "guide/trajectory",
      ],
    },
    "cli",
    {
      type: "category",
      label: "Python API",
      items: [
        "python-api/index",
        "python-api/core",
        "python-api/classifier",
        "python-api/analysis",
        "python-api/xtc",
        {
          type: "category",
          label: "Integrations",
          items: [
            "python-api/integrations/index",
            "python-api/integrations/gemmi",
            "python-api/integrations/biopython",
            "python-api/integrations/biotite",
            "python-api/integrations/mdtraj",
            "python-api/integrations/mdanalysis",
          ],
        },
      ],
    },
    {
      type: "category",
      label: "Zig API",
      items: ["zig-api/autodoc"],
    },
    {
      type: "category",
      label: "Benchmarks",
      items: [
        "benchmarks/single-file",
        "benchmarks/batch",
        "benchmarks/md",
        "benchmarks/validation",
      ],
    },
  ],
};

export default sidebars;
