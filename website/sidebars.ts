import type { SidebarsConfig } from "@docusaurus/plugin-content-docs";

const sidebars: SidebarsConfig = {
  docs: [
    "index",
    "getting-started",
    {
      type: "category",
      label: "User Guide",
      items: [
        "guide/choosing-tool",
        "guide/algorithms",
        "guide/classifiers",
        "guide/trajectory",
      ],
    },
    {
      type: "category",
      label: "CLI Reference",
      items: [
        "cli/commands",
        "cli/input",
        "cli/output",
      ],
    },
  ],
  pythonApi: [
    "python-api/index",
    "python-api/core",
    "python-api/classifier",
    "python-api/analysis",
    "python-api/xtc",
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
  ],
  benchmarks: [
    "benchmarks/single-file",
    "benchmarks/batch",
    "benchmarks/md",
    "benchmarks/validation",
  ],
};

export default sidebars;
