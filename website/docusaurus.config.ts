import type { Config } from "@docusaurus/types";
import type * as Preset from "@docusaurus/preset-classic";

const config: Config = {
  title: "zsasa",
  tagline: "High-performance SASA calculation in Zig",
  favicon: "img/favicon.ico",

  url: "https://n283t.github.io",
  baseUrl: "/zsasa/",

  organizationName: "N283T",
  projectName: "zsasa",

  onBrokenLinks: "throw",
  onBrokenMarkdownLinks: "warn",
  onBrokenAnchors: "warn",

  markdown: {
    format: "detect",
  },

  i18n: {
    defaultLocale: "en",
    locales: ["en"],
  },

  presets: [
    [
      "classic",
      {
        docs: {
          routeBasePath: "/",
          sidebarPath: "./sidebars.ts",
          editUrl: "https://github.com/N283T/zsasa/tree/main/website/",
        },
        pages: false,
        blog: false,
        theme: {
          customCss: "./src/css/custom.css",
        },
      } satisfies Preset.Options,
    ],
  ],

  themeConfig: {
    navbar: {
      title: "zsasa",
      items: [
        {
          type: "docSidebar",
          sidebarId: "docs",
          position: "left",
          label: "Docs",
        },
        {
          href: "https://github.com/N283T/zsasa",
          label: "GitHub",
          position: "right",
        },
      ],
    },
    footer: {
      style: "dark",
      links: [
        {
          title: "Docs",
          items: [
            { label: "Getting Started", to: "/getting-started" },
            { label: "Python API", to: "/python-api" },
            { label: "CLI Reference", to: "/cli" },
          ],
        },
        {
          title: "More",
          items: [
            { label: "GitHub", href: "https://github.com/N283T/zsasa" },
            { label: "PyPI", href: "https://pypi.org/project/zsasa/" },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} zsasa contributors. MIT License.`,
    },
    prism: {
      additionalLanguages: ["bash", "python", "toml", "csv"],
    },
    colorMode: {
      defaultMode: "light",
      respectPrefersColorScheme: true,
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
