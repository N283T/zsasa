import type { ReactNode } from "react";
import clsx from "clsx";
import Link from "@docusaurus/Link";
import useDocusaurusContext from "@docusaurus/useDocusaurusContext";
import Layout from "@theme/Layout";
import Heading from "@theme/Heading";

import styles from "./index.module.css";

function HomepageHeader() {
  const { siteConfig } = useDocusaurusContext();
  return (
    <header className={clsx("hero hero--primary", styles.heroBanner)}>
      <div className="container">
        <Heading as="h1" className="hero__title">
          {siteConfig.title}
        </Heading>
        <p className="hero__subtitle">{siteConfig.tagline}</p>
        <div className={styles.buttons}>
          <Link
            className="button button--secondary button--lg"
            to="/docs/getting-started"
          >
            Getting Started
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/guide/choosing-tool"
          >
            User Guide
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/cli/commands"
          >
            CLI Reference
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/python-api"
          >
            Python API
          </Link>
          <Link
            className="button button--secondary button--lg"
            to="/docs/benchmarks/single-file"
          >
            Benchmarks
          </Link>
        </div>
      </div>
    </header>
  );
}

export default function Home(): ReactNode {
  const { siteConfig } = useDocusaurusContext();
  return (
    <Layout title="Home" description={siteConfig.tagline}>
      <HomepageHeader />
      <main>
        <section style={{ padding: "2rem 0" }}>
          <div className="container">
            <div className="row">
              <div className="col col--4">
                <Heading as="h2">Fast</Heading>
                <p>
                  Up to 3x faster than FreeSASA C with f64 precision. SIMD
                  optimization and multi-threading for maximum throughput.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h2">Two Algorithms</Heading>
                <p>
                  Shrake-Rupley (fast, recommended) and Lee-Richards (precise).
                  Bitmask LUT optimization for additional speedup.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h2">Multiple Formats</Heading>
                <p>
                  mmCIF, PDB, and JSON input. XTC and DCD trajectory support
                  with automatic unit conversion.
                </p>
              </div>
            </div>
            <div className="row" style={{ marginTop: "2rem" }}>
              <div className="col col--4">
                <Heading as="h2">Python Bindings</Heading>
                <p>
                  NumPy integration with pre-built wheels for Linux, macOS, and
                  Windows. Gemmi, BioPython, Biotite, MDTraj, and MDAnalysis support.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h2">Analysis Tools</Heading>
                <p>
                  Per-residue aggregation, RSA calculation, and polar/nonpolar
                  classification. Three built-in atom classifiers.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h2">Cross-Platform</Heading>
                <p>
                  Linux, macOS, and Windows. CLI binary, Python package on PyPI,
                  and native Zig library with interactive autodoc.
                </p>
              </div>
            </div>
          </div>
        </section>
      </main>
    </Layout>
  );
}
