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
        {/* Install section */}
        <section className={styles.install}>
          <div className="container">
            <div className="row">
              <div className="col col--6">
                <Heading as="h3">Python</Heading>
                <pre className={styles.installCode}>
                  <code>pip install zsasa</code>
                </pre>
              </div>
              <div className="col col--6">
                <Heading as="h3">CLI</Heading>
                <pre className={styles.installCode}>
                  <code>curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh</code>
                </pre>
              </div>
            </div>
          </div>
        </section>

        {/* Highlight section */}
        <section className={styles.highlight}>
          <div className="container">
            <div className="row">
              <div className="col col--6">
                <Heading as="h2">Up to 3x Faster than FreeSASA C</Heading>
                <p>
                  SIMD-optimized Shrake-Rupley with multi-threading. Bitmask LUT
                  optimization for additional speedup. Consistent performance
                  with no pathological slowdowns.
                </p>
                <Link to="/docs/comparison">
                  See how zsasa compares →
                </Link>
              </div>
              <div className="col col--6">
                <Heading as="h2">Zero Dependencies</Heading>
                <p>
                  Pure Zig with no external libraries. Single{" "}
                  <code>zig build</code> command. Python wheels available on PyPI
                  for Linux, macOS, and Windows.
                </p>
                <Link to="/docs/getting-started">
                  Install in 30 seconds →
                </Link>
              </div>
            </div>
          </div>
        </section>

        {/* Feature grid */}
        <section style={{ padding: "2rem 0" }}>
          <div className="container">
            <div className="row">
              <div className="col col--4">
                <Heading as="h3">Two Algorithms</Heading>
                <p>
                  Shrake-Rupley (fast, recommended) and Lee-Richards (precise).
                  Selectable f64/f32 precision.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h3">Multiple Formats</Heading>
                <p>
                  mmCIF, PDB, and JSON input. XTC and DCD trajectory support
                  with automatic unit conversion.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h3">Python Bindings</Heading>
                <p>
                  NumPy integration with Gemmi, BioPython, Biotite, MDTraj, and
                  MDAnalysis support.
                </p>
              </div>
            </div>
            <div className="row" style={{ marginTop: "1.5rem" }}>
              <div className="col col--4">
                <Heading as="h3">Analysis Tools</Heading>
                <p>
                  Per-residue aggregation, RSA calculation, and polar/nonpolar
                  classification with three built-in classifiers.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h3">Batch & Trajectory</Heading>
                <p>
                  Native directory batch processing and MD trajectory analysis.
                  Proteome-scale datasets in seconds.
                </p>
              </div>
              <div className="col col--4">
                <Heading as="h3">Cross-Platform</Heading>
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
