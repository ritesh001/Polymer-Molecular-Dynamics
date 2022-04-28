import React from "react";
import clsx from "clsx";
import Layout from "@theme/Layout";
import Link from "@docusaurus/Link";
import useDocusaurusContext from "@docusaurus/useDocusaurusContext";
import styles from "./index.module.css";
import HomepageFeatures from "@site/src/components/HomepageFeatures";
import CodeBlock from "@theme/CodeBlock";

function HomepageHeader() {
  const { siteConfig } = useDocusaurusContext();
  return (
    <header className={clsx("hero hero--primary", styles.heroBanner)}>
      <div className="container">
        <div
          className="row"
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-around",
            gap: "2em",
          }}
        >
          <div style={{ textAlign: "left", width: "35em", maxWidth: "95%" }}>
            <h1 className="hero__title">{siteConfig.title}</h1>
            <p className="hero__subtitle">{siteConfig.tagline}</p>
            {/* <div className={styles.buttons}>
              <div style={{ width: 240 }}>
                <CodeBlock className="language-bash">
                  {"pip install pmd"}
                </CodeBlock>
              </div>
            </div> */}
            <div style={{ display: "flex", gap: 30 }}>
              <Link
                className="button button--secondary button--lg"
                to="/docs/intro"
                style={{ width: 250 }}
              >
                Get Started
              </Link>
              <Link
                className="button button--secondary button--lg"
                to="/docs/intro"
                style={{
                  width: 250,
                  backgroundColor: "#003057",
                  borderColor: "#003057",
                  color: "white",
                }}
              >
                Examples
              </Link>
            </div>
          </div>
          <div style={{ textAlign: "left" }}>
            <CodeBlock className="language-js" title="example-mkinput.py">
              {"import pmd\n" +
                "\n" +
                "smiles = ['*CC*', '*CC(*)CC','*CC(*)c1ccccc1']\n" +
                "\n" +
                "for s in smiles:\n" +
                "    lmp = pmd.Lammps(s, force_field='gaff2')\n" +
                "    lmp.add_procedure(pmd.Minimization())\n" +
                "    lmp.add_procedure(pmd.Equilibration())\n" +
                "    lmp.add_procedure(pmd.TgMeasurement())\n" +
                "    lmp.write_input(output_dir=s)\n" +
                "\n" +
                "    job = pmd.Job(jobname=s,\n" +
                "                  project='Your-project-ID',\n" +
                "                  nodes=2, ppn=24,\n" +
                "                  walltime='48:00:00')\n" +
                "    job.write_pbs(output_dir=s)"}
            </CodeBlock>
          </div>
        </div>
      </div>
    </header>
  );
}

export default function Home() {
  const { siteConfig } = useDocusaurusContext();
  return (
    <Layout
      title={`${siteConfig.title}`}
      description="Description will go into a meta tag in <head />"
    >
      <HomepageHeader />
      <main>
        <HomepageFeatures />
      </main>
    </Layout>
  );
}
