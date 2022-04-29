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
          <div style={{ textAlign: "left", width: "35em", maxWidth: "100%" }}>
            <h1 className="hero__title">{siteConfig.title}</h1>
            <p className="hero__subtitle">{siteConfig.tagline}</p>
            {/* <div className={styles.buttons}>
              <div style={{ width: 240 }}>
                <CodeBlock className="language-bash">
                  {"pip install pmd"}
                </CodeBlock>
              </div>
            </div> */}
            <div style={{ display: "flex", gap: 20, flexWrap: "wrap" }}>
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
          <div style={{ textAlign: "left", maxWidth: "100%" }}>
            <CodeBlock className="language-python" title="example-mkinput.py">
              {"import pmd\n" +
                "\n" +
                "for smiles in ['*CC*', '*CC(*)CC','*CC(*)c1ccccc1']:\n" +
                "    # Define system specs and make the data file\n" +
                "    s = pmd.System(smiles=smiles, force_field='opls',\n" +
                "                   density=0.8, natoms_total=5000,\n" +
                "                   natoms_per_chain=150)\n" +
                "    s.write_data(output_dir=smiles)\n" +
                "\n" +
                "    # Customize LAMMPS simulation and make the input file\n" +
                "    lmp = pmd.Lammps(s)\n" +
                "    lmp.add_procedure(pmd.Minimization())\n" +
                "    lmp.add_procedure(pmd.Equilibration())\n" +
                "    lmp.add_procedure(pmd.TgMeasurement())\n" +
                "    lmp.write_input(output_dir=smiles)\n" +
                "\n" +
                "    # Create job scheduler file\n" +
                "    job = pmd.Job(jobname=smiles, project='Your-pid',\n" +
                "                  nodes=2, ppn=24, walltime='48:00:00')\n" +
                "    job.write_pbs(output_dir=smiles)"}
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
      description="Python toolkit for molecular dynamics prediction of polymer properties"
    >
      <HomepageHeader />
      <main>
        <HomepageFeatures />
      </main>
    </Layout>
  );
}
