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
                "    # Define polymer and system specs\n" +
                "    syst = pmd.System(smiles=smiles, force_field=pmd.OPLS(),\n" +
                "                      density=0.8, natoms_total=10000,\n" +
                "                      natoms_per_chain=150)\n" +
                "    \n" +
                "    # Customize LAMMPS simulation\n" +
                "    lmp = pmd.Lammps(read_data_from=syst, procedures=[\n" +
                "                     pmd.Minimization(min_style='cg'),\n" +
                "                     pmd.Equilibration(Teq=600, Tmax=800),\n" +
                "                     pmd.TgMeasurement(Tinit=600, Tfinal=200)])\n" +
                "    \n" +
                "    # Create job scheduler settings\n" +
                "    job = pmd.Torque(run_lammps=lmp, jobname=smiles,\n" +
                "                     project='Your-pid', nodes=2, ppn=24,\n" +
                "                     walltime='48:00:00')\n" +
                "    \n" +
                "    # Generate all necessary simulation files\n" +
                "    run = pmd.Pmd(system=syst, lammps=lmp, job=job)\n" +
                "    run.create(output_dir=smiles, save_config=True)"}
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
