// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require("prism-react-renderer/themes/github");
const darkCodeTheme = require("prism-react-renderer/themes/dracula");

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: "High-Throughput Polymer Molecular Dynamics Simulations",
  tagline:
    "Python toolkit and guides for in-silico prediction of polymer properties",
  url: "https://high-throughput-pmd.netlify.app/",
  baseUrl: "/",
  onBrokenLinks: "throw",
  onBrokenMarkdownLinks: "warn",
  favicon: "img/favicon.ico",
  organizationName: "Ramprasad Group",
  projectName: "High-Throughput Polymer MD Simulations",

  presets: [
    [
      "classic",
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          path: "guide",
          routeBasePath: "guide",
          sidebarPath: require.resolve("./sidebars.js"),
          sidebarCollapsed: false,
          editUrl:
            "https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/website",
        },
        // blog: {
        //   showReadingTime: true,
        //   editUrl:
        //     "https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/website",
        // },
        theme: {
          customCss: require.resolve("./src/css/custom.css"),
        },
      }),
    ],
  ],

  plugins: [
    [
      "content-docs",
      /** @type {import('@docusaurus/plugin-content-docs').Options} */
      ({
        id: "api",
        path: "api",
        routeBasePath: "api",
        editUrl:
          "https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/website",
        editCurrentVersion: true,
        sidebarPath: require.resolve("./sidebars.js"),
        sidebarCollapsed: false,
        showLastUpdateAuthor: true,
        showLastUpdateTime: true,
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: "High Throughput Polymer MD Simulations",
        logo: {
          alt: "My Site Logo",
          src: "img/logo.svg",
        },
        items: [
          {
            type: "doc",
            docId: "intro",
            to: "/guide",
            position: "left",
            label: "Guide",
          },
          {
            to: "/api/overview",
            label: "API Docs",
            position: "left",
          },
          {
            href: "https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations",
            label: "GitHub",
            position: "right",
          },
        ],
      },
      footer: {
        style: "dark",
        links: [
          {
            title: "Guide",
            items: [
              {
                label: "Guide",
                to: "/guide/intro",
              },
            ],
          },
          {
            title: "Community",
            items: [
              {
                label: "Stack Overflow",
                href: "https://stackoverflow.com/questions/tagged/docusaurus",
              },
              {
                label: "Discord",
                href: "https://discordapp.com/invite/docusaurus",
              },
              {
                label: "Twitter",
                href: "https://twitter.com/docusaurus",
              },
            ],
          },
          {
            title: "More",
            items: [
              {
                label: "API Docs",
                to: "/api/overview",
              },
              {
                label: "GitHub",
                href: "https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations",
              },
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} Ramprasad Group. Built with Docusaurus.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
