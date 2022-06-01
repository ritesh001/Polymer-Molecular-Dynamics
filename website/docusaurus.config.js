// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require("prism-react-renderer/themes/github");
const darkCodeTheme = require("prism-react-renderer/themes/dracula");

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: "Polymer Molecular Dynamics",
  tagline:
    "Python toolkit and guides for molecular dynamics prediction of polymer properties",
  url: "https://polymer-molecular-dynamics.netlify.app/",
  baseUrl: "/",
  onBrokenLinks: "throw",
  onBrokenMarkdownLinks: "warn",
  favicon: "img/favicon.ico",
  organizationName: "Ramprasad Group",
  projectName: "Polymer Molecular Dynamics",

  presets: [
    [
      "classic",
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          path: "docs",
          routeBasePath: "docs",
          sidebarPath: require.resolve("./sidebars.js"),
          sidebarCollapsed: false,
          showLastUpdateAuthor: true,
          showLastUpdateTime: true,
          editUrl:
            "https://github.com/Ramprasad-Group/Polymer-Molecular-Dynamics/tree/main/website",
        },
        // blog: {
        //   showReadingTime: true,
        //   editUrl:
        //     "https://github.com/Ramprasad-Group/Polymer-Molecular-Dynamics/tree/main/website",
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
          "https://github.com/Ramprasad-Group/Polymer-Molecular-Dynamics/tree/main/website",
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
      image: "img/metaImg.png",
      navbar: {
        title: "Polymer Molecular Dynamics",
        logo: {
          alt: "Logo",
          src: "img/logo.svg",
        },
        items: [
          {
            type: "doc",
            docId: "intro",
            to: "/docs",
            position: "left",
            label: "Docs",
          },
          {
            to: "/api/overview",
            label: "API References",
            position: "left",
          },
          {
            href: "https://github.com/Ramprasad-Group/Polymer-Molecular-Dynamics",
            position: "right",
            className: "header-github-link",
            "aria-label": "GitHub repository",
          },
        ],
      },
      footer: {
        style: "dark",
        links: [
          {
            title: "Docs",
            items: [
              {
                label: "Introduction",
                to: "/docs/intro",
              },
              {
                label: "Getting Started",
                to: "/docs/getting-started/installation",
              },
            ],
          },
          {
            title: "Developer",
            items: [
              {
                label: "Kuan-Hsuan (Kevin) Shen",
                href: "https://github.com/kevinshen56714",
              },
            ],
          },
          {
            title: "More",
            items: [
              {
                label: "API References",
                to: "/api/overview",
              },
              {
                label: "GitHub",
                href: "https://github.com/Ramprasad-Group/Polymer-Molecular-Dynamics",
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
