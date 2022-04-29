# Website

This website is built using [Docusaurus 2](https://docusaurus.io/), a modern static website generator.

### Installation

```
$ yarn
```

### Local Development

```
$ yarn start
```

This command starts a local development server and opens up a browser window. Most changes are reflected live without having to restart the server.

### Build

```
$ yarn build
```

This command generates static content into the `build` directory and can be served using any static contents hosting service.

### Deployment

Using SSH:

```
$ USE_SSH=true yarn deploy
```

Not using SSH:

```
$ GIT_USER=<Your GitHub username> yarn deploy
```

If you are using GitHub pages for hosting, this command is a convenient way to build the website and push to the `gh-pages` branch.

# API doc

The API documentation is generated using [pydoc-markdown](https://github.com/NiklasRosenstein/pydoc-markdown)

### Installation

```
$ pip install pydoc-markdown[novella]
```

### Build

```
$ pydoc-markdown
```

Run `python cleanup_api.py` after `pydoc-markdown` to replace `####` with `###` in all the generated api markdown files. This is because pydoc-markdown make methods as h4 in markdown files, but I think methods deserve to be h3.
