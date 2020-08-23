# VORONOI Documentation

* To view Markdown documentation, navigate to `content/docs`.
* To view compiled (HTML) documentation, visit https://lanl.github.io/voronoi, or follow the build steps below.

# Building VORONOI Documentation

Static site generation (Markdown-to-HTML) is done through the [Hugo static site generator](https://gohugo.io).

* To start a Hugo server with VORONOI documentation, run: `hugo server -w`.
* To compile documentation with Hugo, run: `hugo -D -t kube -v`

Both of these commands must be run in the `docs/` directory.