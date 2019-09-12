# VORONOI Documentation

Full VORONOI documentation is provided in both HTML and Markdown formats,
with HTML being the preferred method.

* To view compiled documentation, open `index.html` in this directory.
* To view Markdown documentation, navigate to `hugo/content/docs`.

Static site generation (Markdown-to-HTML) is done through the [Hugo static site generator](https://gohugo.io).

* To start a Hugo server with VORONOI documentation, run: `hugo server -w`.
* To compile documentation with Hugo, run: `hugo -D -t kube -v`

Both of these commands must be run in the `hugo/` directory.