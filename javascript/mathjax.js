window.MathJax = {
  tex: {
    tags: 'all',
    inlineMath: [["$","$"],["\\(", "\\)"]],
    displayMath: [["$$","$$"],["\\[", "\\]"]],
    processEscapes: true,
    processEnvironments: true
  },
  options: {
    ignoreHtmlClass: ".*|",
    processHtmlClass: "arithmatex"
  }
};

document$.subscribe(() => {
  MathJax.startup.output.clearCache()
  MathJax.typesetClear()
  MathJax.texReset()
  MathJax.typesetPromise()
})
