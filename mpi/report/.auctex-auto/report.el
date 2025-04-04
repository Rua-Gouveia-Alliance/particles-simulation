;; -*- lexical-binding: t; -*-

(TeX-add-style-hook
 "report"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "paper=a4" "fontsize=11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("fontenc" "T1") ("fourier" "") ("babel" "english") ("microtype" "protrusion=true" "expansion=true") ("amsmath" "fleqn") ("amsfonts" "") ("amsthm" "") ("graphicx" "pdftex") ("url" "") ("sectsty" "") ("fancyhdr" "") ("hyperref" "") ("listings" "") ("pgfplots" "")))
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "fontenc"
    "fourier"
    "babel"
    "microtype"
    "amsmath"
    "amsfonts"
    "amsthm"
    "graphicx"
    "url"
    "sectsty"
    "fancyhdr"
    "hyperref"
    "listings"
    "pgfplots"))
 :latex)

