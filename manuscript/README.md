Build the manuscript with latexmk (already part of TeXLive):

```
latexmk            # build main.pdf (runs pdflatex + biber as needed)
latexmk -pvc       # live preview: rebuild on every save
latexmk -c         # remove auxiliary files (keep main.pdf)
latexmk -C         # also remove main.pdf
```

Build configuration lives in `latexmkrc`.
