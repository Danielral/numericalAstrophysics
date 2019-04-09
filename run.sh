echo "Assignment1"

if [ ! -d "plots" ]; then
  mkdir plots
fi

pdflatex template.tex
bibtex template.aux
pdflatex template.tex
pdflatex template.tex
