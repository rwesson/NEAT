#!/bin/bash

pdflatex NEAT_paper_I.tex
bibtex NEAT_paper_I
pdflatex NEAT_paper_I.tex
pdflatex NEAT_paper_I.tex

#pdflatex NEAT_paper_II.tex
#bibtex NEAT_paper_II
#pdflatex NEAT_paper_II.tex
#pdflatex NEAT_paper_II.tex
