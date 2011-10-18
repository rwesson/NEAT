#!/bin/bash

pdflatex NEAT_paper2.tex
bibtex NEAT_paper2
pdflatex NEAT_paper2.tex
pdflatex NEAT_paper2.tex
