#!/bin/bash

pdflatex NEAT_paper.tex
bibtex NEAT_paper
pdflatex NEAT_paper.tex
pdflatex NEAT_paper.tex
