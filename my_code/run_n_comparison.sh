#!/bin/bash

# Run n_comparison.m from command line without GUI,
# preventing figures from displaying.

nohup cat n_comparison.m | matlab -nojvm -nodisplay -nosplash -nodesktop &> ../logs/n_comparison.log &

#-r "set(0,'DefaultFigureVisible','off')" 
