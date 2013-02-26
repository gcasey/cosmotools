#!/bin/bash
## This script is used within a HACC output directory to generate
## a file consisting of all the analysis time-steps. The script
## assumes that the simulation has dumped haloparticle tags out
dir=`pwd`
ls *.haloparticletags | cut -d'.' -f2 | sort -n > $dir/timesteps.dat
