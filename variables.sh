#!/bin/bash

# the home directory of git repository
dir_home=""

# the number of available threads
n_parallel=

# The name of input file (fq.gz)
input=""

# Location of minimap reference of TAIR10.
#ref=/path/to/TAIR10.mmi
ref=""

# Location of masksam.awk. 
awkloc=resources/masksam.awk

export dir_home
export n_parallel
export input
export ref
export awkloc
