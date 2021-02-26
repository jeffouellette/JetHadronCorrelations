#! /bin/bash

declare -a vardirs=("Nominal" "JetES5PercUpVar" "JetES5PercDownVar" "JetES2PercUpVar" "JetES2PercDownVar")

for vardir in ${vardirs[@]}; do
  condor_submit vardir=${vardir} RunCorrelator.job
done
  
