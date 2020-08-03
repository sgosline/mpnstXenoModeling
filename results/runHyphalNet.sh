#!/bin/bash


cpath=$PWD/../../hyphalnet/examples/pdxDrugResponse/testDrugResponse.py
echo $cpath
args='--graph=../../hyphalnet/data/igraphPPI.pkl'
python3 $cpath $args

Rscript ../../hyphalnet/bin/netStats.R --distFile='panPDXDistances.csv' --commFile='mpnstPDXmuts_communityStats.csv' --synapseProj='syn21984813'
