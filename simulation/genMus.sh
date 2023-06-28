#/bin/bash 
source /opt/ilcsoft/muonc/init_ilcsoft.sh

ddsim --steeringFile $1 --inputFiles $2 --outputFile $3 --numberOfEvents $4
