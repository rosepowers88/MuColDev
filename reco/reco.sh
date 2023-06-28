#bin/bash
source /opt/ilcsoft/muonc/init_ilcsoft.sh
Marlin --global.LCIOInputFiles=$1 --global.MaxRecordNumber=$2 --Output_REC.LCIOOutputFile=$3 $4
