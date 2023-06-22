#/bin/bash
source /opt/ilcsoft/muonc/init_ilcsoft.sh
export MARLIN_DLL=$(realpath build/libAnaProcessors.so):${MARLIN_DLL}
Marlin --global.LCIOInputFiles=$(realpath $1) --GenHistMaker.PDG=$2 --global.MaxRecordNumber=$3 --MyAIDAProcessor.FileName=$4 $5


