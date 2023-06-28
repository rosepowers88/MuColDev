#/bin/bash
source /opt/ilcsoft/muonc/init_ilcsoft.sh
export MARLIN_DLL="/work/rosep8/TauRecoDev/analysis/build/libAnaProcessors.so"
Marlin --global.LCIOInputFiles=$1 --global.MaxRecordNumber=$2 --GenHistMaker.PDG=$3 --MyAIDAProcessor.FileName=$4 $5


