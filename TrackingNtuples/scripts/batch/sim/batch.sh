#! /bin/bash

source /afs/cern.ch/work/s/smehta/public/tsc-cern/trackingproject/heptrkx-cmssw/cmsset_default.sh
set -o errexit
set -o nounset

SANDBOX=$PWD
jobid=$1
nevts=$2
batch_path=/afs/cern.ch/work/s/smehta/public/tsc-cern/trackingsim
cd /afs/cern.ch/work/s/smehta/public/tsc-cern/trackingproject/CMSSW_10_2_5_Patatrack/src
eval `scramv1 runtime -sh`
cd $SANDBOX

cmsDriver.py TTbar_13TeV_TuneCUETP8M1_cfi  --conditions auto:phase1_2018_realistic -n $nevts --era Run2_2018 --eventcontent FEVTDEBUG --relval 9000,50 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic25ns13TeVEarly2018Collision --geometry DB:Extended --fileout file:TT_GEN_SIM_$jobid.root  > gen_sim_$jobid.log  2>&1

cmsDriver.py step2  --conditions auto:phase1_2018_realistic -s DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@relval2018 --datatier GEN-SIM-DIGI-RAW -n -1 --geometry DB:Extended --era Run2_2018 --eventcontent FEVTDEBUGHLT --filein  file:TT_GEN_SIM_$jobid.root  --fileout file:TT_GEN_SIM_DIGI_RAW_$jobid.root  > gen_sim_digi_raw_$jobid.log 2>&1

cmsRun $batch_path/step3_RAW2DIGI_RECO_VALIDATION_DQM.py inputFiles=TT_GEN_SIM_DIGI_RAW_$jobid.root outputFile=TT_NTUPLE_$jobid.root

rm step3*.root
mv *.root /eos/home-s/smehta/trackingsim/.
