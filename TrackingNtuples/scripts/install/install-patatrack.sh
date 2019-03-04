#!/usr/bin/bash

# Create and run in a new directory
mkdir trackingproject && cd trackingproject

# Set the environment variables
export VO_CMS_SW_DIR=$PWD && export SCRAM_ARCH=slc7_amd64_gcc700 && export LANG=C

wget http://cmsrep.cern.ch/cmssw/repos/bootstrap.sh -O $VO_CMS_SW_DIR/bootstrap.sh
chmod a+x $VO_CMS_SW_DIR/bootstrap.sh

$VO_CMS_SW_DIR/bootstrap.sh -a slc7_amd64_gcc700 -r cms -path $VO_CMS_SW_DIR setup

# Patatrack may already have moved several versions ahead since development is underway (currently on 10_5_0)
$VO_CMS_SW_DIR/common/cmspkg -a slc7_amd64_gcc700 install -y cms+cmssw+CMSSW_10_2_5

# If this next command doesn't work it's probably because of shitty formatting, refer https://github.com/cms-patatrack/patatrack-website/blob/8a184fccb571b3cebab4dd5ad9d4790e500c916f/wiki/PatatrackReleases.md
# This is the old (10_2_5) release that I am on, as I mentioned - they have already released another version 10_2_6

$VO_CMS_SW_DIR/common/cmspkg -a slc7_amd64_gcc700 -- rpm --prefix=$VO_CMS_SW_DIR --upgrade --nodeps --replacepkgs \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+cms-git-tools+180901.0-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+cmssw+CMSSW_10_2_5_Patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+cmssw-tool-conf+44.0-patatrack3-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+coral+CORAL_2_3_21-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+coral-tool-conf+2.1-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+distcc-gcc-toolfile+2.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+gcc-toolfile+13.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+icc-gcc-toolfile+3.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/cms+llvm-gcc-toolfile+13.0-patatrack4-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+cub+1.8.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+cuda+9.2.148-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+cuda-api-wrappers+20180504-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+cuda-gdb-wrapper+1.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+eigen+1f44b667dd9aeeb153284b15fc7fe159d2c09329-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+gbl+V02-01-03-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+llvm+7.0.0-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+lwtnn+2.4-patatrack3-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+professor2+2.2.1-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+py2-dxr+1.0-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+py2-dxr-toolfile+1.0-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+py2-tensorflow+1.6.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+py2-tensorflow-toolfile+1.0-patatrack-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+python_tools+1.0-patatrack2-1-1.slc7_amd64_gcc700.rpm \
https://cern.ch/fwyzard/patatrack/rpms/CMSSW_10_2_X/slc7_amd64_gcc700/external+tensorflow+1.6.0-patatrack-1-1.slc7_amd64_gcc700.rpm


# Now follow the instructions at https://github.com/cms-patatrack/patatrack-website/blob/a5f3fca866b928c2f88ae1cc0c75151ea83a5a91/wiki/PatatrackDevelopment.md
# This will set up your working area for patatrack 10_2_5
# The instructions below are the same as the ones at the link above.
cd ..

scram list CMSSW_10_2_5
cmsrel CMSSW_10_2_5_Patatrack
cd CMSSW_10_2_5_Patatrack/src
cmsenv

git cms-init --upstream-only || true
# you will see the error
#     fatal: 'CMSSW_10_2_5_Patatrack' is not a commit and a branch 'from-CMSSW_10_2_5_Patatrack' cannot be created from it
# it is expected, just follow the rest of the instructions

# add the Patatrack remote and branches
git cms-remote add cms-patatrack
git checkout CMSSW_10_2_5_Patatrack -b CMSSW_10_2_X_Patatrack
git branch -u cms-patatrack/CMSSW_10_2_X_Patatrack
git checkout CMSSW_10_2_5_Patatrack -b from-CMSSW_10_2_5_Patatrack

# enable the developer's repository
git cms-init

# Clone this project TrackingNtuples into the 'src' directory, scram build it, and you should be ready to go
# cd CMSSW_10_2_5_Patatrack/src

git clone https://github.com/SwapneelM/TrackingNtuples.git
cd TrackingNtuples/
scram b -j 4

