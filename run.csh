#!/bin/csh

set workspace=`pwd`

set runs=`echo "Normal\nUncXsecPlus\nUncXsecMinus\nUncPUPlus\nUncPUMinus\nUncJESPlus\nUncJESMinus\nUncJERPlus\nUncJERMinus\nUncQsquare\nUncTrigPlus\nUncTrigMinus\nUncPhoIDPlus\nUncPhoIDMinus\nUncLepIDPlus\nUncLepIDMinus\nUncTopPtPlus\nUncTopPtMinus"`

#source /cvmfs/cms.cern.ch/cmsset_default.csh;setenv SCRAM_ARCH slc5_amd64_gcc462;cd /wk3/cmsdata/Production_CMSSW536/Production_Data_DCSONLY_CMSSW53X_AOD/LUMI_SCRIPT/CMSSW_5_3_2_patch4/src;cmsenv;cd -
#source /cvmfs/cms.cern.ch/cmsset_default.csh;setenv SCRAM_ARCH slc5_amd64_gcc462;cd /wk3/cmsdata/Production_CMSSW536/Production_Data_DCSONLY_CMSSW53X_AOD/LUMI_SCRIPT/CMSSW_5_3_2_patch4/src;cmsenv;cd -
source /cvmfs/cms.cern.ch/cmsset_default.csh;setenv SCRAM_ARCH slc5_amd64_gcc462;cd /wk1/ymtzeng/BumpHuntingLimit/103014/MCMC/LatestVersion/CMSSW_5_3_13/src/;cmsenv;cd -
source /cvmfs/cms.cern.ch/cmsset_default.csh;setenv SCRAM_ARCH slc5_amd64_gcc462;cd /wk1/ymtzeng/BumpHuntingLimit/103014/MCMC/LatestVersion/CMSSW_5_3_13/src/;cmsenv;cd -

set time_=`date +%y%m%d%H%M%S`

foreach run($runs)
    echo "Running "$run
    cd ${workspace}/${run}
    root -l -b -q proj_anyregions_shape_loose_norunCR.cc++ >& log_${run}_${time_}_`hostname` &
end
