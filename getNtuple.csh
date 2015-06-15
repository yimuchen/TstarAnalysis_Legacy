#!/bin/csh

set runs=`echo "Normal\nUncXsecPlus\nUncXsecMinus\nUncPUPlus\nUncPUMinus\nUncJESPlus\nUncJESMinus\nUncJERPlus\nUncJERMinus\nUncQsquare\nUncTrigPlus\nUncTrigMinus\nUncPhoIDPlus\nUncPhoIDMinus\nUncLepIDPlus\nUncLepIDMinus\nUncTopPtPlus\nUncTopPtMinus"`

source /cvmfs/cms.cern.ch/cmsset_default.csh;setenv SCRAM_ARCH slc5_amd64_gcc462;cd /wk3/cmsdata/Production_CMSSW536/Production_Data_DCSONLY_CMSSW53X_AOD/LUMI_SCRIPT/CMSSW_5_3_2_patch4/src;cmsenv;cd -

set time_=`date +%y%m%d%H%M%S`

if ( -e .note_${time_} ) then
rm .note_${time_}
endif
touch .note_${time_}

echo "hadd SumNtuples.root " | awk '{printf("%s %s",$1,$2)}' >> .note_${time_}

foreach run($runs)

echo $run | awk '{printf(" %s/MCTemplatesTree.root ",$1)}' >> .note_${time_}

end

cat .note_${time_} 
cat .note_${time_} | csh
rm .note_${time_}
echo ""
echo "Done"
