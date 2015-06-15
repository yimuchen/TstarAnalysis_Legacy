#!/bin/csh

set version_=v19

set runs=`echo "Normal UncXsecPlus UncXsecMinus UncPUPlus UncPUMinus UncJESPlus UncJESMinus UncJERPlus UncJERMinus UncTrigPlus UncTrigMinus UncPhoIDPlus UncPhoIDMinus UncLepIDPlus UncLepIDMinus UncTopPtPlus UncTopPtMinus UncQsquarePlus UncQsquareMinus UncMatchingPlus UncMatchingMinus"`

foreach run($runs)
    if ( -e $run ) then
        echo "Copying TemplateFit_${version_}.C "$run
        sed "s/int RunStatus_ = Normal/int RunStatus_ = ${run}/g" TemplateFit_${version_}.C | sed 's/bool AppliedSysOnShape =/bool AppliedSysOnShape = false;\/\//g' | sed 's/bool SkipFit =/bool SkipFit = true;\/\//g' >& .tmp__
        mv .tmp__ $run/TemplateFit_${version_}.C
        echo "Running TemplateFit_${version_}.C in $run"
	cd $run; root -l -b -q TemplateFit_${version_}.C++ >& log_${run} ; cd ..
        echo "Finish TemplateFit_${version_}.C in $run"
    else
        echo "[WARNING] No such folder of "${run}
    endif
end
