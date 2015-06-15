#!/bin/csh

set runs=`echo "Normal\nUncXsecPlus\nUncXsecMinus\nUncPUPlus\nUncPUMinus\nUncJESPlus\nUncJESMinus\nUncJERPlus\nUncJERMinus\nUncQsquare\nUncTrigPlus\nUncTrigMinus\nUncPhoIDPlus\nUncPhoIDMinus\nUncLepIDPlus\nUncLepIDMinus\nUncTopPtPlus\nUncTopPtMinus"`

foreach run($runs)
    echo "Prepare "$run
    if ( -e $run ) then
        echo "Remove "$run
        rm -r $run
    endif
    if ( ! -e $run ) then
        echo "Copying "$run
        cp -a template $run
        sed "s/int RunStatus_ = Normal/int RunStatus_ = ${run}/g" template/proj_anyregions_shape_loose_norunCR.cc >& .tmp__
        mv .tmp__ $run/proj_anyregions_shape_loose_norunCR.cc
    endif
end
