#!/bin/csh

set runs=`echo "Normal UncXsecPlus UncXsecMinus UncPUPlus UncPUMinus UncJESPlus UncJESMinus UncJERPlus UncJERMinus UncTrigPlus UncTrigMinus UncPhoIDPlus UncPhoIDMinus UncLepIDPlus UncLepIDMinus UncTopPtPlus UncTopPtMinus UncQsquarePlus UncQsquareMinus UncMatchingPlus UncMatchingMinus"`


foreach run($runs)
    echo "Prepare "$run
    if ( -e $run ) then
        echo "Remove "$run
        rm -r $run
    endif
    if ( ! -e $run ) then
        echo "Preparing "$run
        mkdir $run
        cd $run
        ln -s ../MCTemplatesTree.root .
        mkdir -p FittingOutput/UNC
#        mkdir -p Systematics/sys/
#        touch Systematics/sys/sysforshape.root
#ln -s /Users/panda/files_IBM/ubuntu/CMS/WORKSPACE/T2TopGamma/Analysis/T2TopPlusPhoton/TemplateFit/DEV/data/102014/interface .
	ln -s /Users/panda/files_IBM/ubuntu/CMS/WORKSPACE/T2TopGamma/Analysis/T2TopPlusPhoton/TemplateFit/DEV/data/122414/TgammaAnalysis/interface .
        cd ..
    endif
end
