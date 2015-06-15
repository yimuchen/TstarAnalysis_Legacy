#include "TTree.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "interface/format_mini.h"
//#include "interface/format.h"
#include "interface/MCTruth.h"
#include "interface/TriggerBooking.h"
#include "interface/MCTruthChannel.h"
//#include "interface/HitFitInfoBranches.h"
#include "interface/DPHI.h"
#include "interface/DR.h"
#include "interface/JER.h"
//#include "interface/samples_tr_mu.h"
#include "interface/samples_tr_1lepton4jets.h"
#include "interface/ScalePhotonID.h"
#include "interface/ScaleLeptonID.h"
//#include "interface/samples_tr.h"

int RunStatus_ = Normal;   
/*
enum RunStatus{
    Normal,
    UncXsecPlus,
    UncXsecMinus,
    UncPUPlus,
    UncPUMinus,
    UncJESPlus,
    UncJESMinus,
    UncJERPlus,
    UncJETMinus,
    UncQsquare,
    RunStatusSize
};
*/
#include "interface/DRAWSTACK.h"
#include "interface/HistMerge.h"
#include "interface/checkEvt.h"
#include "interface/RecoLooseLeptonSelection.h"
#include "interface/RecoLeptonSelection.h"
#include "interface/RecoLeptonForCleaningSelection.h"
#include "interface/RecoPhotonSelection.h"
#include "interface/RecoPhotonForJetCleaningSelection.h"
#include "interface/RecoLoosePhotonSelection.h"
#include "interface/RecoPhotonSelectionCheckingDR.h"
#include "interface/RecoJetSelection.h"
#include "interface/SolutionOfWNeutrino.h"
#include "interface/ReduceTree.h"
#include "interface/PUReweighting.h"
#include "interface/TopTriggerEfficiencyProvider.cc"
//#include "interface/TMVAClassification_CutsD.class.C"
#include "TChain.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "interface/RecoZbosonFordiPhoton.h"
#include "interface/RecoTstar.h"
#include "interface/RecoFakeTstar.h"
#include "interface/dMRecoTstar.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "stdlib.h"
#include <fstream>
#include <map>
#include <vector>
#include <time.h>
#include <math.h>
#include "interface/RecoZboson.h"
#include "interface/RecoMoflljj.h"
#include "interface/RecoMZb.h"
#include "interface/PhotonMCTruth.h"
#include "interface/RecoNoCutPhotonSelection.h"
#include "interface/CorrectionOnTstar.h"
#include "interface/MainSysUnc.h"



void analyze(int sample_idx,bool IfConsider)
{
    // PDF information
    FILE *pdffileCR;
    char buffer[128];
    sprintf(buffer,"PDF/%s_noselection.txt",SAMPLE_PDFnoSelection[sample_idx].tag);
    pdffileCR = fopen(buffer,"w");

    if(IfConsider){
        TChain *root = new TChain("bprimeKit/root");
        root->Add(SAMPLE_PDFnoSelection[sample_idx].filename);

        EvtInfoBranches EvtInfo;
        EvtInfo.Register(root);
        ReduceTree(root);


        int nevents_total = root->GetEntries();
        for(int entry=0;entry<root->GetEntries();entry++) {
            if ((entry%10000) == 0)	printf("Loading event #%d of %d.\n",entry,nevents_total);
            root->GetEntry(entry);
            fprintf(pdffileCR,"%i %i %f %f %f %f %f %f\n",
                    EvtInfo.PDFid1,
                    EvtInfo.PDFid2,
                    EvtInfo.PDFx1,
                    EvtInfo.PDFx2,
                    EvtInfo.PDFscale,
                    EvtInfo.PDFv1,
                    EvtInfo.PDFv2,
                    1.0
                   );
        }
        root->Delete();
    }else{
        fprintf(pdffileCR,"0 0 0 0 0 0 0 0\n");
    }
    fclose(pdffileCR);
}

void dump_pdf_withoutSelection(){
    int time1 = time(NULL);
    int time2 = 0;
    int time3 = 0;

    bool IfConsider = true;
	for(int i=0;i<samples_order_size;i++) {
        IfConsider = true;
		std::cout<<"file : "<<SAMPLE_PDFnoSelection[i].tag<<std::endl;
        if(i==Data||(i>=TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1 &&
                    i<=QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1))
            IfConsider = false;
		analyze(i, IfConsider);
	}
	time2=time(NULL); time3=time2-time1; cout<<"---- Spending time : " << (time3/60.0) <<" (mins) ----"<<endl;
	std::cout<<"done "<<std::endl;
}
