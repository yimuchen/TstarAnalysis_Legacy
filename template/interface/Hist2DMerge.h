#ifndef Hist2DMerge_H
#define Hist2DMerge_H

#include "samples_tr_1lepton4jets.h"

// Need to modify

// hmerge[0,1,2,3,4] = [TTbar, singleTop, boson(s), ttX , QCD]
void Hist2DMerge(TH2D *hists[], TH2D *hmerge[]){

    for(int i=TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
            i<=TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) hmerge[0]->Add(hists[i]);

    for(int i=T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
            i<=Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) hmerge[1]->Add(hists[i]);

    for(int i=DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1;
            i<=ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) hmerge[2]->Add(hists[i]);

    for(int i=TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1;
            i<=TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) hmerge[3]->Add(hists[i]);

    for(int i=QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2;
            i<=QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) hmerge[4]->Add(hists[i]);

}


#endif
