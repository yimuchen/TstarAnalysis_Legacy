#ifndef ScalePhotonID_H
#define ScalePhotonID_H

#include "ConstantNumbers.h"
#include "TFile.h"
#include "TH2F.h"
#include <iostream>

TFile *FileScalePhoton = new TFile("interface/Photon_ID_CSEV_SF_Jan22rereco_Full2012_S10_MC_V01.root");
TH2F *PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
TH2F *PhotonIDSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonIDSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01");
TH2F *PhotonIDSF_TightWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonIDSF_TightWP_Jan22rereco_Full2012_S10_MC_V01");
TH2F *PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01");
TH2F *PhotonCSEVSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonCSEVSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01");
TH2F *PhotonCSEVSF_TightWP_Jan22rereco_Full2012_S10_MC_V01 = (TH2F*) FileScalePhoton->Get("PhotonCSEVSF_TightWP_Jan22rereco_Full2012_S10_MC_V01");

float ScalePhotonID(float pT, float eta, Photon_WP WP, int MODE){
    /*
       https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonID2012#Data_MC_Efficiency_scale_factors
       */

    if (WP==P_LOOSE){
        int binNumber = PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01->FindBin(pT, fabs(eta));
        return 
            (PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) + 
             (MODE-1.)*PhotonIDSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber) )* 
            (PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) +
             (MODE-1.)*PhotonCSEVSF_LooseWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber));
    }
    if (WP==P_MEDIUM){
        int binNumber = PhotonIDSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01->FindBin(pT, fabs(eta));
        return 
            (PhotonIDSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) + 
             (MODE-1.)*PhotonIDSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber) )* 
            (PhotonCSEVSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) +
             (MODE-1.)*PhotonCSEVSF_MediumWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber));
    }
    if (WP==P_TIGHT){
        int binNumber = PhotonIDSF_TightWP_Jan22rereco_Full2012_S10_MC_V01->FindBin(pT, fabs(eta));
        return 
            (PhotonIDSF_TightWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) + 
             (MODE-1.)*PhotonIDSF_TightWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber) )* 
            (PhotonCSEVSF_TightWP_Jan22rereco_Full2012_S10_MC_V01->GetBinContent(binNumber) +
             (MODE-1.)*PhotonCSEVSF_TightWP_Jan22rereco_Full2012_S10_MC_V01->GetBinError(binNumber));
    }
    std::cout<<"[WARNING] No corresponding scale factor for photon ID"<<std::endl;
    return 1.0;
}

#endif
