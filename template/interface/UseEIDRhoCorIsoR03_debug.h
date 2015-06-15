#ifndef UseEIDRhoCorIsoR03_H
#define UseEIDRhoCorIsoR03_H

#include "format.h"
#include "ElectronEffectiveArea.h"
#include <iostream>

enum WorkingPoint {
    VETO,
    LOOSE,
    MEDIUM,
    TIGHT
};

bool UseEIDRhoCorIsoR03(LepInfoBranches LepInfo,EvtInfoBranches EvtInfo,int nl, int workingPoint)
{
    // From http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/src/EGammaCutBasedEleId.cc?view=markup
    bool useEID_=false;
    double pt = LepInfo.Pt[nl];

    float cut_dEtaIn[2]         = {999.9, 999.9};
    float cut_dPhiIn[2]         = {999.9, 999.9};
    float cut_sigmaIEtaIEta[2]  = {999.9, 999.9};
    float cut_hoe[2]            = {999.9, 999.9};
    float cut_ooemoop[2]        = {999.9, 999.9};
    float cut_d0vtx[2]          = {999.9, 999.9};
    float cut_dzvtx[2]          = {999.9, 999.9};
    float cut_iso[2]            = {999.9, 999.9};
    bool cut_vtxFit[2]          = {false, false};
    unsigned int cut_mHits[2]   = {999, 999};

    if (workingPoint == VETO) {
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.010;
        cut_dPhiIn[0]        = 0.800; cut_dPhiIn[1]        = 0.700;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.150; cut_hoe[1]           = 999.9;
        cut_ooemoop[0]       = 999.9; cut_ooemoop[1]       = 999.9;
        cut_d0vtx[0]         = 0.040; cut_d0vtx[1]         = 0.040;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = false; cut_vtxFit[1]        = false;
        cut_mHits[0]         = 999  ; cut_mHits[1]         = 999;
        cut_iso[0]           = 0.150; cut_iso[1]           = 0.150;
    } 
    else if (workingPoint == LOOSE) {
        cut_dEtaIn[0]        = 0.007; cut_dEtaIn[1]        = 0.009;
        cut_dPhiIn[0]        = 0.150; cut_dPhiIn[1]        = 0.100;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.200; cut_dzvtx[1]         = 0.200;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    } 
    else if (workingPoint == MEDIUM) {
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.007;
        cut_dPhiIn[0]        = 0.060; cut_dPhiIn[1]        = 0.030;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 1    ; cut_mHits[1]         = 1;
        if (pt >= 20.0) {
            cut_iso[0] = 0.150; cut_iso[1] = 0.150;
        }
        else {
            cut_iso[0] = 0.150; cut_iso[1] = 0.100;
        }
    } 
    else if (workingPoint == TIGHT) {
        cut_dEtaIn[0]        = 0.004; cut_dEtaIn[1]        = 0.005;
        cut_dPhiIn[0]        = 0.030; cut_dPhiIn[1]        = 0.020;
        cut_sigmaIEtaIEta[0] = 0.010; cut_sigmaIEtaIEta[1] = 0.030;
        cut_hoe[0]           = 0.120; cut_hoe[1]           = 0.100;
        cut_ooemoop[0]       = 0.050; cut_ooemoop[1]       = 0.050;
        cut_d0vtx[0]         = 0.020; cut_d0vtx[1]         = 0.020;
        cut_dzvtx[0]         = 0.100; cut_dzvtx[1]         = 0.100;
        cut_vtxFit[0]        = true ; cut_vtxFit[1]        = true;
        cut_mHits[0]         = 0    ; cut_mHits[1]         = 0;
        if (pt >= 20.0) {
            cut_iso[0] = 0.100; cut_iso[1] = 0.100;
        }
        else {
            cut_iso[0] = 0.100; cut_iso[1] = 0.070;
        }
    } 
    else {
        std::cout << "[TestWP] Undefined working point" << std::endl;
    }

    if(LepInfo.LeptonType[nl]!=11){
        std::cout<<"[Warning] it's not a electron type"<<std::endl;
        return useEID_;
    }


    int idx=-1;
    if((fabs(LepInfo.Eta[nl])<=1.4442)){
        idx = 0;
    //}else if(((fabs(LepInfo.Eta[nl])<2.4)&&(fabs(LepInfo.Eta[nl])>=1.566))){
    }else if(((fabs(LepInfo.Eta[nl])<2.5)&&(fabs(LepInfo.Eta[nl])>=1.566))){    // 11252013 modify
        idx = 1;
    }

    if(idx == -1 ) return useEID_;
    // effective area for isolation
    ElectronEffectiveArea::ElectronEffectiveAreaTarget EATarget = ElectronEffectiveArea::kEleEAFall11MC; 
    float AEffR03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, LepInfo.Eta[nl], EATarget); 
    if(!EvtInfo.McFlag){
        EATarget = ElectronEffectiveArea::kEleEAData2012;
        AEffR03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, LepInfo.Eta[nl], EATarget);
        //AEffR03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleNeutralHadronIso03, LepInfo.Eta[nl], EATarget) + ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaIso03, LepInfo.Eta[nl], EATarget) ;
    }

    // apply to neutrals 
    double rhoPrime = std::max((double)EvtInfo.RhoPU[0], 0.0);
    double iso_n = std::max(LepInfo.NeutralHadronIsoR03[nl] + LepInfo.PhotonIsoR03[nl] - rhoPrime * AEffR03, 0.0);

    double iso_ch = LepInfo.ChargedHadronIsoR03[nl];

    // compute final isolation 
    double iso = (iso_n + iso_ch) / pt;

    double dEtaIn = LepInfo.EldeltaEta[nl];
    double dPhiIn = LepInfo.EldeltaPhi[nl];
    double sigmaIEtaIEta = LepInfo.ElsigmaIetaIeta[nl];
    double hoe = LepInfo.ElHadoverEm[nl];
    float  ooemoop       = (1.0/LepInfo.ElEcalE[nl] - LepInfo.ElEoverP[nl]/LepInfo.ElEcalE[nl]);
    float d0vtx = LepInfo.ElTrackDxy_PV[nl];
    float dzvtx = LepInfo.ElTrackDz[nl];
    bool vtxFitConversion = LepInfo.ElhasConv[nl];
    float mHits = LepInfo.NumberOfExpectedInnerHits[nl];

    if( true ) {
        printf("[Debug ElecID] : ID-%i dEtaIn-%1.4f dPhiIn-%1.4f hoe-%1.4f ooemoop-%1.6f d0vtx-%1.6f dzvtx-%1.6f vtxFitConversion-%1.6f mHits-%i iso-%1.6d\n",
                workingPoint,
                dEtaIn,
                dPhiIn,
                sigmaIEtaIEta,
                hoe,
                ooemoop,
                d0vtx,
                dzvtx,
                vtxFitConversion,
                mHits,
                iso
                );
        std::cout<<"[iso] : "<<(float)iso<<std::endl;
    }

    // test cuts
    if (fabs(dEtaIn) < cut_dEtaIn[idx])             
    if (fabs(dPhiIn) < cut_dPhiIn[idx])             
    if (sigmaIEtaIEta < cut_sigmaIEtaIEta[idx])     
    if (hoe < cut_hoe[idx])                         
    if (fabs(ooemoop) < cut_ooemoop[idx])           
    if (fabs(d0vtx) < cut_d0vtx[idx])               
    if (fabs(dzvtx) < cut_dzvtx[idx])               
    if (!cut_vtxFit[idx] || !vtxFitConversion)      
    if (mHits <= cut_mHits[idx])                    
    if ((float)iso < cut_iso[idx])                         
        useEID_=true;

    return useEID_;
}
#endif
