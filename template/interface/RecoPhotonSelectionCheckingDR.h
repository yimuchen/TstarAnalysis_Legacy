#ifndef RecoPhotonSelectionCheckingDR_H
#define RecoPhotonSelectionCheckingDR_H

#include "format.h"
#include "DPHI.h"
#include "EffectiveAreaPhoton.h"

void RecoPhotonSelectionCheckingDR(LepInfoBranches LepInfo, JetInfoBranches JetInfo, PhotonInfoBranches PhotonInfo, 
        int NMuons, int M_Index[], int NElectrons, int E_Index[], int &NJets, int *J_Index, int &NPhotons, int *P_Index,float RhoPU){

    //RhoPU = 0;  // turn off Rho correction
    double LPDR = 0.3;
    for(int nl=0;nl<PhotonInfo.Size;nl++){
        // Medium Photon
        if(     
                (fabs(PhotonInfo.Pt[nl])>30) &&
                ((// Barrel
                  fabs(PhotonInfo.Eta[nl])<=1.4442 &&
		  PhotonInfo.passelectronveto[nl] &&
		  PhotonInfo.hadTowOverEm[nl] <=0.05 &&
                  PhotonInfo.SigmaIetaIeta[nl]<0.011 &&
                  (max(PhotonInfo.phoPFChIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],0)*RhoPU,(float)0.) < 0.7) &&
                  (max(PhotonInfo.phoPFNeuIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],1)*RhoPU,(float)0.)<0.4+0.04*PhotonInfo.Pt[nl])
		  &&
                  (max(PhotonInfo.phoPFPhoIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],2)*RhoPU,(float)0.)<0.5+0.005*PhotonInfo.Pt[nl])
                  )
                 ||
                 (// Endcap
		  (fabs(PhotonInfo.Eta[nl])<2.1)&&(fabs(PhotonInfo.Eta[nl])>=1.566) &&
		  PhotonInfo.SigmaIetaIeta[nl]<0.031 &&
		  PhotonInfo.passelectronveto[nl] &&
		  PhotonInfo.hadTowOverEm[nl] <=0.05 &&
                  (max(PhotonInfo.phoPFChIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],0)*RhoPU,(float)0.) < 0.5) &&
                  (max(PhotonInfo.phoPFNeuIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],1)*RhoPU,(float)0.)<1.5+0.04*PhotonInfo.Pt[nl])
		  &&
                  (max(PhotonInfo.phoPFPhoIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],2)*RhoPU,(float)0.)<1.0+0.005*PhotonInfo.Pt[nl])
                 ))
          ){

            P_Index[NPhotons] = nl;
            NPhotons++;
        }
    }   

}

#endif
