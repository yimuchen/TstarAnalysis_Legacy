#ifndef RecoPhotonSelection_H
#define RecoPhotonSelection_H

#include "format.h"
#include "DPHI.h"
#include "EffectiveAreaPhoton.h"
#include "TLorentzVector.h"

void RecoPhotonSelection(LepInfoBranches LepInfo, PhotonInfoBranches PhotonInfo, 
        int NMuons, int M_Index[], int NElectrons, int E_Index[], int &NPhotons, int *P_Index,float RhoPU){

    //RhoPU = 0;  // turn off Rho correction
    double LPDR = 0.3;
    for(int nl=0;nl<PhotonInfo.Size;nl++){
        /*

           BARREL   Loose (90%) Medium (80%)    Tight (70%)
           Conversion safe electron veto     Yes     Yes     Yes
           Single tower H/E  0.05    0.05    0.05
           0.012   0.011   0.011
           Rho corrected PF charged hadron isolation   2.6     1.5     0.7
           Rho corrected PF neutral hadron isolation   3.5 + 0.04*pho_Pt   1.0 + 0.04*pho_Pt   0.4 + 0.04*pho_Pt
           Rho corrected PF photon isolation   1.3 + 0.005*pho_Pt  0.7 + 0.005*pho_Pt  0.5 + 0.005*pho_Pt

           ENDCAPS    Loose (85%) Medium (75%)    Tight (65%)
           Conversion safe electron veto   Yes     Yes     Yes
           Single tower H/E    0.05    0.05    0.05
           0.034   0.033   0.031
           Rho corrected PF charged hadron isolation   2.3     1.2     0.5
           Rho corrected PF neutral hadron isolation   2.9 + 0.04*pho_Pt   1.5 + 0.04*pho_Pt   1.5 + 0.04*pho_Pt
           Rho corrected PF photon isolation   -   1.0 + 0.005*pho_Pt  1.0 + 0.005*pho_Pt
           */
        // Tight Photon
        if(     
                (fabs(PhotonInfo.Pt[nl])>25) &&
                //(fabs(PhotonInfo.Pt[nl])>20) &&
                //(fabs(PhotonInfo.Pt[nl])>30) &&
                ((// Barrel
                  fabs(PhotonInfo.Eta[nl])<=1.4442 &&
                  PhotonInfo.passelectronveto[nl] &&
                  PhotonInfo.hadTowOverEm[nl] <=0.05 &&
                  PhotonInfo.SigmaIetaIeta[nl]<0.012 &&
                  (max(PhotonInfo.phoPFChIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],0)*RhoPU,(float)0.) < 2.6) &&
                  (max(PhotonInfo.phoPFNeuIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],1)*RhoPU,(float)0.)<3.5+0.04*PhotonInfo.Pt[nl])
		  &&
                  (max(PhotonInfo.phoPFPhoIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],2)*RhoPU,(float)0.)<1.3+0.005*PhotonInfo.Pt[nl])
                  )
                 ||
                 (// Endcap
                  (fabs(PhotonInfo.Eta[nl])<2.1)&&(fabs(PhotonInfo.Eta[nl])>=1.566) &&
                  PhotonInfo.SigmaIetaIeta[nl]<0.034 &&
                  PhotonInfo.passelectronveto[nl] &&
                  PhotonInfo.hadTowOverEm[nl] <=0.05 &&
                  (max(PhotonInfo.phoPFChIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],0)*RhoPU,(float)0.) < 2.3) &&
                  (max(PhotonInfo.phoPFNeuIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],1)*RhoPU,(float)0.)<2.9+0.04*PhotonInfo.Pt[nl])
		  //&&
          //        (max(PhotonInfo.phoPFPhoIsoDR03[nl]-EffectiveAreaPhoton(PhotonInfo.Eta[nl],2)*RhoPU,(float)0.)<1.0+0.005*PhotonInfo.Pt[nl])
                 ))
          ){
            // Photon cleaning
            int dup_flag = 0;
            for(int nm=0;nm<NMuons;nm++){
                double dphi       = DPHI(PhotonInfo.Phi[nl],LepInfo.Phi[M_Index[nm]]);
                double deta       = fabs(PhotonInfo.Eta[nl]-LepInfo.Eta[M_Index[nm]]);
                double dR         = pow(dphi*dphi+deta*deta,0.5);
                if (dR<=LPDR) {
                    dup_flag = 1;
                    break;
                }   
            }   
            if(dup_flag) continue;
            dup_flag = 0;
            for(int ne=0;ne<NElectrons;ne++){
                TLorentzVector electronTemp;
                electronTemp.SetPxPyPzE(LepInfo.Px[E_Index[ne]],
                        LepInfo.Py[E_Index[ne]],
                        LepInfo.Pz[E_Index[ne]],
                        LepInfo.Energy[E_Index[ne]]);

                double dphi       = DPHI(PhotonInfo.Phi[nl],electronTemp.Phi());
                double deta       = fabs(PhotonInfo.Eta[nl]-electronTemp.Eta());
                double dR         = pow(dphi*dphi+deta*deta,0.5);
                if (dR<=LPDR) {
                    dup_flag = 1;
                    break;
                }   
            }   
            if(dup_flag) continue;

            P_Index[NPhotons] = nl;
            NPhotons++;
        }
    }   

}

#endif
