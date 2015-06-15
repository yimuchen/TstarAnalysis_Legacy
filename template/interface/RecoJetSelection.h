#ifndef RecoJetSelection_H
#define RecoJetSelection_H

#include "format.h"
#include "DPHI.h"
#include "TLorentzVector.h"

void RecoJetSelection(LepInfoBranches LepInfo, JetInfoBranches JetInfo, PhotonInfoBranches PhotonInfo, 
        int NMuons, int M_Index[], int NElectrons, int E_Index[], int &NJets, int *J_Index, int NPhotons,int  P_Index[]){

    double LJDR = 0.5;
    for(int nj=0;nj<JetInfo.Size;nj++){
        // Jet cleaning
        int dup_flag = 0;
        for(int nm=0;nm<NMuons;nm++){
            double dphi       = DPHI(JetInfo.Phi[nj],LepInfo.Phi[M_Index[nm]]);
            double deta       = fabs(JetInfo.Eta[nj]-LepInfo.Eta[M_Index[nm]]);
            double dR         = pow(dphi*dphi+deta*deta,0.5);
            if (dR<=LJDR) {
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

            double dphi       = DPHI(JetInfo.Phi[nj],electronTemp.Phi());
            double deta       = fabs(JetInfo.Eta[nj]-electronTemp.Eta());
            double dR         = pow(dphi*dphi+deta*deta,0.5);
            if (dR<=LJDR) {
                dup_flag = 1;
                break;
            }
        }
        if(dup_flag) continue;
        dup_flag = 0;
        for(int np=0;np<NPhotons;np++){
            double dphi       = DPHI(JetInfo.Phi[nj],PhotonInfo.Phi[P_Index[np]]);
            double deta       = fabs(JetInfo.Eta[nj]-PhotonInfo.Eta[P_Index[np]]);
            double dR         = pow(dphi*dphi+deta*deta,0.5);
            if (dR<=LJDR) {
                dup_flag = 1;
                break;
            }
        }
        if(dup_flag) continue;

        // Jet selection
        if((JetInfo.JetIDLOOSE[nj]==1)&&
                (fabs(JetInfo.Pt[nj])>20) &&
                //(fabs(JetInfo.Pt[nj])>25) &&
                //(fabs(JetInfo.Pt[nj])>30) &&
                (fabs(JetInfo.Eta[nj])<2.5)
          ) {
            J_Index[NJets] = nj;
            NJets++;
        }
    }   

}

#endif
