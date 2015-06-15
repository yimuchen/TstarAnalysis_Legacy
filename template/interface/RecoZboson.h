#ifndef RECOZBOSON_H
#define RECOZBOSON_H

#include "format.h"
#include "ConstantNumbers.h"

double RecoZboson(LepInfoBranches LepInfo, int NMuons, int M_Index[], int NElectrons, int E_Index[]){
    double Zmass = 0.0;
    // Lepton charge & flavor
    bool Lcf = false;
    const int N_lep_pair = 30;
    int lep_pair[N_lep_pair][2];
    for(int vv=0;vv<N_lep_pair;vv++) for(int lp=0;lp<2;lp++) lep_pair[vv][lp] = 0; 
    int N_lp = 0; 

    if(NElectrons>=2){
        for(int l1=0;l1<NElectrons;l1++){
            for(int l2=l1+1;l2<NElectrons;l2++){
                if(LepInfo.Charge[E_Index[l1]]==0) continue;
                if(LepInfo.Charge[E_Index[l2]]==0) continue;
                bool lcf = (LepInfo.Charge[E_Index[l1]]+LepInfo.Charge[E_Index[l2]])==0;
                if(lcf) {Lcf = lcf; lep_pair[N_lp][0]=E_Index[l1]; lep_pair[N_lp][1]=E_Index[l2]; N_lp++;}
            }    
        }    
    }    

    if(NMuons>=2){
        for(int l1=0;l1<NMuons;l1++){
            for(int l2=l1+1;l2<NMuons;l2++){
                if(LepInfo.Charge[M_Index[l1]]==0) continue;
                if(LepInfo.Charge[M_Index[l2]]==0) continue;
                bool lcf = (LepInfo.Charge[M_Index[l1]]+LepInfo.Charge[M_Index[l2]])==0;
                if(lcf) {Lcf = lcf; lep_pair[N_lp][0]=M_Index[l1]; lep_pair[N_lp][1]=M_Index[l2]; N_lp++;}
            }    
        }    
    }    

    // Z selection
    int     Z_lep_cand[2]   = {-1,-1};
    float   Zwin_df     = 1000.;
    float   Z_pt_max    = 0.;
    TLorentzVector Z_cand;
    if(N_lp!=0&&Lcf){
        for(int lp=0;lp<N_lp;lp++){
            float lepmass = ELECTRON_MASS;
            if(LepInfo.LeptonType[lep_pair[lp][0]]==13) lepmass = MUON_MASS;

            TLorentzVector LEP_Z[3];
            LEP_Z[0].SetPtEtaPhiM(LepInfo.Pt[lep_pair[lp][0]],
                    LepInfo.Eta[lep_pair[lp][0]],
                    LepInfo.Phi[lep_pair[lp][0]],
                    lepmass);
            LEP_Z[1].SetPtEtaPhiM(LepInfo.Pt[lep_pair[lp][1]],
                    LepInfo.Eta[lep_pair[lp][1]],
                    LepInfo.Phi[lep_pair[lp][1]],
                    lepmass);
            LEP_Z[2] = LEP_Z[0]+LEP_Z[1];
            float Zwin_df_t = fabs(LEP_Z[2].Mag()-Z_MASS);
            if(Zwin_df_t<Zwin_df){
                Zwin_df     = Zwin_df_t;
                Zmass       = LEP_Z[2].Mag();
                Z_cand      = LEP_Z[2];
                Z_lep_cand[0]   = lep_pair[lp][0];
                Z_lep_cand[1]   = lep_pair[lp][1];
                Z_pt_max = LEP_Z[2].Pt();
            }
        }
    }

    return Zmass;
}

#endif
