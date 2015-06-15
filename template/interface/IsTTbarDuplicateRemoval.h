#ifndef IsTTbarDuplicateRemoval_H
#define IsTTbarDuplicateRemoval_H

#include "format.h"

bool IsTTbarDuplicateRemoval(GenInfoBranches GenInfo){
    // ttgamma overlap removal in ttbar samples AN-13-195
    /* /TTGamma_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_RD1_START53_V7N-v1/AODSIM
       (https://hypernews.cern.ch/HyperNews/CMS/get/generators/2369.html)
       13  = pta       ! minimum pt for the photons

       0.5 = drjj    ! min distance between jets
       0.5   = drbb    ! min distance between b's
       0.3 = draa    ! min distance between gammas
       0.5   = drbj    ! min distance between b and jet
       0.3 = draj    ! min distance between gamma and jet
       0.5 = drjl    ! min distance between jet and lepton
       0.3   = drab    ! min distance between gamma and b
       0.5   = drbl    ! min distance between b and lepton
       0.3 = dral    ! min distance between gamma and lepton
     */
    /* /TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola (MCDBId : 6721)
       0  = pta       ! minimum pt for the photons

       0.001 = drjj    ! min distance between jets
       0.001   = drbb    ! min distance between b's
       0   = draa    ! min distance between gammas
       0.001   = drbj    ! min distance between b and jet
       0   = draj    ! min distance between gamma and jet
       0   = drjl    ! min distance between jet and lepton
       0   = drab    ! min distance between gamma and b
       0   = drbl    ! min distance between b and lepton
       0   = dral    ! min distance between gamma and lepton
     */

    // Step1 : find all of photons
    int Ntop = 0; 
    int GenTopIdx[2] = {-1, -1}; 
    int Npho = 0; 
    int GenPhoIdx[GenInfo.Size];
    int Njet = 0; // including gluon, light quark, and heavy quark (b, c)
    int GenJetIdx[GenInfo.Size];
    int Nlep = 0;
    int GenLepIdx[GenInfo.Size];
    for(int g=0;g<GenInfo.Size;g++){
        GenPhoIdx[g] = -1;  // initialization
        GenJetIdx[g] = -1;  // initialization
        GenLepIdx[g] = -1;  // initialization
        // pick top quark
        if(GenInfo.Status[g]==3 && fabs(GenInfo.PdgID[g]) == 6){
            if(Ntop<2){
                GenTopIdx[Ntop] = g;
            }else{
                printf("[ERROR] more than 2 tops\n");
            }
            Ntop++;
        }
        // pick photon
        if(fabs(GenInfo.PdgID[g]) == 22){
            GenPhoIdx[Npho] = g;
            Npho++;
        }
        // pick jet
        if((fabs(GenInfo.PdgID[g]) <= 5 && GenInfo.PdgID[g]!=0 ) || GenInfo.PdgID[g] == 21){
            GenJetIdx[Njet] = g;
            Njet++;
        }

        // pick lep
        if((fabs(GenInfo.PdgID[g]) == 11 || fabs(GenInfo.PdgID[g]) == 13 || fabs(GenInfo.PdgID[g]) == 15  )){
            GenLepIdx[Nlep] = g;
            Nlep++;
        }
    }

    // Step2 : ensure if photon is associated with top
    // & Step3 : if ttbar event is in a phase space of ttG (Et > 10GeV and drb > 0)
    bool IsTTbarDuplicate = false;
    for(int g=0;g<Npho;g++){
        for(int t=0;t<Ntop;t++){
            double dr_ = DR(GenInfo.Eta[GenPhoIdx[g]],GenInfo.Phi[GenPhoIdx[g]],
                    GenInfo.Eta[GenTopIdx[t]],GenInfo.Phi[GenTopIdx[t]]);
            double dpT_ = fabs(GenInfo.Pt[GenPhoIdx[g]]-GenInfo.Pt[GenTopIdx[t]])/GenInfo.Pt[GenTopIdx[t]];
            if(  dr_ < 0.5 && GenInfo.Pt[GenPhoIdx[g]] >13  ){
                IsTTbarDuplicate = true;
            }
        }
    }
    bool IsdRjj = false;    // 0.5
    bool IsdRjl = false;    // 0.5
    bool IsdRaa = false;    // 0.3
    bool IsdRaj = false;    // 0.3
    bool IsdRal = false;    // 0.3

    // 1). IsdRjj
    for(int j1 = 0; j1 < Njet; j1++){
        for(int j2 = 0; j2 < Njet; j2++){
            if(j1==j2) continue;
            double dr_ = DR(GenInfo.Eta[GenJetIdx[j1]],GenInfo.Phi[GenJetIdx[j1]],
                    GenInfo.Eta[GenJetIdx[j2]],GenInfo.Phi[GenJetIdx[j2]]);
            if(dr_>0.5){
                IsdRjj = true;
                break;
            }
        }
        if(IsdRjj) break;
    }
    // 2). IsdRjl
    for(int j1 = 0; j1 < Njet; j1++){
        for(int j2 = 0; j2 < Nlep; j2++){
            double dr_ = DR(GenInfo.Eta[GenJetIdx[j1]],GenInfo.Phi[GenJetIdx[j1]],
                    GenInfo.Eta[GenLepIdx[j2]],GenInfo.Phi[GenLepIdx[j2]]);
            if(dr_>0.5){
                IsdRjl = true;
                break;
            }
        }
        if(IsdRjl) break;
    }
    // 3). IsdRaa
    for(int j1 = 0; j1 < Npho; j1++){
        for(int j2 = 0; j2 < Npho; j2++){
            if(j1==j2) continue;
            double dr_ = DR(GenInfo.Eta[GenPhoIdx[j1]],GenInfo.Phi[GenPhoIdx[j1]],
                    GenInfo.Eta[GenPhoIdx[j2]],GenInfo.Phi[GenPhoIdx[j2]]);
            if(dr_>0.3){
                IsdRaa = true;
                break;
            }
        }
        if(IsdRaa) break;
    }
    // 4). IsdRaj
    for(int j1 = 0; j1 < Njet; j1++){
        for(int j2 = 0; j2 < Npho; j2++){
            double dr_ = DR(GenInfo.Eta[GenJetIdx[j1]],GenInfo.Phi[GenJetIdx[j1]],
                    GenInfo.Eta[GenPhoIdx[j2]],GenInfo.Phi[GenPhoIdx[j2]]);
            if(dr_>0.3){
                IsdRaj = true;
                break;
            }
        }
        if(IsdRaj) break;
    }
    // 5). IsdRal
    for(int j1 = 0; j1 < Nlep; j1++){
        for(int j2 = 0; j2 < Npho; j2++){
            double dr_ = DR(GenInfo.Eta[GenLepIdx[j1]],GenInfo.Phi[GenLepIdx[j1]],
                    GenInfo.Eta[GenPhoIdx[j2]],GenInfo.Phi[GenPhoIdx[j2]]);
            if(dr_>0.3){
                IsdRal = true;
                break;
            }
        }
        if(IsdRal) break;
    }

    return (IsTTbarDuplicate && IsdRjj && IsdRjl && IsdRaa && IsdRaj && IsdRal);
}

#endif
