#ifndef RecoTstar_H
#define RecoTstar_H

#include "format.h"
#include "ConstantNumbers.h"
#include "ObjectMCTruth.h"

bool UsingBtag = false;
//double massRes[5] = {20.,40.,40.,52.6,52.6};   // [W, Th, Tl, Tsh, Tsl]
double massRes[5] = {9.296,16.49,22.43,31.59,31.63};   // [W, Th, Tl, Tsh, Tsl]

//double RecoTstar(
double *RecoTstar(
        EvtInfoBranches EvtInfo,
        GenInfoBranches GenInfo,
        LepInfoBranches LepInfo, 
        JetInfoBranches JetInfo,
        PhotonInfoBranches PhotonInfo,
        int NLeptons, int L_Index[], 
        int NJets, int J_Index[], 
        int NPhotons, int P_Index[], 
        int MCTruth[], // [0] : photon position1; [1] : photon position2; [2] : lepton
        int Nu_mode){ // neutrino's solution mode [0,1,2] = [chi2, smaller,larger]

	double Mtstar_[2] = {-1,-1};	// [hadronic, leptonic]

    if(NJets>6) NJets=6;
    if(NJets>=4&&NPhotons==2&&NLeptons==1){

        float nuz[2]={-100000000,-100000000};
        double chi2_=100000000;

        int indx_ = L_Index[0];
        TLorentzVector vlep,nu;
        vlep.SetPxPyPzE(LepInfo.Px[indx_],
                LepInfo.Py[indx_],
                LepInfo.Pz[indx_],
                LepInfo.Energy[indx_]);
        SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

        int nuIdxStart = 0;
        int nuIdxEnd   = 2;

        if(Nu_mode!=0){
            nuIdxStart = Nu_mode-1;
            nuIdxEnd = Nu_mode;
        }

        for(int nuIdx = nuIdxStart ; nuIdx < nuIdxEnd;nuIdx++){
            nu.SetPxPyPzE(
                    EvtInfo.PFMET*cos(EvtInfo.PFMETPhi),
                    EvtInfo.PFMET*sin(EvtInfo.PFMETPhi),
                    nuz[nuIdx],
                    sqrt(EvtInfo.PFMET*EvtInfo.PFMET+nuz[nuIdx]*nuz[nuIdx]));

            // [j1,j2] for Wjj, [j3,j4] for bjets coming from [Topjjb, Toplvb]
            for(int j1=0;j1<NJets;j1++)
                for(int j2=j1+1;j2<NJets;j2++)
                    for(int j3=0;j3<NJets;j3++){
                        if (j3==j1||j3==j2) continue;
                        for(int j4=0;j4<NJets;j4++){
                            if (j4==j1||j4==j2||j4==j3) continue;
                            if(UsingBtag)
                                if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                            JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

                            TLorentzVector Jets_[4],gamma_[2],Wboson[2],topquark[2], tstar1[2], tstar2[2];
                            // for 4 jets
                            int indexTemp_[4] = {j1,j2,j3,j4};
                            for(int fillindx_=0;fillindx_<4;fillindx_++){
                                Jets_[fillindx_].SetPtEtaPhiM(
                                        JetInfo.Pt[J_Index[indexTemp_[fillindx_]]],
                                        JetInfo.Eta[J_Index[indexTemp_[fillindx_]]],
                                        JetInfo.Phi[J_Index[indexTemp_[fillindx_]]],
                                        0);
                            }
                            Wboson[0] = Jets_[0] + Jets_[1];
                            Wboson[1] = vlep + nu;  // leptonically decay

                            topquark[0] = Wboson[0] + Jets_[2];
                            topquark[1] = Wboson[1] + Jets_[3]; // leptonically decay

                            gamma_[0].SetPtEtaPhiM(
                                    PhotonInfo.Pt[P_Index[0]],
                                    PhotonInfo.Eta[P_Index[0]],
                                    PhotonInfo.Phi[P_Index[0]],
                                    0);
                            gamma_[1].SetPtEtaPhiM(
                                    PhotonInfo.Pt[P_Index[1]],
                                    PhotonInfo.Eta[P_Index[1]],
                                    PhotonInfo.Phi[P_Index[1]],
                                    0);
                            // first combination of top + gamma
                            tstar1[0] = topquark[0] + gamma_[0];
                            tstar1[1] = topquark[1] + gamma_[1];    // leptonically decay
                            // second combination of top + gamma
                            tstar2[0] = topquark[0] + gamma_[1];
                            tstar2[1] = topquark[1] + gamma_[0];    // leptonically decay

                            //double massRes[5] = {20.,40.,40.,52.6,52.6};  // [W, Th, Tl, Tsh, Tsl]

                            double chi2_1 =
                                pow(Wboson[0].Mag()-W_MASS,2)/pow(massRes[0],2)+
                                pow(topquark[0].Mag()-top_MASS,2)/pow(massRes[1],2)+
                                pow(topquark[1].Mag()-top_MASS,2)/pow(massRes[2],2)+
                                pow(tstar1[0].Mag()-tstar1[1].Mag(),2)/pow(sqrt(massRes[3]*massRes[3]+massRes[4]*massRes[4]),2);

                            double chi2_2 =
                                pow(Wboson[0].Mag()-W_MASS,2)/pow(massRes[0],2)+
                                pow(topquark[0].Mag()-top_MASS,2)/pow(massRes[1],2)+
                                pow(topquark[1].Mag()-top_MASS,2)/pow(massRes[2],2)+
                                pow(tstar2[0].Mag()-tstar2[1].Mag(),2)/pow(sqrt(massRes[3]*massRes[3]+massRes[4]*massRes[4]),2);

                            if (chi2_1 < chi2_){
                                chi2_ = chi2_1;
                                Mtstar_[0] = tstar1[0].Mag();
                                Mtstar_[1] = tstar1[1].Mag();
                                MCTruth[0] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                MCTruth[1] = ObjectMCTruth(P_Index[1],GenInfo,PhotonInfo);  // leptonically decay
                                MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                            }
                            if (chi2_2 < chi2_){
                                chi2_ = chi2_2;
                                Mtstar_[0] = tstar2[0].Mag();
                                Mtstar_[1] = tstar2[1].Mag();
                                MCTruth[1] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);  // leptonically decay
                                MCTruth[0] = ObjectMCTruth(P_Index[1],GenInfo,PhotonInfo);
                                MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                            }

                        }
                    }
        }
    }

    return Mtstar_;

}

#endif
