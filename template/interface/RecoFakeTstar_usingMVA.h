#ifndef RecoFakeTstar_H
#define RecoFakeTstar_H

#include "format.h"
#include "ConstantNumbers.h"
#include "ObjectMCTruth.h"
#include "TMVAClassification_MLPBNN.class.C"
#include <string>
#include <vector>


/*
   1). fake photon (r) from jet
   - 1r + 1l + >=5j
   - 0r + 1l + >=6j
   2). fake r from lepton
   - 1r + 2l + >=4j
   - 0r + 3l + >=4j
   3). fake r from lepton and jet
   - 0r + 2l + >=5j
   4). fake lepton from jet
   - 2r + 0l + >=5j
   5). fake lepton from r
   - 3r + 0l + >=4j

 */


double RecoFakeTstar(
        EvtInfoBranches EvtInfo,
        GenInfoBranches GenInfo,
        LepInfoBranches LepInfo, 
        JetInfoBranches JetInfo,
        PhotonInfoBranches PhotonInfo,
        int NLeptons, int L_Index[], 
        int NJets, int J_Index[], 
        int NPhotons, int P_Index[],
        int MCTruth[], // [0] : photon position1; [1] : photon position2; [2] : lepton
        int category){ 

    // Get MLPBNN value to pick a possible jet faking to a photon
    std::vector<std::string> theInputVars;
    theInputVars.push_back("CHF");
    theInputVars.push_back("NHF");
    theInputVars.push_back("CHFandNHF");
    theInputVars.push_back("CHForder");
    theInputVars.push_back("NHForder");
    theInputVars.push_back("CHFandNHForder");
    ReadMLPBNN ReadMLPBNN_(theInputVars);
    // use jets with the first Nx_ maximum of MLPBNN for fake photons
    float MLPBNNcut = 0.2;
    const int Nx_ = 7;
    int FakePhotonCandidatesFromJets[Nx_] ;  // index in J_Index 
    float MLPBNNValues[Nx_] ;  // index in J_Index 
    for(int inx = 0;inx<Nx_;inx++){
        FakePhotonCandidatesFromJets[inx] = -1.;
        float MaximumOfMLPBNN= -1000.;
        for(int j1=0;j1<NJets;j1++){    // selected jet
            bool IsDuplicate = false;
            for(int jnx = 0; jnx < inx; jnx++){
                if(j1== FakePhotonCandidatesFromJets[jnx])
                    IsDuplicate =true;
            }
            if(IsDuplicate) continue;

            int Order_=0;
            int OrderCHF_=0;
            int OrderNHF_=0;
            for(int j2=0;j2<JetInfo.Size;j2++){ // reco jet used to find an order
                if(j2==J_Index[j1]) continue;
                if((JetInfo.CHF[j2] + JetInfo.NHF[j2]) < (JetInfo.CHF[J_Index[j1]] + JetInfo.NHF[J_Index[j1]]))
                    Order_++;
                if((JetInfo.CHF[j2] ) < (JetInfo.CHF[J_Index[j1]] ))
                    OrderCHF_++;
                if((JetInfo.NHF[j2] ) < (JetInfo.NHF[J_Index[j1]] ))
                    OrderNHF_++;
            }   

            std::vector<double> inputValues;
            inputValues.push_back(JetInfo.CHF[J_Index[j1]]); // CHF
            inputValues.push_back(JetInfo.NHF[J_Index[j1]]); // NHF
            inputValues.push_back(JetInfo.CHF[J_Index[j1]]+JetInfo.NHF[J_Index[j1]]); // CHF+NHF
            inputValues.push_back(OrderCHF_); // CHF order
            inputValues.push_back(OrderNHF_); // NHF order
            inputValues.push_back(Order_); // CHF+NHF order

            double mlpbnn_ = ReadMLPBNN_.GetMvaValue(inputValues);
            if(mlpbnn_ > MaximumOfMLPBNN){
                FakePhotonCandidatesFromJets[inx] = j1; 
                MLPBNNValues[inx] = (float)mlpbnn_; 
                MaximumOfMLPBNN = (float)mlpbnn_;
            }   
        }   
    }

    if(Nx_>2)
        for(int inx = 2;inx<Nx_;inx++)
            if((MLPBNNValues[inx])<MLPBNNcut)
                FakePhotonCandidatesFromJets[inx] = FakePhotonCandidatesFromJets[inx-1];

    if(NJets>8) NJets = 8;  // reduce some combination

    if(category==0){
        if(!(NJets>=5&&NPhotons==1&&NLeptons==1)) return -1; 

        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; double Mtstar_ = -100000000;

        int indx_ = L_Index[0];
        TLorentzVector vlep,nu;
        vlep.SetPxPyPzE(LepInfo.Px[indx_],
                LepInfo.Py[indx_],
                LepInfo.Pz[indx_],
                LepInfo.Energy[indx_]);
        SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

        for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                            for(int j5=0;j5<NJets;j5++){
                                if (j5==j1||j5==j2||j5==j3||j5==j4) continue;
                                //if(UsingBtag)
                                //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                Wboson[1] = vlep + nu;

                                topquark[0] = Wboson[0] + Jets_[2];
                                topquark[1] = Wboson[1] + Jets_[3];

                                gamma_[0].SetPtEtaPhiM(
                                        PhotonInfo.Pt[P_Index[0]],
                                        PhotonInfo.Eta[P_Index[0]],
                                        PhotonInfo.Phi[P_Index[0]],
                                        0);

                                bool IsJ5in = false;
                                for(int inx = 0;inx<Nx_;inx++) if(j5==FakePhotonCandidatesFromJets[inx]) IsJ5in = true;
                                if(!IsJ5in) continue;

                                gamma_[1].SetPtEtaPhiM(
                                        JetInfo.Pt[J_Index[j5]]*(JetInfo.CEF[J_Index[j5]]+JetInfo.NEF[J_Index[j5]]),
                                        //JetInfo.Pt[J_Index[j5]],
                                        JetInfo.Eta[J_Index[j5]],
                                        JetInfo.Phi[J_Index[j5]],
                                        0);
                                // first combination of top + gamma
                                tstar1[0] = topquark[0] + gamma_[0];
                                tstar1[1] = topquark[1] + gamma_[1];
                                // second combination of top + gamma
                                tstar2[0] = topquark[0] + gamma_[1];
                                tstar2[1] = topquark[1] + gamma_[0];

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
                                    Mtstar_ = tstar1[1].Mag();
                                    MCTruth[0] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[1] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                                if (chi2_2 < chi2_){
                                    chi2_ = chi2_2;
                                    Mtstar_ = tstar2[1].Mag();
                                    MCTruth[1] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[0] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                            }

                        }
                    }
        }

        return Mtstar_;
    }else if(category==1){
        if(!(NJets>=6&&NLeptons==1) ) return -1;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; double Mtstar_ = -100000000;
        int indx_ = L_Index[0];
        TLorentzVector vlep,nu;
        vlep.SetPxPyPzE(LepInfo.Px[indx_],
                LepInfo.Py[indx_],
                LepInfo.Pz[indx_],
                LepInfo.Energy[indx_]);
        SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

        for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                            for(int j5=0;j5<NJets;j5++){
                                if (j5==j1||j5==j2||j5==j3||j5==j4) continue;
                                for(int j6=0;j6<NJets;j6++){
                                    if (j6==j1||j6==j2||j6==j3||j6==j4||j6==j5) continue;
                                    //if(UsingBtag)
                                    //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                    //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                    Wboson[1] = vlep + nu;

                                    topquark[0] = Wboson[0] + Jets_[2];
                                    topquark[1] = Wboson[1] + Jets_[3];

                                    bool IsJ5in = false;
                                    for(int inx = 0;inx<Nx_;inx++) if(j5==FakePhotonCandidatesFromJets[inx]) IsJ5in = true;
                                    if(!IsJ5in) continue;
                                    gamma_[0].SetPtEtaPhiM(
                                            JetInfo.Pt[J_Index[j5]]*(JetInfo.CEF[J_Index[j5]]+JetInfo.NEF[J_Index[j5]]),
                                            //JetInfo.Pt[J_Index[j5]],
                                            JetInfo.Eta[J_Index[j5]],
                                            JetInfo.Phi[J_Index[j5]],
                                            0);
                                    bool IsJ6in = false;
                                    for(int inx = 0;inx<Nx_;inx++) if(j6==FakePhotonCandidatesFromJets[inx]) IsJ6in = true;
                                    if(!IsJ6in) continue;
                                    gamma_[1].SetPtEtaPhiM(
                                            JetInfo.Pt[J_Index[j6]]*(JetInfo.CEF[J_Index[j6]]+JetInfo.NEF[J_Index[j6]]),
                                            //JetInfo.Pt[J_Index[j6]],
                                            JetInfo.Eta[J_Index[j6]],
                                            JetInfo.Phi[J_Index[j6]],
                                            0);
                                    // first combination of top + gamma
                                    tstar1[0] = topquark[0] + gamma_[0];
                                    tstar1[1] = topquark[1] + gamma_[1];
                                    // second combination of top + gamma
                                    tstar2[0] = topquark[0] + gamma_[1];
                                    tstar2[1] = topquark[1] + gamma_[0];

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
                                        Mtstar_ = tstar1[1].Mag();
                                        MCTruth[0] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                        MCTruth[1] = ObjectMCTruth(J_Index[j6],GenInfo,JetInfo);
                                        MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                    }
                                    if (chi2_2 < chi2_){
                                        chi2_ = chi2_2;
                                        Mtstar_ = tstar2[1].Mag();
                                        MCTruth[1] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                        MCTruth[0] = ObjectMCTruth(J_Index[j6],GenInfo,JetInfo);
                                        MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                    }
                                }

                            }
                        }
                    }
        }
        return Mtstar_;
    }else if (category==2){
        if(!(NJets>=4&&NPhotons==1&&NLeptons==2)) return -1;
        double Mtstar_ = -100000000;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; 
        for(int lep_idx=0;lep_idx<NLeptons;lep_idx++){
            int indx_ = L_Index[lep_idx];
            TLorentzVector vlep,nu;
            vlep.SetPxPyPzE(LepInfo.Px[indx_],
                    LepInfo.Py[indx_],
                    LepInfo.Pz[indx_],
                    LepInfo.Energy[indx_]);
            SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

            for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                                //if(UsingBtag)
                                //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                Wboson[1] = vlep + nu;

                                topquark[0] = Wboson[0] + Jets_[2];
                                topquark[1] = Wboson[1] + Jets_[3];

                                gamma_[0].SetPtEtaPhiM(
                                        PhotonInfo.Pt[P_Index[0]],
                                        PhotonInfo.Eta[P_Index[0]],
                                        PhotonInfo.Phi[P_Index[0]],
                                        0);
                                gamma_[1].SetPtEtaPhiM(
                                        LepInfo.Pt[L_Index[NLeptons-lep_idx-1]],
                                        LepInfo.Eta[L_Index[NLeptons-lep_idx-1]],
                                        LepInfo.Phi[L_Index[NLeptons-lep_idx-1]],
                                        0);
                                // first combination of top + gamma
                                tstar1[0] = topquark[0] + gamma_[0];
                                tstar1[1] = topquark[1] + gamma_[1];
                                // second combination of top + gamma
                                tstar2[0] = topquark[0] + gamma_[1];
                                tstar2[1] = topquark[1] + gamma_[0];

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
                                    Mtstar_ = tstar1[1].Mag();
                                    MCTruth[0] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[1] = ObjectMCTruth(L_Index[NLeptons-lep_idx-1],GenInfo,LepInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                                if (chi2_2 < chi2_){
                                    chi2_ = chi2_2;
                                    Mtstar_ = tstar2[1].Mag();
                                    MCTruth[1] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[0] = ObjectMCTruth(L_Index[NLeptons-lep_idx-1],GenInfo,LepInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                            }

                        }
            }
        }
        return Mtstar_;
    }else if(category==3){
        if(!(NJets>=4&&NLeptons==3)) return -1;
        double Mtstar_ = -100000000;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; 
        for(int lep_idx=0;lep_idx<NLeptons;lep_idx++){
            int indx_ = L_Index[lep_idx];
            TLorentzVector vlep,nu;
            vlep.SetPxPyPzE(LepInfo.Px[indx_],
                    LepInfo.Py[indx_],
                    LepInfo.Pz[indx_],
                    LepInfo.Energy[indx_]);
            SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

            for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                                //if(UsingBtag)
                                //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                Wboson[1] = vlep + nu;

                                topquark[0] = Wboson[0] + Jets_[2];
                                topquark[1] = Wboson[1] + Jets_[3];

                                int fake_photon_from_lep_idx[2]={-1,-1}; 
                                if(lep_idx==0){
                                    fake_photon_from_lep_idx[0] = 1;
                                    fake_photon_from_lep_idx[1] = 2;
                                }else if(lep_idx==1){
                                    fake_photon_from_lep_idx[0] = 0;
                                    fake_photon_from_lep_idx[1] = 2;
                                }else if(lep_idx==2){
                                    fake_photon_from_lep_idx[0] = 0;
                                    fake_photon_from_lep_idx[1] = 1;
                                }

                                gamma_[0].SetPtEtaPhiM(
                                        LepInfo.Pt[L_Index[fake_photon_from_lep_idx[0]]],
                                        LepInfo.Eta[L_Index[fake_photon_from_lep_idx[0]]],
                                        LepInfo.Phi[L_Index[fake_photon_from_lep_idx[0]]],
                                        0);
                                gamma_[1].SetPtEtaPhiM(
                                        LepInfo.Pt[L_Index[fake_photon_from_lep_idx[1]]],
                                        LepInfo.Eta[L_Index[fake_photon_from_lep_idx[1]]],
                                        LepInfo.Phi[L_Index[fake_photon_from_lep_idx[1]]],
                                        0);
                                // first combination of top + gamma
                                tstar1[0] = topquark[0] + gamma_[0];
                                tstar1[1] = topquark[1] + gamma_[1];
                                // second combination of top + gamma
                                tstar2[0] = topquark[0] + gamma_[1];
                                tstar2[1] = topquark[1] + gamma_[0];

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
                                    Mtstar_ = tstar1[1].Mag();
                                    MCTruth[0] = ObjectMCTruth(L_Index[fake_photon_from_lep_idx[0]],GenInfo,LepInfo);
                                    MCTruth[1] = ObjectMCTruth(L_Index[fake_photon_from_lep_idx[1]],GenInfo,LepInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                                if (chi2_2 < chi2_){
                                    chi2_ = chi2_2;
                                    Mtstar_ = tstar2[1].Mag();
                                    MCTruth[1] = ObjectMCTruth(L_Index[fake_photon_from_lep_idx[0]],GenInfo,LepInfo);
                                    MCTruth[0] = ObjectMCTruth(L_Index[fake_photon_from_lep_idx[1]],GenInfo,LepInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                            }

                        }
            }
        }
        return Mtstar_;
    }else if(category==4){
        if(!(NJets>=5&&NLeptons==2)) return -1;
        double Mtstar_ = -100000000;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; 
        for(int lep_idx=0;lep_idx<NLeptons;lep_idx++){
            int indx_ = L_Index[lep_idx];
            TLorentzVector vlep,nu;
            vlep.SetPxPyPzE(LepInfo.Px[indx_],
                    LepInfo.Py[indx_],
                    LepInfo.Pz[indx_],
                    LepInfo.Energy[indx_]);
            SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

            for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                                for(int j5=0;j5<NJets;j5++){
                                    if (j5==j1||j5==j2||j5==j3 || j5==j4) continue;
                                    //if(UsingBtag)
                                    //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                    //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                    Wboson[1] = vlep + nu;

                                    topquark[0] = Wboson[0] + Jets_[2];
                                    topquark[1] = Wboson[1] + Jets_[3];

                                    bool IsJ5in = false;
                                    for(int inx = 0;inx<Nx_;inx++) if(j5==FakePhotonCandidatesFromJets[inx]) IsJ5in = true;
                                    if(!IsJ5in) continue;
                                    gamma_[0].SetPtEtaPhiM(
                                            JetInfo.Pt[J_Index[j5]]*(JetInfo.CEF[J_Index[j5]]+JetInfo.NEF[J_Index[j5]]),
                                            //JetInfo.Pt[J_Index[j5]],
                                            JetInfo.Eta[J_Index[j5]],
                                            JetInfo.Phi[J_Index[j5]],
                                            0);
                                    gamma_[1].SetPtEtaPhiM(
                                            LepInfo.Pt[L_Index[NLeptons-lep_idx-1]],
                                            LepInfo.Eta[L_Index[NLeptons-lep_idx-1]],
                                            LepInfo.Phi[L_Index[NLeptons-lep_idx-1]],
                                            0);
                                    // first combination of top + gamma
                                    tstar1[0] = topquark[0] + gamma_[0];
                                    tstar1[1] = topquark[1] + gamma_[1];
                                    // second combination of top + gamma
                                    tstar2[0] = topquark[0] + gamma_[1];
                                    tstar2[1] = topquark[1] + gamma_[0];

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
                                        Mtstar_ = tstar1[1].Mag();
                                        MCTruth[0] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                        MCTruth[1] = ObjectMCTruth(L_Index[NLeptons-lep_idx-1],GenInfo,LepInfo);
                                        MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                    }
                                    if (chi2_2 < chi2_){
                                        chi2_ = chi2_2;
                                        Mtstar_ = tstar2[1].Mag();
                                        MCTruth[1] = ObjectMCTruth(J_Index[j5],GenInfo,JetInfo);
                                        MCTruth[0] = ObjectMCTruth(L_Index[NLeptons-lep_idx-1],GenInfo,LepInfo);
                                        MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                    }
                                }
                            }

                        }
            }
        }
        return Mtstar_;
    }else if(category==5){
        if(!(NJets>=5&&NPhotons==2&&NLeptons==0)) return -1;
        double Mtstar_ = -100000000;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; 
        for(int lep_idx=0;lep_idx<NJets;lep_idx++){
            int indx_ = J_Index[lep_idx];
            TLorentzVector vlep,nu;
            vlep.SetPxPyPzE(
                    JetInfo.Px[indx_],
                    JetInfo.Py[indx_],
                    JetInfo.Pz[indx_],
                    JetInfo.Energy[indx_]);
            SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

            for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                                if (j1==lep_idx || j2==lep_idx || j3==lep_idx || j4 == lep_idx) continue;
                                //if(UsingBtag)
                                //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                Wboson[1] = vlep + nu;

                                topquark[0] = Wboson[0] + Jets_[2];
                                topquark[1] = Wboson[1] + Jets_[3];

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
                                tstar1[1] = topquark[1] + gamma_[1];
                                // second combination of top + gamma
                                tstar2[0] = topquark[0] + gamma_[1];
                                tstar2[1] = topquark[1] + gamma_[0];

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
                                    Mtstar_ = tstar1[1].Mag();
                                    MCTruth[0] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[1] = ObjectMCTruth(P_Index[1],GenInfo,PhotonInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                                if (chi2_2 < chi2_){
                                    chi2_ = chi2_2;
                                    Mtstar_ = tstar2[1].Mag();
                                    MCTruth[1] = ObjectMCTruth(P_Index[0],GenInfo,PhotonInfo);
                                    MCTruth[0] = ObjectMCTruth(P_Index[1],GenInfo,PhotonInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                            }

                        }
            }
        }
        return Mtstar_;
    }else if(category==6){
        if(!(NJets>=4&&NPhotons==3&&NLeptons==0)) return -1;
        double Mtstar_ = -100000000;
        float nuz[2]={-100000000,-100000000}; double chi2_=100000000; 
        for(int lep_idx=0;lep_idx<NPhotons;lep_idx++){
            int indx_ = P_Index[lep_idx];
            TLorentzVector vlep,nu;
            vlep.SetPtEtaPhiM(
                    PhotonInfo.Pt[indx_],
                    PhotonInfo.Eta[indx_],
                    PhotonInfo.Phi[indx_],
                    0);
            SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz[0],nuz[1]);

            for(int nuIdx = 0; nuIdx < 2;nuIdx++){
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
                                //if(UsingBtag)
                                //    if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 &&
                                //                JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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
                                Wboson[1] = vlep + nu;

                                topquark[0] = Wboson[0] + Jets_[2];
                                topquark[1] = Wboson[1] + Jets_[3];

                                int photon_idx[2]={-1,-1}; 
                                if(lep_idx==0){
                                    photon_idx[0] = 1;
                                    photon_idx[1] = 2;
                                }else if(lep_idx==1){
                                    photon_idx[0] = 0;
                                    photon_idx[1] = 2;
                                }else if(lep_idx==2){
                                    photon_idx[0] = 0;
                                    photon_idx[1] = 1;
                                }

                                gamma_[0].SetPtEtaPhiM(
                                        PhotonInfo.Pt[P_Index[photon_idx[0]]],
                                        PhotonInfo.Eta[P_Index[photon_idx[0]]],
                                        PhotonInfo.Phi[P_Index[photon_idx[0]]],
                                        0);
                                gamma_[1].SetPtEtaPhiM(
                                        PhotonInfo.Pt[P_Index[photon_idx[1]]],
                                        PhotonInfo.Eta[P_Index[photon_idx[1]]],
                                        PhotonInfo.Phi[P_Index[photon_idx[1]]],
                                        0);
                                // first combination of top + gamma
                                tstar1[0] = topquark[0] + gamma_[0];
                                tstar1[1] = topquark[1] + gamma_[1];
                                // second combination of top + gamma
                                tstar2[0] = topquark[0] + gamma_[1];
                                tstar2[1] = topquark[1] + gamma_[0];

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
                                    Mtstar_ = tstar1[1].Mag();
                                    MCTruth[0] = ObjectMCTruth(P_Index[photon_idx[0]],GenInfo,PhotonInfo);
                                    MCTruth[1] = ObjectMCTruth(P_Index[photon_idx[1]],GenInfo,PhotonInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                                if (chi2_2 < chi2_){
                                    chi2_ = chi2_2;
                                    Mtstar_ = tstar2[1].Mag();
                                    MCTruth[1] = ObjectMCTruth(P_Index[photon_idx[0]],GenInfo,PhotonInfo);
                                    MCTruth[0] = ObjectMCTruth(P_Index[photon_idx[1]],GenInfo,PhotonInfo);
                                    MCTruth[2] = ObjectMCTruth(indx_,GenInfo,LepInfo);
                                }
                            }

                        }
            }
        }
        return Mtstar_;
    }else{
        return -1;
    }

}

#endif
