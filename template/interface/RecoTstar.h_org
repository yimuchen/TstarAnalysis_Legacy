#ifndef RecoTstar_H
#define RecoTstar_H

#include "format.h"
#include "ConstantNumbers.h"

double RecoTstar(
		EvtInfoBranches EvtInfo,
		LepInfoBranches LepInfo, 
		JetInfoBranches JetInfo,
		PhotonInfoBranches PhotonInfo,
		int NLeptons, int L_Index[], 
		int NJets, int J_Index[], 
		int NPhotons, int P_Index[]){


	if(NJets>=4&&NPhotons==2&&NLeptons==1){

		float nuz1=-100000000;
		float nuz2=-100000000;

		int indx_ = L_Index[0];
		TLorentzVector vlep,nu;
		//vlep.SetPtEtaPhiM(LepInfo.Pt[indx_],
				//LepInfo.Eta[indx_],
				//LepInfo.Phi[indx_],
				//0);
		vlep.SetPxPyPzE(LepInfo.Px[indx_],
				LepInfo.Py[indx_],
				LepInfo.Pz[indx_],
				LepInfo.Energy[indx_]);
		SolutionOfWNeutrino(vlep,EvtInfo.PFMET,EvtInfo.PFMETPhi,80.398,nuz1,nuz2);
		nu.SetPxPyPzE(
				EvtInfo.PFMET*cos(EvtInfo.PFMETPhi),
				EvtInfo.PFMET*sin(EvtInfo.PFMETPhi),
				nuz2,
				sqrt(EvtInfo.PFMET*EvtInfo.PFMET+nuz2*nuz2));
		double chi2_=100000000;
		double Mtstar_ = -100000000;

		// [j1,j2] for Wjj, [j3,j4] for bjets coming from [Topjjb, Toplvb]
		for(int j1=0;j1<NJets;j1++)
			for(int j2=j1+1;j2<NJets;j2++)
				for(int j3=0;j3<NJets;j3++){
					if (j3==j1||j3==j2) continue;
					for(int j4=0;j4<NJets;j4++){
						if (j4==j1||j4==j2||j4==j3) continue;
						//if ((JetInfo.TrackCountHiPurBJetTags[J_Index[j3]]<1.93 ||
						//            JetInfo.TrackCountHiPurBJetTags[J_Index[j4]]<1.93) ) continue;

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

						double massRes[3] = {20.,40.,60.};

						double chi2_1 =
							pow(Wboson[0].Mag()-W_MASS,2)/pow(massRes[0],2)+
							pow(topquark[0].Mag()-top_MASS,2)/pow(massRes[1],2)+
							pow(topquark[1].Mag()-top_MASS,2)/pow(massRes[1],2)+
							pow(tstar1[0].Mag()-tstar1[1].Mag(),2)/pow(massRes[2],2);

						double chi2_2 =
							pow(Wboson[0].Mag()-W_MASS,2)/pow(massRes[0],2)+
							pow(topquark[0].Mag()-top_MASS,2)/pow(massRes[1],2)+
							pow(topquark[1].Mag()-top_MASS,2)/pow(massRes[1],2)+
							pow(tstar2[0].Mag()-tstar2[1].Mag(),2)/pow(massRes[2],2);

						if (chi2_1 < chi2_){
							chi2_ = chi2_1;
							Mtstar_ = tstar1[1].Mag();
						}
						if (chi2_2 < chi2_){
							chi2_ = chi2_2;
							Mtstar_ = tstar2[1].Mag();
						}

					}
				}

		return Mtstar_;
	}else{
		return -1;
	}

}

#endif
