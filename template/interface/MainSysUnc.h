#ifndef MainSysUnc_H
#define MainSysUnc_H

#include "samples_tr_1lepton4jets.h"

bool IncludeDYMatchScaleUnc = false;    // false for loose and signal region, but true for control region (due to limited stats.)

// Xsec + QScale + Matching
TH1F *MainSysUnc(TH1F *hists[], int rebin_){    /*Relative uncertainty*/

    // from TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1 ~ samples_order_size
    printf("\n");

    int STAND_INDEX = TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    int STAND_DY_INDEX = DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1;
    // 1). XSEC
    TH1F *hMainSysUnc = new TH1F("hMainSysUnctmp","",hists[STAND_INDEX]->GetNbinsX(),hists[STAND_INDEX]->GetXaxis()->GetXmin(),
            hists[STAND_INDEX]->GetXaxis()->GetXmax());

    TH1F *hstd_ = new TH1F("hstd_","",hists[STAND_INDEX]->GetNbinsX(),hists[STAND_INDEX]->GetXaxis()->GetXmin(),
            hists[STAND_INDEX]->GetXaxis()->GetXmax());

    TH1F *hstdTTbar_ = (TH1F*) hists[STAND_INDEX]->Clone();
    hstdTTbar_->Rebin( rebin_ );

    TH1F *hstdDY_ = (TH1F*) hists[STAND_DY_INDEX]->Clone();
    hstdDY_->Rebin( rebin_ );

    TH1F *hvaried_ = new TH1F("hvaried_","",hists[STAND_INDEX]->GetNbinsX(),hists[STAND_INDEX]->GetXaxis()->GetXmin(),
            hists[STAND_INDEX]->GetXaxis()->GetXmax());
    for(int ibin = 1;ibin <= hvaried_->GetNbinsX(); ibin++){
        hvaried_->SetBinContent(ibin,0.);
    }

    for(int i=TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
            i<samples_order_size;i++){
        hstd_->Add(hists[i]);
        for(int ibin = 1;ibin <= hists[i]->GetNbinsX(); ibin++){
            hvaried_->SetBinContent(ibin, hvaried_->GetBinContent(ibin) + (1.+SAMPLE[i].unc)*hists[i]->GetBinContent(ibin));
        }
    }

    hMainSysUnc->Rebin( rebin_ );
    hvaried_->Rebin( rebin_ );
    hstd_->Rebin( rebin_ );
    for(int ibin = 1;ibin <= hMainSysUnc->GetNbinsX(); ibin++){
        if(hstd_->GetBinContent(ibin)!=0){
            hMainSysUnc->SetBinContent(ibin,fabs(hvaried_->GetBinContent(ibin) - hstd_->GetBinContent(ibin))/hstd_->GetBinContent(ibin));
        }else{
            hMainSysUnc->SetBinContent(ibin,0.);
        }
        printf("[Xsec uncerntainty %ith] %f\n", ibin, hMainSysUnc->GetBinContent(ibin));
    }

    // 2). Qscale
    // ttbar
    float uncmax_ = 0.;
    int IndxM = TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    int IndxP = TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    TH1F *hQscaleM = (TH1F*) hists[IndxM]->Clone();
    TH1F *hQscaleP = (TH1F*) hists[IndxP]->Clone();
    hQscaleM->Rebin( rebin_ );
    hQscaleP->Rebin( rebin_ );
    for(int ibin = 1;ibin <= hMainSysUnc->GetNbinsX(); ibin++){
        uncmax_ = 0.;
        float std_ = hstdTTbar_->GetBinContent(ibin);
        if(std_!=0){
            if(uncmax_< fabs(std_ - hQscaleM->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hQscaleM->GetBinContent(ibin))/hstd_->GetBinContent(ibin);

            if(uncmax_< fabs(std_ - hQscaleP->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hQscaleP->GetBinContent(ibin))/hstd_->GetBinContent(ibin);
        }
        hMainSysUnc->SetBinContent(ibin, sqrt( hMainSysUnc->GetBinContent(ibin)*hMainSysUnc->GetBinContent(ibin) + 
                    uncmax_*uncmax_));
    }
    // DY 
    uncmax_ = 0.;
    IndxM = DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1 ;
    IndxP = DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    TH1F *hQscaleDYM = (TH1F*) hists[IndxM]->Clone();
    TH1F *hQscaleDYP = (TH1F*) hists[IndxP]->Clone();
    hQscaleDYM->Rebin( rebin_ );
    hQscaleDYP->Rebin( rebin_ );
    for(int ibin = 1;ibin <= hMainSysUnc->GetNbinsX(); ibin++){
        uncmax_ = 0.;
        float std_ = hstdDY_->GetBinContent(ibin);
        if(std_!=0){
            if(uncmax_< fabs(std_ - hQscaleDYM->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hQscaleDYM->GetBinContent(ibin))/hstd_->GetBinContent(ibin);

            if(uncmax_< fabs(std_ - hQscaleDYP->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hQscaleDYP->GetBinContent(ibin))/hstd_->GetBinContent(ibin);
        }
        if(IncludeDYMatchScaleUnc)
            hMainSysUnc->SetBinContent(ibin, sqrt( hMainSysUnc->GetBinContent(ibin)*hMainSysUnc->GetBinContent(ibin) + 
                        uncmax_*uncmax_));
    }

    // 3). Matching
    // ttbar
    IndxM = TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    IndxP = TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    TH1F *hMatchingM = (TH1F*) hists[IndxM]->Clone();
    TH1F *hMatchingP = (TH1F*) hists[IndxP]->Clone();
    hMatchingM->Rebin( rebin_ );
    hMatchingP->Rebin( rebin_ );
    for(int ibin = 1;ibin <= hMainSysUnc->GetNbinsX(); ibin++){
        uncmax_ = 0.;
        float std_ = hstdTTbar_->GetBinContent(ibin);
        if(std_!=0){
            if(uncmax_< fabs(std_ - hMatchingM->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hMatchingM->GetBinContent(ibin))/hstd_->GetBinContent(ibin);

            if(uncmax_< fabs(std_ - hMatchingP->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hMatchingP->GetBinContent(ibin))/hstd_->GetBinContent(ibin);
        }
        hMainSysUnc->SetBinContent(ibin, sqrt( hMainSysUnc->GetBinContent(ibin)*hMainSysUnc->GetBinContent(ibin) + 
                    uncmax_*uncmax_));
    }
    // DY
    IndxM = DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1;
    IndxP = DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    TH1F *hMatchingDYM = (TH1F*) hists[IndxM]->Clone();
    TH1F *hMatchingDYP = (TH1F*) hists[IndxP]->Clone();
    hMatchingDYM->Rebin( rebin_ );
    hMatchingDYP->Rebin( rebin_ );
    for(int ibin = 1;ibin <= hMainSysUnc->GetNbinsX(); ibin++){
        uncmax_ = 0.;
        float std_ = hstdDY_->GetBinContent(ibin);
        if(std_!=0){
            if(uncmax_< fabs(std_ - hMatchingDYM->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hMatchingDYM->GetBinContent(ibin))/hstd_->GetBinContent(ibin);

            if(uncmax_< fabs(std_ - hMatchingDYP->GetBinContent(ibin))/hstd_->GetBinContent(ibin))
                uncmax_ = fabs(std_ - hMatchingDYP->GetBinContent(ibin))/hstd_->GetBinContent(ibin);
        }
        if(IncludeDYMatchScaleUnc)
            hMainSysUnc->SetBinContent(ibin, sqrt( hMainSysUnc->GetBinContent(ibin)*hMainSysUnc->GetBinContent(ibin) + 
                        uncmax_*uncmax_));
    }


    return hMainSysUnc;
}

#endif
