#include <string>
#include <iostream>
#include "TH1F.h"
#include "TFile.h"

bool IsBinCorrected = false;
bool IsBinByBinOrVariedShape = false;
const int N = 6;
string UncNames[N] = {
    "Mergecate00",
    "Mergecate01",
    "Mergecate02",
    "Mergecate03",
    "SignalAndControlRegionscate0",
    "SignalAndControlRegionscate4"
};


TFile *file = new TFile("../../FittingOutput/output_theta_1.0.root");
char buffer[128];
int mass = 950;

TH1F *GetUnc(string unc_name);//{
void SysOnBkgEstimateBinByBinShape(){
    TH1F *htotalunc;
    for(int i=0;i<N;i++){
	TH1F *htmp = GetUnc(UncNames[i]);
	if(i==0){
	    htotalunc = (TH1F*) htmp->Clone();
	}else{
	    for(int ibin=1;ibin <= htmp->GetNbinsX(); ibin++){
		float qradsum = sqrt(htotalunc->GetBinContent(ibin) * htotalunc->GetBinContent(ibin) + 
			htmp->GetBinContent(ibin) * htmp->GetBinContent(ibin)	);
		htotalunc->SetBinContent(ibin, qradsum);
	    }
	}
	delete htmp;
    }
    //totalunc = sqrt(totalunc);
    //std::cout<<"----------------------------------"<<std::endl;
    //std::cout<<"Total unc : "<<totalunc<<std::endl;
    htotalunc->SetName("htotalunc950");
    htotalunc->Draw();
    htotalunc->SaveAs("htotalunc950.root");
}

TH1F *GetUnc(string unc_name){

    sprintf(buffer,"Tgamma%iGeV__pdf_bkg", mass);
    TH1F *hstd = (TH1F*) file->Get(buffer);
    sprintf(buffer,"Tgamma%iGeV__pdf_bkg__%s__plus", mass,unc_name.c_str());
    TH1F *hplus = (TH1F*) file->Get(buffer);
    sprintf(buffer,"Tgamma%iGeV__pdf_bkg__%s__minus", mass,unc_name.c_str());
    TH1F *hminus = (TH1F*) file->Get(buffer);

    sprintf(buffer,"Tgamma%iGeV__pdf_bkg__%s", mass,unc_name.c_str());
    TH1F *hGetUnc = new TH1F(buffer,"",hstd->GetNbinsX(), hstd->GetXaxis()->GetXmin(), hstd->GetXaxis()->GetXmax());

    double variation_plus = 0.;
    double variation_minus = 0.;
    double Vstd = hstd->Integral();
    if(IsBinByBinOrVariedShape){
	for(int ibin = 1; ibin <= hstd->GetNbinsX(); ibin++){
	    if(!IsBinCorrected){
		if(fabs(hstd->GetBinContent(ibin) - hplus->GetBinContent(ibin)) > 
			fabs(hstd->GetBinContent(ibin) - hminus->GetBinContent(ibin)))
		    hGetUnc->SetBinContent(ibin, 
			    fabs(hstd->GetBinContent(ibin) - hplus->GetBinContent(ibin))/hstd->GetBinContent(ibin));
		else
		    hGetUnc->SetBinContent(ibin, 
			    fabs(hstd->GetBinContent(ibin) - hminus->GetBinContent(ibin))/hstd->GetBinContent(ibin));
	    }
	}
    }else{
	if(fabs(hplus->Integral() - hstd->Integral()) > fabs(hminus->Integral() - hstd->Integral())){
	    for(int ibin = 1; ibin <= hstd->GetNbinsX(); ibin++){
		hGetUnc->SetBinContent(ibin, 
			fabs(hstd->GetBinContent(ibin) - hplus->GetBinContent(ibin))/hstd->GetBinContent(ibin));
	    }
	}else{
	    for(int ibin = 1; ibin <= hstd->GetNbinsX(); ibin++){
		hGetUnc->SetBinContent(ibin, 
			fabs(hstd->GetBinContent(ibin) - hminus->GetBinContent(ibin))/hstd->GetBinContent(ibin));
	    }
	}
    }
    delete hstd;
    delete hplus;
    delete hminus;

    return hGetUnc;
}
