#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

float mass = 950.;
char buffer[128];
float GetRelativeUncInTotal(TH1F *hSTD, TH1F *hUNC);
void compare(){

    gROOT->ProcessLine(".L interface/setTDRStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetErrorX(0.5);

    TCanvas *c1 = new TCanvas("c1","",640,640);
    c1->cd();
    c1->cd(1)->SetLogy(1);
    // sys
    sprintf(buffer,"htotalunc%i.root",(int)mass);
    TFile *file_sys = new TFile(buffer);
    sprintf(buffer,"htotalunc%i",(int)mass);
    TH1F *htotalunc950 = (TH1F*) file_sys->Get(buffer);
    htotalunc950->SetLineColor(2);
    htotalunc950->SetLineWidth(2);
    htotalunc950->GetYaxis()->SetRangeUser(0.04,5);
    htotalunc950->GetYaxis()->SetTitle("#DeltaN/N");
    htotalunc950->GetXaxis()->SetTitle("M_{t#gamma} (GeV/c^2)");
    htotalunc950->Draw();

    // stat.
    //TFile *file = new TFile("/Users/panda/files_IBM/ubuntu/CMS/WORKSPACE/T2TopGamma/Analysis/T2TopPlusPhoton/TemplateFit/DEV/data/102014/FittingOutput/output_theta_1.0.root");
    TFile *file = new TFile("../../FittingOutput/output_theta_1.0.root");

    TH1F *cate0 = (TH1F*) file->Get("data_cate-0");
    TH1F *unc = new TH1F("unc","",cate0->GetNbinsX(), cate0->GetXaxis()->GetXmin(), cate0->GetXaxis()->GetXmax());

    sprintf(buffer,"Tgamma%iGeV__pdf_bkg", mass);
    TH1F *hstd = (TH1F*) file->Get(buffer);

    double error = 0.; 
    double content_ = 0.; 
    for(int ibin = 1;ibin<= cate0->GetNbinsX();ibin++){
	content_ = cate0->IntegralAndError(ibin,ibin,error);
	unc->SetBinContent(ibin,error/content_); 
    }   
    unc->SetLineWidth(2);
    unc->Draw("same");

    TLegend *Unclegend_nm = new TLegend(0.577-0.4,0.305+0.45,0.5,0.9);
    Unclegend_nm->SetBorderSize(0);
    Unclegend_nm->SetFillColor(0);
    Unclegend_nm->SetFillStyle(0);
    Unclegend_nm->SetNColumns(1);
    Unclegend_nm->SetTextSize(0.04);
    Unclegend_nm->SetTextSizePixels(25);

    Unclegend_nm->AddEntry(htotalunc950,"syst. unc.","l");
    Unclegend_nm->AddEntry(unc,"stat. unc.","l");
    Unclegend_nm->Draw();

    c1->SaveAs("sysVSstat.pdf");
    printf("Total uncertainty bin-by-bin : sys(%1.2f%%) + stat(%1.2f%%)\n",
	100*GetRelativeUncInTotal(hstd, htotalunc950),
	100*GetRelativeUncInTotal(hstd, unc));
	    
}

float GetRelativeUncInTotal(TH1F *hSTD, TH1F *hUNC){
    if(hSTD->GetNbinsX()!=hUNC->GetNbinsX()){
	printf("[ERROR] Inconsistent Nbins : hSTD(%i) != hUNC(%i)\n", hSTD->GetNbinsX(), hUNC->GetNbinsX());
	return 0.;
    }
    float absUnc = 0.;
    for(int ibin=1; ibin<=hSTD->GetNbinsX(); ibin++){
	absUnc += (hSTD->GetBinContent(ibin) * hUNC->GetBinContent(ibin)) * 
	    (hSTD->GetBinContent(ibin) * hUNC->GetBinContent(ibin));
    }
    return sqrt(absUnc)/hSTD->Integral();
}
