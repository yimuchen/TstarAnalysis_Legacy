#include <math.h>
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "TLatex.h"
#include "TMath.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooAbsCollection.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "interface/samples_tr_1lepton4jets.h"
#include <iostream>
#include "RooWorkspace.h"
#include "TPad.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
//#include "interface/CMSStyle2015/CMS_lumi.h"

int RunStatus_ = Normal;
/*
   enum RunStatus{
    Normal,
    UncXsecPlus,
    UncXsecMinus,
    UncPUPlus,
    UncPUMinus,
    UncJESPlus,
    UncJESMinus,
    UncJERPlus,
    UncJERMinus,
    UncQsquare,     // for generating
    UncTrigPlus,
    UncTrigMinus,
    UncTopPtPlus,
    UncTopPtMinus,
    UncPhoIDPlus,
    UncPhoIDMinus,
    UncLepIDPlus,
    UncLepIDMinus,
    // splite for getting uncertainties on shapes
    UncQsquarePlus,
    UncQsquareMinus,
    UncMatchingPlus,
    UncMatchingMinus,
    RunStatusSize
   };

   string RunStatusNames[RunStatusSize] = {
   */
string RunStatusNamesOnlyForTheta[RunStatusSize] = { 
    "", 
    "Xsec__plus",
    "Xsec__minus",
    "PU__plus",
    "PU__minus",
    "JES__plus",
    "JES__minus",
    "JER__plus",
    "JER__minus",
    "Qsquare",   // for generating 
    "Trig__plus",
    "Trig__minus",
    "TopPt__plus",
    "TopPt__minus",
    "PhoID__plus",
    "PhoID__minus",
    "LepID__plus",
    "LepID__minus",
    // splite for getting uncertainties on shapes
    "Qsquare__plus",
    "Qsquare__minus",
    "Matching__plus",
    "Matching__minus"
};

using namespace RooFit;

float BR = 1.0;
/*
   IsIndividualBkgShape
   	false : normalization factor obtained using varied individual shapes + total varied shape (current, v)
   	true  : normalization factor + varied individual shapes
   */
bool IsIndividualBkgShape = false;	// "true" mode only has 0.005% difference on bkg estimation 
bool IsForThetaOrHiggsCombined = true; // true : output_theta.root for theta; false : output.root for HiggsComb
bool usingBiasStudy=false;
bool UsingLogyScaleInMass=true;
bool write_template=true;
bool usingSignalInjection=true;
bool usingCate1insteadofCate3 = false;
bool usingCate1insteadofCate3OnlyForSR = false;
bool usingExtrapolateForCRMCShape = true;
bool IsPhysicsProcessConsidered = false;
bool AppliedSysOnShape = true;
bool IsMoneyPlot = false;
bool IsMoneyPlotReduced = false;
bool SkipFit = false;
bool OnlyLoadTemplates = false;
int Index_Cate1or3 = 1;
int statistics = 10;
const int numberbins = 10;
float Xregions[2] = {50,1550};
//TLatex * latexCMS = new TLatex(0.48,0.96,"2012 CMS 19.7fb^{-1} #sqrt{s}=8TeV");

int _signal_start = Tgamma600GeV;
int _sample_start = (int)TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
float ranges[2] = {0.0,1.0};
const int NCate_ = 5;
int counters[NCate_] = {0,0,0,0,0};
int unc_on_icate_mode_tag[2] = {-1,-1};  // [icate, mode]   : icate (0~2 : cate0, cate3, cate4) and mode (0~2 : +, -, normal)
double purity[3][2] = { // [cate0,cate3, cate4][match,total]
    {0,0},
    {0,0},
    {0,0}
};
double purityError[3][2] = { // [cate0,cate3, cate4][match,total]
    {0,0},
    {0,0},
    {0,0}
};
double fractions_[Tgamma1500GeV-Tgamma500GeV+1][3];
double NBgkWithUncOnMerge[Tgamma1500GeV-Tgamma500GeV+1][3][3];  //[sample][template][down, normal, up]
double NBgkWithUncOnSRandCR[Tgamma1500GeV-Tgamma500GeV+1][3][3];  //[sample][template][down, normal, up]
double NBgkWithUncOnDATAandMC[Tgamma1500GeV-Tgamma500GeV+1][3][3];  //[sample][template][down, normal, up]

const int Nsamples_ =(int)samples_order_size - 
TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
int colors[Nsamples_] = {
    kRed-2,
    kBlue-2,
    kYellow-2,
    kOrange+7,
    kGreen-2,
    kViolet-2,
    kPink-2,
    kSpring-2,
    kRed+2,
    kBlue+2,
    kYellow+2,
    kOrange+2,
    kGreen+2,
    kViolet+2,
    kPink+2,
    kSpring+2,
    kRed-4,
    kBlue-4,
    kYellow-4,
    kOrange-4,
    kGreen-4,
    kViolet-4,
    kPink-4,
    kSpring-4,
    kRed-6,
    kBlue-6,
    kBlack
};

string CateTitles[NCate_] = {
    "cate0 : [photon + lepton]",
    "cate1 : [photon + jet]",
    "cate2 : [2 leptons]",
    "cate3 : [lepton + jet]",
    "cate4 : [2 jets]"
};
string outputfile[NCate_] = {
    "cateM",
    "cate1",
    "cate2",
    "cate3",
    "cate4"
};
string outputfileB[NCate_] = {
    "cate0",
    "cate1",
    "cate2",
    "cate3",
    "cate4"
};
string legendName[NCate_] = {
    "Cat M",
    "Cat 1",
    "Cat 2",
    "Cat 3",
    "Cat 4"
};
string legendNameB[NCate_] = {
    "Cat 0",
    "Cat 1",
    "Cat 2",
    "Cat 3",
    "Cat 4"
};
string outputfileCR[NCate_] = {
    "Control Regions#0~3",//"1 #gamma + 2 leptons #geq 4 jets",
    "1 #gamma + 1 leptons #geq 5 jets",
    "0 #gamma + 3 leptons #geq 4 jets",
    "0 #gamma + 2 leptons #geq 5 jets",
    "0 #gamma + 1 leptons #geq 6 jets"
};
char buffer[128];
char buffer2[128];
const int Ntmp = 53;
string patchNames[Ntmp] = {
    // TTbar group
        // sub-group 1
        "TTJets_HadronicMGDecays",
        "TTJets_SemiLeptMGDecays",
        "TTJets_FullLeptMGDecays",    // TTJets_FullLeptMGDecays for cate 4
        // sub-group 2
        "TTJets_MSDecays",
        // sub-group 3
        "TTJets_powheg",
        // sub-group 4
        "TT_700_1000_powheg",
        "TT_1000_inf_powheg",     // Mtt cut is a bias for cate4
        // sub-group 5
        "TTJets",
    // Z group
        // sub-group 1
        "DYJetsToLL_PtZ-100",     // DYJetsToLL_PtZ-100   for cate0
        "DYJetsToLL_PtZ-50To70",  // DYJetsToLL_PtZ-50To70 for cate 4
        "DYJetsToLL_PtZ-70To100",
        // sub-group 2
        "DYToEE_M-20",
        "DYToMuMu_M-20",
        // sub-group 3
        "DYJetsToLL_TuneZ2_M-50",
        "DYJetsToLL_TuneZ2_M-10To50",
        // sub-group 4
        "ZGToLLG",    // ZGToLLG  for cate0

    // Di-bosons group
        // sub-group 1
        "ZZJetsTo2L2Q",
        // sub-group 2
        "ZZTo2e2mu",
        "ZZTo2e2tau",
        "ZZTo2mu2tau",
        "ZZTo4e",
        "ZZTo4mu",
        "ZZTo4tau",
        // sub-group 3
        "DiPhotonJets",
        "WGToLNuG",
        "WJetsToLNu",
        "WW",         
        "WWGJets",
        "WWWJets",    // for cate4
        "WZ",         
        "ZZ",
    // Higgs group
    "ggToH2ZZTo4l",
    "ggToH2ZZTo2l2q",
    "VBFToH2ZZTo2l2q",
    "TTbarH_HToZZTo4L",
    "VBFToH2ZZTo4l",
    "WHToH2ZZTo4l",
    "ZHToH2ZZTo4l",

    // single top group
    "T_TuneZ2_s-channel",
    "T_TuneZ2_t-channel",
    "T_TuneZ2_tW-channel",  // T_TuneZ2_tW-channel for cate 4
    "Tbar_TuneZ2_s-channel",
    "Tbar_TuneZ2_t-channel",  
    "Tbar_TuneZ2_tW-channel",

    // gamma group
    "GJets_HT-40To100",
    "GJets_HT-100To200",
    "GJets_HT-200To400",
    "GJets_HT-400ToInf",

    // ttX group
    "ttW",
    "ttWW",
    "ttG",    // ttG  for cate0
    "ttH",
    "ttZ"
};

TFile *file = new TFile("MCTemplatesTree.root"); // much better
TFile *filepatch = new TFile("MCTemplatesTree.root");

TFile *output_;

TH1F *h_templates[samples_order_size][NCate_]; //[][photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_templatesunc[samples_order_size][NCate_]; //[][photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_templatesInSR[samples_order_size][NCate_]; //[][photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_templatesPatchInSR[Ntmp][NCate_]; //[][photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplates[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplatesPure[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplatesWOSysUnc[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplatesInSR[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplatesUnc[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MergeTemplatesUncInSR[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
TH1F *h_MCinSignalRegion;
TH1F *h_MCinSignalRegionCate4AndRest;
TH1F *h_MCinSignalRegion2Photons;
TH1F *h_TemplatesExtractedFromData[4]; // [cate0,cate1/cate3, cate4, signal region]   changed at 04/11/14
TH1F *h_TemplatesExtractedFromDataTemp[5]; // [cate0, cate1, cate2, cate3, cate4]   changed at 08/30/14
double fitparas[5][2];
TF1  *f_TemplatesExtractedFromDataTemp[5]; // [cate0, cate1, cate2, cate3, cate4]   changed at 08/30/14
TF1  *f_templatesPatchInSR[5]; // [cate0, cate1, cate2, cate3, cate4]   changed at 08/30/14
TH1F *h_SignalTemplateInSR[Tgamma1500GeV-Tgamma500GeV+1]; 

TCanvas *canvas[NCate_];
TCanvas *canvasUncOnDPP[NCate_];
TCanvas *canvasCompSRandCR[NCate_];
TCanvas *canvasCompDATAandMC[NCate_];
TCanvas *CheckTemplateCanvas = new TCanvas("CheckTemplateCanvas","",640,640);
TCanvas *SRCanvas[Tgamma1500GeV-Tgamma500GeV+1];
TCanvas *ClosureSRCanvas[Tgamma1500GeV-Tgamma500GeV+1];
TCanvas *MergeCanvas;
TCanvas *biasCanvas;
TCanvas *SignalInjectionCanvas;
TRandom3 *rnd = new TRandom3(123);

TLegend *legend_nm[NCate_];
TLegend *legend_nm_ComparisonSRandCR[NCate_];
TLegend *legend_nm_ComparisonDATAandMC[NCate_];
TLegend *Mergelegend_nm;

void randomTemplate(TH1F *hpercentagetemplate[],float percentage/*0~1*/);
void BiasStudy(TH1F *hdata, TH1F* hcate0, TH1F* hcate3, TH1F* hcate4, TH1F* hsig, TH1F* hbias[], int nToy);
void BiasStudy(TH1F *hdata, TH1F* hcate0, TH1F* hcate3, TH1F* hcate4, TH1F* hsig, TH1F* hbias, TCanvas* TCanvasTemp);
float DrawFit(int sidx_/*signal index*/,TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[],TH1F *h_MCinSignalRegion, RooWorkspace *wspace,
        bool mode/*0:data, 1:MC*/);
void DrawFitUsingMCShape(int sidx_/*signal index*/,TH1F *h_SignalTemplateInSR[],
	        TH1F *h_MCinSignalRegion, bool mode/*0:data, 1:MC*/);
void DrawBkgOnlyFit( TH1F *h_TemplatesExtractedFromData[] );
void SignalInjectionFit(int sidx_/*signal index*/,TH1F *h_MCinSignalRegion,TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[], TH1F *bias, int nToy);
TGraphErrors *SignalInjectionFitBias(int sidx_/*signal index*/,TH1F *h_MCinSignalRegion, TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[], int nToy);
void DrawGausFitResult(TF1 *func);
void loadTemplates();
void uncOnDifferentPhysicsProcesses();
void uncOnDataMCComparison();
void uncOnMerge();
void uncOnSignalAndControlRegions();
void patchForSR();
void Scale(TH1F *hist, TH1F *histtmp/*unc w.r.t.*/, TH1F* hunc, int mode);   // mode : 0(Up) 1(Down) 2(Normal)
void FitOnlyForDataWithUncShape(int sidx_, TH1F *hCate0, TH1F *hCate3,TH1F *hCate4, float *fracs);
void DrawRATIO(TH1F *numeritor, TH1F *denominator);
void DrawPULL(TH1F *numeritor, TH1F *denominator,string pullname);
void PrintUncTable();
void TrimChar(char *OrignalChar_, char *TrimChar_);
void drawCMSLumi();
TH1F *hshape(int M_, int iUncName);

TH1F *htotalunc950;	// only for money plot
// Systematics/SysOnBkgEstimate/htotalunc950.root

void TemplateFit_v21(){

    gROOT->ProcessLine(".L interface/setTDRStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetErrorX(0.5);


    bool WARNING_totalunc950 = false;
    if(AppliedSysOnShape){
	TFile *fhtotalunc950 = new TFile("Systematics/SysOnBkgEstimate/htotalunc950.root");
	if(!fhtotalunc950) {
	    printf("[Error] There is no Systematics/SysOnBkgEstimate/htotalunc950.root\n");
	    WARNING_totalunc950 = true;
	    htotalunc950 = new TH1F("htotalunc950","",numberbins,Xregions[0],Xregions[1]);
	    for(int ibin = 1; ibin<=numberbins; ibin++)
		htotalunc950->SetBinContent(ibin, 0.);
	}else{
	    htotalunc950 = (TH1F*) fhtotalunc950->Get("htotalunc950");
	}
    }

    MergeCanvas = new TCanvas("MergeCanvas","",640,640);
    biasCanvas = new TCanvas("biasCanvas","",640,640);
    SignalInjectionCanvas = new TCanvas("SignalInjectionCanvas","",640,640);

    // output file for theta limit or HiggsComine 
    sprintf(buffer,"FittingOutput/output.root");
    if(IsForThetaOrHiggsCombined) sprintf(buffer,"FittingOutput/output_theta_%1.1f.root",BR);
    output_ = new TFile(buffer,"recreate");

    if(!usingCate1insteadofCate3)
        Index_Cate1or3 = 3;
    int TemplateCateNumbers[3] = {0,(int)Index_Cate1or3,4};
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

    loadTemplates();

    TLegend *legend_Comparison2Photon = new TLegend(0.5077,0.705,0.9,0.927);
    TCanvas *Comparison2Photon = new TCanvas("Comparison2Photon","",640,640);
    Comparison2Photon->cd();
    Comparison2Photon->cd(1)->SetLogy(1);
    TH1F * h_MergeTemplates__0 = (TH1F*) h_MergeTemplates[0]->Clone();

    h_MergeTemplates__0->Scale(1./h_MergeTemplates__0->Integral());
    h_MCinSignalRegion2Photons->Scale(1./h_MCinSignalRegion2Photons->Integral());
    h_MergeTemplates__0->GetYaxis()->SetTitle("Arbitrary Unit");
    h_MergeTemplates__0->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
    h_MergeTemplates__0->Draw("error");
    h_MCinSignalRegion2Photons->Draw("same");
    legend_Comparison2Photon->AddEntry(h_MergeTemplates__0,"Cat M(CR)","l");
    legend_Comparison2Photon->AddEntry(h_MCinSignalRegion2Photons,"2 photons Bkg.(SR)","l");
    legend_Comparison2Photon->SetBorderSize(0);
    legend_Comparison2Photon->SetFillColor(0);
    legend_Comparison2Photon->SetFillStyle(0);
    legend_Comparison2Photon->SetNColumns(1);
    legend_Comparison2Photon->SetTextSize(0.04);
    legend_Comparison2Photon->SetTextSizePixels(25);
    legend_Comparison2Photon->Draw();
    Comparison2Photon->SaveAs("FittingOutput/Comparison2Photon.pdf");


    if(OnlyLoadTemplates)
	return;
    // Draw plots in 5 cates
    for(int i=0;i<NCate_;i++){
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/%s.pdf",outputfileB[i].c_str());
	else
	    sprintf(buffer,"FittingOutput/%s_%1.1f.pdf",outputfileB[i].c_str(),BR);
        canvas[i]->SaveAs(buffer);
    }

    // Comparison among different-categories templates
    MergeCanvas->cd();
    for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
        h_MergeTemplates[idxCate_]->Scale(1./h_MergeTemplates[idxCate_]->Integral());
	h_MergeTemplatesWOSysUnc[idxCate_] = (TH1F*) h_MergeTemplatesPure[idxCate_]->Clone();	// with syst.+stat.
        h_MergeTemplatesWOSysUnc[idxCate_]->Scale(1./h_MergeTemplatesWOSysUnc[idxCate_]->Integral());
        if(idxCate_==0){
            h_MergeTemplatesWOSysUnc[idxCate_]->GetYaxis()->SetRangeUser(ranges[0],ranges[1]);
            h_MergeTemplatesWOSysUnc[idxCate_]->GetYaxis()->SetTitle("Arbitrary Unit");
            h_MergeTemplatesWOSysUnc[idxCate_]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
            h_MergeTemplatesWOSysUnc[idxCate_]->Draw("hist,error");
        }else{
            h_MergeTemplatesWOSysUnc[idxCate_]->Draw("hist,error,same");
        }
        sprintf(buffer,"%s",legendNameB[idxCate_].c_str());
        Mergelegend_nm->AddEntry(h_MergeTemplatesWOSysUnc[idxCate_],buffer,"l");
        Mergelegend_nm->Draw();
    }
    if(!usingExtrapolateForCRMCShape){
        h_MergeTemplates[0]->Add(h_MergeTemplates[1]);
        h_MergeTemplates[0]->Add(h_MergeTemplates[2]);
        h_MergeTemplates[0]->Add(h_MergeTemplates[3]);
        h_MergeTemplates[0]->Scale(1./h_MergeTemplates[0]->Integral());
    }

    // (Cate0+Cate2)  + (Cate1+Cate3)
    //h_MergeTemplates[0]->Add(h_MergeTemplates[2]);
    //h_MergeTemplates[0]->Scale(1./h_MergeTemplates[0]->Integral());

    //h_MergeTemplates[1]->Add(h_MergeTemplates[3]);
    //h_MergeTemplates[1]->Scale(1./h_MergeTemplates[1]->Integral());

    //h_MergeTemplates[0]->Add(h_MergeTemplates[1]);
    //h_MergeTemplates[0]->Scale(1./h_MergeTemplates[0]->Integral());

    if(BR==1.0){
	MergeCanvas->SaveAs("FittingOutput/Combined_cates.pdf");
	MergeCanvas->SaveAs("FittingOutput/Combined_cates.root");

	h_MergeTemplatesWOSysUnc[0]->GetYaxis()->SetRangeUser(0.000001,1000);
	MergeCanvas->cd(1)->SetLogy(1);
	Mergelegend_nm->SetX1NDC(0.577);
	Mergelegend_nm->SetX2NDC(0.90);
	Mergelegend_nm->SetY1NDC(0.555);
	Mergelegend_nm->SetY2NDC(0.87);
	Mergelegend_nm->Draw();
	MergeCanvas->Modified();
	MergeCanvas->Update();
	MergeCanvas->SaveAs("FittingOutput/Combined_catesInLog.pdf");
    }else{
	sprintf(buffer,"FittingOutput/Combined_cates_%1.1f.pdf",BR);
	MergeCanvas->SaveAs(buffer);
	sprintf(buffer,"FittingOutput/Combined_cates_%1.1f.root",BR);
	MergeCanvas->SaveAs(buffer);
    }

    // add more statistics of ttbar and DY in SR
    patchForSR();
    // Comparison between SR and CR
    for(int i=0;i<=4;i++)
        if(h_MergeTemplatesInSR[i]->Integral()>0)
            h_MergeTemplatesInSR[i]->Scale(1./h_MergeTemplatesInSR[i]->Integral());

    h_MergeTemplatesInSR[0]->Add(h_MergeTemplatesInSR[1]);
    h_MergeTemplatesInSR[0]->Add(h_MergeTemplatesInSR[2]);
    h_MergeTemplatesInSR[0]->Add(h_MergeTemplatesInSR[3]);

    for(int i=0;i<NCate_;i++){
        output_->cd();
        sprintf(buffer,"%s_AfterMerging",h_MergeTemplatesInSR[i]->GetName());
        h_MergeTemplatesInSR[i]->SetName(buffer);
        h_MergeTemplatesInSR[i]->Write();
    }

    printf("[DEBUG CATE4] first bin : %f\n",h_MergeTemplatesInSR[4]->GetBinContent(1));

    for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
        int idxCate_SR = idxCate_;
        if(usingCate1insteadofCate3OnlyForSR && idxCate_==3)
            idxCate_SR = 1;

        float padtop_ = 0.15;
        float padbottom_ = padtop_ + 0.043;
        TPad *tp2 = new TPad("tp2","",0,0,1,padbottom_);
        TPad *tp1 = new TPad("tp1","",0,padtop_,1,1);

        canvasCompSRandCR[idxCate_]->cd();

        tp1->Draw();
        tp1->cd();
        tp1->SetBottomMargin(0.05);
        tp1->SetTopMargin(0.06);

        h_MergeTemplatesInSR[idxCate_SR]->Scale(1./h_MergeTemplatesInSR[idxCate_SR]->Integral());

        double min_ = 1.;
        double mintemp_ = 1.;
        for(int ibin = 1;ibin <= h_MergeTemplates[idxCate_]->GetNbinsX();ibin++){
            mintemp_ = h_MergeTemplates[idxCate_]->GetBinContent(ibin);
            if(min_>mintemp_ && mintemp_>0.00001) 
                min_ = mintemp_;
        }
        for(int ibin = 1;ibin <= h_MergeTemplatesInSR[idxCate_SR]->GetNbinsX();ibin++){
            mintemp_ = h_MergeTemplatesInSR[idxCate_SR]->GetBinContent(ibin);
            if(min_>mintemp_ && mintemp_>0.00001) 
                min_ = mintemp_;
        }
        std::cout<<"[JackyCate"<< idxCate_<<"] : min_ "<< min_<<std::endl;

        h_MergeTemplates[idxCate_]->GetYaxis()->SetRangeUser(0.1*min_,1.0);
        //h_MergeTemplates[idxCate_]->GetYaxis()->SetRangeUser(ranges[0],ranges[1]);
        h_MergeTemplates[idxCate_]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
        h_MergeTemplates[idxCate_]->GetYaxis()->SetTitle("Arbitrary Unit");
        h_MergeTemplates[idxCate_]->Draw();

        h_MergeTemplatesInSR[idxCate_SR]->Draw("same");

        sprintf(buffer,"%s(CRs)",legendName[idxCate_].c_str());
        legend_nm_ComparisonSRandCR[idxCate_]->AddEntry(h_MergeTemplates[idxCate_],buffer,"l");
        sprintf(buffer,"%s(SR)",legendName[idxCate_SR].c_str());
        legend_nm_ComparisonSRandCR[idxCate_]->AddEntry(h_MergeTemplatesInSR[idxCate_SR],buffer,"l");
        legend_nm_ComparisonSRandCR[idxCate_]->Draw();

        tp1->SetLogy(1);

        canvasCompSRandCR[idxCate_]->cd();
        tp2->Draw();
        tp2->cd();
        tp2->SetTopMargin(0.0);
        tp2->SetBottomMargin(0.4);
        //DrawRATIO(h_MergeTemplatesInSR[idxCate_SR],h_MergeTemplates[idxCate_]);
        //DrawRATIO(h_MergeTemplates[idxCate_],h_MergeTemplatesInSR[idxCate_SR]);
        //DrawPULL(h_MergeTemplates[idxCate_],h_MergeTemplatesInSR[idxCate_SR]);
        DrawPULL(h_MergeTemplatesInSR[idxCate_SR],h_MergeTemplates[idxCate_],"#frac{SR-CRs}{#Delta_{CRs}}");

	if(BR==1.0){
	    sprintf(buffer,"FittingOutput/Comparison_between_SRandCR_%s.pdf",outputfileB[idxCate_].c_str());
	    canvasCompSRandCR[idxCate_]->SaveAs(buffer);

	    sprintf(buffer,"FittingOutput/Comparison_between_SRandCR_%s.root",outputfileB[idxCate_].c_str());
	    canvasCompSRandCR[idxCate_]->SaveAs(buffer);
	}else{
	    sprintf(buffer,"FittingOutput/Comparison_between_SRandCR_%s_%1.1f.pdf",outputfileB[idxCate_].c_str(),BR);
	    canvasCompSRandCR[idxCate_]->SaveAs(buffer);

	    sprintf(buffer,"FittingOutput/Comparison_between_SRandCR_%s_%1.1f.root",outputfileB[idxCate_].c_str(),BR);
	    canvasCompSRandCR[idxCate_]->SaveAs(buffer);
	}
    }

    std::cout<<"Purity for cate0 : "<<purity[0][1]/purity[0][0]*100.<<"+/-"<< 
       (purity[0][1]/purity[0][0]*(1-purity[0][1]/purity[0][0]))/purity[0][0]*100. <<"%"<<std::endl;
    std::cout<<"Purity for cate3 : "<<purity[1][1]/purity[1][0]*100.<<"+/-"<<
        (purity[1][1]/purity[1][0]*(1-purity[1][1]/purity[1][0]))/purity[1][0]*100.<<"%"<<std::endl;
    std::cout<<"Purity for cate4 : "<<purity[2][1]/purity[2][0]*100.<<"+/-"<<
        (purity[2][1]/purity[2][0]*(1-purity[2][1]/purity[2][0]))/purity[2][0]*100.<<"%"<<std::endl;


    // Templates : cate0(data), cate1/cate3(data), cate4(data)
    CheckTemplateCanvas->Divide(1,3);
    CheckTemplateCanvas->cd(1);
    TH1F *h_TemplatesExtractedFromData__0 = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
    h_TemplatesExtractedFromData__0->Scale(1./h_TemplatesExtractedFromData__0->Integral());
    h_TemplatesExtractedFromData__0->GetYaxis()->SetRangeUser(0,1);
    h_TemplatesExtractedFromData__0->Draw();
    h_MergeTemplates[0]->Draw("same");

    CheckTemplateCanvas->cd(2);
    TH1F *h_TemplatesExtractedFromData__1 = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
    h_TemplatesExtractedFromData__1->Scale(1./h_TemplatesExtractedFromData__1->Integral());
    h_TemplatesExtractedFromData__1->GetYaxis()->SetRangeUser(0,1);
    h_TemplatesExtractedFromData__1->Draw();
    h_MergeTemplates[Index_Cate1or3]->Draw("same");

    CheckTemplateCanvas->cd(3);
    TH1F *h_TemplatesExtractedFromData__2 = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
    h_TemplatesExtractedFromData__2->Scale(1./h_TemplatesExtractedFromData__2->Integral());
    h_TemplatesExtractedFromData__2->GetYaxis()->SetRangeUser(0,1);
    h_TemplatesExtractedFromData__2->Draw();
    h_MergeTemplates[4]->Draw("same");
    if(BR==1.0)
	CheckTemplateCanvas->SaveAs("FittingOutput/comparison_dataMC.pdf");
    else{
	sprintf(buffer,"FittingOutput/comparison_dataMC_%1.1f.pdf",BR);
	CheckTemplateCanvas->SaveAs(buffer);
    }

    for(int idxCate_ = 0;idxCate_<3;idxCate_++){

        float padtop_ = 0.15;
        float padbottom_ = padtop_ + 0.043;
        TPad *tp2 = new TPad("tp2","",0,0,1,padbottom_);
        TPad *tp1 = new TPad("tp1","",0,padtop_,1,1);

        int idxTempCate_ = TemplateCateNumbers[idxCate_];
        canvasCompDATAandMC[idxTempCate_]->cd();

        tp1->Draw();
        tp1->cd();
        tp1->SetBottomMargin(0.05);
        tp1->SetTopMargin(0.06);

        h_MergeTemplates[idxTempCate_]->GetYaxis()->SetRangeUser(ranges[0],ranges[1]);
        h_MergeTemplates[idxTempCate_]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
        h_MergeTemplates[idxTempCate_]->GetYaxis()->SetTitle("Arbitrary Unit");

        double min_ = 1.;
        double mintemp_ = 1.;
        //for(int ibin = 1;ibin <= h_MergeTemplates[idxTempCate_]->GetNbinsX();ibin++){
        for(int ibin = 1;ibin < h_MergeTemplates[idxTempCate_]->GetNbinsX();ibin++){
            mintemp_ = h_MergeTemplates[idxTempCate_]->GetBinContent(ibin);
            if(min_>mintemp_ && mintemp_>0.000001) 
                min_ = mintemp_;
        }
        for(int ibin = 1;ibin <= h_TemplatesExtractedFromData[idxCate_]->GetNbinsX();ibin++){
            mintemp_ = h_TemplatesExtractedFromData[idxCate_]->GetBinContent(ibin);
            if(min_>mintemp_ && mintemp_>0.000001) 
                min_ = mintemp_;
        }


        h_MergeTemplates[idxTempCate_]->GetYaxis()->SetRangeUser(0.9*min_,1.0);
        h_MergeTemplates[idxTempCate_]->Draw();

        TH1F *h_TemplatesExtractedFromData__ = (TH1F*) h_TemplatesExtractedFromData[idxCate_]->Clone();
        h_TemplatesExtractedFromData__->Sumw2();
        h_TemplatesExtractedFromData__->Scale(1./h_TemplatesExtractedFromData__->Integral());
        h_TemplatesExtractedFromData__->Draw("same");

        legend_nm_ComparisonDATAandMC[idxTempCate_]->AddEntry(h_TemplatesExtractedFromData__,"Data","l");
        sprintf(buffer,"MC (%s)",legendName[idxTempCate_].c_str());
        legend_nm_ComparisonDATAandMC[idxTempCate_]->AddEntry(h_MergeTemplates[idxTempCate_],buffer,"l");
        legend_nm_ComparisonDATAandMC[idxTempCate_]->Draw();

        sprintf(buffer,"(%s)",outputfileCR[idxTempCate_].c_str());
        TLatex *latex_ = new TLatex(0.587,0.76,buffer);
        latex_->SetTextSize(0.03);
        latex_->SetTextColor(kBlack);
        latex_->SetNDC();
        latex_->Draw();

        tp1->SetLogy(1);

        canvasCompDATAandMC[idxTempCate_]->cd();
        tp2->Draw();
        tp2->cd();
        tp2->SetTopMargin(0.0);
        tp2->SetBottomMargin(0.4);
        //DrawRATIO(h_TemplatesExtractedFromData__, h_MergeTemplates[idxTempCate_]);
        //DrawRATIO(h_MergeTemplates[idxTempCate_], h_TemplatesExtractedFromData__);
        DrawPULL(h_TemplatesExtractedFromData__, h_MergeTemplates[idxTempCate_],"#frac{data-MC}{#Delta_{MC}}");

	if(BR==1.0){
	    sprintf(buffer,"FittingOutput/Comparison_between_DATAandMC_%s.pdf",outputfileB[idxTempCate_].c_str());
	    canvasCompDATAandMC[idxTempCate_]->SaveAs(buffer);

	    sprintf(buffer,"FittingOutput/Comparison_between_DATAandMC_%s.root",outputfileB[idxTempCate_].c_str());
	    canvasCompDATAandMC[idxTempCate_]->SaveAs(buffer);
	}else{
	    sprintf(buffer,"FittingOutput/Comparison_between_DATAandMC_%s_%1.1f.pdf",outputfileB[idxTempCate_].c_str(),BR);
	    canvasCompDATAandMC[idxTempCate_]->SaveAs(buffer);

	    sprintf(buffer,"FittingOutput/Comparison_between_DATAandMC_%s_%1.1f.root",outputfileB[idxTempCate_].c_str(),BR);
	    canvasCompDATAandMC[idxTempCate_]->SaveAs(buffer);
	}

        delete h_TemplatesExtractedFromData__;
        delete latex_;
        h_MergeTemplates[idxTempCate_]->SetTitle("");
        delete tp1;
        delete tp2;
    }

    if(SkipFit)
	return;

    /*
    Templates : 
        cate0(data)   :   h_TemplatesExtractedFromData[0] 
        cate1/3(data) :   h_TemplatesExtractedFromData[1]
        cate4(data)   :   h_TemplatesExtractedFromData[2]
        SR(data)      :   h_TemplatesExtractedFromData[3]
        SR(MC)        :   h_MCinSignalRegion
        SR(S 950)     :   h_SignalTemplateInSR
    */
    // Fit on MC (Closure test) using MC shape
    for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
	ClosureSRCanvas[sidx_]->cd();
	DrawFitUsingMCShape(sidx_/*signal index*/,h_SignalTemplateInSR,
		h_MCinSignalRegion, 1/*0:data, 1:MC*/);
	if(UsingLogyScaleInMass)
	    ClosureSRCanvas[sidx_]->cd()->SetLogy(1);
	ClosureSRCanvas[sidx_]->Update();
	if(BR==1.0){
	    sprintf(buffer,"FittingOutput/fitWith%s_closure_test_MC_shape.pdf",SAMPLE[Tgamma500GeV+sidx_].tag);
	}else{
	    sprintf(buffer,"FittingOutput/fitWith%s_closure_test_MC_shape_%1.1f.pdf",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
	}
	ClosureSRCanvas[sidx_]->SaveAs(buffer);
    }

    // Fit on MC (Closure test)
    for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
        ClosureSRCanvas[sidx_]->cd();
        float fitfbkg = DrawFit(sidx_,h_SignalTemplateInSR, h_TemplatesExtractedFromData,h_MCinSignalRegion,wspace,1);
        if(UsingLogyScaleInMass)
            ClosureSRCanvas[sidx_]->cd()->SetLogy(1);
        ClosureSRCanvas[sidx_]->Update();
	if(BR==1.0){
	    sprintf(buffer,"FittingOutput/fitWith%s_closure_test.pdf",SAMPLE[Tgamma500GeV+sidx_].tag);
	}else{
	    sprintf(buffer,"FittingOutput/fitWith%s_closure_test_%1.1f.pdf",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
	}
	ClosureSRCanvas[sidx_]->SaveAs(buffer);
    }
    // Fit on data
    for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
        SRCanvas[sidx_]->cd();
        SRCanvas[sidx_]->Update();
	IsMoneyPlot = false;
	IsMoneyPlotReduced = false;
        float fitfbkg = DrawFit(sidx_/*signal index*/,h_SignalTemplateInSR, h_TemplatesExtractedFromData,h_MCinSignalRegion,wspace,0);
        if(UsingLogyScaleInMass)
            SRCanvas[sidx_]->cd()->SetLogy(1);
        SRCanvas[sidx_]->Update();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/fitWith%s.pdf",SAMPLE[Tgamma500GeV+sidx_].tag);
	else
	    sprintf(buffer,"FittingOutput/fitWith%s_%1.1f.pdf",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
        SRCanvas[sidx_]->SaveAs(buffer);

        if(sidx_==(Tgamma950GeV - Tgamma500GeV)){
	    if(BR==1.0)
		sprintf(buffer,"FittingOutput/fitWith%s.root",SAMPLE[Tgamma500GeV+sidx_].tag);
	    else
		sprintf(buffer,"FittingOutput/fitWith%s_%1.1f.root",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
            SRCanvas[sidx_]->SaveAs(buffer);
        }
    }
    // Fit on data (money plot)
    for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
        SRCanvas[sidx_]->cd();
        SRCanvas[sidx_]->Update();
	IsMoneyPlot = true;
	IsMoneyPlotReduced = false;
        float fitfbkg = DrawFit(sidx_/*signal index*/,h_SignalTemplateInSR, h_TemplatesExtractedFromData,h_MCinSignalRegion,wspace,0);
	IsMoneyPlot = false;
	IsMoneyPlotReduced = false;
        if(UsingLogyScaleInMass)
            SRCanvas[sidx_]->cd()->SetLogy(1);
        SRCanvas[sidx_]->Update();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/fitWith%s_moneyplot.pdf",SAMPLE[Tgamma500GeV+sidx_].tag);
	else
	    sprintf(buffer,"FittingOutput/fitWith%s_%1.1f_moneyplot.pdf",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
        SRCanvas[sidx_]->SaveAs(buffer);

        if(sidx_==(Tgamma950GeV - Tgamma500GeV)){
	    if(BR==1.0)
		sprintf(buffer,"FittingOutput/fitWith%s_moneyplot.root",SAMPLE[Tgamma500GeV+sidx_].tag);
	    else
		sprintf(buffer,"FittingOutput/fitWith%s_%1.1f_moneyplot.root",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
            SRCanvas[sidx_]->SaveAs(buffer);
        }
    }
    // Fit on data (money plot --> reduced)
    for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
        SRCanvas[sidx_]->cd();
        SRCanvas[sidx_]->Update();
	IsMoneyPlot = true;
	IsMoneyPlotReduced = true;
        float fitfbkg = DrawFit(sidx_/*signal index*/,h_SignalTemplateInSR, h_TemplatesExtractedFromData,h_MCinSignalRegion,wspace,0);
	IsMoneyPlot = false;
	IsMoneyPlotReduced = false;
        if(UsingLogyScaleInMass)
            SRCanvas[sidx_]->cd()->SetLogy(1);
        SRCanvas[sidx_]->Update();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/fitWith%s_moneyplot_reduced.pdf",SAMPLE[Tgamma500GeV+sidx_].tag);
	else
	    sprintf(buffer,"FittingOutput/fitWith%s_%1.1f_moneyplot_reduced.pdf",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
        SRCanvas[sidx_]->SaveAs(buffer);

	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/fitWith%s_moneyplot_reduced.png",SAMPLE[Tgamma500GeV+sidx_].tag);
	else
	    sprintf(buffer,"FittingOutput/fitWith%s_%1.1f_moneyplot_reduced.png",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
        SRCanvas[sidx_]->SaveAs(buffer);

        if(sidx_==(Tgamma950GeV - Tgamma500GeV)){
	    if(BR==1.0)
		sprintf(buffer,"FittingOutput/fitWith%s_moneyplot_reduced.root",SAMPLE[Tgamma500GeV+sidx_].tag);
	    else
		sprintf(buffer,"FittingOutput/fitWith%s_%1.1f_moneyplot_reduced.root",SAMPLE[Tgamma500GeV+sidx_].tag,BR);
            SRCanvas[sidx_]->SaveAs(buffer);
        }
    }
    if(write_template){
	if(BR==1.0)
	    wspace->writeToFile("FittingOutput/workspace.root");
	else{
	    sprintf(buffer,"FittingOutput/workspace_%1.1f.root",BR);
	    wspace->writeToFile(buffer);
	}
    }

    TCanvas *ConlyBkgFit = new TCanvas("ConlyBkgFit","",640,640);
    ConlyBkgFit->cd();
    ConlyBkgFit->Update();
    if(UsingLogyScaleInMass)
	ConlyBkgFit->cd()->SetLogy(1);
    ConlyBkgFit->Update();
    DrawBkgOnlyFit( h_TemplatesExtractedFromData );
    if(BR==1.0){
	ConlyBkgFit->SaveAs("FittingOutput/BkgOnlyFit.pdf");
	ConlyBkgFit->SaveAs("FittingOutput/BkgOnlyFit.root");
    }else{
	sprintf(buffer,"FittingOutput/BkgOnlyFit_%1.1f.pdf",BR);
	ConlyBkgFit->SaveAs(buffer);
	sprintf(buffer,"FittingOutput/BkgOnlyFit_%1.1f.root",BR);
	ConlyBkgFit->SaveAs(buffer);
    }

    // Bias Study
    if(usingBiasStudy){
        biasCanvas->cd();
        TH1F* hbias[5]; // [bias total/fbkg0/fbkg3/fbkg4, chi2] 
        for(int i=0;i<4;i++){
            sprintf(buffer,"hbais_%i",i);
            if(i==0){
                hbias[i] = new TH1F(buffer,"",25,-1,1);
            }else{
                hbias[i] = new TH1F(buffer,"",100,-3,3);
            }
        }
        hbias[4] = new TH1F("hchi2","",100,0,100);
        BiasStudy(
                //h_MCinSignalRegion,
                h_TemplatesExtractedFromData[3],
                h_TemplatesExtractedFromData[0], 
                h_TemplatesExtractedFromData[1], 
                h_TemplatesExtractedFromData[2], 
                h_SignalTemplateInSR[Tgamma950GeV-Tgamma500GeV], hbias, 1000);
        hbias[0]->GetXaxis()->SetTitle("(Normal(fbkg) - Toy(fbkg) ) / Normal(fbkg)");
        hbias[0]->Fit("gaus");
        hbias[0]->Draw();
        biasCanvas->Update();
	if(BR==1.0){
	    biasCanvas->SaveAs("FittingOutput/bias_950.root");
	    biasCanvas->SaveAs("FittingOutput/bias_950.pdf");
	}else{
	    sprintf(buffer,"FittingOutput/bias_950_%1.1f.root",BR);
	    biasCanvas->SaveAs(buffer);
	    sprintf(buffer,"FittingOutput/bias_950_%1.1f.pdf",BR);
	    biasCanvas->SaveAs(buffer);
	}

        TF1 *func = new TF1("func","gaus",-0.5,0.5);
        int fracIdex[3] = {0,3,4};

        for(int fracIdex_ =0;fracIdex_<3;fracIdex_++){
            func->SetParameter(0,hbias[fracIdex_+1]->GetMaximum());
            func->SetParameter(1,hbias[fracIdex_+1]->GetMean());
            func->SetParameter(2,3*hbias[fracIdex_+1]->GetRMS());
            sprintf(buffer,"(Normal(fbgk%i) - Toy(fbkg%i) ) / #sigma (fbgk%i)",fracIdex[fracIdex_],fracIdex[fracIdex_],fracIdex[fracIdex_]);

            hbias[fracIdex_+1]->GetYaxis()->SetRangeUser(hbias[fracIdex_+1]->GetMinimum(),1.2*hbias[fracIdex_+1]->GetMaximum());
            hbias[fracIdex_+1]->GetXaxis()->SetTitle(buffer);
            hbias[fracIdex_+1]->Fit("func","M","",-3*hbias[fracIdex_+1]->GetRMS(),3*hbias[fracIdex_+1]->GetRMS());
            hbias[fracIdex_+1]->Draw();
            DrawGausFitResult(func);

            biasCanvas->Update();
	    if(BR==1.0)
		sprintf(buffer,"FittingOutput/bias_fbkg%i_950.root",fracIdex[fracIdex_]);
	    else
		sprintf(buffer,"FittingOutput/bias_fbkg%i_950_%1.1f.root",fracIdex[fracIdex_],BR);
            biasCanvas->SaveAs(buffer);
	    if(BR==1.0)
		sprintf(buffer,"FittingOutput/bias_fbkg%i_950.pdf",fracIdex[fracIdex_]);
	    else
		sprintf(buffer,"FittingOutput/bias_fbkg%i_950_%1.1f.pdf",fracIdex[fracIdex_],BR);
            biasCanvas->SaveAs(buffer);
        }

        hbias[4]->GetXaxis()->SetTitle("Chi2 [bias study]");
        hbias[4]->Draw();
        biasCanvas->Update();
	if(BR==1.0)
	    biasCanvas->SaveAs("FittingOutput/bias_chi2.pdf");
	else{
	    sprintf(buffer,"FittingOutput/bias_chi2_%1.1f.pdf",BR);
	    biasCanvas->SaveAs(buffer);
	}
    }

    if (usingSignalInjection){
        gROOT->ProcessLine(".L interface/setTDRStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetErrorX(0.5);
        SignalInjectionCanvas->cd();
        SignalInjectionCanvas->cd()->SetLogy(1);
        TGraphErrors *hbias = 
            //SignalInjectionFitBias(Tgamma750GeV-Tgamma500GeV,h_MCinSignalRegion, 
            SignalInjectionFitBias(Tgamma950GeV-Tgamma500GeV,h_MCinSignalRegion, 
                h_SignalTemplateInSR,h_TemplatesExtractedFromData, 50);
	if(BR==1.0){
	    SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection950GeV.pdf");
	    //SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection750GeV.pdf");
	}else{
	    sprintf(buffer,"FittingOutput/SignalInjection950GeV_%1.1f.pdf",BR);
	    SignalInjectionCanvas->SaveAs(buffer);
	    //SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection750GeV.pdf");
	}

        /*
        TF1 *func = new TF1("func","gaus",-1,1);
        hbias->SetTitle("Pull");
        hbias->GetXaxis()->SetTitle("(fitting_Sig_yields - inject_yields)/ #sigma_{fitting_Sig_yields}");
        hbias->Draw();
        //hbias->Fit("func");
        hbias->Fit("func","M","",-3*hbias->GetRMS(),3*hbias->GetRMS());
        DrawGausFitResult(func);
        */
        SignalInjectionCanvas->cd()->SetLogy(0);
        hbias->SetMarkerStyle(22);
        hbias->SetFillStyle(3002);
        hbias->SetFillColor(kBlue-7);
        hbias->SetTitle("");
        hbias->GetXaxis()->SetTitle("Injection [#times Signal]");
        hbias->GetYaxis()->SetTitle("Fitting Result [#times Signal]");
        SignalInjectionCanvas->SetGridx(1);
        SignalInjectionCanvas->SetGridy(1);
        hbias->Draw("APL, error3");

	if(BR==1.0){
	    SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection950GeV_bias.pdf");
	    SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection950GeV_bias.root");
	}else{
	    sprintf(buffer,"FittingOutput/SignalInjection950GeV_bias_%1.1f.pdf",BR);
	    SignalInjectionCanvas->SaveAs(buffer);
	    sprintf(buffer,"FittingOutput/SignalInjection950GeV_bias_%1.1f.root",BR);
	    SignalInjectionCanvas->SaveAs(buffer);
	}

        //SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection750GeV_bias.pdf");
        //SignalInjectionCanvas->SaveAs("FittingOutput/SignalInjection750GeV_bias.root");
    }

    /* uncertainties on templates :
       - different physics processes 
            h_templates[samples_order_size][NCate_]; //[][photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
            h_MergeTemplates[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
        
       - reduction of templates 
            h_MergeTemplates[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
                cate0 <-> cate2
                cate1 <-> cate3

       - data and MC
            h_MergeTemplates[NCate_]; //[photon + lepton, photon + jet, 2 leptons, lepton + jet, 2 jets]
            h_TemplatesExtractedFromData[4]; // [cate0,cate1/3, cate4, signal region]
       */
    // 1). uncertainty on different physics processes
    //uncOnDifferentPhysicsProcesses();

    // 2). uncertainty on signal and control regions
    //uncOnSignalAndControlRegions();

    //// 3). uncertainty on Merge
    //uncOnMerge();

    //// 4). uncertainty on DataMCComparison
    //uncOnDataMCComparison();

    //// For uncertainty table
    //PrintUncTable();

    if(WARNING_totalunc950)
	printf("[WARNING MESSAGE] There is no Systematics/SysOnBkgEstimate/htotalunc950.root\n");
}

void BiasStudy(TH1F *hdata, TH1F* hcate0, TH1F* hcate3, TH1F* hcate4, TH1F* hsig, TH1F* hbias[], int nToy){
        

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data","dataset with x",x,hdata);

        RooDataHist dCate0("dCate0","dCate0set with x",x,hcate0);
        RooHistPdf Cate0("Cate0","Cate0",x,dCate0,0) ;

        RooDataHist dCate3("dCate3","dCate3set with x",x,hcate3);
        RooHistPdf Cate3("Cate3","Cate3",x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,hcate4);
        RooHistPdf Cate4("Cate4","Cate4",x,dCate4,0) ;

        RooDataHist dSignal("dSignal","dSignalset with x",x,hsig);
        RooHistPdf Signal("Signal","Signal",x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.0,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.0,1) ; 
        RooRealVar fbkg4("fbkg4","fraction",0.5,0.0,1) ; 
        // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
        RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
        model.fitTo(data,Minimizer("Minuit2", "minimize"));

        //model.fitTo(data,InitialHesse(kTRUE),Hesse(kTRUE));

        //model.chi2FitTo(data,Extended(kTRUE));
        //model.fitTo(data,Extended(kTRUE));

        float Norm_fbkg = 1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal()) ;

        float Norm_fbkg0= fbkg0.getVal(); 
        float Norm_fbkg3= (1.-fbkg0.getVal())*(fbkg3.getVal()); 
        float Norm_fbkg4= (1.-fbkg0.getVal())*(1.-fbkg3.getVal())*fbkg4.getVal(); 

        float Norm_fbkg0_error= fbkg0.getError(); 
        float Norm_fbkg3_error= sqrt( fbkg3.getError()*fbkg3.getError() + 
                fbkg3.getVal()*fbkg3.getVal()*fbkg0.getError()*fbkg0.getError() +
                fbkg0.getVal()*fbkg0.getVal()*fbkg3.getError()*fbkg3.getError() 
                );
        float Norm_fbkg4_error=  sqrt( 
                fbkg4.getError()*fbkg4.getError() +
                (fbkg0.getVal()*fbkg0.getVal()+1)*(fbkg4.getVal()*fbkg4.getVal()*fbkg3.getError()*fbkg3.getError() +
                    fbkg3.getVal()*fbkg3.getVal()*fbkg4.getError()*fbkg4.getError() ) +
                (fbkg4.getVal()*fbkg4.getVal()*fbkg0.getError()*fbkg0.getError() + 
                 fbkg0.getVal()*fbkg0.getVal()*fbkg4.getError()*fbkg4.getError() ) +
                fbkg3.getVal()*fbkg3.getVal()*fbkg4.getVal()*fbkg4.getVal()*fbkg0.getError()*fbkg0.getError()
                );

        std::cout<<"[Jacky test] Norm_fbkg0 : "<<Norm_fbkg0<<std::endl; 

        for(int itoy=0;itoy<nToy;itoy++){
            /*
               cate 0 <-> hcate0
               cate 3 <-> hcate3
               cate 4 <-> hcate4
               randomTemplate(hcate_toy,0.67);
               */

            TH1F *hcate_toy[3]; 
            hcate_toy[0] = new TH1F("temp_cate0","",numberbins,Xregions[0],Xregions[1]);
            hcate_toy[1] = new TH1F("temp_cate3","",numberbins,Xregions[0],Xregions[1]);
            hcate_toy[2] = new TH1F("temp_cate4","",numberbins,Xregions[0],Xregions[1]);

            for(int i=1;i<=hcate0->GetNbinsX();i++)
                hcate_toy[0]->Fill(hcate0->GetBinCenter(i),
                        (float)rnd->Poisson((double)hcate0->GetBinContent(i)));  
                        //(float)rnd->Gaus((double)hcate0->GetBinContent(i), (double)hcate0->GetBinError(i)));  
            for(int i=1;i<=hcate3->GetNbinsX();i++)
                hcate_toy[1]->Fill(hcate3->GetBinCenter(i),
                        (float)rnd->Poisson((double)hcate3->GetBinContent(i)));  
                        //(float)rnd->Gaus((double)hcate3->GetBinContent(i), (double)hcate3->GetBinError(i)));  
            for(int i=1;i<=hcate4->GetNbinsX();i++)
                hcate_toy[2]->Fill(hcate4->GetBinCenter(i),
                        (float)rnd->Poisson((double)hcate4->GetBinContent(i)));  
                        //(float)rnd->Gaus((double)hcate4->GetBinContent(i), (double)hcate4->GetBinError(i)));  

            RooDataHist dCate0_toy("dCate0_toy","dCate0set with x",x,hcate_toy[0]);
            RooHistPdf Cate0_toy("Cate0_toy","Cate0_toy",x,dCate0_toy,0) ;

            RooDataHist dCate3_toy("dCate3_toy","dCate3set with x",x,hcate_toy[1]);
            RooHistPdf Cate3_toy("Cate3_toy","Cate3_toy",x,dCate3_toy,0) ;

            RooDataHist dCate4_toy("dCate4_toy","dCate4set with x",x,hcate_toy[2]);
            RooHistPdf Cate4_toy("Cate4_toy","Cate4_toy",x,dCate4_toy,0) ;

            RooDataHist data_toy("data_toy","dataset with x",x,hdata);
            // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            RooAddPdf model_toy("model_toy","model_toy",
                    RooArgList(Cate0_toy,Cate3_toy,Cate4_toy,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);

            model_toy.fitTo(data_toy,Minimizer("Minuit2", "minimize"));
            //model_toy.fitTo(data_toy,InitialHesse(kTRUE),Hesse(kTRUE));
            //model_toy.chi2FitTo(data_toy,Extended(kTRUE));
            //model_toy.fitTo(data_toy,Extended(kTRUE));

            hbias[0]->Fill((Norm_fbkg- (1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal())) )/Norm_fbkg);

            hbias[1]->Fill((Norm_fbkg0-fbkg0.getVal())/Norm_fbkg0_error);
            hbias[2]->Fill((Norm_fbkg3-(1-fbkg0.getVal())*(fbkg3.getVal()))/Norm_fbkg3_error);
            hbias[3]->Fill((Norm_fbkg4-(1-fbkg0.getVal())*(1-fbkg3.getVal())*fbkg4.getVal())/Norm_fbkg4_error);

            //hbias[1]->Fill((Norm_fbkg0-fbkg0.getVal())/Norm_fbkg0);
            //hbias[2]->Fill((Norm_fbkg3-(1-fbkg0.getVal())*(fbkg3.getVal()))/Norm_fbkg3);
            //hbias[3]->Fill((Norm_fbkg4-(1-fbkg0.getVal())*(1-fbkg3.getVal())*fbkg4.getVal())/Norm_fbkg4);

            std::cout<<"[Jacky test] fbkg0.getVal() : "<<fbkg0.getVal()<<std::endl; 

            hcate_toy[0]->Delete();
            hcate_toy[1]->Delete();
            hcate_toy[2]->Delete();

            RooPlot* xframe = x.frame() ; 
            //xframe->SetTitle("In signal region");
            xframe->SetTitle("");
            data_toy.plotOn(xframe);
            model.plotOn(xframe,LineColor(kBlue-2));
            model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
            model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
            model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
            model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
            xframe->Draw() ;
	    drawCMSLumi();

            hbias[4]->Fill(xframe->chiSquare("model","data_toy",3));
        }
}

void BiasStudy(TH1F *hdata, TH1F* hcate0, TH1F* hcate3, TH1F* hcate4, TH1F* hsig, TH1F* hbias, TCanvas* TCanvasTemp){
        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;

        RooDataHist data("data","dataset with x",x,hdata);

        RooDataHist dCate0("dCate0","dCate0set with x",x,hcate0);
        RooHistPdf Cate0("Cate0","Cate0",x,dCate0,0) ;

        RooDataHist dCate3("dCate3","dCate3set with x",x,hcate3);
        RooHistPdf Cate3("Cate3","Cate3",x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,hcate4);
        RooHistPdf Cate4("Cate4","Cate4",x,dCate4,0) ;

        RooDataHist dSignal("dSignal","dSignalset with x",x,hsig);
        RooHistPdf Signal("Signal","Signal",x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0,1) ; 
        RooRealVar fbkg4("fbkg4","fraction",0.5,0,1) ; 
        // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
        RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
        model.fitTo(data,Minimizer("Minuit2", "minimize"));

        float Norm_fbkg = 1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal()) ;

        TCanvasTemp->Divide(9,9);
        int nToy = 1000;
        int Ncount_ = 0;
        for(int itoy=0;itoy<nToy;itoy++){

            TH1F * hdata_toy = new TH1F("hdata_toy","",
                    hdata->GetNbinsX(),hdata->GetXaxis()->GetXmin(),hdata->GetXaxis()->GetXmax());
            for(int i=1;i<=hdata->GetNbinsX();i++)
                hdata_toy->Fill(hdata_toy->GetBinCenter(i),(float)rnd->Poisson((double)hdata->GetBinContent(i)));  

            RooDataHist data_toy("data_toy","dataset with x",x,hdata_toy);
            // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            RooAddPdf model_toy("model_toy","model_toy",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
            model.fitTo(data_toy,Minimizer("Minuit2", "minimize"));
            hbias->Fill((Norm_fbkg- (1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal())) )/Norm_fbkg);
            delete hdata_toy;

            //if(!((Norm_fbkg- (1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal())) )/Norm_fbkg<-0.08)) continue;
            if(Ncount_>80) continue;

            std::cout<<"Toy data try : "<<Ncount_<<" / "<<itoy<<std::endl;

            TCanvasTemp->cd(Ncount_+1);
            RooPlot* xframe = x.frame() ; 
            //xframe->SetTitle("In signal region");
            xframe->SetTitle("");
            data_toy.plotOn(xframe);
            model.plotOn(xframe,LineColor(kBlue-2));
            model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
            model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
            model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
            model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
            xframe->Draw() ;
	    drawCMSLumi();


            TLegend *legend_nmFit = new TLegend(0.5077,0.605,0.9,0.87);
            TH1F *htmp[6];
            //string nametemp[6] = {"toy data","Fit model","Cate0(data)","Cate3(data)","Cate4(data)","Signal"};
            string nametemp[6] = {"toy data","Fit model","Cat M(data)","Cat 3(data)","Cat 4(data)","Signal"};
            int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
            int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
            for(int i=0;i<6;i++){
                sprintf(buffer,"%s",nametemp[i].c_str());
                htmp[i] = new TH1F(buffer,"",10,0,10);
                htmp[i]->SetLineColor(Lcolors[i]);
                htmp[i]->SetLineStyle(Lstyle[i]);
                htmp[i]->SetLineWidth(3);

                if(i==0){
                    htmp[i]->SetMarkerStyle(20);
                    htmp[i]->SetLineWidth(1);
                    legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
                }else{
                    legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
                }   
            }   

            legend_nmFit->SetBorderSize(0);
            legend_nmFit->SetFillColor(0);
            legend_nmFit->SetFillStyle(0);
            legend_nmFit->SetNColumns(1);
            legend_nmFit->SetTextSize(0.04);
            legend_nmFit->SetTextSizePixels(25);
            legend_nmFit->Draw();

            // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            TLatex *latexComp = new TLatex(0.52,0.51,"Model = ");
            latexComp->SetTextSize(0.03);
            latexComp->SetTextColor(kBlue-2);
            latexComp->SetNDC();
            latexComp->Draw();

            //sprintf(buffer,"#splitline{%1.1ExCate0 +}{%1.1ExCate3 +}",
            sprintf(buffer,"#splitline{%1.1ExCat M +}{%1.1ExCat 3 +}",
                    fbkg0.getVal(),
                    (1-fbkg0.getVal())*fbkg3.getVal()
                   );
            latexComp = new TLatex(0.63,0.50,buffer);
            latexComp->SetTextSize(0.03);
            latexComp->SetTextColor(kBlue-2);
            latexComp->SetNDC();
            latexComp->Draw();

            sprintf(buffer,"#splitline{%1.1ExCat 4 + }{%1.1ExSig}",
                    (1 - fbkg0.getVal())*(1 - fbkg3.getVal())*fbkg4.getVal(),
                    (1 - fbkg0.getVal())*(1 - fbkg3.getVal())*(1-fbkg4.getVal())
                   );
            latexComp = new TLatex(0.63,0.44,buffer);
            latexComp->SetTextSize(0.03);
            latexComp->SetTextColor(kBlue-2);
            latexComp->SetNDC();
            latexComp->Draw();
            Ncount_++;
        }
}

/*
Templates : 
cate0(data)   :   h_TemplatesExtractedFromData[0] 
cate1/3(data) :   h_TemplatesExtractedFromData[1]
cate4(data)   :   h_TemplatesExtractedFromData[2]
SR(data)      :   h_TemplatesExtractedFromData[3]
SR(MC)        :   h_MCinSignalRegion
SR(S 950)     :   h_SignalTemplateInSR
 */
float DrawFit(int sidx_/*signal index*/,TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[],TH1F *h_MCinSignalRegion, RooWorkspace *wspace,
        bool mode/*0:data, 1:MC*/){

        char TrimChar_[128];
	TH1F *hdata_toy;
	TH1F *hsig_temp = (TH1F*) h_SignalTemplateInSR[sidx_]->Clone();
	TH1F *hbkg_cate0 = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
	TH1F *hbkg_cate4 = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
	sprintf(buffer,"hbkg_cate0_%f",(float)rnd->Uniform(0,1));
	hbkg_cate0->SetName(buffer);
	hbkg_cate0->Scale(1./hbkg_cate0->Integral());
	sprintf(buffer,"hbkg_cate4_%f",(float)rnd->Uniform(0,1));
	hbkg_cate4->SetName(buffer);
	hbkg_cate4->Scale(1./hbkg_cate4->Integral());
	TH1F *htotal_bkg = (TH1F*) hbkg_cate0->Clone(); // will be scaled according to fitting result

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data_obs","data_obs",x,h_TemplatesExtractedFromData[3]);
        RooDataHist MC("MC","MC with x",x,h_MCinSignalRegion);

        RooDataHist dCate0("dCate0","dCate0set with x",x,h_TemplatesExtractedFromData[0]);
        sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

        char buffer13[256];
        sprintf(buffer13,"dCate3set with x");
        sprintf(buffer,"dCate3");
        if(usingCate1insteadofCate3){
            sprintf(buffer,"dCate1");
            sprintf(buffer13,"dCate1set with x");
        }

        //RooDataHist dCate3("dCate3","dCate3set with x",x,h_TemplatesExtractedFromData[1]);
        RooDataHist dCate3(buffer,buffer13,x,h_TemplatesExtractedFromData[1]);
        sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        if(usingCate1insteadofCate3){
            sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        }
        RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,h_TemplatesExtractedFromData[2]);
        sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

        sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
        RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 
        //RooRealVar fbkg4("fbkg4","fraction",0.5,0.00,1) ; 
            // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) 
            // RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
        // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) 
        RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);

        if(mode==0){
            hdata_toy = (TH1F*) h_TemplatesExtractedFromData[3]->Clone();
        }else{
            hdata_toy = (TH1F*) h_MCinSignalRegion->Clone();
        }

	float xmin_ = -0.000000000001;
	float xmax_ = 1.5*hdata_toy->Integral();
	RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
	RooRealVar NSign("NSign","N(Sig)",2.,xmin_,xmax_) ;

        //RooRealVar Nbkgall("Nbkgall","yield",0.5,0.00,hdata_toy->Integral()) ; 
        //RooRealVar NSign("NSign","yield",0.5,0.00,hdata_toy->Integral()) ;
        //RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));


        //model.chi2FitTo(data,Extended(kTRUE));
        if(mode==0){
            //model.fitTo(data,Minimizer("Minuit2", "minimize"));
            model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            // Although we used "extended maximum likelihood fit", somethimes it cannot work well. 
            while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata_toy->Integral() - 1.) > 0.01){
                printf("[In while loop (DrawFit-DATA)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata_toy->Integral(), 
                        Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
                model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            } 
        }else{
            //model.fitTo(MC,Minimizer("Minuit2", "minimize"));
            model.fitTo(MC,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata_toy->Integral() - 1.) > 0.01){
                printf("[In while loop (DrawFit-MC)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata_toy->Integral(), 
                        Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
                model.fitTo(MC,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            } 
        }

        // http://roofit-scanner-y1s-emu.googlecode.com/svn/branches/EPS-HEP_Cuts/background/One_Range/old/Fit_Multi_Exp_Var.C
        RooArgSet obs(x);
        RooArgSet* flparams = (RooArgSet*) model.getParameters(obs)->selectByAttrib("Constant",kFALSE) ; 

        RooPlot* xframe = x.frame() ; 
        RooPlot* xframe_tmp = x.frame() ; 
        //xframe->SetTitle("In signal region");
        xframe->SetTitle("");
        xframe_tmp->SetTitle("");
        float ymax_ = 1.5*h_TemplatesExtractedFromData[3]->GetMaximum();
        if(mode==0){
            data.plotOn(xframe);
            data.plotOn(xframe_tmp);
        }else{
            ymax_ = 1.5*h_MCinSignalRegion->GetMaximum();
            MC.plotOn(xframe);
        }
	if(!IsMoneyPlot){
	    model.plotOn(xframe,LineColor(kBlue-2));
	    model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
	    //model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
	    model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
	    model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
	    //xframe->GetYaxis()->SetRangeUser(0.,xframe->GetYaxis()->GetXmax());
	    xframe->GetYaxis()->SetRangeUser(0.0,ymax_);
	    if(UsingLogyScaleInMass)
		xframe->GetYaxis()->SetRangeUser(0.001,pow(10,2.0*log10(ymax_)));
	    xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
	    xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
	    xframe->Draw() ;
	    drawCMSLumi();
	}else{
	    if(!IsMoneyPlotReduced){
		hsig_temp->SetFillColor(kBlue-4);
		hsig_temp->SetLineColor(kBlue-4);
		//hsig_temp->SetFillStyle(3001);
		hsig_temp->SetFillStyle(3004);
		model.plotOn(xframe,LineColor(kBlue-2));
		model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
		model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
		model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
		xframe->GetYaxis()->SetRangeUser(0.0,ymax_);
		if(UsingLogyScaleInMass)
		    xframe->GetYaxis()->SetRangeUser(0.001,pow(10,2.5*log10(ymax_)));
		xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
		xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
		xframe->Draw() ;
		hsig_temp->Draw("hist,same");
		drawCMSLumi();
	    }else{
		hsig_temp->SetFillColor(kBlue-4);
		hsig_temp->SetLineColor(kBlue-4);
		//hsig_temp->SetFillStyle(3001);
		hsig_temp->SetFillStyle(3004);
		//model.plotOn(xframe,LineColor(kBlue-2));
		model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
		model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
		model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
		xframe_tmp->GetYaxis()->SetRangeUser(0.0,ymax_);
		if(UsingLogyScaleInMass)
		    xframe_tmp->GetYaxis()->SetRangeUser(0.001,pow(10,2.5*log10(ymax_)));
		xframe_tmp->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
		xframe_tmp->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
		xframe_tmp->Draw() ;
		hsig_temp->GetYaxis()->SetRangeUser(0.001,pow(10,2.5*log10(ymax_)));
		hsig_temp->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
		hsig_temp->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
		//hsig_temp->Draw("hist");
		hsig_temp->Draw("hist,same");
		//model.plotOn(xframe,LineColor(kBlue-2));
		//xframe->Draw("same") ;
		drawCMSLumi();
	    }
	}


        // polish the plot here
	float fbkg = Nbkgall.getVal()/(Nbkgall.getVal()+NSign.getVal());
	if(!IsMoneyPlot){
	    float offsetLegend = 0.0;
	    if(UsingLogyScaleInMass)
		offsetLegend=0.06;
	    TLegend *legend_nmFit = new TLegend(0.5077,0.605+offsetLegend,0.9,0.87+offsetLegend);
	    TH1F *htmp[6];
	    TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
	    sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
	    string nametemp[6] = {"data","Fit model","Cat M (data)","Cat 3 (data)","Cat 4 (data)",buffer};
	    if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";
	    if(mode==1){
		nametemp[0] = "MC bkg.";
	    }
	    int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
	    int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
	    for(int i=0;i<6;i++){
		if(i==3) continue;
		sprintf(buffer,"%s",nametemp[i].c_str());
		htmp[i] = new TH1F(buffer,"",10,0,10);
		htmp[i]->SetLineColor(Lcolors[i]);
		htmp[i]->SetLineStyle(Lstyle[i]);
		htmp[i]->SetLineWidth(3);

		if(i==0){
		    htmp[i]->SetMarkerStyle(20);
		    htmp[i]->SetLineWidth(1);
		    legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
		}else{
		    legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
		}   
	    }   

	    legend_nmFit->SetBorderSize(0);
	    legend_nmFit->SetFillColor(0);
	    legend_nmFit->SetFillStyle(0);
	    legend_nmFit->SetNColumns(1);
	    legend_nmFit->SetTextSize(0.04);
	    legend_nmFit->SetTextSizePixels(25);
	    legend_nmFit->Draw();

	    // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
	    // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) + (1- fbkg) * Signal
	    //         Nbkg                                                                   NSign
	    float TextOffsetX = 0.01;
	    float TextOffsetY = -0.01;
	    TLatex *latexComp;
	    if(UsingLogyScaleInMass)
		latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
	    else
		latexComp = new TLatex(0.52,0.51,"Fit model = ");
	    latexComp->SetTextSize(0.03);
	    latexComp->SetTextColor(kBlue-2);
	    latexComp->SetNDC();
	    latexComp->Draw();


	    sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4 +}",
		    fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
	    if(usingCate1insteadofCate3)
		sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
			fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
	    latexComp = new TLatex(0.63,0.50,buffer);
	    if(UsingLogyScaleInMass)
		latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
	    latexComp->SetTextSize(0.03);
	    latexComp->SetTextColor(kBlue-2);
	    latexComp->SetNDC();
	    latexComp->Draw();

	    sprintf(buffer,"%1.2fxSig",
		    //fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
		    1.-fbkg
		   );
	    latexComp = new TLatex(0.63,0.46,buffer);
	    if(UsingLogyScaleInMass)
		latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
	    latexComp->SetTextSize(0.03);
	    latexComp->SetTextColor(kBlue-2);
	    latexComp->SetNDC();
	    latexComp->Draw();

	    sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
		    Nbkgall.getVal(),
		    Nbkgall.getError(),
		    NSign.getVal(),
		    NSign.getError()//, hdata_toy->Integral()
		   );
	    latexComp = new TLatex(0.52,0.36,buffer);
	    if(UsingLogyScaleInMass)
		latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
	    latexComp->SetTextSize(0.03);
	    latexComp->SetTextColor(kBlue-2);
	    latexComp->SetNDC();
	    latexComp->Draw();

	    sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
	    latexComp = new TLatex(0.52,0.56,buffer);
	    if(UsingLogyScaleInMass)
		latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
	    latexComp->SetTextSize(0.03);
	    latexComp->SetTextColor(kBlack);
	    latexComp->SetNDC();
	    latexComp->Draw();
	}else{
	    if(!IsMoneyPlotReduced){
		float offsetLegend = 0.0;
		if(UsingLogyScaleInMass)
		    offsetLegend=0.06;
		float OffX = 0.02;
		float ft_ = 0.005;
		// data+fit
		TLegend *legend_nmFit = new TLegend(0.5077-OffX,0.605-0.05+0.225+offsetLegend,0.9-OffX,0.87+offsetLegend);
		// Total bkg
		TLegend *legend_nmFit5 = new TLegend(0.5077,0.605-0.05+0.18+ft_+offsetLegend,0.9,0.78-ft_+offsetLegend); 
		// cateM+4
		TLegend *legend_nmFit2 = new TLegend(0.5077+4*OffX,0.605-0.05+0.09+offsetLegend,0.9+4*OffX,0.735+offsetLegend);

		// t*
		TLegend *legend_nmFit3 = new TLegend(0.5077,0.605-0.05+0.045+offsetLegend,0.9,0.645+offsetLegend);      
		// t* MC
		TLegend *legend_nmFit4 = new TLegend(0.5077-OffX,0.605-0.05+ft_+offsetLegend,0.9-OffX,0.6-ft_+offsetLegend);    

		const int Nhtmp = 8;
		TH1F *htmp[Nhtmp];
		TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
		sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
		sprintf(buffer2,"t* (m_{t*}= %1.0d GeV/c^{2}, MC)",atoi(TrimChar_));
		string nametemp[Nhtmp] = {"data","Fit model","Cat M (data)","Cat 3 (data)","Cat 4 (data)",buffer, buffer2,
		    "Total Bkg & unc.(sys.)"};
		if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";
		if(mode==1){
		    nametemp[0] = "MC bkg.";
		}
		int Lcolors[Nhtmp] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2, kBlue-4, 800-7};
		int Lstyle[Nhtmp] = {1,1,kDashed,kDashed,kDashed,kDashed, 1, 1};
		for(int i=0;i<Nhtmp;i++){
		    if(i==3) continue;
		    sprintf(buffer,"%s",nametemp[i].c_str());
		    htmp[i] = new TH1F(buffer,"",10,0,10);
		    htmp[i]->SetLineColor(Lcolors[i]);
		    htmp[i]->SetLineStyle(Lstyle[i]);
		    htmp[i]->SetLineWidth(3);

		    if(i==0){
			htmp[i]->SetMarkerStyle(20);
			htmp[i]->SetLineWidth(1);
			legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
		    }else if(i==6){
			htmp[i]->SetFillColor(kBlue-4);
			htmp[i]->SetLineColor(kBlue-4);
			//htmp[i]->SetFillStyle(3001);
			htmp[i]->SetFillStyle(3004);
			legend_nmFit4->AddEntry(htmp[i],nametemp[i].c_str(),"f");
		    }else if(i==7){
			htmp[i]->SetMarkerSize(0.001);
			htmp[i]->SetFillStyle(3005);
			htmp[i]->SetFillColor(800-7);
			legend_nmFit5->AddEntry(htmp[i],nametemp[i].c_str(),"lf");
		    }else{
			if(i==1)
			    legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
			if(i==2||i==4)
			    legend_nmFit2->AddEntry(htmp[i],nametemp[i].c_str(),"l");
			if(i==5)
			    legend_nmFit3->AddEntry(htmp[i],nametemp[i].c_str(),"l");
		    }
		}
		float TextSize_ = 0.035;
		TLatex *latexComp;

		legend_nmFit->SetBorderSize(0);
		legend_nmFit->SetFillColor(0);
		legend_nmFit->SetFillStyle(0);
		legend_nmFit->SetNColumns(1);
		legend_nmFit->SetTextSize(TextSize_);
		legend_nmFit->SetTextSizePixels(25);
		legend_nmFit->Draw();

		legend_nmFit2->SetBorderSize(0);
		legend_nmFit2->SetFillColor(0);
		legend_nmFit2->SetFillStyle(0);
		legend_nmFit2->SetNColumns(1);
		legend_nmFit2->SetTextSize(TextSize_);
		legend_nmFit2->SetTextSizePixels(25);
		legend_nmFit2->Draw();

		legend_nmFit3->SetBorderSize(0);
		legend_nmFit3->SetFillColor(0);
		legend_nmFit3->SetFillStyle(0);
		legend_nmFit3->SetNColumns(1);
		legend_nmFit3->SetTextSize(TextSize_);
		legend_nmFit3->SetTextSizePixels(25);
		legend_nmFit3->Draw();

		legend_nmFit4->SetBorderSize(0);
		legend_nmFit4->SetFillColor(0);
		legend_nmFit4->SetFillStyle(0);
		legend_nmFit4->SetNColumns(1);
		legend_nmFit4->SetTextSize(TextSize_);
		legend_nmFit4->SetTextSizePixels(25);
		legend_nmFit4->Draw();

		legend_nmFit5->SetBorderSize(0);
		legend_nmFit5->SetFillColor(0);
		legend_nmFit5->SetFillStyle(0);
		legend_nmFit5->SetNColumns(1);
		legend_nmFit5->SetTextSize(TextSize_);
		legend_nmFit5->SetTextSizePixels(25);
		legend_nmFit5->Draw();

		htotal_bkg->Scale(Nbkgall.getVal()*fbkg0.getVal());
		hbkg_cate4->Scale(Nbkgall.getVal()*(1-fbkg0.getVal()));
		htotal_bkg->Add(hbkg_cate4);
		htotal_bkg->SetLineColor(800-7);

		sprintf(buffer,"unc_%f",(float)rnd->Uniform(0,1));
		TH1F *unc = new TH1F(buffer,"",htotal_bkg->GetNbinsX(),
			htotal_bkg->GetXaxis()->GetXmin(), htotal_bkg->GetXaxis()->GetXmax());
		for(int ibin = 1; ibin <= htotal_bkg->GetNbinsX(); ibin++){
		    if(AppliedSysOnShape){
			unc->SetBinError(ibin, htotal_bkg->GetBinContent(ibin) * htotalunc950->GetBinContent(ibin));
		    }else{
			unc->SetBinError(ibin, htotal_bkg->GetBinContent(ibin) );
		    }
		    unc->SetBinContent(ibin, htotal_bkg->GetBinContent(ibin) );
		}
		htotal_bkg->Draw("hist,same");

		unc->SetMarkerSize(0.001);
		unc->SetFillStyle(3005);
		unc->SetFillColor(800-7);
		unc->Draw("E2,same");

		hbkg_cate0->SetLineColor(Lcolors[2]);
		hbkg_cate0->SetLineStyle(Lstyle[2]);
		hbkg_cate0->SetLineWidth(3);
		hbkg_cate0->Scale(Nbkgall.getVal()*fbkg0.getVal());
		hbkg_cate0->Draw("hist,same");

		hbkg_cate4->SetLineColor(Lcolors[4]);
		hbkg_cate4->SetLineStyle(Lstyle[4]);
		hbkg_cate4->SetLineWidth(3);
		hbkg_cate4->Scale(1./hbkg_cate4->Integral());
		hbkg_cate4->Scale(Nbkgall.getVal()*(1-fbkg0.getVal()));
		hbkg_cate4->Draw("hist,same");

		// model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
		// model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) + (1- fbkg) * Signal
		//         Nbkg                                                                   NSign
		float TextOffsetX = 0.01;
		float TextOffsetY = -0.01;
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
		else
		    latexComp = new TLatex(0.52,0.51,"Fit model = ");
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlue-2);
		latexComp->SetNDC();
		latexComp->Draw();


		sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4 +}",
			fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
		if(usingCate1insteadofCate3)
		    sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
			    fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
		latexComp = new TLatex(0.63,0.50,buffer);
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlue-2);
		latexComp->SetNDC();
		latexComp->Draw();

		sprintf(buffer,"%1.2fxSig",
			//fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
			1.-fbkg
		       );
		latexComp = new TLatex(0.63,0.46,buffer);
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlue-2);
		latexComp->SetNDC();
		latexComp->Draw();

		sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
			Nbkgall.getVal(),
			Nbkgall.getError(),
			NSign.getVal(),
			NSign.getError()//, hdata_toy->Integral()
		       );
		latexComp = new TLatex(0.52,0.36,buffer);
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlue-2);
		latexComp->SetNDC();
		latexComp->Draw();

		sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
		latexComp = new TLatex(0.52,0.56,buffer);
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlack);
		latexComp->SetNDC();
		latexComp->Draw();
	    }else{
		float offsetLegend = 0.0;
		if(UsingLogyScaleInMass)
		    offsetLegend=0.06;
		float OffX = 0.02;
		float ft_ = 0.005;
		// data+fit
		TLegend *legend_nmFit = new TLegend(0.5077-OffX,0.605-0.05+0.225+offsetLegend,0.9-OffX,0.87+offsetLegend);
		// Total bkg
		TLegend *legend_nmFit5 = new TLegend(0.5077-OffX,0.605-0.05+0.18+ft_+offsetLegend,0.9-OffX,0.78-ft_+offsetLegend); 
		// t* MC
		TLegend *legend_nmFit4 = new TLegend(0.5077-OffX,0.605-0.05+0.09+offsetLegend,0.9-OffX,0.735+offsetLegend);    

		const int Nhtmp = 3;
		TH1F *htmp[Nhtmp];
		TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
		sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
		sprintf(buffer2,"t* (m_{t*}= %1.0d GeV/c^{2}, MC)",atoi(TrimChar_));
		string nametemp[Nhtmp] = {"data", buffer2, "Total Bkg & unc.(sys.)"};
		int Lcolors[Nhtmp] = {1, kBlue-4, 800-7};
		int Lstyle[Nhtmp] = {1, 1, 1};
		for(int i=0;i<Nhtmp;i++){
		    sprintf(buffer,"%s",nametemp[i].c_str());
		    htmp[i] = new TH1F(buffer,"",10,0,10);
		    htmp[i]->SetLineColor(Lcolors[i]);
		    htmp[i]->SetLineStyle(Lstyle[i]);
		    htmp[i]->SetLineWidth(3);

		    if(i==0){
			htmp[i]->SetMarkerStyle(20);
			htmp[i]->SetLineWidth(1);
			legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
		    }else if(i==1){
			htmp[i]->SetFillColor(kBlue-4);
			htmp[i]->SetLineColor(kBlue-4);
			//htmp[i]->SetFillStyle(3001);
			htmp[i]->SetFillStyle(3004);
			legend_nmFit4->AddEntry(htmp[i],nametemp[i].c_str(),"f");
		    }else if(i==2){
			htmp[i]->SetMarkerSize(0.001);
			htmp[i]->SetFillStyle(3005);
			htmp[i]->SetFillColor(800-7);
			legend_nmFit5->AddEntry(htmp[i],nametemp[i].c_str(),"lf");
		    }
		}
		float TextSize_ = 0.035;
		TLatex *latexComp;

		legend_nmFit->SetBorderSize(0);
		legend_nmFit->SetFillColor(0);
		legend_nmFit->SetFillStyle(0);
		legend_nmFit->SetNColumns(1);
		legend_nmFit->SetTextSize(TextSize_);
		legend_nmFit->SetTextSizePixels(25);
		legend_nmFit->Draw();

		legend_nmFit4->SetBorderSize(0);
		legend_nmFit4->SetFillColor(0);
		legend_nmFit4->SetFillStyle(0);
		legend_nmFit4->SetNColumns(1);
		legend_nmFit4->SetTextSize(TextSize_);
		legend_nmFit4->SetTextSizePixels(25);
		legend_nmFit4->Draw();

		legend_nmFit5->SetBorderSize(0);
		legend_nmFit5->SetFillColor(0);
		legend_nmFit5->SetFillStyle(0);
		legend_nmFit5->SetNColumns(1);
		legend_nmFit5->SetTextSize(TextSize_);
		legend_nmFit5->SetTextSizePixels(25);
		legend_nmFit5->Draw();

		htotal_bkg->Scale(Nbkgall.getVal()*fbkg0.getVal());
		hbkg_cate4->Scale(Nbkgall.getVal()*(1-fbkg0.getVal()));
		htotal_bkg->Add(hbkg_cate4);
		htotal_bkg->SetLineColor(800-7);

		sprintf(buffer,"unc_%f",(float)rnd->Uniform(0,1));
		TH1F *unc = new TH1F(buffer,"",htotal_bkg->GetNbinsX(),
			htotal_bkg->GetXaxis()->GetXmin(), htotal_bkg->GetXaxis()->GetXmax());
		for(int ibin = 1; ibin <= htotal_bkg->GetNbinsX(); ibin++){
		    if(AppliedSysOnShape){
			unc->SetBinError(ibin, htotal_bkg->GetBinContent(ibin) * htotalunc950->GetBinContent(ibin));
		    }else{
			unc->SetBinError(ibin, htotal_bkg->GetBinContent(ibin) );
		    }
		    unc->SetBinContent(ibin, htotal_bkg->GetBinContent(ibin) );
		}
		htotal_bkg->Draw("hist,same");

		unc->SetMarkerSize(0.001);
		unc->SetFillStyle(3005);
		unc->SetFillColor(800-7);
		unc->Draw("E2,same");

		//hbkg_cate0->SetLineColor(Lcolors[2]);
		//hbkg_cate0->SetLineStyle(Lstyle[2]);
		//hbkg_cate0->SetLineWidth(3);
		//hbkg_cate0->Scale(Nbkgall.getVal()*fbkg0.getVal());
		//hbkg_cate0->Draw("hist,same");

		//hbkg_cate4->SetLineColor(Lcolors[4]);
		//hbkg_cate4->SetLineStyle(Lstyle[4]);
		//hbkg_cate4->SetLineWidth(3);
		//hbkg_cate4->Scale(1./hbkg_cate4->Integral());
		//hbkg_cate4->Scale(Nbkgall.getVal()*(1-fbkg0.getVal()));
		//hbkg_cate4->Draw("hist,same");

		// model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
		// model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) + (1- fbkg) * Signal
		//         Nbkg                                                                   NSign
		float TextOffsetX = 0.01;
		float TextOffsetY = -0.01;
		//if(UsingLogyScaleInMass)
		//    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
		//else
		//    latexComp = new TLatex(0.52,0.51,"Fit model = ");
		//latexComp->SetTextSize(0.03);
		//latexComp->SetTextColor(kBlue-2);
		//latexComp->SetNDC();
		//latexComp->Draw();


		//sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4 +}",
		//	fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
		//if(usingCate1insteadofCate3)
		//    sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
		//	    fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
		//latexComp = new TLatex(0.63,0.50,buffer);
		//if(UsingLogyScaleInMass)
		//    latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
		//latexComp->SetTextSize(0.03);
		//latexComp->SetTextColor(kBlue-2);
		//latexComp->SetNDC();
		//latexComp->Draw();

		//sprintf(buffer,"%1.2fxSig",
		//	//fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
		//	1.-fbkg
		//       );
		//latexComp = new TLatex(0.63,0.46,buffer);
		//if(UsingLogyScaleInMass)
		//    latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
		//latexComp->SetTextSize(0.03);
		//latexComp->SetTextColor(kBlue-2);
		//latexComp->SetNDC();
		//latexComp->Draw();

		//sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
		//	Nbkgall.getVal(),
		//	Nbkgall.getError(),
		//	NSign.getVal(),
		//	NSign.getError()//, hdata_toy->Integral()
		//       );
		//latexComp = new TLatex(0.52,0.36,buffer);
		//if(UsingLogyScaleInMass)
		//    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
		//latexComp->SetTextSize(0.03);
		//latexComp->SetTextColor(kBlue-2);
		//latexComp->SetNDC();
		//latexComp->Draw();

		sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
		//sprintf(buffer,"#chi^{2} / ndf = %1.2f/ %1.0f",(float)xframe->chiSquare(),(float)flparams->getSize());
		latexComp = new TLatex(0.52,0.56,buffer);
		if(UsingLogyScaleInMass)
		    latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
		latexComp->SetTextSize(0.03);
		latexComp->SetTextColor(kBlack);
		latexComp->SetNDC();
		latexComp->Draw();
	    }
	}


        /*
        RooRealVar Nbkg0("Nbkg0","yield",0.5,0.00,hdata_toy->Integral()) ; 
        RooRealVar Nbkg3("Nbkg3","yield",0.5,0.00,hdata_toy->Integral()) ; 
        RooRealVar Nbkg4("Nbkg4","yield",0.5,0.00,hdata_toy->Integral()) ; 
        RooRealVar NSign("NSign","yield",0.5,0.00,hdata_toy->Integral()) ; 
        RooAddPdf model_toy2("model_toy2","model_toy2",
                RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(Nbkg0,Nbkg3,Nbkg4,NSign));     

        if(mode==0){
            //model_toy2.fitTo(data,SumW2Error(kTRUE),Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE));
            model_toy2.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),Hesse(kTRUE));
        }else{
            //model_toy2.fitTo(MC,SumW2Error(kTRUE),Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE));
            model_toy2.fitTo(MC,SumW2Error(kFALSE),Extended(kTRUE),Hesse(kTRUE));
        }
        std::cout<<"[Fit mode "<< mode<<" for "<<SAMPLE[Tgamma500GeV+sidx_].tag<<"] : ( "<<
            Nbkg0.getVal()<<" +/- "<<Nbkg0.getError()<<" , "<<
            Nbkg3.getVal()<<" +/- "<<Nbkg3.getError()<<" , "<<
            Nbkg4.getVal()<<" +/- "<<Nbkg4.getError()<<" , "
            <<" ) : "<< NSign.getVal() <<" +/- "<< NSign.getError() <<" ( "<<NSign.getErrorLo()<< " + "<<
            NSign.getErrorHi()<<")"<<std::endl;
        */


        // Store template for limit calculation
        if(mode==0&&write_template){
            //float frac_total_BG = 1. - (1-fbkg0.getVal())*(1-fbkg3.getVal())*(1- fbkg4.getVal());
            float frac_total_BG = fbkg;
            sprintf(buffer,"fbkg0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
            //RooRealVar fbkg0_temp(buffer,"fraction",fbkg0.getVal()/frac_total_BG) ; 
            RooRealVar fbkg0_temp(buffer,"fraction",fbkg0.getVal()) ; 
            sprintf(buffer,"fbkg3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
            RooRealVar fbkg3_temp(buffer,"fraction",fbkg3.getVal()) ; 
                    //(((1-fbkg0.getVal())*fbkg3.getVal())/frac_total_BG)/(1. - fbkg0.getVal()/frac_total_BG ) ) ; 
            // fbkg4_temp = 1 - (fbkg0_temp+fbkg3_temp), when we use RooAddPdf 
            sprintf(buffer,"pdf_bkg_%s",SAMPLE[Tgamma500GeV+sidx_].tag);
            //old: pdf_bkg = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 ))
            //RooAddPdf pdf_bkg(buffer,buffer,RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0_temp,fbkg3_temp));

            // pdf_bkg = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3) * Cate4 )
            //RooAddPdf pdf_bkg(buffer,buffer,RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0_temp,fbkg3_temp),kTRUE);
            RooAddPdf pdf_bkg(buffer,buffer,RooArgList(Cate0,Cate4),RooArgList(fbkg0_temp),kTRUE);
            RooArgSet CoefNorm = model.getCoefNormalization() ;
            pdf_bkg.fixCoefNormalization(CoefNorm);
            wspace->import(Signal);
            wspace->import(pdf_bkg);

            TH1F *signal_dummy = (TH1F*) h_SignalTemplateInSR[sidx_]->Clone();
            sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);

            if(IsForThetaOrHiggsCombined){
                TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
                sprintf(buffer,"%s__signal%1.0i", SAMPLE[Tgamma500GeV+sidx_].tag,atoi(TrimChar_));
            }
            signal_dummy->SetName(buffer);
            if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV)){
                signal_dummy->Write();
		for(int iunc=0;iunc<RunStatusSize;iunc++){
		    if(iunc==Normal) continue;
		    if(iunc==UncXsecPlus) continue;
		    if(iunc==UncXsecMinus) continue;
		    if(iunc==UncQsquare) continue;
		    if(iunc==UncTopPtPlus) continue;
		    if(iunc==UncTopPtMinus) continue;
		    if(iunc>UncLepIDMinus) continue;
		    TH1F *signal_dummy_unc = hshape(atoi(TrimChar_),iunc);
		    sprintf(buffer,"%s__signal%1.0i__%s", SAMPLE[Tgamma500GeV+sidx_].tag,atoi(TrimChar_),
			    RunStatusNamesOnlyForTheta[iunc].c_str());
		    signal_dummy_unc->SetName(buffer);
		    signal_dummy_unc->Write();
		    delete signal_dummy_unc;
		}
	    }
            delete signal_dummy;

            //old : model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
            //         Nbkg                                                                   NSign
            float total_data = h_TemplatesExtractedFromData[3]->Integral();
            TH1F *cate0_dummy = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
            //cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fbkg0.getVal());
            cate0_dummy->Scale(1./cate0_dummy->Integral()*Nbkgall.getVal()*fbkg0.getVal());

            TH1F *cate3_dummy = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
            //cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data* (1. - fbkg0.getVal())* fbkg3.getVal());
            cate3_dummy->Scale(1./cate3_dummy->Integral()*Nbkgall.getVal()* (1. - fbkg0.getVal())* fbkg3.getVal());

            TH1F *cate4_dummy = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
            //cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data* (1.-fbkg0.getVal())*(1.-fbkg3.getVal())*fbkg4.getVal());
            cate4_dummy->Scale(1./cate4_dummy->Integral()*Nbkgall.getVal()* (1. - fbkg0.getVal()));

            sprintf(buffer,"pdf_bkg_%s",SAMPLE[Tgamma500GeV+sidx_].tag);
            if(IsForThetaOrHiggsCombined)
                sprintf(buffer,"%s__pdf_bkg",SAMPLE[Tgamma500GeV+sidx_].tag);
            TH1F *htotal_bg = (TH1F*) cate0_dummy->Clone();
            htotal_bg->SetName(buffer);
            //htotal_bg->Add(cate3_dummy);
            htotal_bg->Add(cate4_dummy);
	    if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV)){
		if(!IsIndividualBkgShape){
		    htotal_bg->Write();
		}else{
		    sprintf(buffer,"%s__pdf_bkgM",
			    SAMPLE[Tgamma500GeV+sidx_].tag);
		    cate0_dummy->SetName(buffer);
		    sprintf(buffer,"%s__pdf_bkg4",
			    SAMPLE[Tgamma500GeV+sidx_].tag);
		    cate4_dummy->SetName(buffer);
		    cate0_dummy->Write();
		    cate4_dummy->Write();
		}
	    }
            delete htotal_bg;

            fractions_[sidx_][0] = fbkg0.getVal();
            fractions_[sidx_][1] = fbkg3.getVal();
            //fractions_[sidx_][2] = fbkg4.getVal();


            if(sidx_==0){
                wspace->import(data);
                TH1F *hdata = (TH1F*) h_TemplatesExtractedFromData[3]->Clone();
                hdata->SetName("data_obs");
                sprintf(buffer,"%s__DATA",SAMPLE[Tgamma500GeV+sidx_].tag);
                if(IsForThetaOrHiggsCombined)
                    hdata->SetName(buffer);
                if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV))
                    hdata->Write();
                delete hdata;
            }else if(IsForThetaOrHiggsCombined){
                TH1F *hdata = (TH1F*) h_TemplatesExtractedFromData[3]->Clone();
                sprintf(buffer,"%s__DATA",SAMPLE[Tgamma500GeV+sidx_].tag);
                hdata->SetName(buffer);
                if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV))
                    hdata->Write();
                delete hdata;
            }
        }

        //return (float)(1. - (1-fbkg0.getVal())*(1.-fbkg3.getVal())*(1.- fbkg4.getVal()));
        return fbkg;
}

void DrawFitUsingMCShape(int sidx_/*signal index*/,TH1F *h_SignalTemplateInSR[],
        TH1F *h_MCinSignalRegion, bool mode/*0:data, 1:MC*/){

        char TrimChar_[128];
        TH1F *hdata_toy;
	TH1F *hsig_temp = (TH1F*) h_SignalTemplateInSR[sidx_]->Clone();
	TH1F *hbkg_cate0 = (TH1F*) h_MergeTemplatesWOSysUnc[0]->Clone();
	TH1F *hbkg_cate1 = (TH1F*) h_MergeTemplatesWOSysUnc[1]->Clone();
	TH1F *hbkg_cate2 = (TH1F*) h_MergeTemplatesWOSysUnc[2]->Clone();
	TH1F *hbkg_cate3 = (TH1F*) h_MergeTemplatesWOSysUnc[3]->Clone();
	TH1F *hbkg_cate4 = (TH1F*) h_MergeTemplatesWOSysUnc[4]->Clone();
	sprintf(buffer,"hbkg_cate0_%f",(float)rnd->Uniform(0,1));
	hbkg_cate0->SetName(buffer);
	hbkg_cate0->Scale(1./hbkg_cate0->Integral());

	sprintf(buffer,"hbkg_cate1_%f",(float)rnd->Uniform(0,1));
	hbkg_cate1->SetName(buffer);
	hbkg_cate1->Scale(1./hbkg_cate1->Integral());

	sprintf(buffer,"hbkg_cate2_%f",(float)rnd->Uniform(0,1));
	hbkg_cate2->SetName(buffer);
	hbkg_cate2->Scale(1./hbkg_cate2->Integral());

	sprintf(buffer,"hbkg_cate3_%f",(float)rnd->Uniform(0,1));
	hbkg_cate3->SetName(buffer);
	hbkg_cate3->Scale(1./hbkg_cate3->Integral());

	sprintf(buffer,"hbkg_cate4_%f",(float)rnd->Uniform(0,1));
	hbkg_cate4->SetName(buffer);
	hbkg_cate4->Scale(1./hbkg_cate4->Integral());

	TH1F *htotal_bkg = (TH1F*) hbkg_cate0->Clone();	// will be scaled according to fitting result

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data_obs","data_obs",x,h_TemplatesExtractedFromData[3]);
        RooDataHist MC("MC","MC with x",x,h_MCinSignalRegion);

        RooDataHist dCate0("dCate0","dCate0set with x",x,h_MergeTemplatesWOSysUnc[0]);
        sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

        RooDataHist dCate1("dCate1","dCate1set with x",x,h_MergeTemplatesWOSysUnc[1]);
        sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate1(buffer,buffer,x,dCate1,0) ;

        RooDataHist dCate2("dCate2","dCate2set with x",x,h_MergeTemplatesWOSysUnc[2]);
        sprintf(buffer,"Cate2_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate2(buffer,buffer,x,dCate2,0) ;

        RooDataHist dCate3("dCate3","dCate3set with x",x,h_MergeTemplatesWOSysUnc[3]);
        sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,h_MergeTemplatesWOSysUnc[4]);
        sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

        sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
        RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg1("fbkg1","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg2("fbkg2","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 
        //RooRealVar fbkg4("fbkg4","fraction",0.5,0.00,1) ; 
            // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) 
            // RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
        // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) 
        //RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);
        RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate1,Cate2,Cate3,Cate4),RooArgList(fbkg0,fbkg1,fbkg2,fbkg3),kTRUE);

        if(mode==0){
            hdata_toy = (TH1F*) h_TemplatesExtractedFromData[3]->Clone();
        }else{
            hdata_toy = (TH1F*) h_MCinSignalRegion->Clone();
        }
	float xmin_ = -0.000000000001;
	float xmax_ = 1.5*hdata_toy->Integral();
	RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
	RooRealVar NSign("NSign","N(Sig)",2.,xmin_,xmax_) ;

	//RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
	/*
	   model = fbkg ( 
	   fbkg0 * Cate0 + (1-fbkg0) (
	   fbkg1* Cate1 + (1-fbkg1)( 
	   fbkg2* Cate2 + (1-fbkg2)*( 
	   fbkg3* Cate3 + (1-fbkg3)* Cate4)
	   )
	   ) 
	   ) + (1- fbkg) * Signal
	 */
        RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));


        //model.chi2FitTo(data,Extended(kTRUE));
        if(mode==0){
            //model.fitTo(data,Minimizer("Minuit2", "minimize"));
            model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            // Although we used "extended maximum likelihood fit", somethimes it cannot work well. 
            while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata_toy->Integral() - 1.) > 0.01){
                printf("[In while loop (DrawFit-DATA)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata_toy->Integral(), 
                        Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
                model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            } 
        }else{
            //model.fitTo(MC,Minimizer("Minuit2", "minimize"));
            model.fitTo(MC,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata_toy->Integral() - 1.) > 0.01){
                printf("[In while loop (DrawFit-MC)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata_toy->Integral(), 
                        Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
                model.fitTo(MC,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
            } 
        }

        // http://roofit-scanner-y1s-emu.googlecode.com/svn/branches/EPS-HEP_Cuts/background/One_Range/old/Fit_Multi_Exp_Var.C
        RooArgSet obs(x);
        RooArgSet* flparams = (RooArgSet*) model.getParameters(obs)->selectByAttrib("Constant",kFALSE) ; 

        RooPlot* xframe = x.frame() ; 
        xframe->SetTitle("");
        float ymax_ = 1.5*h_TemplatesExtractedFromData[3]->GetMaximum();
        if(mode==0){
            data.plotOn(xframe);
        }else{
            ymax_ = 1.5*h_MCinSignalRegion->GetMaximum();
            MC.plotOn(xframe);
        }
	hsig_temp->SetFillColor(kBlue-4);
	hsig_temp->SetLineColor(kBlue-4);
	//hsig_temp->SetFillStyle(3001);
	hsig_temp->SetFillStyle(3004);
        model.plotOn(xframe,LineColor(kBlue-2));
	// colors[idxCate_]
        model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(colors[0]));
        model.plotOn(xframe,Components(Cate1),LineStyle(kDashed),LineColor(colors[1]));
        model.plotOn(xframe,Components(Cate2),LineStyle(kDashed),LineColor(colors[2]));
        model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(colors[3]));
        model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(colors[4]));
        model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
        xframe->GetYaxis()->SetRangeUser(0.0,ymax_);
        if(UsingLogyScaleInMass)
            xframe->GetYaxis()->SetRangeUser(0.00001,pow(10,2.5*log10(ymax_)));
        xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
        xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
        xframe->Draw() ;
	//hsig_temp->Draw("hist,same");
        //xframe->Draw("same") ;
	drawCMSLumi();


        // polish the plot here
        float offsetLegend = 0.0;
        if(UsingLogyScaleInMass)
            offsetLegend=0.06;
	float OffX = 0.02;
	float ft_ = 0.005;
        TLegend *legend_nmFit = new TLegend(0.5077-OffX,0.605+offsetLegend,0.9-OffX,0.87+offsetLegend);// data+fit

	const int Nhtmp = 8;
        TH1F *htmp[Nhtmp];
        TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
        sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
        string nametemp[Nhtmp] = {"MC Bkg.","Fit model"
	    ,"Cat 0(MC)"
	    ,"Cat 1(MC)"
	    ,"Cat 2(MC)"
	    ,"Cat 3(MC)"
	    ,"Cat 4(MC)"
	    ,buffer
	    };
        int Lcolors[Nhtmp] = {1,kBlue-2,colors[0],colors[1],colors[2],colors[3],colors[4], 432+2};
        int Lstyle[Nhtmp] = {1,1,kDashed,kDashed,kDashed,kDashed, kDashed, kDashed};
        for(int i=0;i<Nhtmp;i++){
            //if(i==3) continue;
            sprintf(buffer,"%s",nametemp[i].c_str());
            htmp[i] = new TH1F(buffer,"",10,0,10);
            htmp[i]->SetLineColor(Lcolors[i]);
            htmp[i]->SetLineStyle(Lstyle[i]);
            htmp[i]->SetLineWidth(3);

            if(i==0){
                htmp[i]->SetMarkerStyle(20);
                htmp[i]->SetLineWidth(1);
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
            }else{
		legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
            }   
        }   

	float TextSize_ = 0.035;
        TLatex *latexComp;

        legend_nmFit->SetBorderSize(0);
        legend_nmFit->SetFillColor(0);
        legend_nmFit->SetFillStyle(0);
        legend_nmFit->SetNColumns(1);
        legend_nmFit->SetTextSize(TextSize_);
        legend_nmFit->SetTextSizePixels(25);
        legend_nmFit->Draw();


        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) + (1- fbkg) * Signal
        //         Nbkg                                                                   NSign
        float TextOffsetX = 0.01;
        float TextOffsetY = -0.01;
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
        else
            latexComp = new TLatex(0.52,0.51,"Fit model = ");
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        float fbkg = Nbkgall.getVal()/(Nbkgall.getVal()+NSign.getVal());

        sprintf(buffer,"#splitline{%1.2fxCat 0 +}{%1.2fxCat 1 +}",
                fbkg*fbkg0.getVal(), 
		fbkg*(1-fbkg0.getVal())*fbkg1.getVal()
		);
        latexComp = new TLatex(0.63,0.50,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#splitline{%1.2fxCat 2 +}{%1.2fxCat 3 +}",
		fbkg*(1-fbkg0.getVal())*(1-fbkg1.getVal())*fbkg2.getVal(),
		fbkg*(1-fbkg0.getVal())*(1-fbkg1.getVal())*(1-fbkg2.getVal())*fbkg3.getVal()
		);
        latexComp = new TLatex(0.63,0.50,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005-0.053+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#splitline{%1.2fxCat 4 +}{%1.2fxSig}",
		fbkg*(1-fbkg0.getVal())*(1-fbkg1.getVal())*(1-fbkg2.getVal())*(1-fbkg3.getVal()),
                1.-fbkg
                );
        latexComp = new TLatex(0.63,0.46,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005-0.06+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
        latexComp = new TLatex(0.52,0.56,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlack);
        latexComp->SetNDC();
        latexComp->Draw();
}


void DrawBkgOnlyFit( TH1F *h_TemplatesExtractedFromData[] ){
    	int sidx_ = 0; // dummy

        char TrimChar_[128];
        TH1F *hdata_toy;

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data_obs","data_obs",x,h_TemplatesExtractedFromData[3]);
        RooDataHist MC("MC","MC with x",x,h_MCinSignalRegion);

        RooDataHist dCate0("dCate0","dCate0set with x",x,h_TemplatesExtractedFromData[0]);
        sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

        char buffer13[256];
        sprintf(buffer13,"dCate3set with x");
        sprintf(buffer,"dCate3");
        if(usingCate1insteadofCate3){
            sprintf(buffer,"dCate1");
            sprintf(buffer13,"dCate1set with x");
        }

        //RooDataHist dCate3("dCate3","dCate3set with x",x,h_TemplatesExtractedFromData[1]);
        RooDataHist dCate3(buffer,buffer13,x,h_TemplatesExtractedFromData[1]);
        sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        if(usingCate1insteadofCate3){
            sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        }
        RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,h_TemplatesExtractedFromData[2]);
        sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

        //sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        //RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
        //RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 
        //RooRealVar fbkg4("fbkg4","fraction",0.5,0.00,1) ; 
            // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) 
            // RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
        // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) 
        RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);

	hdata_toy = (TH1F*) h_TemplatesExtractedFromData[3]->Clone();

	float xmin_ = -0.000000000001;
	float xmax_ = 1.5*hdata_toy->Integral();
	RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
        //RooRealVar Nbkgall("Nbkgall","yield",0.5,0.00,hdata_toy->Integral()) ; 
        //RooRealVar NSign("NSign","yield",0.5,0.00,hdata_toy->Integral()) ;
        //RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        //RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));
        RooAddPdf model("model","model",RooArgList(bkgmodel),RooArgList(Nbkgall));

	//model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE),Range(50,650));
	model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE),Range(50,500));

        // http://roofit-scanner-y1s-emu.googlecode.com/svn/branches/EPS-HEP_Cuts/background/One_Range/old/Fit_Multi_Exp_Var.C
        RooArgSet obs(x);
        RooArgSet* flparams = (RooArgSet*) model.getParameters(obs)->selectByAttrib("Constant",kFALSE) ; 

        RooPlot* xframe = x.frame() ; 
        //xframe->SetTitle("In signal region");
        xframe->SetTitle("");
        float ymax_ = 1.5*h_TemplatesExtractedFromData[3]->GetMaximum();
	data.plotOn(xframe);
        model.plotOn(xframe,LineColor(kBlue-2));
        model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
        //model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
        model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
        //model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
        //xframe->GetYaxis()->SetRangeUser(0.,xframe->GetYaxis()->GetXmax());
        xframe->GetYaxis()->SetRangeUser(0.0,ymax_);
        if(UsingLogyScaleInMass)
            xframe->GetYaxis()->SetRangeUser(0.001,pow(10,2.0*log10(ymax_)));
        xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
        xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
        xframe->Draw() ;
	drawCMSLumi();

        // polish the plot here
        float offsetLegend = 0.0;
        if(UsingLogyScaleInMass)
            offsetLegend=0.06;
        TLegend *legend_nmFit = new TLegend(0.5077+0.1,0.605+offsetLegend,0.9,0.87+offsetLegend);
        TH1F *htmp[6];
        TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
        sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
        string nametemp[6] = {"data","Fit model","Cat M (data)","Cat 3 (data)","Cat 4 (data)",buffer};
        if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";
        int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
        int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
        for(int i=0;i<6;i++){
            if(i==3) continue;
            if(i==5) continue;
            sprintf(buffer,"%s",nametemp[i].c_str());
            htmp[i] = new TH1F(buffer,"",10,0,10);
            htmp[i]->SetLineColor(Lcolors[i]);
            htmp[i]->SetLineStyle(Lstyle[i]);
            htmp[i]->SetLineWidth(3);

            if(i==0){
                htmp[i]->SetMarkerStyle(20);
                htmp[i]->SetLineWidth(1);
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
            }else{
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
            }   
        }   

        legend_nmFit->SetBorderSize(0);
        legend_nmFit->SetFillColor(0);
        legend_nmFit->SetFillStyle(0);
        legend_nmFit->SetNColumns(1);
        legend_nmFit->SetTextSize(0.04);
        legend_nmFit->SetTextSizePixels(25);
        legend_nmFit->Draw();

        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) * Cate4 ) + (1- fbkg) * Signal
        //         Nbkg                                                                   NSign
        float TextOffsetX = 0.01;
        float TextOffsetY = -0.01;
        TLatex *latexComp;
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
        else
            latexComp = new TLatex(0.52,0.51,"Fit model = ");
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        float fbkg = Nbkgall.getVal()/(Nbkgall.getVal());

        sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4}",
                fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
        if(usingCate1insteadofCate3)
            sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
                fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
        latexComp = new TLatex(0.63,0.50,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        //sprintf(buffer,"%1.2fxSig",
        //        1.-fbkg
        //        );
        //latexComp = new TLatex(0.63,0.46,buffer);
        //if(UsingLogyScaleInMass)
        //    latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
        //latexComp->SetTextSize(0.03);
        //latexComp->SetTextColor(kBlue-2);
        //latexComp->SetNDC();
        //latexComp->Draw();

	//sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
	sprintf(buffer,"#splitline{N(Fitted-Bkg + Extrapolation) : }{                                %1.1f#pm%1.1f}",
		Nbkgall.getVal()*h_TemplatesExtractedFromData[0]->Integral()/h_TemplatesExtractedFromData[0]->Integral(1,3),
		Nbkgall.getError()*h_TemplatesExtractedFromData[0]->Integral()/h_TemplatesExtractedFromData[0]->Integral(1,3)
                );
        latexComp = new TLatex(0.52,0.36,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        //sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
        sprintf(buffer,"#chi^{2} / ndf = %1.2f / %1.2f",xframe->chiSquare(), (float)flparams->getSize());
        latexComp = new TLatex(0.52,0.56,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlack);
        latexComp->SetNDC();
        //latexComp->Draw();

}

void randomTemplate(TH1F *hpercentagetemplate[],float percentage/*0~1*/){
    if(RunStatus_>=UncQsquarePlus){
        sprintf(buffer,"tree_data%s",
                RunStatusNames[UncQsquare].c_str());
    }else{
        sprintf(buffer,"tree_data%s",
                RunStatusNames[RunStatus_].c_str());
    }
    TTree *MCTemplatesTree = (TTree*)file->Get(buffer);

    float weight_template=-1;
    float Mass_template=-1;
    float category_template=-1;
    MCTemplatesTree->SetBranchAddress("weight",&weight_template);
    MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
    MCTemplatesTree->SetBranchAddress("category",&category_template);

    hpercentagetemplate[0] = new TH1F("temp_cate0","",numberbins,Xregions[0],Xregions[1]);
    sprintf(buffer,"temp_cate3");
    if(usingCate1insteadofCate3) 
        sprintf(buffer,"temp_cate1");
    hpercentagetemplate[1] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
    hpercentagetemplate[2] = new TH1F("temp_cate4","",numberbins,Xregions[0],Xregions[1]);
    for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
        MCTemplatesTree->GetEntry(entry);
        if ( weight_template == -1 ) continue;
        if ( Mass_template == -1 ) continue;
            /*
               cate 0 <-> region 2
               cate 3 <-> region 4
               cate 4 <-> region 1
               */
        if ( category_template == 2 ){ 
            if(rnd->Uniform(0,1)<percentage)
                hpercentagetemplate[0]->Fill(Mass_template,weight_template);
        }
        if ( category_template == 4 ){ 
            if(rnd->Uniform(0,1)<percentage)
                hpercentagetemplate[1]->Fill(Mass_template,weight_template);
        }
        if ( category_template == 1 ){ 
            if(rnd->Uniform(0,1)<percentage)
                hpercentagetemplate[2]->Fill(Mass_template,weight_template);
        }
    }
}

/*
Templates : 
cate0(data)   :   h_TemplatesExtractedFromData[0] 
cate1/3(data) :   h_TemplatesExtractedFromData[1]
cate4(data)   :   h_TemplatesExtractedFromData[2]
SR(data)      :   h_TemplatesExtractedFromData[3]
SR(MC)        :   h_MCinSignalRegion
SR(S 950)     :   h_SignalTemplateInSR
 */
TGraphErrors *SignalInjectionFitBias(int sidx_/*signal index*/,TH1F *h_MCinSignalRegion, TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[], int nToy){

        char TrimChar_[128];

        TH1F *hdata = (TH1F*) h_MCinSignalRegion->Clone();
        hdata->Add(h_SignalTemplateInSR[sidx_]);

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data_obs","data_obs",x,hdata);

        RooDataHist dCate0("dCate0","dCate0set with x",x,h_TemplatesExtractedFromData[0]);
        sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

        char buffer13[256];
        sprintf(buffer,"dCate3");
        sprintf(buffer13,"dCate3set with x");
        if(usingCate1insteadofCate3){
            sprintf(buffer,"dCate1");
            sprintf(buffer13,"dCate1set with x");
        }

        RooDataHist dCate3(buffer,buffer13,x,h_TemplatesExtractedFromData[1]);
        sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        if(usingCate1insteadofCate3)
            sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,h_TemplatesExtractedFromData[2]);
        sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

        sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
        RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 

        // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) 
        //RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
        RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);

	float xmin_ = -0.000000000001;
	float xmax_ = 1.5*hdata->Integral();
	RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
	RooRealVar NSign("NSign","N(Sig)",2.,xmin_,xmax_) ;
        //RooRealVar Nbkgall("Nbkgall","yield",0.5,0.00,hdata->Integral()) ; 
        //RooRealVar NSign("NSign","yield",0.5,0.00,hdata->Integral()) ;
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));

        model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
        while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata->Integral() - 1.) > 0.01){
            printf("[In while loop (Injection)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata->Integral(), 
                    Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
            model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
        } 

        RooArgSet obs(x);
        RooArgSet* flparams = (RooArgSet*) model.getParameters(obs)->selectByAttrib("Constant",kFALSE) ;

        RooPlot* xframe = x.frame() ; 
        //xframe->SetTitle("In signal region");
        xframe->SetTitle("");
        float ymax_ = 1.5*hdata->GetMaximum();
        data.plotOn(xframe);
        model.plotOn(xframe,LineColor(kBlue-2));
        model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
        //model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
        model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
        model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
        xframe->GetYaxis()->SetRangeUser(0.,ymax_);
        if(UsingLogyScaleInMass)
            xframe->GetYaxis()->SetRangeUser(0.001,pow(10,2.0*log10(ymax_)));
        xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
        xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
        xframe->Draw() ;
	drawCMSLumi();

        // polish the plot here
        float offsetLegend = 0.0; 
        if(UsingLogyScaleInMass)
            offsetLegend=0.06;
        TLegend *legend_nmFit = new TLegend(0.5077,0.605+offsetLegend,0.9,0.87+offsetLegend);
        TH1F *htmp[6];
        TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
        sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
        string nametemp[6] = {"Bkg.+Sig.(simulation)","Fit model","Cat M(data)","Cat 3(data)","Cat 4(data)",buffer};
        if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";
        int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
        int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
        for(int i=0;i<6;i++){
            if(i==3) continue;
            sprintf(buffer,"%s",nametemp[i].c_str());
            htmp[i] = new TH1F(buffer,"",10,0,10);
            htmp[i]->SetLineColor(Lcolors[i]);
            htmp[i]->SetLineStyle(Lstyle[i]);
            htmp[i]->SetLineWidth(3);

            if(i==0){
                htmp[i]->SetMarkerStyle(20);
                htmp[i]->SetLineWidth(1);
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
            }else{
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
            }   
        }   

        legend_nmFit->SetBorderSize(0);
        legend_nmFit->SetFillColor(0);
        legend_nmFit->SetFillStyle(0);
        legend_nmFit->SetNColumns(1);
        legend_nmFit->SetTextSize(0.04);
        legend_nmFit->SetTextSizePixels(25);
        legend_nmFit->Draw();

        //old : model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        //         Nbkg                                                                   NSign
        float TextOffsetX = 0.01;
        float TextOffsetY = -0.01;
        TLatex *latexComp;
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
        else
            latexComp = new TLatex(0.52,0.51,"Fit model = ");
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        float fbkg = Nbkgall.getVal()/(Nbkgall.getVal()+NSign.getVal());
        sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4 +}",
                fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
        if(usingCate1insteadofCate3)
            sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
                    fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
        latexComp = new TLatex(0.63,0.50,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"%1.2fxSig",
                //fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
                1.-fbkg
               );
        latexComp = new TLatex(0.63,0.46,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
                Nbkgall.getVal(),
                Nbkgall.getError(),
                NSign.getVal(),
                NSign.getError()//, hdata_toy->Integral()
               );
        latexComp = new TLatex(0.52,0.36,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
        latexComp = new TLatex(0.52,0.56,buffer);
        if(UsingLogyScaleInMass)
            latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlack);
        latexComp->SetNDC();
        latexComp->Draw();

        float Injections[nToy];
        float InjectionsError[nToy];
        float Fitted[nToy];
        float FittedError[nToy];
        for(int itoy=0;itoy<nToy;itoy++){
            Injections[itoy] = 0.;
            InjectionsError[itoy] = 0.;
            Fitted[itoy] = 0.;
            FittedError[itoy] = 0.;
        }

        float step_ = 0.2;
        for(int itoy=0;itoy<nToy;itoy++){

            TH1F * hdata_toy = new TH1F("hdata_toy","",
                    h_MCinSignalRegion->GetNbinsX(),
                    h_MCinSignalRegion->GetXaxis()->GetXmin(),
                    h_MCinSignalRegion->GetXaxis()->GetXmax());

            TH1F * hSig_toy = new TH1F("hSig_toy","",
                    h_SignalTemplateInSR[sidx_]->GetNbinsX(),
                    h_SignalTemplateInSR[sidx_]->GetXaxis()->GetXmin(),
                    h_SignalTemplateInSR[sidx_]->GetXaxis()->GetXmax());

            for(int i=1;i<=h_MCinSignalRegion->GetNbinsX();i++)
                hdata_toy->Fill(h_MCinSignalRegion->GetBinCenter(i),
                        ((double)h_MCinSignalRegion->GetBinContent(i)));  
                        //(float)rnd->Poisson((double)h_MCinSignalRegion->GetBinContent(i)));  

            for(int i=1;i<=h_SignalTemplateInSR[sidx_]->GetNbinsX();i++)
                hSig_toy->Fill(h_SignalTemplateInSR[sidx_]->GetBinCenter(i),
                        step_*(itoy+1.)*((double)h_SignalTemplateInSR[sidx_]->GetBinContent(i)));  
                        //(float)rnd->Poisson((double)h_SignalTemplateInSR[sidx_]->GetBinContent(i)));  

            Injections[itoy] = step_*(itoy+1.);

            hdata_toy->Add(hSig_toy);
            float Total_yields = hdata_toy->Integral();
            float inject_yields = h_SignalTemplateInSR[sidx_]->Integral();

            RooDataHist data_toy("data_toy","dataset with x",x,hdata_toy);
            // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            //RooAddPdf model_toy("model_toy","model_toy",
            //        RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);            

            //RooAddPdf bkgmodel_toy("bkgmodel_toy","bkgmodel_toy",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
            RooAddPdf bkgmodel_toy("bkgmodel_toy","bkgmodel_toy",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);
            // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
            //         Nbkg                                                                   NSign
	    float xmin_ = -0.000000000001;
	    float xmax_ = 1.5*hdata_toy->Integral();
	    RooRealVar Nbkgall_toy("Nbkgall_toy","N(bkg)",33.,xmin_,xmax_) ;
	    RooRealVar NSign_toy("NSign_toy","N(Sig)",2.,xmin_,xmax_) ;
            //RooRealVar Nbkgall_toy("Nbkgall_toy","yield",0.5,0.00,hdata_toy->Integral()) ; 
            //RooRealVar NSign_toy("NSign_toy","yield",0.5,0.00,hdata_toy->Integral()) ;
            RooAddPdf model_toy("model_toy","model_toy",RooArgList(bkgmodel_toy,Signal),RooArgList(Nbkgall_toy,NSign_toy));

            //model_toy.fitTo(data_toy,Minimizer("Minuit2", "minimize"));
            model_toy.fitTo(data_toy,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));

            float fitting_Sig_yields = NSign_toy.getVal();
            float fitting_Sig_yields_error = NSign_toy.getError();
            Fitted[itoy] = fitting_Sig_yields/inject_yields;
            FittedError[itoy] = fitting_Sig_yields_error/inject_yields;

            std::cout<<"[Signal Injection "<< inject_yields <<" ] : "<<fitting_Sig_yields<<" +/- "<<fitting_Sig_yields_error<<std::endl;

            hdata_toy->Delete();
            hSig_toy->Delete();
        }

        TGraphErrors *bias = new TGraphErrors(nToy, Injections, Fitted, InjectionsError, FittedError );

        return bias;
}
void SignalInjectionFit(int sidx_/*signal index*/,TH1F *h_MCinSignalRegion, TH1F *h_SignalTemplateInSR[],
        TH1F *h_TemplatesExtractedFromData[], TH1F *bias, int nToy){

        TH1F *hdata = (TH1F*) h_MCinSignalRegion->Clone();
        hdata->Add(h_SignalTemplateInSR[sidx_]);

        RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
        RooDataHist data("data_obs","data_obs",x,hdata);

        RooDataHist dCate0("dCate0","dCate0set with x",x,h_TemplatesExtractedFromData[0]);
        sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

        char buffer13[256];
        sprintf(buffer,"dCate3");
        sprintf(buffer13,"dCate3set with x");
        if(usingCate1insteadofCate3){
            sprintf(buffer,"dCate1");
            sprintf(buffer13,"dCate1set with x");
        }

        RooDataHist dCate3(buffer,buffer13,x,h_TemplatesExtractedFromData[1]);
        sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        if(usingCate1insteadofCate3)
            sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

        RooDataHist dCate4("dCate4","dCate4set with x",x,h_TemplatesExtractedFromData[2]);
        sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

        sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
        RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
        RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

        RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
        RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 
        //RooRealVar fbkg4("fbkg4","fraction",0.5,0.00,1) ; 
        //RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);

        // bkgmodel = ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) 
        RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
	float xmin_ = -0.000000000001;
	float xmax_ = 1.5*hdata->Integral();
	RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
	RooRealVar NSign("NSign","N(Sig)",2.,xmin_,xmax_) ;
        //RooRealVar Nbkgall("Nbkgall","yield",0.5,0.00,hdata->Integral()) ; 
        //RooRealVar NSign("NSign","yield",0.5,0.00,hdata->Integral()) ;
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));

        //model.fitTo(data,Minimizer("Minuit2", "minimize"));
        model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
        while (fabs((Nbkgall.getVal()+NSign.getVal())/hdata->Integral() - 1.) > 0.01){
            printf("[In while loop (Injection)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", hdata->Integral(), 
                    Nbkgall.getVal(), NSign.getVal(), SAMPLE[Tgamma500GeV+sidx_].tag );   
            model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
        } 

        RooPlot* xframe = x.frame() ; 
        //xframe->SetTitle("In signal region");
        xframe->SetTitle("");
        float ymax_ = 1.5*hdata->GetMaximum();
        data.plotOn(xframe);
        model.plotOn(xframe,LineColor(kBlue-2));
        model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
        model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
        model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
        model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
        xframe->GetYaxis()->SetRangeUser(0.,ymax_);
        xframe->Draw() ;
	drawCMSLumi();


        TLegend *legend_nmFit = new TLegend(0.5077,0.605,0.9,0.87);
        TH1F *htmp[6];
        //string nametemp[6] = {"MC+signal injection","Fit model","Cate0(data)","Cate3(data)","Cate4(data)",SAMPLE[Tgamma500GeV+sidx_].tag};
        string nametemp[6] = {"MC bkg.+signal injection","Fit model","Cat M(data)","Cat 3(data)","Cat 4(data)",SAMPLE[Tgamma500GeV+sidx_].tag};
        if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";
        int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
        int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
        for(int i=0;i<6;i++){
            sprintf(buffer,"%s",nametemp[i].c_str());
            htmp[i] = new TH1F(buffer,"",10,0,10);
            htmp[i]->SetLineColor(Lcolors[i]);
            htmp[i]->SetLineStyle(Lstyle[i]);
            htmp[i]->SetLineWidth(3);

            if(i==0){
                htmp[i]->SetMarkerStyle(20);
                htmp[i]->SetLineWidth(1);
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
            }else{
                legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
            }   
        }   

        legend_nmFit->SetBorderSize(0);
        legend_nmFit->SetFillColor(0);
        legend_nmFit->SetFillStyle(0);
        legend_nmFit->SetNColumns(1);
        legend_nmFit->SetTextSize(0.04);
        legend_nmFit->SetTextSizePixels(25);
        legend_nmFit->Draw();

        //old : model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
        // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
        //         Nbkg                                                                   NSign
        TLatex *latexComp = new TLatex(0.52,0.51,"Model = ");
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        float fbkg = Nbkgall.getVal()/(Nbkgall.getVal()+NSign.getVal());
        //sprintf(buffer,"#splitline{%1.1ExCate0 +}{%1.1ExCate3 +}",
        sprintf(buffer,"#splitline{%1.1ExCat M +}{%1.1ExCat 3 +}",
                fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());
        if(usingCate1insteadofCate3)
            //sprintf(buffer,"#splitline{%1.1ExCate0 +}{%1.1ExCate1 +}",
            sprintf(buffer,"#splitline{%1.1ExCat M +}{%1.1ExCat 1 +}",
                    fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());

        latexComp = new TLatex(0.63,0.50,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#splitline{%1.1ExCat 4 + }{%1.1ExSig}",
                fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
                1.-fbkg
               );
        latexComp = new TLatex(0.63,0.44,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#splitline{Fitted-Bkg : %1.1f#pm%1.1f}{Fitted-Sig : %1.1f#pm%1.1f }",
                Nbkgall.getVal(),
                Nbkgall.getError(),
                NSign.getVal(),
                NSign.getError()//, hdata_toy->Integral()
               );
        latexComp = new TLatex(0.52,0.36,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        //sprintf(buffer,"#chi^{2}/ndf = %g",xframe->chiSquare(3));
        sprintf(buffer,"#chi^{2} = %g",xframe->chiSquare());
        latexComp = new TLatex(0.52,0.56,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlack);
        latexComp->SetNDC();
        latexComp->Draw();

        for(int itoy=0;itoy<nToy;itoy++){

            TH1F * hdata_toy = new TH1F("hdata_toy","",
                    h_MCinSignalRegion->GetNbinsX(),
                    h_MCinSignalRegion->GetXaxis()->GetXmin(),
                    h_MCinSignalRegion->GetXaxis()->GetXmax());

            TH1F * hSig_toy = new TH1F("hSig_toy","",
                    h_SignalTemplateInSR[sidx_]->GetNbinsX(),
                    h_SignalTemplateInSR[sidx_]->GetXaxis()->GetXmin(),
                    h_SignalTemplateInSR[sidx_]->GetXaxis()->GetXmax());

            for(int i=1;i<=h_MCinSignalRegion->GetNbinsX();i++)
                hdata_toy->Fill(h_MCinSignalRegion->GetBinCenter(i),
                        (float)rnd->Poisson((double)h_MCinSignalRegion->GetBinContent(i)));  

            for(int i=1;i<=h_SignalTemplateInSR[sidx_]->GetNbinsX();i++)
                hSig_toy->Fill(h_SignalTemplateInSR[sidx_]->GetBinCenter(i),
                        (float)rnd->Poisson((double)h_SignalTemplateInSR[sidx_]->GetBinContent(i)));  

            hdata_toy->Add(hSig_toy);
            float Total_yields = hdata_toy->Integral();
            int inject_yields = h_SignalTemplateInSR[sidx_]->Integral();

            RooDataHist data_toy("data_toy","dataset with x",x,hdata_toy);
            // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
            //RooAddPdf model_toy("model_toy","model_toy",
            //        RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);            

            RooAddPdf bkgmodel_toy("bkgmodel_toy","bkgmodel_toy",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
            // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
            //         Nbkg                                                                   NSign
            RooAddPdf model_toy("model_toy","model_toy",RooArgList(bkgmodel_toy,Signal),RooArgList(Nbkgall,NSign));

            //model_toy.fitTo(data_toy,Minimizer("Minuit2", "minimize"));
            model_toy.fitTo(data_toy,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));

            //float fitting_Sig_yields = Total_yields * (1-fbkg0.getVal())*(1.-fbkg3.getVal())*(1.- fbkg4.getVal());
            float fitting_Sig_yields = NSign.getVal();

            /*
               float newfbkg[3] = {
               (1-fbkg0.getVal()),
               (1-fbkg3.getVal()),
               (1-fbkg4.getVal())};
               float newfbkgerr[3] = {
               fbkg0.getError(),
               fbkg3.getError(),
               fbkg4.getError()
               };

               float ab = newfbkg[0] * newfbkg[1];
               float aberr = ab*sqrt(pow(newfbkgerr[0]/newfbkg[0],2) + pow(newfbkgerr[1]/newfbkg[1],2));
               float abc = ab*newfbkg[2];
               float abcerr= abc*sqrt( pow(aberr/ab,2) + pow( newfbkgerr[2]/newfbkg[2],2) );
             */

            //float fitting_Sig_yields_error = Total_yields * abcerr;
            float fitting_Sig_yields_error = NSign.getError();

            //bias->Fill((fitting_Sig_yields-inject_yields)/inject_yields);
            //bias->Fill((fitting_Sig_yields - inject_yields)/fitting_Sig_yields_error);
            std::cout<<"[Signal Injection "<< inject_yields <<" ] : "<<fitting_Sig_yields<<" +/- "<<fitting_Sig_yields_error<<std::endl;

            /*
            RooRealVar Nbkg0("Nbkg0","yield",0.5,0.00,hdata_toy->Integral()) ; 
            RooRealVar Nbkg3("Nbkg3","yield",0.5,0.00,hdata_toy->Integral()) ; 
            RooRealVar Nbkg4("Nbkg4","yield",0.5,0.00,hdata_toy->Integral()) ; 
            RooRealVar NSign("NSign","yield",0.5,0.00,hdata_toy->Integral()) ; 
            RooAddPdf model_toy2("model_toy2","model_toy2",
                    RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(Nbkg0,Nbkg3,Nbkg4,NSign));            
            //model_toy2.fitTo(data_toy,SumW2Error(kTRUE),Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE));
            model_toy2.fitTo(data_toy,SumW2Error(kFALSE),Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE));
            //model_toy2.fitTo(data_toy,Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE));
            //model_toy2.fitTo(data_toy,Extended(kTRUE),Hesse(kTRUE),Minos(kTRUE),Minimizer("Minuit2", "minimize"));
            //model_toy2.fitTo(data_toy,Extended(kTRUE),Minimizer("Minuit2", "minimize"));
            std::cout<<"[2nd Signal Injection "<< (int)inject_yields <<" / "<<Total_yields
                <<" ] : "<< NSign.getVal() <<" +/- "<< NSign.getError() <<" ( "<<NSign.getErrorLo()<< " + "<<
                NSign.getErrorHi()<<")"<<std::endl;
            */
            if(NSign.getVal()<inject_yields){
                //bias->Fill((- NSign.getVal() + (int)inject_yields)/NSign.getErrorHi());
                bias->Fill((- NSign.getVal() + inject_yields)/NSign.getErrorHi());
            }else{
                //bias->Fill((NSign.getVal() - (int)inject_yields)/NSign.getErrorLo());
                bias->Fill((NSign.getVal() - inject_yields)/NSign.getErrorLo());
            }
            //bias->Fill((NSign.getVal() - (int)inject_yields)/NSign.getError());

            hdata_toy->Delete();
            hSig_toy->Delete();
        }

}

void DrawGausFitResult(TF1 *func){
        float mean_value = func->GetParameter(1);
        float sigma_value = func->GetParameter(2);

        sprintf(buffer,"mean : %1.2f #pm %1.2f",mean_value,func->GetParError(1));
        TLatex *latexComp = new TLatex(0.63,0.8,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();

        sprintf(buffer,"#sigma : %1.2f #pm %1.2f",sigma_value,func->GetParError(2));
        latexComp = new TLatex(0.63,0.75,buffer);
        latexComp->SetTextSize(0.03);
        latexComp->SetTextColor(kBlue-2);
        latexComp->SetNDC();
        latexComp->Draw();
}

void loadTemplates(){
    // Load systematic uncertainties on background shapes
    if(AppliedSysOnShape){
        TFile *fileSysUnc = new TFile("Systematics/sys/sysforshape.root");
        if(!fileSysUnc) {
            printf("[Error] There is no Systematics/sys/sysforshape.root\n");
            //exit(1);
        }
        for(int icate=0;icate<NCate_;icate++){
            sprintf(buffer,"SumSys_cate%i",icate);
            h_MergeTemplatesUnc[icate] = (TH1F*) fileSysUnc->Get(buffer);

            sprintf(buffer,"SumSys_cate%i in SR",icate);
            h_MergeTemplatesUncInSR[icate] = (TH1F*) fileSysUnc->Get(buffer);
        }
    }
    output_->cd();

    for(int i=0;i<Tgamma1500GeV-Tgamma500GeV+1;i++){
        for(int j=0;j<3;j++){
            NBgkWithUncOnMerge[i][0][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnMerge[i][1][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnMerge[i][2][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnSRandCR[i][0][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnSRandCR[i][1][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnSRandCR[i][2][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnDATAandMC[i][0][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnDATAandMC[i][1][j] = 0;  //[sample][template][down, normal, up]
            NBgkWithUncOnDATAandMC[i][2][j] = 0;  //[sample][template][down, normal, up]
        }
    }

    for(int i=0;i<=Tgamma1500GeV-Tgamma500GeV;i++){
        sprintf(buffer,"SRCanvas(%s)",SAMPLE[Tgamma500GeV+i].tag);
        SRCanvas[i] = new TCanvas(buffer,"",640,640);

        sprintf(buffer,"ClosureSRCanvas(%s)",SAMPLE[Tgamma500GeV+i].tag);
        ClosureSRCanvas[i] = new TCanvas(buffer,"",640,640);

        for(int j=0;j<3;j++)
            fractions_[i][j] = -1;
    }
    for(int i=0;i<NCate_;i++){
        sprintf(buffer,"c_%i",i);
        canvas[i] = new TCanvas(buffer,"",640,640);
        sprintf(buffer,"c_CompSRandCR_%i",i);
        canvasCompSRandCR[i] = new TCanvas(buffer,"",640,640);
        sprintf(buffer,"c_CompDATAandMC_%i",i);
        canvasCompDATAandMC[i] = new TCanvas(buffer,"",640,640);
    }

    for(int i=0;i<NCate_;i++){
        legend_nm[i] = new TLegend(0.577,0.305,0.9,0.87);
        legend_nm[i]->SetBorderSize(0);
        legend_nm[i]->SetFillColor(0);
        legend_nm[i]->SetFillStyle(0);
        legend_nm[i]->SetNColumns(1);
        legend_nm[i]->SetTextSize(0.04);
        legend_nm[i]->SetTextSizePixels(25);

        legend_nm_ComparisonSRandCR[i] = new TLegend(0.577,0.305+0.3,0.9,0.87);
        legend_nm_ComparisonSRandCR[i]->SetBorderSize(0);
        legend_nm_ComparisonSRandCR[i]->SetFillColor(0);
        legend_nm_ComparisonSRandCR[i]->SetFillStyle(0);
        legend_nm_ComparisonSRandCR[i]->SetNColumns(1);
        legend_nm_ComparisonSRandCR[i]->SetTextSize(0.04);
        legend_nm_ComparisonSRandCR[i]->SetTextSizePixels(25);

        legend_nm_ComparisonDATAandMC[i] = new TLegend(0.577,0.305+0.3,0.9,0.87);
        legend_nm_ComparisonDATAandMC[i]->SetBorderSize(0);
        legend_nm_ComparisonDATAandMC[i]->SetFillColor(0);
        legend_nm_ComparisonDATAandMC[i]->SetFillStyle(0);
        legend_nm_ComparisonDATAandMC[i]->SetNColumns(1);
        legend_nm_ComparisonDATAandMC[i]->SetTextSize(0.04);
        legend_nm_ComparisonDATAandMC[i]->SetTextSizePixels(25);
    }
    Mergelegend_nm = new TLegend(0.577,0.305,0.9,0.87);
    Mergelegend_nm->SetBorderSize(0);
    Mergelegend_nm->SetFillColor(0);
    Mergelegend_nm->SetFillStyle(0);
    Mergelegend_nm->SetNColumns(1);
    Mergelegend_nm->SetTextSize(0.04);
    Mergelegend_nm->SetTextSizePixels(25);

    for(int idxCate_ =0 ;idxCate_<NCate_;idxCate_++){
        sprintf(buffer,"Merge_cate-%i",idxCate_);
        h_MergeTemplates[idxCate_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
        h_MergeTemplates[idxCate_]->Sumw2();
        h_MergeTemplates[idxCate_]->SetLineWidth(2);
        h_MergeTemplates[idxCate_]->SetLineColor(colors[idxCate_]);
        h_MergeTemplates[idxCate_]->SetMarkerStyle(1);

        sprintf(buffer,"Merge_cate-%i in SR",idxCate_);
        h_MergeTemplatesInSR[idxCate_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
        h_MergeTemplatesInSR[idxCate_]->Sumw2();
        h_MergeTemplatesInSR[idxCate_]->SetLineWidth(2);
        h_MergeTemplatesInSR[idxCate_]->SetLineColor(1);
        h_MergeTemplatesInSR[idxCate_]->SetMarkerStyle(1);
    }
    h_MCinSignalRegion = new TH1F("MC in signal region","",numberbins,Xregions[0],Xregions[1]);
    h_MCinSignalRegion->Sumw2();
    h_MCinSignalRegionCate4AndRest = new TH1F("h_MCinSignalRegionCate4AndRest","",numberbins,Xregions[0],Xregions[1]);
    h_MCinSignalRegionCate4AndRest->Sumw2();

    h_MCinSignalRegion2Photons = new TH1F("h_MCinSignalRegion2Photons","",numberbins,Xregions[0],Xregions[1]);
    h_MCinSignalRegion2Photons->Sumw2();

    const int NSRMC = 13;
    double SRMCYields[NSRMC] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    double SRMCErrors[NSRMC] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    string SRMCCateName[NSRMC] = {
                "isolated-photon + 1 fake photon from lepton",
                "isolated-photon + 1 fake photon from jet ",
                "fake photons from leptons",
                "fake photon from lepton + 1 fake photon from jet ",
                "fake photons from jets",
		"2 prompt photon",
		"1 prompt photon + 1 non-prompt photon (ISR/FSR/devay in flight)",
		"2 non-prompt photon (ISR/FSR/devay in flight)",
		"1 prompt photon + 1 unknown",
		"1 non-prompt photon (ISR/FSR/devay in flight) + 1 unknown",
                "1 fake photon from lepton + 1 unknown",
                "1 fake photon from jet + 1 unknown",
		"unknown"
    };

    for(int i=0; i<samples_order_size-_sample_start; i++){

        /*
           UncQsquarePlus,
           UncQsquareMinus,
           UncMatchingPlus,
           UncMatchingMinus,
         */

        if(RunStatus_>=UncQsquarePlus){
            sprintf(buffer,"tree_%s%s",
                    SAMPLE[i+_sample_start].tag,RunStatusNames[UncQsquare].c_str());

            if(RunStatus_==UncQsquarePlus){
                if(i+_sample_start==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(i+_sample_start==DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncQsquareMinus){
                if(i+_sample_start==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(i+_sample_start==DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncMatchingPlus){
                if(i+_sample_start==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                //if(i+_sample_start==DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1)
                //    sprintf(buffer,"tree_%s%s",
                //            SAMPLE[DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                //            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncMatchingMinus){
                if(i+_sample_start==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1)
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                //if(i+_sample_start==DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1)
                //    sprintf(buffer,"tree_%s%s",
                //            SAMPLE[DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                //            RunStatusNames[UncQsquare].c_str());
            }
        }else{
            sprintf(buffer,"tree_%s%s",
                    SAMPLE[i+_sample_start].tag,RunStatusNames[RunStatus_].c_str());
        }
        printf("[Loading] %s\n",buffer);
        //if(i==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1-_sample_start)
        //    sprintf(buffer,"tree_TTJets_powheg");
        TTree *MCTemplatesTree = (TTree*)file->Get(buffer);

        for(int idxCate_ =0 ;idxCate_<NCate_;idxCate_++){
            sprintf(buffer,"%s_cate-%i",
                    SAMPLE[i+_sample_start].tag,
                    idxCate_);
            h_templates[i][idxCate_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
            h_templates[i][idxCate_]->Sumw2();
            h_templates[i][idxCate_]->SetLineWidth(2);

            sprintf(buffer,"%s_cate-%i in SR",
                    SAMPLE[i+_sample_start].tag,
                    idxCate_);
            h_templatesInSR[i][idxCate_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
            h_templatesInSR[i][idxCate_]->Sumw2();
            h_templatesInSR[i][idxCate_]->SetLineWidth(2);
        }

        int sample_ID=-1;
        int Nr_template=-1;
        int Nj_template=-1;
        int Nl_template=-1;
        int Nlloose_template=-1;
        int Nrloose_template=-1;
        float weight_template=-1;
        float Mass_template=-1;
        float category_template=-1;
        int MCTruth1 = -1; 
        int MCTruth2 = -1; 
        int MCTruth3 = -1; 

        MCTemplatesTree->SetBranchAddress("sampleID",&sample_ID);
        MCTemplatesTree->SetBranchAddress("weight",&weight_template);
        MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
        MCTemplatesTree->SetBranchAddress("category",&category_template);
        MCTemplatesTree->SetBranchAddress("Nr",&Nr_template);
        MCTemplatesTree->SetBranchAddress("Nj",&Nj_template);
        MCTemplatesTree->SetBranchAddress("Nl",&Nl_template);
        MCTemplatesTree->SetBranchAddress("Nl_loose",&Nlloose_template);
        MCTemplatesTree->SetBranchAddress("Nr_loose",&Nrloose_template);
        MCTemplatesTree->SetBranchAddress("MCTruth1",&MCTruth1);
        MCTemplatesTree->SetBranchAddress("MCTruth2",&MCTruth2);
        MCTemplatesTree->SetBranchAddress("MCTruth3",&MCTruth3);


        for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
            MCTemplatesTree->GetEntry(entry);
            if ( weight_template == -1 ) continue;
            if ( Mass_template == -1 ) continue;

            // purity for cate4 in control region of 1 lep + >= 6jets
            if(category_template == 1&&MCTruth1>=0&&MCTruth2>=0){
            //if(category_template == 1){
                purity[2][0]+=weight_template;
                purityError[2][0]+=weight_template*weight_template;
                if( (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) ){
                    purity[2][1]+=weight_template;
                    purityError[2][1]+=weight_template*weight_template;
                }
            }
            //weight_template = 1; // no weighting factor
            /*
               - MC Truth table : 
                -2: matched, but not considered
                -1: unknown 
                1 : W 
                2 : Z 
                3 : from b or B meson
                4 : from c or D meson/bayron
                5 : matched to a parton (q,g)
                x6 : photon
                    0x6 : prompt photon
                    1x6: photon from decay in flight
                    2x6: photon from ISR 
                    3x6: photon from FSR 
                7 : tau 
                8 : muon
                9 : electron

               - 5 categories + 1 signal region
                0 : 1 isolated-photon + 1 fake photon from lepton
                1 : 1 isolated-photon + 1 fake photon from jet 
                2 : 2 fake photons from leptons
                3 : 1 fake photon from lepton + 1 fake photon from jet 
                4 : 2 fake photons from jets
                SR: 2 isolated-photons

                - 6 control regions (using category_template due to history):category_template should be named as "control region"
                    0 : NPhotons==1&&NLeptons==1&&NLeptonsV==1&&NJets>=5 
                    1 : NPhotons==0&&NLeptons==1&&NLeptonsV==1&&NJets>=6
                    2 : NPhotons==1&&NLeptons==2&&NLeptonsV==2&&NJets>=4
                    3 : NPhotons==0&&NLeptons>=3&&NLeptonsV==3&&NJets>=4
                    4 : NPhotons==0&&NLeptons==2&&NLeptonsV==2&&NJets>=5
                    10: NPhotons>=2&&NLeptons>=1&&NJets>=4

               */

            if(category_template == 10){
                h_MCinSignalRegion->Fill(Mass_template,weight_template);
	    }

            if(!(MCTruth3==1||MCTruth3==2||MCTruth3==7||MCTruth3==8||MCTruth3==9)) continue;

            if( (MCTruth1%6==0&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2%6==0&&(MCTruth1>=7&&MCTruth1<=9)) ){
                if(category_template!=10)
                h_templates[i][0]->Fill(Mass_template,weight_template);

                if(category_template==10)
                    h_templatesInSR[i][0]->Fill(Mass_template,weight_template);
                if(category_template==10)
		    SRMCYields[0] += weight_template;
                if(category_template==10)
		    SRMCErrors[0] += weight_template*weight_template;
                if(category_template==10)
		    printf("[MC components in SR : Cate0] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
            }
            if( (MCTruth1%6==0&&(MCTruth2>=3&&MCTruth2<=5)) || (MCTruth2%6==0&&(MCTruth1>=3&&MCTruth1<=5)) ){
                if(category_template!=10)
                h_templates[i][1]->Fill(Mass_template,weight_template);

                if(category_template==10)
                    h_templatesInSR[i][1]->Fill(Mass_template,weight_template);
                if(category_template==10)
		    SRMCYields[1] += weight_template;
                if(category_template==10)
		    SRMCErrors[1] += weight_template*weight_template;
                if(category_template==10)
		    printf("[MC components in SR : Cate1] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
            }
            if( (MCTruth2>=7&&MCTruth2<=9) && (MCTruth1>=7&&MCTruth1<=9) ){
                if(category_template!=10)
                h_templates[i][2]->Fill(Mass_template,weight_template);

                if(category_template==10)
                    h_templatesInSR[i][2]->Fill(Mass_template,weight_template);
                if(category_template==10)
		    SRMCYields[2] += weight_template;
                if(category_template==10)
		    SRMCErrors[2] += weight_template*weight_template;
                if(category_template==10)
		    printf("[MC components in SR : Cate2] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
            }
            if( ((MCTruth1>=7&&MCTruth1<=9)&&(MCTruth2>=3&&MCTruth2<=5)) || 
                    ((MCTruth2>=7&&MCTruth2<=9)&&(MCTruth1>=3&&MCTruth1<=5)) ){
                if(category_template!=10)
                h_templates[i][3]->Fill(Mass_template,weight_template);

                if(category_template==10)
                    h_templatesInSR[i][3]->Fill(Mass_template,weight_template);
                if(category_template==10)
		    SRMCYields[3] += weight_template;
                if(category_template==10)
		    SRMCErrors[3] += weight_template*weight_template;
                if(category_template==10)
		    printf("[MC components in SR : Cate3] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
            }
            if( (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) ){
                if(category_template!=10)
                h_templates[i][4]->Fill(Mass_template,weight_template);

                if(category_template==10)
                    h_templatesInSR[i][4]->Fill(Mass_template,weight_template);
                if(category_template==10)
		    SRMCYields[4] += weight_template;
                if(category_template==10)
		    SRMCErrors[4] += weight_template*weight_template;
                if(category_template==10)
		    printf("[MC components in SR : Cate4] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		if(category_template==10)
		    h_MCinSignalRegionCate4AndRest->Fill(Mass_template,weight_template);
            }

	    if(category_template==10){
		if(MCTruth1==0&&MCTruth2==0){
		    SRMCYields[5] += weight_template;
		    SRMCErrors[5] += weight_template*weight_template;
		    printf("[MC components in SR : Cate5] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		    h_MCinSignalRegionCate4AndRest->Fill(Mass_template,weight_template);
		    h_MCinSignalRegion2Photons->Fill(Mass_template,weight_template);
		}

		if(
			((MCTruth1==0&&MCTruth2%6==0&&MCTruth2/6!=0) ||
			(MCTruth2==0&&MCTruth1%6==0&&MCTruth1/6!=0) )
		  ){
		    SRMCYields[6] += weight_template;
		    SRMCErrors[6] += weight_template*weight_template;
		    printf("[MC components in SR : Cate6] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		    h_MCinSignalRegionCate4AndRest->Fill(Mass_template,weight_template);
		    h_MCinSignalRegion2Photons->Fill(Mass_template,weight_template);
		}

		if(
			((MCTruth2%6==0&&MCTruth2/6!=0) &&
			(MCTruth1%6==0&&MCTruth1/6!=0) )
		  ){
		    SRMCYields[7] += weight_template;
		    SRMCErrors[7] += weight_template*weight_template;
		    printf("[MC components in SR : Cate7] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		    h_MCinSignalRegionCate4AndRest->Fill(Mass_template,weight_template);
		    h_MCinSignalRegion2Photons->Fill(Mass_template,weight_template);
		}

		if(
			((MCTruth1==0&&MCTruth2==-1) ||
			(MCTruth2==0&&MCTruth1==-1) )
		  ){
		    SRMCYields[8] += weight_template;
		    SRMCErrors[8] += weight_template*weight_template;
		    printf("[MC components in SR : Cate8] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		}

		if(
			((MCTruth1==-1&&MCTruth2%6==0&&MCTruth2/6!=0) ||
			(MCTruth2==-1&&MCTruth1%6==0&&MCTruth1/6!=0) )
		  ){
		    SRMCYields[9] += weight_template;
		    SRMCErrors[9] += weight_template*weight_template;
		    printf("[MC components in SR : Cate9] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		}

		if(
			((MCTruth1==-1&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2==-1&&(MCTruth1>=7&&MCTruth1<=9)) )
			) {
		    SRMCYields[10] += weight_template;
		    SRMCErrors[10] += weight_template*weight_template;
		    printf("[MC components in SR : Cate10] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		}

		if(
			((MCTruth1==-1&&(MCTruth2>=3&&MCTruth2<=5)) || (MCTruth2==-1&&(MCTruth1>=3&&MCTruth1<=5)) )
			) {
		    SRMCYields[11] += weight_template;
		    SRMCErrors[11] += weight_template*weight_template;
		    printf("[MC components in SR : Cate11] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		}
	    }
	    // no matched
            if(! 
		    ((MCTruth1%6==0&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2%6==0&&(MCTruth1>=7&&MCTruth1<=9)) || 
		    (MCTruth1%6==0&&(MCTruth2>=3&&MCTruth2<=5)) || (MCTruth2%6==0&&(MCTruth1>=3&&MCTruth1<=5)) ||
		    (MCTruth2>=7&&MCTruth2<=9) && (MCTruth1>=7&&MCTruth1<=9) ||
		    ((MCTruth1>=7&&MCTruth1<=9)&&(MCTruth2>=3&&MCTruth2<=5))||
		    	((MCTruth2>=7&&MCTruth2<=9)&&(MCTruth1>=3&&MCTruth1<=5)) ||
			(MCTruth1==0&&MCTruth2==0)||
			((MCTruth1==0&&MCTruth2%6==0&&MCTruth2/6!=0) ||
			(MCTruth2==0&&MCTruth1%6==0&&MCTruth1/6!=0) ) ||
			((MCTruth2%6==0&&MCTruth2/6!=0) &&
			(MCTruth1%6==0&&MCTruth1/6!=0) ) ||
			((MCTruth1==0&&MCTruth2==-1) ||
			(MCTruth2==0&&MCTruth1==-1) ) ||
			((MCTruth1==-1&&MCTruth2%6==0&&MCTruth2/6!=0) ||
			(MCTruth2==-1&&MCTruth1%6==0&&MCTruth1/6!=0) ) ||
			((MCTruth1==-1&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2==-1&&(MCTruth1>=7&&MCTruth1<=9)) ) ||
			((MCTruth1==-1&&(MCTruth2>=3&&MCTruth2<=5)) || (MCTruth2==-1&&(MCTruth1>=3&&MCTruth1<=5)) ) ||
		    (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) )
	      )
                if(category_template==10){
		    SRMCYields[12] += weight_template;
		    SRMCErrors[12] += weight_template*weight_template;
		    printf("[MC components in SR : Cate12] %s : %f MCTruth-Photon(%i,%i)\n",
			    SAMPLE[i+_sample_start].tag, weight_template, MCTruth1, MCTruth2 );
		}

            // purity for cate0 in control region of 1photon + 2leptons +>= 4jets
            if(category_template == 2&&MCTruth1>=0&&MCTruth2>=0){
                purity[0][0]+=weight_template;
                purityError[0][0]+=weight_template*weight_template;
                if( (MCTruth1%6==0&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2%6==0&&(MCTruth1>=7&&MCTruth1<=9)) ){
                    purity[0][1]+=weight_template;
                    purityError[0][1]+=weight_template*weight_template;
                }
            }
            // purity for cate3 in control region of 2leptons +>= 5jets
            if(category_template == 4&&MCTruth1>=0&&MCTruth2>=0){
                purity[1][0]+=weight_template;
                purityError[1][0]+=weight_template*weight_template;
                if( ((MCTruth1>=7&&MCTruth1<=9)&&(MCTruth2>=3&&MCTruth2<=5)) || 
                        ((MCTruth2>=7&&MCTruth2<=9)&&(MCTruth1>=3&&MCTruth1<=5)) ){
                    purity[1][1]+=weight_template;
                    purityError[1][1]+=weight_template*weight_template;
                }
            }
            // purity for cate4 in control region of 1 lep + >= 6jets
            //if(category_template == 1&&MCTruth1>=0&&MCTruth2>=0){
            //    purity[2][0]+=weight_template;
            //    purityError[2][0]+=weight_template*weight_template;
            //    if( (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) ){
            //        purity[2][1]+=weight_template;
            //        purityError[2][1]+=weight_template*weight_template;
            //    }
            //}
        }
	
        for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
            //if(h_templates[i][idxCate_]->GetEntries()<statistics) continue;
            canvas[idxCate_]->cd();
            if(IsPhysicsProcessConsidered)
                h_templates[i][idxCate_]->Scale(1./h_templates[i][idxCate_]->Integral()); // normalization is excluded 
            if(idxCate_==0)
                h_templates[i][idxCate_]->Write();  // for check
            h_MergeTemplates[idxCate_]->Add(h_templates[i][idxCate_]);  
            h_templates[i][idxCate_]->SetLineColor(colors[i]);
            if(counters[idxCate_]==0){
                h_templates[i][idxCate_]->GetYaxis()->SetRangeUser(ranges[0],ranges[1]);
                h_templates[i][idxCate_]->SetTitle(CateTitles[idxCate_].c_str());
                h_templates[i][idxCate_]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                h_templates[i][idxCate_]->Draw();
            }else{
                h_templates[i][idxCate_]->Draw("same");
            }
            sprintf(buffer,"%s",
                    SAMPLE[i+_sample_start].tag);
            legend_nm[idxCate_]->AddEntry(h_templates[i][idxCate_],buffer,"l");
            legend_nm[idxCate_]->Draw();
            counters[idxCate_]++;
        }
        for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
            //h_templatesInSR[i][idxCate_]->Scale(1./h_templatesInSR[i][idxCate_]->Integral());
            //if(h_templatesInSR[i][idxCate_]->Integral()>0)
            //    h_MergeTemplatesInSR[idxCate_]->Add(h_templatesInSR[i][idxCate_]);  // normalization is excluded
        }

    }

    for(int icate=0;icate<NSRMC;icate++){
	double error_ = 0.;
	double integral_ = h_MCinSignalRegion->IntegralAndError(
		1/*first bin*/, h_MCinSignalRegion->GetNbinsX(),error_, "");
	double deltaAoverNa = sqrt(SRMCErrors[icate])/SRMCYields[icate];
	double deltaBoverNb = error_/integral_;
	printf("[MC components in SR] Cate%i(%s) : %1.1f +/- %1.1f %%\n",
		icate,
		SRMCCateName[icate].c_str(),
		100.*SRMCYields[icate]/h_MCinSignalRegion->Integral(),
		100.*SRMCYields[icate]/h_MCinSignalRegion->Integral() * 
		sqrt(deltaAoverNa*deltaAoverNa + deltaBoverNb*deltaBoverNb)
		);
    }

    // with Alex's comment for bin 8, 9
    // WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2  - _sample_start
    double wjet_error = 0.; // take the last one
    for(int ibin = 1; ibin <= h_MergeTemplates[4]->GetNbinsX(); ibin++){
        if(h_templates[WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2-_sample_start][4]->GetBinError(ibin) != 0)
    	wjet_error = h_templates[WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2-_sample_start][4]->GetBinError(ibin) ;
    }    
    h_MergeTemplates[4]->SetBinError( 8,  sqrt(pow(h_MergeTemplates[4]->GetBinError(8),2) + wjet_error*wjet_error));
    h_MergeTemplates[4]->SetBinError( 9,  sqrt(pow(h_MergeTemplates[4]->GetBinError(9),2) + wjet_error*wjet_error));
    h_MergeTemplates[4]->SetBinError( 10,  sqrt(pow(h_MergeTemplates[4]->GetBinError(10),2) + wjet_error*wjet_error));
    
    h_MCinSignalRegion->Write();
    h_MCinSignalRegionCate4AndRest->Write();
    h_MCinSignalRegion2Photons->Write();
    // Extract templates from Data
    if(true){
        if(RunStatus_>=UncQsquarePlus){
            sprintf(buffer,"tree_data%s",
                    RunStatusNames[UncQsquare].c_str());
        }else{
            sprintf(buffer,"tree_data%s",
                    RunStatusNames[RunStatus_].c_str());
        }
        TTree *MCTemplatesTree = (TTree*)file->Get(buffer);
        h_TemplatesExtractedFromData[0] = new TH1F("data_cate-0","",numberbins,Xregions[0],Xregions[1]);
        h_TemplatesExtractedFromData[1] = new TH1F("data_cate-3","",numberbins,Xregions[0],Xregions[1]);
        h_TemplatesExtractedFromData[2] = new TH1F("data_cate-4","",numberbins,Xregions[0],Xregions[1]);
        h_TemplatesExtractedFromData[3] = new TH1F("data (in signal region)","",numberbins,Xregions[0],Xregions[1]);

        for(int i=0;i<4;i++){
            h_TemplatesExtractedFromData[i]->Sumw2();
            h_TemplatesExtractedFromData[i]->SetLineWidth(2);
            h_TemplatesExtractedFromData[i]->SetLineColor(1);
            h_TemplatesExtractedFromData[i]->SetMarkerStyle(21);
            h_TemplatesExtractedFromData[i]->SetMarkerSize(1.0);
        }

        for(int i=0;i<5;i++){
            sprintf(buffer,"data_cate_%i_tmp",i);
            h_TemplatesExtractedFromDataTemp[i] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
            h_TemplatesExtractedFromDataTemp[i]->Sumw2();
            h_TemplatesExtractedFromDataTemp[i]->SetLineWidth(2);
            h_TemplatesExtractedFromDataTemp[i]->SetLineColor(1);
            h_TemplatesExtractedFromDataTemp[i]->SetMarkerStyle(21);
            h_TemplatesExtractedFromDataTemp[i]->SetMarkerSize(1.0);

            sprintf(buffer,"fit_data_cate_%i_tmp",i);
            f_TemplatesExtractedFromDataTemp[i] = new TF1(buffer,"exp([0]+[1]*x)");
            f_TemplatesExtractedFromDataTemp[i]->SetParameters(3.34434e+00,-1.30296e-02);
        }

        int sample_ID=-1;
        int Nr_template=-1;
        int Nj_template=-1;
        int Nl_template=-1;
        int Nlloose_template=-1;
        int Nrloose_template=-1;
        float weight_template=-1;
        float Mass_template=-1;
        float category_template=-1;
        int MCTruth1 = -1; 
        int MCTruth2 = -1; 

        MCTemplatesTree->SetBranchAddress("sampleID",&sample_ID);
        MCTemplatesTree->SetBranchAddress("weight",&weight_template);
        MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
        MCTemplatesTree->SetBranchAddress("category",&category_template);
        MCTemplatesTree->SetBranchAddress("Nr",&Nr_template);
        MCTemplatesTree->SetBranchAddress("Nj",&Nj_template);
        MCTemplatesTree->SetBranchAddress("Nl",&Nl_template);
        MCTemplatesTree->SetBranchAddress("Nl_loose",&Nlloose_template);
        MCTemplatesTree->SetBranchAddress("Nr_loose",&Nrloose_template);
        MCTemplatesTree->SetBranchAddress("MCTruth1",&MCTruth1);
        MCTemplatesTree->SetBranchAddress("MCTruth2",&MCTruth2);

        for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
            MCTemplatesTree->GetEntry(entry);
            if ( weight_template == -1 ) continue;
            if ( Mass_template == -1 ) continue;

            // shape for cate0 in control region of 1photon + 1lepton +>= 4jets
            if(category_template == 2){
                h_TemplatesExtractedFromData[0]->Fill(Mass_template);
                h_TemplatesExtractedFromDataTemp[0]->Fill(Mass_template);
            }
            // shape for cate1 in control region of 1photon + 1leptons +>= 5jets
            if(category_template == 0){
                if(usingCate1insteadofCate3)
                    h_TemplatesExtractedFromData[1]->Fill(Mass_template);
                h_TemplatesExtractedFromDataTemp[1]->Fill(Mass_template);
            }
            // shape for cate2 in control region of 0photon + 3leptons +>= 4jets
            if(category_template == 3){
                h_TemplatesExtractedFromDataTemp[2]->Fill(Mass_template);
            }
            // shape for cate3 in control region of 2leptons +>= 5jets
            if(category_template == 4){
                if(!usingCate1insteadofCate3)
                    h_TemplatesExtractedFromData[1]->Fill(Mass_template);
                h_TemplatesExtractedFromDataTemp[3]->Fill(Mass_template);
            }
            // shape for cate4 in control region of 1 lep + >= 6jets
            if(category_template == 1){
                h_TemplatesExtractedFromData[2]->Fill(Mass_template);
                h_TemplatesExtractedFromDataTemp[4]->Fill(Mass_template);
            }
            if(category_template == 10){
                h_TemplatesExtractedFromData[3]->Fill(Mass_template);
            }
        }
    }
    h_TemplatesExtractedFromData[0]->Scale(1./h_TemplatesExtractedFromData[0]->Integral());

    // trim cate2
    if(true){
	sprintf(buffer,"fit_data_cate_%i_tmp",2);
	std::cout<<"[first fit]"<<std::endl;
	h_TemplatesExtractedFromDataTemp[2]->Fit(buffer,"IE","",350,1500);
	h_TemplatesExtractedFromDataTemp[2]->Fit(buffer,"IEMLL","",350,1500);
	double values__[2] = {
	    f_TemplatesExtractedFromDataTemp[2]->Integral(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(6),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(6) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(6)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(6),
	    f_TemplatesExtractedFromDataTemp[2]->IntegralError(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(6),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(6) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(6)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(6)
	};
	std::cout<<"[second fit]"<<std::endl;
	h_TemplatesExtractedFromDataTemp[2]->Fit(buffer,"IE","",500,1500);
	h_TemplatesExtractedFromDataTemp[2]->Fit(buffer,"IEMLL","",500,1500);

	h_TemplatesExtractedFromDataTemp[2]->SetBinContent(7,
		f_TemplatesExtractedFromDataTemp[2]->Integral(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(7),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(7) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(7)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(7));
	h_TemplatesExtractedFromDataTemp[2]->SetBinError(7,
		f_TemplatesExtractedFromDataTemp[2]->IntegralError(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(7),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(7) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(7)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(7));


	//h_TemplatesExtractedFromDataTemp[2]->Fit(buffer,"MEI","",500,1500);
	h_TemplatesExtractedFromDataTemp[2]->SetBinContent(8,
		f_TemplatesExtractedFromDataTemp[2]->Integral(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(8),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(8) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(8)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(8));
	h_TemplatesExtractedFromDataTemp[2]->SetBinError(8,
		f_TemplatesExtractedFromDataTemp[2]->IntegralError(
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(8),
		    h_TemplatesExtractedFromDataTemp[2]->GetBinLowEdge(8) +
		    h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(8)
		    )/h_TemplatesExtractedFromDataTemp[2]->GetBinWidth(8));

	h_TemplatesExtractedFromDataTemp[2]->SetBinContent(6,values__[0]);
	h_TemplatesExtractedFromDataTemp[2]->SetBinError(6,values__[1]);
    }

    for(int i=0;i<5;i++){
        h_TemplatesExtractedFromDataTemp[i]->Scale(1./h_TemplatesExtractedFromDataTemp[i]->Integral());
    }
    //h_TemplatesExtractedFromData[0]->Add(h_TemplatesExtractedFromDataTemp[1]);
    //h_TemplatesExtractedFromData[0]->Add(h_TemplatesExtractedFromDataTemp[2]);
    //h_TemplatesExtractedFromData[0]->Add(h_TemplatesExtractedFromDataTemp[3]);

    double mergebincontent = 0.;
    double mergebincontentError = 0.;
    int Nmergebincontent = 0.;
    for(int i=1;i<=h_TemplatesExtractedFromData[0]->GetNbinsX();i++){
        mergebincontent = 0.;
        mergebincontentError = 0.;
        Nmergebincontent = 0.;
        for(int j=0;j<4;j++){
            if(h_TemplatesExtractedFromDataTemp[j]->GetBinContent(i) > 0){
                Nmergebincontent++;
                mergebincontent += h_TemplatesExtractedFromDataTemp[j]->GetBinContent(i);
                mergebincontentError += pow(h_TemplatesExtractedFromDataTemp[j]->GetBinError(i),2);
            }else{
                //if(false)
                if(i==6||i==7||i==8){
               // if(i==6||i==7){
                //if(i==6||i==8){
                    sprintf(buffer,"fit_data_cate_%i_tmp",j);
                    //h_TemplatesExtractedFromDataTemp[j]->Fit(buffer,"LL,M,R","",200,1550);
                    //h_TemplatesExtractedFromDataTemp[j]->Fit(buffer,"","",500,1550);
                    //h_TemplatesExtractedFromDataTemp[j]->Fit(buffer,"","",500,950);
                    if(i==6){
                        h_TemplatesExtractedFromDataTemp[j]->Fit(buffer,"","",350,1500);
                    }else{
                        h_TemplatesExtractedFromDataTemp[j]->Fit(buffer,"","",500,1500);
                    }
                    f_TemplatesExtractedFromDataTemp[j]->GetParameters(fitparas[j]);
                    Nmergebincontent++;
                    mergebincontent += 
                        f_TemplatesExtractedFromDataTemp[j]->Integral(
                            h_TemplatesExtractedFromDataTemp[j]->GetBinLowEdge(i),
                            h_TemplatesExtractedFromDataTemp[j]->GetBinLowEdge(i) + 
                            h_TemplatesExtractedFromDataTemp[j]->GetBinWidth(i)
                            //fitparas[j]
                            )/h_TemplatesExtractedFromDataTemp[j]->GetBinWidth(i);
                    mergebincontentError += 
                        pow(
                                f_TemplatesExtractedFromDataTemp[j]->IntegralError(
                                    h_TemplatesExtractedFromDataTemp[j]->GetBinLowEdge(i),
                                    h_TemplatesExtractedFromDataTemp[j]->GetBinLowEdge(i) + 
                                    h_TemplatesExtractedFromDataTemp[j]->GetBinWidth(i)
                                    //fitparas[j]
                                    )/h_TemplatesExtractedFromDataTemp[j]->GetBinWidth(i),
                                2);
                }
            }
        }
        printf("[DEBUG MERGE] %i bin : %f / %i (%f)\n", i, mergebincontent/Nmergebincontent, Nmergebincontent, 
                pow(mergebincontentError,0.5)/Nmergebincontent);
        if(Nmergebincontent!=0){
            h_TemplatesExtractedFromData[0]->SetBinContent(i, mergebincontent/Nmergebincontent);
            h_TemplatesExtractedFromData[0]->SetBinError(i, pow(mergebincontentError,0.5)/Nmergebincontent);
        }else{
            h_TemplatesExtractedFromData[0]->SetBinContent(i, 0.);
        }

    }
    for(int i=0;i<5;i++){
        h_TemplatesExtractedFromDataTemp[i]->Write();
        if(i!=4)
            h_TemplatesExtractedFromData[i]->Write();
    }

    // For control-region MC shape
    if(usingExtrapolateForCRMCShape){
        for(int j=0;j<5;j++){
            h_MergeTemplates[j]->Write();
            h_MergeTemplates[j]->Scale(1./h_MergeTemplates[j]->Integral());
            h_MergeTemplatesWOSysUnc[j] = (TH1F*) h_MergeTemplates[j]->Clone();
            if(AppliedSysOnShape)
                for(int ibin = 1;ibin <= h_MergeTemplates[j]->GetXaxis()->GetNbins() ; ibin++){
                    float unc_ = h_MergeTemplates[j]->GetBinError(ibin);    // absolute error
                    unc_ = unc_*unc_ +
                        pow( h_MergeTemplatesUnc[j]->GetBinContent(ibin) * h_MergeTemplates[j]->GetBinContent(ibin) , 2);
                    h_MergeTemplates[j]->SetBinError(ibin, sqrt(unc_));
                }
            h_MergeTemplatesPure[j] = (TH1F*) h_MergeTemplates[j]->Clone();

            sprintf(buffer,"%sPure",h_MergeTemplates[j]->GetName());
            h_MergeTemplatesPure[j]->SetName(buffer);
            h_MergeTemplatesPure[j]->Write();

        }
        for(int i=1;i<=h_MergeTemplates[0]->GetNbinsX();i++){
            mergebincontent = 0.;
            mergebincontentError = 0.;
            Nmergebincontent = 0.;
            for(int j=0;j<4;j++){
                if(i>=4&&h_MergeTemplates[j]->GetBinContent(i)!=0&&
                        (h_MergeTemplates[j]->GetBinContent(i)/h_MergeTemplates[j]->GetBinContent(i-1) < 0.20||
                         h_MergeTemplates[j]->GetBinContent(i)/h_MergeTemplates[j]->GetBinContent(i-1) > 1.0 )){
                    sprintf(buffer,"fit_data_cate_%i_tmp",j);
                    //h_MergeTemplates[j]->Fit(buffer,"NO","",350,1550);
                    //h_MergeTemplates[j]->Fit(buffer,"NO","",500,950);
                    //h_MergeTemplates[j]->Fit(buffer,"NO","",500,950);
                    if(i>=7){
                        h_MergeTemplates[j]->Fit(buffer,"NO","",650,950);
                    }else{
                        h_MergeTemplates[j]->Fit(buffer,"NO","",200,500);
                    }
                    Nmergebincontent++;

                    double bincontent_ = 
                        f_TemplatesExtractedFromDataTemp[j]->Integral(
                                h_MergeTemplates[j]->GetBinLowEdge(i),
                                h_MergeTemplates[j]->GetBinLowEdge(i) + 
                                h_MergeTemplates[j]->GetBinWidth(i)
                                )/h_MergeTemplates[j]->GetBinWidth(i);
                    mergebincontent += bincontent_;
                    double bincontentError_ = 
                        pow(
                                f_TemplatesExtractedFromDataTemp[j]->IntegralError(
                                    h_MergeTemplates[j]->GetBinLowEdge(i),
                                    h_MergeTemplates[j]->GetBinLowEdge(i) + 
                                    h_MergeTemplates[j]->GetBinWidth(i)
                                    )/h_MergeTemplates[j]->GetBinWidth(i),
                                2);
                    mergebincontentError += bincontentError_;
                    // Saving the fitted value into an original template will change the next fitting result. 
                    //h_MergeTemplates[j]->SetBinContent(i, bincontent_);   
                    //h_MergeTemplates[j]->SetBinError(i, sqrt(bincontentError_));
                }else if(h_MergeTemplates[j]->GetBinContent(i) > 0){
                    Nmergebincontent++;
                    mergebincontent += h_MergeTemplates[j]->GetBinContent(i);
                    mergebincontentError += pow(h_MergeTemplates[j]->GetBinError(i),2);
                }
            }
            printf("[DEBUG MERGE] %i bin : %f / %i (%f)\n", i, mergebincontent/Nmergebincontent, Nmergebincontent, 
                    pow(mergebincontentError,0.5)/Nmergebincontent);
            if(Nmergebincontent!=0){
                h_MergeTemplates[0]->SetBinContent(i, mergebincontent/Nmergebincontent);
                h_MergeTemplates[0]->SetBinError(i, pow(mergebincontentError,0.5)/Nmergebincontent);
            //}else{
            //    h_MergeTemplates[0]->SetBinContent(i, 0.);
            }
        }
    }
    for(int j=0;j<5;j++){
        sprintf(buffer,"%sAfterMerging",h_MergeTemplates[j]->GetName());
        h_MergeTemplates[j]->SetName(buffer);
        h_MergeTemplates[j]->Write();
    }
    // combine : cate0+2  +  cate1+3
    //h_TemplatesExtractedFromData[0]->Add(h_TemplatesExtractedFromDataTemp[2]);
    //h_TemplatesExtractedFromData[0]->Scale(1./h_TemplatesExtractedFromData[0]->Integral());

    //h_TemplatesExtractedFromDataTemp[1]->Add(h_TemplatesExtractedFromDataTemp[3]);
    //h_TemplatesExtractedFromDataTemp[1]->Scale(h_TemplatesExtractedFromDataTemp[3]->Integral());

    //h_TemplatesExtractedFromData[0]->Add(h_TemplatesExtractedFromDataTemp[1]);
    //h_TemplatesExtractedFromData[0]->Scale(1./h_TemplatesExtractedFromData[0]->Integral());

    // Signal templates in SR
    for(int sidx_=0;sidx_<=TGluonTgamma1500GeV-Tgamma500GeV;sidx_++){

	printf("signal(%s) in SR\n",SAMPLE[Tgamma500GeV+sidx_].tag);

        if(RunStatus_>=UncQsquarePlus){
            sprintf(buffer,"tree_%s%s",
                    SAMPLE[Tgamma500GeV+sidx_].tag,RunStatusNames[UncQsquare].c_str());
        }else{
            sprintf(buffer,"tree_%s%s",
                    SAMPLE[Tgamma500GeV+sidx_].tag,RunStatusNames[RunStatus_].c_str());
        }
        TTree *MCTemplatesTree = (TTree*)file->Get(buffer);
	if(sidx_< Tgamma1500GeV-Tgamma500GeV+1){
	    sprintf(buffer,"signal(%s) in SR",SAMPLE[Tgamma500GeV+sidx_].tag);
	    h_SignalTemplateInSR[sidx_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);

	    h_SignalTemplateInSR[sidx_]->Sumw2();
	    h_SignalTemplateInSR[sidx_]->SetLineWidth(2);
	}

        int sample_ID=-1;
        int Nr_template=-1;
        int Nj_template=-1;
        int Nl_template=-1;
        int Nlloose_template=-1;
        int Nrloose_template=-1;
        float weight_template=-1;
        float Mass_template=-1;
        float category_template=-1;
        int MCTruth1 = -1; 
        int MCTruth2 = -1; 

        MCTemplatesTree->SetBranchAddress("sampleID",&sample_ID);
        MCTemplatesTree->SetBranchAddress("weight",&weight_template);
        MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
        MCTemplatesTree->SetBranchAddress("category",&category_template);
        MCTemplatesTree->SetBranchAddress("Nr",&Nr_template);
        MCTemplatesTree->SetBranchAddress("Nj",&Nj_template);
        MCTemplatesTree->SetBranchAddress("Nl",&Nl_template);
        MCTemplatesTree->SetBranchAddress("Nl_loose",&Nlloose_template);
        MCTemplatesTree->SetBranchAddress("Nr_loose",&Nrloose_template);
        MCTemplatesTree->SetBranchAddress("MCTruth1",&MCTruth1);
        MCTemplatesTree->SetBranchAddress("MCTruth2",&MCTruth2);


        for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
            MCTemplatesTree->GetEntry(entry);
            if ( weight_template == -1 ) continue;
            if ( Mass_template == -1 ) continue;
            if(category_template == 10){
		if(sidx_< Tgamma1500GeV-Tgamma500GeV +1 ){
		    h_SignalTemplateInSR[sidx_]->Fill(Mass_template,weight_template*BR*BR);
		}else if(sidx_ >= Tgamma1500GeV-Tgamma500GeV +1 && sidx_< TGluon1500GeV-Tgamma500GeV +1){
		    // very low statistics will cause spike, and contribution in BR(10%) is 1/100 of TgTg. Thus, ignore.
		    //h_SignalTemplateInSR[sidx_%(Tgamma1500GeV-Tgamma500GeV +1)]->Fill(Mass_template,weight_template*(1.-BR)*(1.-BR));
		}else if(sidx_ >= TGluon1500GeV-Tgamma500GeV +1 && sidx_< TGluonTgamma1500GeV-Tgamma500GeV +1){
		    h_SignalTemplateInSR[sidx_%(Tgamma1500GeV-Tgamma500GeV +1)]->Fill(Mass_template,weight_template*2.0*BR*(1.-BR));
		}
            }
        }
    }
}

void uncOnDifferentPhysicsProcesses(){
    output_->cd();
    for(int j=0; j<NCate_; j++){
        sprintf(buffer,"c_UncOnDPP%i",j);
        canvasUncOnDPP[j] = new TCanvas(buffer,"",640,640);
        canvasUncOnDPP[j]->cd();
        if(UsingLogyScaleInMass)
            canvasUncOnDPP[j]->cd()->SetLogy(1);
        int Ncounters_ = 0;
        for(int i=0; i<samples_order_size-_sample_start; i++){
            if(!h_templates[i][j]) continue;
            if(h_templates[i][j]->GetEntries()<statistics) continue;

            h_templatesunc[i][j] = (TH1F*) h_templates[i][j]->Clone(); 
            h_templatesunc[i][j]->Divide(h_MergeTemplates[j]);

            if(Ncounters_==0){
                h_templatesunc[i][j]->GetYaxis()->SetRangeUser(-2,5);
                h_templatesunc[i][j]->SetTitle(CateTitles[j].c_str());
                h_templatesunc[i][j]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                h_templatesunc[i][j]->Draw();
            }else{
                h_templatesunc[i][j]->Draw("same");
            }
            legend_nm[j]->Draw();
            Ncounters_++;
        }
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/UncOnDPP_cate%i.pdf",j);
	else
	    sprintf(buffer,"FittingOutput/UncOnDPP_cate%i_%1.1f.pdf",j,BR);
        canvasUncOnDPP[j]->SaveAs(buffer);
    }

    TH1F *h_TemplatesCateUnc[5];
    for(int j=0; j<NCate_; j++){
        sprintf(buffer,"h_TemplatesCateUnc-%i",j);
        h_TemplatesCateUnc[j] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
    }
    for(int j=0; j<NCate_; j++){
        float uncesUp[numberbins];
        float uncesDown[numberbins];
        for(int iu=0;iu<numberbins;iu++){
            uncesUp[iu] = 1.;
            uncesDown[iu] = 1.;
        }
        for(int i=0; i<samples_order_size-_sample_start; i++){
            if(!h_templatesunc[i][j]) continue;
            for(int iu=1;iu<=numberbins;iu++){
                if(h_templatesunc[i][j]->GetBinContent(iu) > uncesUp[iu-1] ) 
                    uncesUp[iu-1] = h_templatesunc[i][j]->GetBinContent(iu);
                if(h_templatesunc[i][j]->GetBinContent(iu) < uncesDown[iu-1] && h_templatesunc[i][j]->GetBinContent(iu)!=0) 
                    uncesDown[iu-1] = h_templatesunc[i][j]->GetBinContent(iu);
            }
        }

        for(int iu=1;iu<=numberbins;iu++){
            if(fabs(uncesUp[iu-1]-1)>fabs(uncesDown[iu-1]-1)){
                //std::cout<<"cate "<<j <<"uncesUp["<<iu<<"] : "<<uncesUp[iu-1]<<std::endl;
                h_TemplatesCateUnc[j]->SetBinContent(iu,uncesUp[iu-1]);
            }else{
                //std::cout<<"cate "<<j <<"uncesDown["<<iu<<"] : "<<uncesDown[iu-1]<<std::endl;
                h_TemplatesCateUnc[j]->SetBinContent(iu,uncesDown[iu-1]);
            }
        }
    }

    string mode_s[3] = {"Up","Down","Normal"};
    string mode_t[3] = {"plus","minus","normal"};
    string cate_s[3] = {"cate0","cate3","cate4"};
    int unc_cannel[3] = {0,3,4};
    if(usingCate1insteadofCate3){
        cate_s[1] = "cate1";
        unc_cannel[1] = 1;
    }

    TCanvas *UncOnShape[3];
    for(int icate = 0;icate<3;icate++){
        sprintf(buffer,"UncOnShapeCate%i",unc_cannel[icate]);
        UncOnShape[icate] = new TCanvas(buffer,"",640,640);
    }
    TH1F *hUncOnShape[3][3];    // [3 cates] [up, down, norm]

    ///*
    for(int imode = 0;imode<3;imode++){
        for(int icate = 0;icate<3;icate++){
            for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
                // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
                float total_data = h_TemplatesExtractedFromData[3]->Integral();
                TH1F *cate0_dummy = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
                if(icate==0)
                    Scale(cate0_dummy, cate0_dummy, h_TemplatesCateUnc[unc_cannel[icate]], imode);

                TH1F *cate3_dummy = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
                if(icate==1)
                    Scale(cate3_dummy, cate3_dummy, h_TemplatesCateUnc[unc_cannel[icate]], imode);

                TH1F *cate4_dummy = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
                if(icate==2)
                    Scale(cate4_dummy, cate4_dummy, h_TemplatesCateUnc[unc_cannel[icate]], imode);

                float fracs[3] = {0,0,0};
                FitOnlyForDataWithUncShape(sidx_, cate0_dummy, cate3_dummy, cate4_dummy, fracs);

                //cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                //cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data* (1. - fracs[0])* fracs[1]);
                //cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data* (1.-fracs[0])*(1.-fracs[1])*fracs[2]);
                cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data*fracs[1]);
                cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracs[2]);

                TH1F *cates_dummy[3];
                cates_dummy[0] = (TH1F*) cate0_dummy->Clone();
                cates_dummy[1] = (TH1F*) cate3_dummy->Clone();
                cates_dummy[2] = (TH1F*) cate4_dummy->Clone();

                if (sidx_==Tgamma950GeV-Tgamma500GeV){
                    // Up and Down, Normal
                    hUncOnShape[icate][imode] = (TH1F*) cates_dummy[icate]->Clone();
                }

                sprintf(buffer,"pdf_bkg_%s_DifferentPhysicsProcesses%s%s",
                        SAMPLE[Tgamma500GeV+sidx_].tag,
                        cate_s[icate].c_str(),
                        mode_s[imode].c_str());
                if(IsForThetaOrHiggsCombined)
                    sprintf(buffer,"%s__pdf_bkg__DifferentPhysicsProcesses%s__%s",
                            SAMPLE[Tgamma500GeV+sidx_].tag,
                            cate_s[icate].c_str(),
                            mode_t[imode].c_str());
                TH1F *htotal_bg = (TH1F*) cate0_dummy->Clone();
                htotal_bg->SetName(buffer);
                htotal_bg->Add(cate3_dummy);
                htotal_bg->Add(cate4_dummy);
                if(!(IsForThetaOrHiggsCombined&&imode==2)&&(!IsForThetaOrHiggsCombined))
                    if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV))
                        htotal_bg->Write();
                delete htotal_bg;
                //delete []cates_dummy;

            }
        }
    }
    //*/
    for(int icate = 0;icate<3;icate++){
        UncOnShape[icate]->cd();

        TLegend *Unclegend_nm = new TLegend(0.577,0.305+0.3,0.9,0.87);
        Unclegend_nm->SetBorderSize(0);
        Unclegend_nm->SetFillColor(0);
        Unclegend_nm->SetFillStyle(0);
        Unclegend_nm->SetNColumns(1);
        Unclegend_nm->SetTextSize(0.04);
        Unclegend_nm->SetTextSizePixels(25);
        for(int imode = 0;imode<3;imode++){
            sprintf(buffer,"%s-%s",
                    cate_s[icate].c_str(),
                    mode_s[imode].c_str());
            hUncOnShape[icate][imode]->SetTitle( "Unc On DifferentPhysicsProcesses" );
            hUncOnShape[icate][imode]->SetLineColor( colors[imode] );
            Unclegend_nm->AddEntry(hUncOnShape[icate][imode],buffer,"l");
            hUncOnShape[icate][imode]->Scale(1./hUncOnShape[icate][imode]->Integral());
            if(imode==0){
                hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0,1);
                if(UsingLogyScaleInMass)
                    hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0.0001,pow(10,2.0*log10(hUncOnShape[icate][imode]->GetMaximum())));
                hUncOnShape[icate][imode]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                hUncOnShape[icate][imode]->GetYaxis()->SetTitle("Arbitrary Unit");
                hUncOnShape[icate][imode]->Draw("hist");
            }else{
                hUncOnShape[icate][imode]->Draw("same,hist");
            }
        }
        if(UsingLogyScaleInMass)
            UncOnShape[icate]->cd()->SetLogy(1);
        Unclegend_nm->Draw();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_DifferentPhysicsProcesses%s.pdf",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_DifferentPhysicsProcesses%s_%1.1f.pdf",cate_s[icate].c_str(),BR);

        UncOnShape[icate]->SaveAs(buffer);
        delete Unclegend_nm;
    }

}

void patchForSR(){
    int cate_SR = 5;    // 5 : loose region x ISO efficiency ; 10 : signal region
    for(int i=0; i<Ntmp; i++){
        if(RunStatus_>=UncQsquarePlus){
            sprintf(buffer,"tree_%s%s",
                    patchNames[i].c_str(),RunStatusNames[UncQsquare].c_str());

            if(RunStatus_==UncQsquarePlus){
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncQsquareMinus){
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncMatchingPlus){
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
            if(RunStatus_==UncMatchingMinus){
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
                if(!strcmp(patchNames[i].c_str(), 
                    SAMPLE[DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1].tag))
                    sprintf(buffer,"tree_%s%s",
                            SAMPLE[DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1].tag,
                            RunStatusNames[UncQsquare].c_str());
            }
        }else{
            sprintf(buffer,"tree_%s%s",
                    patchNames[i].c_str(), RunStatusNames[RunStatus_].c_str());
        }

        TTree *MCTemplatesTree;
        MCTemplatesTree = (TTree*)filepatch->Get(buffer);

        for(int idxCate_ =0 ;idxCate_<NCate_;idxCate_++){
            sprintf(buffer,"%s_cate_%iinSR_patch",
            //sprintf(buffer,"%s_cate-%i in SR",
                    patchNames[i].c_str(),
                    idxCate_);
            h_templatesPatchInSR[i][idxCate_] = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
            h_templatesPatchInSR[i][idxCate_]->Sumw2();
            h_templatesPatchInSR[i][idxCate_]->SetLineWidth(2);
        }

        int sample_ID=-1;
        int Nr_template=-1;
        int Nj_template=-1;
        int Nl_template=-1;
        int Nlloose_template=-1;
        int Nrloose_template=-1;
        float weight_template=-1;
        float Mass_template=-1;
        float category_template=-1;
        int MCTruth1 = -1; 
        int MCTruth2 = -1; 
        int MCTruth3 = -1; 

        MCTemplatesTree->SetBranchAddress("sampleID",&sample_ID);
        MCTemplatesTree->SetBranchAddress("weight",&weight_template);
        MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
        MCTemplatesTree->SetBranchAddress("category",&category_template);
        MCTemplatesTree->SetBranchAddress("Nr",&Nr_template);
        MCTemplatesTree->SetBranchAddress("Nj",&Nj_template);
        MCTemplatesTree->SetBranchAddress("Nl",&Nl_template);
        MCTemplatesTree->SetBranchAddress("Nl_loose",&Nlloose_template);
        MCTemplatesTree->SetBranchAddress("Nr_loose",&Nrloose_template);
        MCTemplatesTree->SetBranchAddress("MCTruth1",&MCTruth1);
        MCTemplatesTree->SetBranchAddress("MCTruth2",&MCTruth2);
        MCTemplatesTree->SetBranchAddress("MCTruth3",&MCTruth3);


        float testSum = 0.;
        for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
            MCTemplatesTree->GetEntry(entry);
            //if ( weight_template == -1 ) continue;
            if ( weight_template <= 0 ) continue;
            if ( Mass_template == -1 ) continue;
            if(!(MCTruth3==1||MCTruth3==2||MCTruth3==7||MCTruth3==8||MCTruth3==9)) continue;
            //if( (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) ){
            if( ((MCTruth1>=7&&MCTruth1<=9)&&(MCTruth2>=3&&MCTruth2<=5)) || 
                    ((MCTruth2>=7&&MCTruth2<=9)&&(MCTruth1>=3&&MCTruth1<=5)) ){
                if(category_template==cate_SR)
                    testSum+=weight_template;
            }
        }

        for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
            MCTemplatesTree->GetEntry(entry);
            if ( weight_template == -1 ) continue;
            if ( Mass_template == -1 ) continue;
            //weight_template = 1.; // For low statistics

            if(!(MCTruth3==1||MCTruth3==2||MCTruth3==7||MCTruth3==8||MCTruth3==9)) continue;

            if( (MCTruth1%6==0&&(MCTruth2>=7&&MCTruth2<=9)) || (MCTruth2%6==0&&(MCTruth1>=7&&MCTruth1<=9)) ){
                if(category_template==cate_SR)
                    h_templatesPatchInSR[i][0]->Fill(Mass_template,weight_template);
            //    if(category_template==cate_SR&&Mass_template>=800&&Mass_template<950)
            //        printf("[DEBUG MODE] %s : %f (%f) \n",patchNames[i].c_str() ,Mass_template,weight_template);
            }
            if( (MCTruth1%6==0&&(MCTruth2>=3&&MCTruth2<=5)) || (MCTruth2%6==0&&(MCTruth1>=3&&MCTruth1<=5)) ){
                if(category_template==cate_SR)
                    h_templatesPatchInSR[i][1]->Fill(Mass_template,weight_template);
            //    if(category_template==cate_SR&&Mass_template>=950&&Mass_template<1100)
            //        printf("[DEBUG MODE] %s : %f (%f) \n",patchNames[i].c_str() ,Mass_template,weight_template);

            }
            if( (MCTruth2>=7&&MCTruth2<=9) && (MCTruth1>=7&&MCTruth1<=9) ){
                if(category_template==cate_SR)
                    h_templatesPatchInSR[i][2]->Fill(Mass_template,weight_template);
            }
            if( ((MCTruth1>=7&&MCTruth1<=9)&&(MCTruth2>=3&&MCTruth2<=5)) || 
                    ((MCTruth2>=7&&MCTruth2<=9)&&(MCTruth1>=3&&MCTruth1<=5)) ){
                if(category_template==cate_SR)
                    h_templatesPatchInSR[i][3]->Fill(Mass_template,weight_template);
                if(category_template==cate_SR&&Mass_template>=650&&Mass_template<800)
                //if(category_template==cate_SR&&Mass_template>=800&&Mass_template<950)
                    printf("[DEBUG MODE] %s : %f (%f) \n",patchNames[i].c_str() ,Mass_template,weight_template/testSum);
            }
            if( (MCTruth2>=3&&MCTruth2<=5) && (MCTruth1>=3&&MCTruth1<=5) ){
                if(category_template==cate_SR)
                    h_templatesPatchInSR[i][4]->Fill(Mass_template,weight_template);
                //if(category_template==cate_SR&&Mass_template>=800&&Mass_template<950)
                //    printf("[DEBUG MODE] %s : %f (%f) \n",patchNames[i].c_str() ,Mass_template,weight_template/testSum);
                //if(category_template==cate_SR&&Mass_template>=50&&Mass_template<200)
                //    printf("[DEBUG MODE] %s : %f (%f) \n",patchNames[i].c_str() ,Mass_template,weight_template/testSum);
            }
        }

        /*
        for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
            // Not use the below one, because of low statistics
            h_templatesPatchInSR[i][idxCate_]->Scale(1./h_templatesPatchInSR[i][idxCate_]->Integral());   
            std::cout<<"[Patch] "<<patchNames[i]<<" : "<<
                h_templatesPatchInSR[i][idxCate_]->GetEntries()<<
                " ( "<< h_templatesPatchInSR[i][idxCate_]->Integral()/h_templatesPatchInSR[i][idxCate_]->GetEntries()  
                <<" ) in cate "<<idxCate_<<std::endl;
            //if(h_templatesPatchInSR[i][idxCate_]->Integral()>0)
            if(h_templatesPatchInSR[i][idxCate_]->GetEntries()>10)
                h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);  // normalization is excluded
            if(h_templatesPatchInSR[i][idxCate_]->GetEntries()>10&&idxCate_==4)
                printf("[DEBUG CATE4] %s (%f)\n", patchNames[i].c_str(), h_templatesPatchInSR[i][idxCate_]->GetBinContent(1));
        }
        */
    }
        for(int idxCate_ = 0;idxCate_<NCate_;idxCate_++){
            sprintf(buffer,"fit_MC_cate_%i_tmp",idxCate_);
            f_templatesPatchInSR[idxCate_] = new TF1(buffer,"exp([0]+[1]*x)");
            f_templatesPatchInSR[idxCate_]->SetParameters(3.34434e+00,-1.30296e-02);
            for(int i=0;i<=Ntmp;i++) 
                std::cout<<"[Patch Number in Cate"<< idxCate_<<" ] "<<patchNames[i]<<" : "<<
                h_templatesPatchInSR[i][idxCate_]->Integral()<<std::endl;
            float weight_area = 1.;
            // TTbar group
                // cope with no events in tail
            for(int idx_=3;idx_<=7;idx_+=4) // TTJets_MSDecays and ttbar
                for(int ibin=4;ibin<=h_templatesPatchInSR[idx_][idxCate_]->GetNbinsX();ibin++){
                    if(
                            h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(ibin-1) > 4  
                            ||
                            h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(ibin-1) < 0.1  
                            //&& h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(ibin-1)!=0
                            ){
                        for(int jbin=ibin;jbin<=h_templatesPatchInSR[idx_][idxCate_]->GetNbinsX();jbin++){
                            sprintf(buffer,"fit_MC_cate_%i_tmp",idxCate_);
                            if(jbin<8){
                                h_templatesPatchInSR[idx_][idxCate_]->Fit(buffer,"","",200,
                                        h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin-1));
                            }else{
                                h_templatesPatchInSR[idx_][idxCate_]->Fit(buffer,"","",500,
                                        h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin-1));
                            }

                            double content_ = f_templatesPatchInSR[idxCate_]->Integral(
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[idx_][idxCate_]->GetBinWidth(jbin);

                            if(fabs(1-h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(jbin)/content_)>0.5&&
                                    h_templatesPatchInSR[idx_][idxCate_]->GetBinContent(jbin)>0){
                                h_templatesPatchInSR[idx_][idxCate_]->SetBinContent(jbin, content_);

                                h_templatesPatchInSR[idx_][idxCate_]->SetBinError(jbin,
                                        f_templatesPatchInSR[idxCate_]->IntegralError(
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[idx_][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[idx_][idxCate_]->GetBinWidth(jbin));

                            }
                        }
                        break;
                    }
                }
                //if(h_templatesPatchInSR[7][idxCate_]->Integral()!=0) // will be a problem when considering Qscale or Match
                    weight_area = h_templatesPatchInSR[7][idxCate_]->Integral();

                for(int i=1;i<=7;i++) 
                    if(i!=2)        // don't use full leptonical ttbar
                        if(i!=6)    // don't use ttbar 1000-Inf
                            //    if(i!=7)// don't use ttbar (due to lower statistics in the tail)
                        {
                            if((i==3&&idxCate_==3)) continue;// don't use MSDecay ttbar (strange shape)
                            h_templatesPatchInSR[0][idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);   
                        }

                if(h_templatesPatchInSR[0][idxCate_]->Integral()!=0) 
                    h_templatesPatchInSR[0][idxCate_]->Scale(1./h_templatesPatchInSR[0][idxCate_]->Integral());   
                h_templatesPatchInSR[0][idxCate_]->Scale(weight_area);
                weight_area = 1.0;

            // Z group
                weight_area = h_templatesPatchInSR[13][idxCate_]->Integral()+h_templatesPatchInSR[14][idxCate_]->Integral();

                for(int i=8;i<=14;i++) 
                    if(i!=13)
                        if(i!=8)   // DY_PtZ-100
                        h_templatesPatchInSR[13][idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);   
                if(h_templatesPatchInSR[13][idxCate_]->Integral()!=0)
                    h_templatesPatchInSR[13][idxCate_]->Scale(1./h_templatesPatchInSR[13][idxCate_]->Integral());   
                h_templatesPatchInSR[13][idxCate_]->Scale(weight_area);   

                weight_area =1.0;

                // sub-group 4 (h_templatesPatchInSR[15][idxCate_])
                //if(h_templatesPatchInSR[15][idxCate_]->Integral()!=0)
                //    h_templatesPatchInSR[15][idxCate_]->Scale(1./h_templatesPatchInSR[15][idxCate_]->Integral());   
                //h_templatesPatchInSR[8][idxCate_]->Add(h_templatesPatchInSR[15][idxCate_]);   
            // Di-bosons group
                // sub-group 1
                //if(h_templatesPatchInSR[16][idxCate_]->Integral()!=0)
                //    h_templatesPatchInSR[16][idxCate_]->Scale(1./h_templatesPatchInSR[16][idxCate_]->Integral());   
                //// sub-group 2
                //for(int i=18;i<=22;i++) 
                //    h_templatesPatchInSR[17][idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);   
                //if(h_templatesPatchInSR[17][idxCate_]->Integral()!=0)
                //    h_templatesPatchInSR[17][idxCate_]->Scale(1./h_templatesPatchInSR[17][idxCate_]->Integral());   
                /*
            // Higgs group
                for(int i=32;i<=37;i++) 
            // single top group
                for(int i=39;i<=43;i++) 
            // gamma group
                for(int i=45;i<=47;i++) 
            // ttX group
                for(int i=49;i<=52;i++) 
                */
                for(int i=24;i<=52;i++) 
                    if(i!=25)   // Wjet 
                    if(i!=50)   // ttG
                    h_templatesPatchInSR[23][idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);   

                h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[0][idxCate_]);  // TTbar

                // cope with no events in tail
                for(int ibin=4;ibin<=h_templatesPatchInSR[13][idxCate_]->GetNbinsX();ibin++){
                    if(
                            h_templatesPatchInSR[13][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[13][idxCate_]->GetBinContent(ibin-1) <0.05  ||
                            h_templatesPatchInSR[13][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[13][idxCate_]->GetBinContent(ibin-1) > 4  
                            ){
                        sprintf(buffer,"fit_MC_cate_%i_tmp",idxCate_);
                        h_templatesPatchInSR[13][idxCate_]->Fit(buffer,"","",200,
                                h_templatesPatchInSR[13][idxCate_]->GetBinLowEdge(ibin-1));
                        for(int jbin=ibin;jbin<=h_templatesPatchInSR[13][idxCate_]->GetNbinsX();jbin++){

                            double content_ = f_templatesPatchInSR[idxCate_]->Integral(
                                            h_templatesPatchInSR[13][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[13][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[13][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[13][idxCate_]->GetBinWidth(jbin);


                            if(
                                    (fabs(1-h_templatesPatchInSR[13][idxCate_]->GetBinContent(jbin)/content_)>0.5&&
                                     h_templatesPatchInSR[13][idxCate_]->GetBinContent(jbin)>0)
                                    ){
                                h_templatesPatchInSR[13][idxCate_]->SetBinContent(jbin,content_ );

                                h_templatesPatchInSR[13][idxCate_]->SetBinError(jbin,
                                        f_templatesPatchInSR[idxCate_]->IntegralError(
                                            h_templatesPatchInSR[13][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[13][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[13][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[13][idxCate_]->GetBinWidth(jbin));

                            }
                        }
                        break;
                    }
                }
                //h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[13][idxCate_]); // Z-boson removed on 09262014 lowstat
                //h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[15][idxCate_]);  // ZG

                // cope with no events in tail
                if(false)
                for(int ibin=4;ibin<=h_templatesPatchInSR[23][idxCate_]->GetNbinsX();ibin++){
                    if(
                            h_templatesPatchInSR[23][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[23][idxCate_]->GetBinContent(ibin-1) <0.05  ||
                            h_templatesPatchInSR[23][idxCate_]->GetBinContent(ibin)
                            / h_templatesPatchInSR[23][idxCate_]->GetBinContent(ibin-1) > 4  
                            ){
                        sprintf(buffer,"fit_MC_cate_%i_tmp",idxCate_);
                        h_templatesPatchInSR[23][idxCate_]->Fit(buffer,"","",200,
                                h_templatesPatchInSR[23][idxCate_]->GetBinLowEdge(ibin));
                        for(int jbin=ibin;jbin<=h_templatesPatchInSR[23][idxCate_]->GetNbinsX();jbin++){

                            double content_ = f_templatesPatchInSR[idxCate_]->Integral(
                                            h_templatesPatchInSR[23][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[23][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[23][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[23][idxCate_]->GetBinWidth(jbin);


                            if(
                                    (fabs(1-h_templatesPatchInSR[23][idxCate_]->GetBinContent(jbin)/content_)>0.5&&
                                     h_templatesPatchInSR[23][idxCate_]->GetBinContent(jbin)>0)
                                    ){
                                h_templatesPatchInSR[23][idxCate_]->SetBinContent(jbin,content_ );

                                h_templatesPatchInSR[23][idxCate_]->SetBinError(jbin,
                                        f_templatesPatchInSR[idxCate_]->IntegralError(
                                            h_templatesPatchInSR[23][idxCate_]->GetBinLowEdge(jbin),
                                            h_templatesPatchInSR[23][idxCate_]->GetBinLowEdge(jbin) + 
                                            h_templatesPatchInSR[23][idxCate_]->GetBinWidth(jbin)
                                            )/h_templatesPatchInSR[23][idxCate_]->GetBinWidth(jbin));

                            }
                        }
                        break;
                    }
                }

                h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[23][idxCate_]); // others 

                h_templatesPatchInSR[0][idxCate_]->Write();
                h_templatesPatchInSR[2][idxCate_]->Write();
                h_templatesPatchInSR[13][idxCate_]->Write();
                h_templatesPatchInSR[25][idxCate_]->Write();
                h_templatesPatchInSR[23][idxCate_]->Write();

                // solve the missing middle point 
                for(int ibin=4;ibin<=h_MergeTemplatesInSR[idxCate_]->GetNbinsX();ibin++){
                    sprintf(buffer,"fit_MC_cate_%i_tmp",idxCate_);
                    h_MergeTemplatesInSR[idxCate_]->Fit(buffer,"NO","",200,
                            h_MergeTemplatesInSR[idxCate_]->GetBinLowEdge(ibin));
                    for(int jbin=ibin;jbin<=h_MergeTemplatesInSR[idxCate_]->GetNbinsX();jbin++){
                        double content_ = f_templatesPatchInSR[idxCate_]->Integral(
                                h_MergeTemplatesInSR[idxCate_]->GetBinLowEdge(jbin),
                                h_MergeTemplatesInSR[idxCate_]->GetBinLowEdge(jbin) + 
                                h_MergeTemplatesInSR[idxCate_]->GetBinWidth(jbin)
                                )/h_MergeTemplatesInSR[idxCate_]->GetBinWidth(jbin);
                        if(
                                h_MergeTemplatesInSR[idxCate_]->GetBinContent(jbin)==0 &&
                                h_MergeTemplatesInSR[idxCate_]->GetBinContent(jbin-1)!=0 &&
                                h_MergeTemplatesInSR[idxCate_]->GetBinContent(jbin+1)!=0
                                ){
                            h_MergeTemplatesInSR[idxCate_]->SetBinContent(jbin,content_ );
                            h_MergeTemplatesInSR[idxCate_]->SetBinError(jbin,
                                    f_templatesPatchInSR[idxCate_]->IntegralError(
                                        h_MergeTemplatesInSR[idxCate_]->GetBinLowEdge(jbin),
                                        h_MergeTemplatesInSR[idxCate_]->GetBinLowEdge(jbin) + 
                                        h_MergeTemplatesInSR[idxCate_]->GetBinWidth(jbin)
                                        )/h_MergeTemplatesInSR[idxCate_]->GetBinWidth(jbin));
                        }
                    }
                    break;
                }

                h_MergeTemplatesInSR[idxCate_]->Write();
                h_MergeTemplatesInSR[idxCate_]->Scale(1./h_MergeTemplatesInSR[idxCate_]->Integral());
                if(AppliedSysOnShape)
                    for(int ibin = 1;ibin <= h_MergeTemplatesInSR[idxCate_]->GetXaxis()->GetNbins() ; ibin++){
                        float unc_ = h_MergeTemplatesInSR[idxCate_]->GetBinError(ibin);    // absolute error
                        unc_ = unc_*unc_ +
                            pow( h_MergeTemplatesUncInSR[idxCate_]->GetBinContent(ibin) * 
                                    h_MergeTemplatesInSR[idxCate_]->GetBinContent(ibin),2);
                        h_MergeTemplatesInSR[idxCate_]->SetBinError(ibin, sqrt(unc_));    // need to check
                    }
                printf("[Patch DEBUG] 0(%f), 13(%f), 15(%f), 23(%f), 38(%f), 44(%f), and 48(%f)\n",
                        h_templatesPatchInSR[0][idxCate_]->Integral(),
                        h_templatesPatchInSR[13][idxCate_]->Integral(),
                        h_templatesPatchInSR[15][idxCate_]->Integral(),
                        h_templatesPatchInSR[23][idxCate_]->Integral(),
                        h_templatesPatchInSR[38][idxCate_]->Integral(),
                        h_templatesPatchInSR[44][idxCate_]->Integral(),
                        h_templatesPatchInSR[48][idxCate_]->Integral()
                        );
            //h_templatesPatchInSR[i][idxCate_]->Scale(1./h_templatesPatchInSR[i][idxCate_]->Integral());   

            //std::cout<<"[Patch] "<<patchNames[i]<<" : "<<
            //    h_templatesPatchInSR[i][idxCate_]->GetEntries()<<
            //    " ( "<< h_templatesPatchInSR[i][idxCate_]->Integral()/h_templatesPatchInSR[i][idxCate_]->GetEntries()  
            //    <<" ) in cate "<<idxCate_<<std::endl;
            //if(h_templatesPatchInSR[i][idxCate_]->Integral()>0)

            //if(h_templatesPatchInSR[i][idxCate_]->GetEntries()>10)
            //    h_MergeTemplatesInSR[idxCate_]->Add(h_templatesPatchInSR[i][idxCate_]);  // normalization is excluded

            //if(h_templatesPatchInSR[i][idxCate_]->GetEntries()>10&&idxCate_==4)
            //    printf("[DEBUG CATE4] %s (%f)\n", patchNames[i].c_str(), h_templatesPatchInSR[i][idxCate_]->GetBinContent(1));
        }
}

void uncOnDataMCComparison(){
    output_->cd();
    TH1F *h_TemplatesCateUnc[3];
    TH1F *h_TemplatesCateTmp[3];
    TH1F *h_TemplatesExtractedFromData_tmp[3];
    h_TemplatesCateUnc[0] = (TH1F*) h_MergeTemplates[0]->Clone();
    h_TemplatesCateUnc[1] = (TH1F*) h_MergeTemplates[Index_Cate1or3]->Clone();
    h_TemplatesCateUnc[2] = (TH1F*) h_MergeTemplates[4]->Clone();
    h_TemplatesCateTmp[0] = (TH1F*) h_MergeTemplates[0]->Clone();
    h_TemplatesCateTmp[1] = (TH1F*) h_MergeTemplates[Index_Cate1or3]->Clone();
    h_TemplatesCateTmp[2] = (TH1F*) h_MergeTemplates[4]->Clone();
    for(int i=0;i<3;i++) {
        h_TemplatesCateUnc[i]->Scale(1./h_TemplatesCateUnc[i]->Integral());
        h_TemplatesCateTmp[i]->Scale(1./h_TemplatesCateTmp[i]->Integral());
        h_TemplatesExtractedFromData_tmp[i] = (TH1F*) h_TemplatesExtractedFromData[i]->Clone();
        h_TemplatesExtractedFromData_tmp[i]->Scale(1./h_TemplatesExtractedFromData_tmp[i]->Integral());
        // Pick the largest uncertainty out of the difference and error
        //h_TemplatesCateUnc[i]->Divide(h_TemplatesExtractedFromData_tmp[i]);
        float maxError = 0.;
        float denominator_ = 1.;
        for(int j=1;j<=h_TemplatesCateUnc[i]->GetNbinsX();j++){
            denominator_ = h_TemplatesExtractedFromData_tmp[i]->GetBinContent(j);
            if(denominator_ == 0.) denominator_ = 1.;
            maxError = h_TemplatesCateUnc[i]->GetBinError(j);
            if(h_TemplatesExtractedFromData_tmp[i]->GetBinError(j) > maxError ) 
                maxError = h_TemplatesExtractedFromData_tmp[i]->GetBinError(j);
            if(fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_TemplatesExtractedFromData_tmp[i]->GetBinContent(j)) > maxError )
                maxError = fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_TemplatesExtractedFromData_tmp[i]->GetBinContent(j));

            h_TemplatesCateUnc[i]->SetBinContent(j, 1.+maxError/denominator_);
        }
    }

    string mode_s[3] = {"Up","Down","Normal"};
    string mode_t[3] = {"plus","minus","normal"};
    string cate_s[3] = {"cate0","cate3","cate4"};
    int unc_cannel[3] = {0,3,4};
    if(usingCate1insteadofCate3){
        cate_s[1] = "cate1";
        unc_cannel[1] = 1;
    }

    TCanvas *UncOnShape[3];
    for(int icate = 0;icate<3;icate++){
        sprintf(buffer,"UncOnDataMCComparisonCate%i",unc_cannel[icate]);
        UncOnShape[icate] = new TCanvas(buffer,"",640,640);
    }
    TH1F *hUncOnShape[3][3];    // [3 cates] [up, down, norm]

    for(int imode = 0;imode<3;imode++){
        for(int icate = 0;icate<3;icate++){
            for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
                // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
                float total_data = h_TemplatesExtractedFromData[3]->Integral();
                TH1F *cate0_dummy = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
                if(icate==0)
                    Scale(cate0_dummy, h_TemplatesCateTmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                TH1F *cate3_dummy = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
                if(icate==1)
                    Scale(cate3_dummy, h_TemplatesCateTmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                TH1F *cate4_dummy = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
                if(icate==2)
                    Scale(cate4_dummy, h_TemplatesCateTmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                float fracs[3] = {0,0,0};
                TCanvas *TCanvas_temp = new TCanvas("TCanvas_temp","",640,640);
                TCanvas_temp->cd();
                if(UsingLogyScaleInMass)
                    TCanvas_temp->cd()->SetLogy(1);
                unc_on_icate_mode_tag[0] = icate;
                unc_on_icate_mode_tag[1] = imode;
                FitOnlyForDataWithUncShape(sidx_, cate0_dummy, cate3_dummy, cate4_dummy, fracs);
		if(BR==1.0)
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_DataMCComparison%s%s.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str());
		else
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_DataMCComparison%s%s_%1.1f.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str(),BR);
                TCanvas_temp->SaveAs(buffer);
                delete TCanvas_temp;

                //cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                //cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data* (1. - fracs[0])* fracs[1]);
                //cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data* (1.-fracs[0])*(1.-fracs[1])*fracs[2]);
                cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data*fracs[1]);
                cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracs[2]);

                TH1F *cates_dummy[3];
                cates_dummy[0] = (TH1F*) cate0_dummy->Clone();
                cates_dummy[1] = (TH1F*) cate3_dummy->Clone();
                cates_dummy[2] = (TH1F*) cate4_dummy->Clone();

                if (sidx_==Tgamma950GeV-Tgamma500GeV){
                    // Up and Down
                    hUncOnShape[icate][imode] = (TH1F*) cates_dummy[icate]->Clone();
                }

                sprintf(buffer,"pdf_bkg_%s_DataMCComparison%s%s",
                        SAMPLE[Tgamma500GeV+sidx_].tag,
                        cate_s[icate].c_str(),
                        mode_s[imode].c_str());
                if(IsForThetaOrHiggsCombined)
                    sprintf(buffer,"%s__pdf_bkg__DataMCComparison%s__%s",
                            SAMPLE[Tgamma500GeV+sidx_].tag,
                            cate_s[icate].c_str(),
                            mode_t[imode].c_str());
                TH1F *htotal_bg = (TH1F*) cate0_dummy->Clone();
                htotal_bg->SetName(buffer);
                htotal_bg->Add(cate3_dummy);
                htotal_bg->Add(cate4_dummy);
                if(!(IsForThetaOrHiggsCombined&&imode==2))
                    if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV))
                        htotal_bg->Write();
                //NBgkWithUncOnDATAandMC[sidx_][icate][imode] = htotal_bg->Integral();  
                NBgkWithUncOnDATAandMC[sidx_][icate][imode] = total_data*(fracs[0]+fracs[1]+fracs[2]);  
                delete htotal_bg;

            }
        }
    }
    for(int icate = 0;icate<3;icate++){
        UncOnShape[icate]->cd();

        TLegend *Unclegend_nm = new TLegend(0.577,0.305+0.3,0.9,0.87);
        Unclegend_nm->SetBorderSize(0);
        Unclegend_nm->SetFillColor(0);
        Unclegend_nm->SetFillStyle(0);
        Unclegend_nm->SetNColumns(1);
        Unclegend_nm->SetTextSize(0.04);
        Unclegend_nm->SetTextSizePixels(25);
        for(int imode = 0;imode<3;imode++){
            sprintf(buffer,"%s%s",
                    cate_s[icate].c_str(),
                    mode_s[imode].c_str());
            hUncOnShape[icate][imode]->SetTitle( "Unc On ComparisonOfDataAndMC" );
            hUncOnShape[icate][imode]->SetLineColor( colors[imode] );
            Unclegend_nm->AddEntry(hUncOnShape[icate][imode],buffer,"l");
            hUncOnShape[icate][imode]->Scale(1./hUncOnShape[icate][imode]->Integral());
            if(imode==0){
                hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0,1);
                if(UsingLogyScaleInMass)
                    hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0.0001,pow(10,2.0*log10(hUncOnShape[icate][imode]->GetMaximum())));
                hUncOnShape[icate][imode]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                hUncOnShape[icate][imode]->GetYaxis()->SetTitle("Arbitrary Unit");

                hUncOnShape[icate][imode]->Draw("hist");
            }else{
                hUncOnShape[icate][imode]->Draw("same,hist");
            }
        }
        if(UsingLogyScaleInMass)
            UncOnShape[icate]->cd()->SetLogy(1);
        Unclegend_nm->Draw();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_ComparisonOfDataAndMC%s.pdf",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_ComparisonOfDataAndMC%s_%1.1f.pdf",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_ComparisonOfDataAndMC%s.root",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_ComparisonOfDataAndMC%s_%1.1f.root",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
        delete Unclegend_nm;
    }
}

void Scale(TH1F *hist, TH1F *histtmp/*unc w.r.t.*/, TH1F* hunc, int mode){   // mode : 0(Up) 1(Down) 2(Normal)
    hist->Scale(1./hist->Integral());
    histtmp->Scale(1./histtmp->Integral());
    float maxUnc = 0.;
    for(int i=1;i<=hist->GetNbinsX();i++)
        if(fabs(1. - fabs(hunc->GetBinContent(i))) > maxUnc && hunc->GetBinContent(i)!=0)
            maxUnc = fabs(1. - fabs(hunc->GetBinContent(i)));

    if(mode==0){
        for(int i=1;i<=hist->GetNbinsX();i++){
            if(hunc->GetBinContent(i)!=0){
                if(hist->GetBinContent(i)!=0){
                    hist->SetBinContent(i,hist->GetBinContent(i) * ( 1. + fabs(1. - fabs(hunc->GetBinContent(i)) )) );
                    hist->SetBinError(i,hist->GetBinError(i) * ( 1. + fabs(1. - fabs(hunc->GetBinContent(i)) )) );
                }else{
                    hist->SetBinContent(i,histtmp->GetBinContent(i) * ( 1. + fabs(1. - fabs(hunc->GetBinContent(i)) )) );
                    hist->SetBinError(i,histtmp->GetBinError(i) * ( 1. + fabs(1. - fabs(hunc->GetBinContent(i)) )) );
                }
            }else{
                hist->SetBinContent(i,histtmp->GetBinContent(i) * ( 1. + maxUnc) );
                hist->SetBinError(i,histtmp->GetBinError(i) * ( 1. + maxUnc) );
            }
        }   
    }else if(mode==1){
        for(int i=1;i<=hist->GetNbinsX();i++){
            if(hunc->GetBinContent(i)!=0){
                if(hist->GetBinContent(i)!=0){
                    hist->SetBinContent(i,hist->GetBinContent(i)*TMath::Max((1.-fabs(1. - fabs(hunc->GetBinContent(i)) ) ),0.));
                    hist->SetBinError(i,hist->GetBinError(i)*TMath::Max(( 1. - fabs(1. - fabs(hunc->GetBinContent(i)) )),0.) );
                }else{
                    hist->SetBinContent(i,histtmp->GetBinContent(i)*TMath::Max((1.-fabs(1. - fabs(hunc->GetBinContent(i)) ) ),0.));
                    hist->SetBinError(i,histtmp->GetBinError(i)*TMath::Max(( 1. - fabs(1. - fabs(hunc->GetBinContent(i)) )),0.) );
                }
            }else{
                hist->SetBinContent(i,histtmp->GetBinContent(i)*TMath::Max(( 1. - maxUnc ),0.) );
                hist->SetBinError(i,histtmp->GetBinError(i)*TMath::Max(( 1. - maxUnc),0.) );
            }
        }   
    }       
}

// (old version) consider cate0 (0<->2) and cate3 (1<->3)
// (new version) consider cate0 (0<->2) and cate0 (0<->1) and cate0 (0<->3)
// (latest version, 10/30/2014) consider cate0 (0<->2) and cate0 (0<->1) and cate0 (0<->3)  h_MergeTemplatesPure
void uncOnMerge(){
    const int NuncOnMerge = 4;
    output_->cd();
    TH1F *h_TemplatesCateUnc[NuncOnMerge];
    TH1F *h_MergeTemplates_tmp[NuncOnMerge];
    int index_usingCate1insteadofCate3[2] = {1,3};
    if(usingCate1insteadofCate3){
        index_usingCate1insteadofCate3[0] = 3;
        index_usingCate1insteadofCate3[1] = 1;
    }
    //h_TemplatesCateUnc[0] = (TH1F*) h_MergeTemplates[2]->Clone();
    //h_TemplatesCateUnc[1] = (TH1F*) h_MergeTemplates[index_usingCate1insteadofCate3[0]]->Clone();

    //h_MergeTemplates_tmp[0] = (TH1F*) h_MergeTemplates[0]->Clone();
    //h_MergeTemplates_tmp[1] = (TH1F*) h_MergeTemplates[index_usingCate1insteadofCate3[1]]->Clone();
    for(int i=0;i<NuncOnMerge;i++){
        h_MergeTemplates_tmp[i] = (TH1F*) h_MergeTemplates[0]->Clone();
        h_TemplatesCateUnc[i] = (TH1F*) h_MergeTemplatesPure[i]->Clone();
    }

    for(int i=0;i<NuncOnMerge;i++){
        h_MergeTemplates_tmp[i]->Scale(1./h_MergeTemplates_tmp[i]->Integral());
        h_TemplatesCateUnc[i]->Scale(1./h_TemplatesCateUnc[i]->Integral());
        // Pick the largest uncertainty out of the difference and error
        float maxError = 0.;
        float denominator_ = 1.;
        for(int j=1;j<=h_TemplatesCateUnc[i]->GetNbinsX();j++){
            denominator_ = h_MergeTemplates_tmp[i]->GetBinContent(j);
            if(denominator_ == 0.) denominator_=1.;
            maxError = h_TemplatesCateUnc[i]->GetBinError(j);
            if(h_MergeTemplates_tmp[i]->GetBinError(j) > maxError ) maxError = h_MergeTemplates_tmp[i]->GetBinError(j);
            if(fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_MergeTemplates_tmp[i]->GetBinContent(j)) > maxError )
                maxError = fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_MergeTemplates_tmp[i]->GetBinContent(j));

            h_TemplatesCateUnc[i]->SetBinContent(j, 1.+maxError/denominator_);
        }
    }


    string mode_s[3] = {"Up","Down","Normal"};
    string mode_t[3] = {"plus","minus","normal"};
    string cate_s[NuncOnMerge] = {"cate00","cate01","cate02","cate03"};
    string cate_sL[NuncOnMerge] = {"cateM0","cateM1","cateM2","cateM3"};
    TCanvas *UncOnShape[NuncOnMerge];
    for(int icate = 0;icate<NuncOnMerge;icate++){
        sprintf(buffer,"UncOnMerge%s",cate_s[icate].c_str());
        UncOnShape[icate] = new TCanvas(buffer,"",640,640);
    }
    TH1F *hUncOnShape[NuncOnMerge][3];    // [cate01, cate02, cate03] [up, down, norm]

    float fracsNorminal[3] = {0,0,0};
    for(int icate = 0;icate<NuncOnMerge;icate++){
	for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
	    fracsNorminal[0] = 0.;
	    fracsNorminal[1] = 0.;
	    fracsNorminal[2] = 0.;
	    for(int imode = 2;imode>=0;imode--){
                // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
                float total_data = h_TemplatesExtractedFromData[3]->Integral();
                TH1F *cate0_dummy = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
                //if(icate==0)
                Scale(cate0_dummy, h_MergeTemplates[icate+1]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                TH1F *cate3_dummy = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
                //if(icate==1)
                //    Scale(cate3_dummy, h_TemplatesCateUnc[icate], imode);

                TH1F *cate4_dummy = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();

                float fracs[3] = {0,0,0};
                TCanvas *TCanvas_temp = new TCanvas("TCanvas_temp","",640,640);
                TCanvas_temp->cd();
                if(UsingLogyScaleInMass)
                    TCanvas_temp->cd()->SetLogy(1);
                //unc_on_icate_mode_tag[0] = icate;
                unc_on_icate_mode_tag[0] = 0;
                unc_on_icate_mode_tag[1] = imode;
                FitOnlyForDataWithUncShape(sidx_, cate0_dummy, cate3_dummy, cate4_dummy, fracs);
		if(BR==1.0)
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_Merge%s%s.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str());
		else
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_Merge%s%s_%1.1f.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str(),BR);
                TCanvas_temp->SaveAs(buffer);
                delete TCanvas_temp;

                //cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                //cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data* (1. - fracs[0])* fracs[1]);
                //cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data* (1.-fracs[0])*(1.-fracs[1])*fracs[2]);
                cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data*fracs[1]);
                cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracs[2]);

		if(imode==2){
		    fracsNorminal[0] = fracs[0];
		    fracsNorminal[1] = fracs[1];
		    fracsNorminal[2] = fracs[2];
		}

                TH1F *cates_dummy[2];
                cates_dummy[0] = (TH1F*) cate0_dummy->Clone();
                cates_dummy[1] = (TH1F*) cate3_dummy->Clone();

                if (sidx_==Tgamma950GeV-Tgamma500GeV){
                    // Up and Down, Normal
                    //hUncOnShape[icate][imode] = (TH1F*) cates_dummy[icate]->Clone();
                    hUncOnShape[icate][imode] = (TH1F*) cates_dummy[0]->Clone();
                }

                sprintf(buffer,"pdf_bkg_%s_Merge%s%s",
                        SAMPLE[Tgamma500GeV+sidx_].tag,
                        cate_s[icate].c_str(),
                        mode_s[imode].c_str());
                if(IsForThetaOrHiggsCombined)
                    sprintf(buffer,"%s__pdf_bkg__Merge%s__%s",
                            SAMPLE[Tgamma500GeV+sidx_].tag,
                            cate_s[icate].c_str(),
                            mode_t[imode].c_str());
                TH1F *htotal_bg = (TH1F*) cate0_dummy->Clone();
                htotal_bg->SetName(buffer);
                htotal_bg->Add(cate3_dummy);
                htotal_bg->Add(cate4_dummy);
                if(!(IsForThetaOrHiggsCombined&&imode==2))
                    if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV)){
			if(!IsIndividualBkgShape){
			    htotal_bg->Write();
			}else{
			    cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracsNorminal[0]);
			    cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracsNorminal[2]);
			    sprintf(buffer,"%s__pdf_bkgM__Merge%s__%s",
				    SAMPLE[Tgamma500GeV+sidx_].tag,
				    cate_s[icate].c_str(),
				    mode_t[imode].c_str());
			    cate0_dummy->SetName(buffer);
			    sprintf(buffer,"%s__pdf_bkg4__Merge%s__%s",
				    SAMPLE[Tgamma500GeV+sidx_].tag,
				    cate_s[icate].c_str(),
				    mode_t[imode].c_str());
			    cate4_dummy->SetName(buffer);
			    cate0_dummy->Write();
			    cate4_dummy->Write();
			}
		    }
                //NBgkWithUncOnMerge[sidx_][icate][imode] = htotal_bg->Integral();  
                NBgkWithUncOnMerge[sidx_][icate][imode] = total_data*(fracs[0]+fracs[1]+fracs[2]);  
                delete htotal_bg;

            }
        }
    }

    for(int icate = 0;icate<NuncOnMerge;icate++){
        UncOnShape[icate]->cd();

        TLegend *Unclegend_nm = new TLegend(0.577,0.305+0.3,0.9,0.87);
        Unclegend_nm->SetBorderSize(0);
        Unclegend_nm->SetFillColor(0);
        Unclegend_nm->SetFillStyle(0);
        Unclegend_nm->SetNColumns(1);
        Unclegend_nm->SetTextSize(0.04);
        Unclegend_nm->SetTextSizePixels(25);
        for(int imode = 0;imode<3;imode++){
            sprintf(buffer,"%s%s",
                    //cate_s[icate].c_str(),
                    cate_sL[icate].c_str(),
                    mode_s[imode].c_str());
            hUncOnShape[icate][imode]->SetTitle( "Unc On Reduction" );
            hUncOnShape[icate][imode]->SetLineColor( colors[imode] );
            Unclegend_nm->AddEntry(hUncOnShape[icate][imode],buffer,"l");
            hUncOnShape[icate][imode]->Scale(1./hUncOnShape[icate][imode]->Integral());
            if(imode==0){
                hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0,1);
                if(UsingLogyScaleInMass)
                    hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0.0001,pow(10,2.0*log10(hUncOnShape[icate][imode]->GetMaximum())));
                hUncOnShape[icate][imode]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                hUncOnShape[icate][imode]->GetYaxis()->SetTitle("Arbitrary Unit");

                hUncOnShape[icate][imode]->Draw("hist");
            }else{
                hUncOnShape[icate][imode]->Draw("same,hist");
            }
        }
        if(UsingLogyScaleInMass)
            UncOnShape[icate]->cd()->SetLogy(1);
        Unclegend_nm->Draw();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_Merge%s.pdf",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_Merge%s_%1.1f.pdf",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_Merge%s.root",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_Merge%s_%1.1f.root",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
        delete Unclegend_nm;
    }
}

void uncOnSignalAndControlRegions(){
    output_->cd();
    TH1F *h_TemplatesCateUnc[3];
    int index_usingCate1insteadofCate3 = 3; if(usingCate1insteadofCate3) index_usingCate1insteadofCate3 = 1;
    int index_usingCate1insteadofCate3_SR = 3; 
    if(usingCate1insteadofCate3OnlyForSR||usingCate1insteadofCate3) index_usingCate1insteadofCate3_SR = 1;
    h_TemplatesCateUnc[0] = (TH1F*) h_MergeTemplatesInSR[0]->Clone();
    h_TemplatesCateUnc[1] = (TH1F*) h_MergeTemplatesInSR[index_usingCate1insteadofCate3_SR]->Clone();
    h_TemplatesCateUnc[2] = (TH1F*) h_MergeTemplatesInSR[4]->Clone();
    TH1F *h_MergeTemplates_tmp[3];
    h_MergeTemplates_tmp[0] = (TH1F*) h_MergeTemplates[0]->Clone();
    h_MergeTemplates_tmp[1] = (TH1F*) h_MergeTemplates[index_usingCate1insteadofCate3]->Clone();
    h_MergeTemplates_tmp[2] = (TH1F*) h_MergeTemplates[4]->Clone();

    for(int i=0;i<3;i++){
        h_MergeTemplates_tmp[i]->Scale(1./h_MergeTemplates_tmp[i]->Integral());
        h_TemplatesCateUnc[i]->Scale(1./h_TemplatesCateUnc[i]->Integral());
        // Pick the largest uncertainty out of the difference and error
        float maxError = 0.;
        float denominator_ = 1.;
        for(int j=1;j<=h_TemplatesCateUnc[i]->GetNbinsX();j++){
            denominator_ = h_MergeTemplates_tmp[i]->GetBinContent(j);
            if(denominator_ == 0.) denominator_ =1.;
            maxError = h_TemplatesCateUnc[i]->GetBinError(j);
            if(h_MergeTemplates_tmp[i]->GetBinError(j) > maxError ) maxError = h_MergeTemplates_tmp[i]->GetBinError(j);
            if(fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_MergeTemplates_tmp[i]->GetBinContent(j)) > maxError )
                maxError = fabs(h_TemplatesCateUnc[i]->GetBinContent(j) - h_MergeTemplates_tmp[i]->GetBinContent(j));

            h_TemplatesCateUnc[i]->SetBinContent(j, 1.+maxError/denominator_);
        }
    }

    string mode_s[3] = {"Up","Down","Normal"};
    string mode_t[3] = {"plus","minus","normal"};
    string cate_s[3] = {"cate0","cate3","cate4"};
    int unc_cannel[3] = {0,3,4};
    if(usingCate1insteadofCate3){
        cate_s[1] = "cate1";
        unc_cannel[1] = 1;
    }

    TCanvas *UncOnShape[3];
    for(int icate = 0;icate<3;icate++){
        sprintf(buffer,"UncOnSignalAndControlRegionsCate%i",unc_cannel[icate]);
        UncOnShape[icate] = new TCanvas(buffer,"",640,640);
    }
    TH1F *hUncOnShape[3][3];    // [3 cates] [up, down, norm]

    float fracsNorminal[3] = {0,0,0};
    for(int icate = 0;icate<3;icate++){
	for(int sidx_=0;sidx_<=Tgamma1500GeV-Tgamma500GeV;sidx_++){
	    fracsNorminal[0] = 0.;
	    fracsNorminal[1] = 0.;
	    fracsNorminal[2] = 0.;
	    for(int imode = 2;imode>=0;imode--){
                // model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
                float total_data = h_TemplatesExtractedFromData[3]->Integral();
                TH1F *cate0_dummy = (TH1F*) h_TemplatesExtractedFromData[0]->Clone();
                if(icate==0)
                    Scale(cate0_dummy, h_MergeTemplates_tmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                TH1F *cate3_dummy = (TH1F*) h_TemplatesExtractedFromData[1]->Clone();
                if(icate==1)
                    Scale(cate3_dummy, h_MergeTemplates_tmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                TH1F *cate4_dummy = (TH1F*) h_TemplatesExtractedFromData[2]->Clone();
                if(icate==2)
                    Scale(cate4_dummy, h_MergeTemplates_tmp[icate]/*unc w.r.t.*/, h_TemplatesCateUnc[icate], imode);

                float fracs[3] = {0,0,0};
                TCanvas *TCanvas_temp = new TCanvas("TCanvas_temp","",640,640);
                TCanvas_temp->cd();
                if(UsingLogyScaleInMass)
                    TCanvas_temp->cd()->SetLogy(1);
                unc_on_icate_mode_tag[0] = icate;
                unc_on_icate_mode_tag[1] = imode;
                FitOnlyForDataWithUncShape(sidx_, cate0_dummy, cate3_dummy, cate4_dummy, fracs);
		if(BR==1.0)
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_SignalAndControlRegions%s%s.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str());
		else
		    sprintf(buffer,"FittingOutput/UNC/pdf_bkg_%s_SignalAndControlRegions%s%s_%1.1f.pdf",
			    SAMPLE[Tgamma500GeV+sidx_].tag,
			    cate_s[icate].c_str(),
			    mode_s[imode].c_str(),BR);
                TCanvas_temp->SaveAs(buffer);
                delete TCanvas_temp;

                //cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                //cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data* (1. - fracs[0])* fracs[1]);
                //cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data* (1.-fracs[0])*(1.-fracs[1])*fracs[2]);
                cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracs[0]);
                cate3_dummy->Scale(1./cate3_dummy->Integral()*total_data*fracs[1]);
                cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracs[2]);

		if(imode==2){
		    fracsNorminal[0] = fracs[0];
		    fracsNorminal[1] = fracs[1];
		    fracsNorminal[2] = fracs[2];
		}

                TH1F *cates_dummy[3];
                cates_dummy[0] = (TH1F*) cate0_dummy->Clone();
                cates_dummy[1] = (TH1F*) cate3_dummy->Clone();
                cates_dummy[2] = (TH1F*) cate4_dummy->Clone();

                if (sidx_==Tgamma950GeV-Tgamma500GeV){
                    // Up and Down, Normal
                    hUncOnShape[icate][imode] = (TH1F*) cates_dummy[icate]->Clone();
                }

                sprintf(buffer,"pdf_bkg_%s_SignalAndControlRegions%s%s",
                        SAMPLE[Tgamma500GeV+sidx_].tag,
                        cate_s[icate].c_str(),
                        mode_s[imode].c_str());
                if(IsForThetaOrHiggsCombined)
                    sprintf(buffer,"%s__pdf_bkg__SignalAndControlRegions%s__%s",
                            SAMPLE[Tgamma500GeV+sidx_].tag,
                            cate_s[icate].c_str(),
                            mode_t[imode].c_str());
                TH1F *htotal_bg = (TH1F*) cate0_dummy->Clone();
                htotal_bg->SetName(buffer);
                htotal_bg->Add(cate3_dummy);
                htotal_bg->Add(cate4_dummy);
                //if(!(IsForThetaOrHiggsCombined&&imode==2))
                if(!(IsForThetaOrHiggsCombined&&(imode==2||icate==1)))
                    if(!(IsForThetaOrHiggsCombined&&sidx_<_signal_start-Tgamma500GeV)){
			if(!IsIndividualBkgShape){
			    htotal_bg->Write();
			}else{
			    cate0_dummy->Scale(1./cate0_dummy->Integral()*total_data*fracsNorminal[0]);
			    cate4_dummy->Scale(1./cate4_dummy->Integral()*total_data*fracsNorminal[2]);
			    sprintf(buffer,"%s__pdf_bkgM__SignalAndControlRegions%s__%s",
				    SAMPLE[Tgamma500GeV+sidx_].tag,
				    cate_s[icate].c_str(),
				    mode_t[imode].c_str());
			    cate0_dummy->SetName(buffer);
			    sprintf(buffer,"%s__pdf_bkg4__SignalAndControlRegions%s__%s",
				    SAMPLE[Tgamma500GeV+sidx_].tag,
				    cate_s[icate].c_str(),
				    mode_t[imode].c_str());
			    cate4_dummy->SetName(buffer);
			    cate0_dummy->Write();
			    cate4_dummy->Write();
			}

		    }
                //NBgkWithUncOnSRandCR[sidx_][icate][imode] = htotal_bg->Integral();  
                NBgkWithUncOnSRandCR[sidx_][icate][imode] = total_data*(fracs[0]+fracs[1]+fracs[2]);  
                delete htotal_bg;

            }
        }
    }
    for(int icate = 0;icate<3;icate++){
        UncOnShape[icate]->cd();

        TLegend *Unclegend_nm = new TLegend(0.577,0.305+0.3,0.9,0.87);
        Unclegend_nm->SetBorderSize(0);
        Unclegend_nm->SetFillColor(0);
        Unclegend_nm->SetFillStyle(0);
        Unclegend_nm->SetNColumns(1);
        Unclegend_nm->SetTextSize(0.04);
        Unclegend_nm->SetTextSizePixels(25);
        for(int imode = 0;imode<3;imode++){
            sprintf(buffer,"%s%s",
                    cate_s[icate].c_str(),
                    mode_s[imode].c_str());
	    if(icate==0)
		sprintf(buffer,"%s%s",
			"cateM",
			mode_s[imode].c_str());
            hUncOnShape[icate][imode]->SetTitle( "Unc On SignalAndControlRegions" );
            hUncOnShape[icate][imode]->SetLineColor( colors[imode] );
            Unclegend_nm->AddEntry(hUncOnShape[icate][imode],buffer,"l");
            hUncOnShape[icate][imode]->Scale(1./hUncOnShape[icate][imode]->Integral());
            if(imode==0){
                hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0,1);
                if(UsingLogyScaleInMass)
                    hUncOnShape[icate][imode]->GetYaxis()->SetRangeUser(0.0001,pow(10,2.0*log10(hUncOnShape[icate][imode]->GetMaximum())));
                hUncOnShape[icate][imode]->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
                hUncOnShape[icate][imode]->GetYaxis()->SetTitle("Arbitrary Unit");

                hUncOnShape[icate][imode]->Draw("hist");
            }else{
                hUncOnShape[icate][imode]->Draw("same,hist");
            }
        }

        if(UsingLogyScaleInMass)
            UncOnShape[icate]->cd()->SetLogy(1);
        Unclegend_nm->Draw();
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_SignalAndControlRegions%s.pdf",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_SignalAndControlRegions%s_%1.1f.pdf",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
	if(BR==1.0)
	    sprintf(buffer,"FittingOutput/Unc_On_SignalAndControlRegions%s.root",cate_s[icate].c_str());
	else
	    sprintf(buffer,"FittingOutput/Unc_On_SignalAndControlRegions%s_%1.1f.root",cate_s[icate].c_str(),BR);
        UncOnShape[icate]->SaveAs(buffer);
        delete Unclegend_nm;
    }
}

void FitOnlyForDataWithUncShape(int sidx_, TH1F *hCate0, TH1F *hCate3,TH1F *hCate4, float *fracs){
    char TrimChar_[128];
    float total_data = h_TemplatesExtractedFromData[3]->Integral();
    RooRealVar x("x","Mass",Xregions[0],Xregions[1]) ;
    RooDataHist data("data_obs","data_obs",x,h_TemplatesExtractedFromData[3]);

    RooDataHist dCate0("dCate0","dCate0set with x",x,hCate0);
    sprintf(buffer,"Cate0_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
    RooHistPdf Cate0(buffer,buffer,x,dCate0,0) ;

    char buffer13[256];
    sprintf(buffer,"dCate3");
    sprintf(buffer13,"dCate3set with x");
    if(usingCate1insteadofCate3){
        sprintf(buffer,"dCate1");
        sprintf(buffer13,"dCate1set with x");
    }

    RooDataHist dCate3(buffer,buffer13,x,hCate3);
    sprintf(buffer,"Cate3_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
    if(usingCate1insteadofCate3)
        sprintf(buffer,"Cate1_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
    RooHistPdf Cate3(buffer,buffer,x,dCate3,0) ;

    RooDataHist dCate4("dCate4","dCate4set with x",x,hCate4);
    sprintf(buffer,"Cate4_%s", SAMPLE[Tgamma500GeV+sidx_].tag);
    RooHistPdf Cate4(buffer,buffer,x,dCate4,0) ;

    sprintf(buffer,"%s", SAMPLE[Tgamma500GeV+sidx_].tag);
    RooDataHist dSignal("dSignal","dSignalset with x",x,h_SignalTemplateInSR[sidx_]);
    RooHistPdf Signal(buffer,buffer,x,dSignal,0) ;

    RooRealVar fbkg0("fbkg0","fraction",0.5,0.00,1) ; 
    RooRealVar fbkg3("fbkg3","fraction",0.5,0.00,1) ; 
    //RooRealVar fbkg4("fbkg4","fraction",0.5,0.00,1) ; 
    //RooAddPdf model("model","model",RooArgList(Cate0,Cate3,Cate4,Signal),RooArgList(fbkg0,fbkg3,fbkg4),kTRUE);

    //RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate3,Cate4),RooArgList(fbkg0,fbkg3),kTRUE);
    RooAddPdf bkgmodel("bkgmodel","bkgmodel",RooArgList(Cate0,Cate4),RooArgList(fbkg0),kTRUE);
    // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
    //         Nbkg                                                                   NSign
    float xmin_ = -0.000000000001;
    float xmax_ = 1.5*total_data;
    RooRealVar Nbkgall("Nbkgall","N(bkg)",33.,xmin_,xmax_) ;
    RooRealVar NSign("NSign","N(Sig)",2.,xmin_,xmax_) ;
    //RooRealVar Nbkgall("Nbkgall","yield",0.5,0.00,total_data) ;
    //RooRealVar NSign("NSign","yield",0.5,0.00,total_data) ;
    RooAddPdf model("model","model",RooArgList(bkgmodel,Signal),RooArgList(Nbkgall,NSign));
    //model.fitTo(data,Minimizer("Minuit2", "minimize"));
    model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
    while (fabs((Nbkgall.getVal()+NSign.getVal())/total_data - 1.) > 0.01){
        printf("[In while loop (UncShape)] %1.0f = (%1.4f + %1.4f) for Mass(%s)\n", total_data, Nbkgall.getVal(), NSign.getVal(),
                SAMPLE[Tgamma500GeV+sidx_].tag
              );   
        model.fitTo(data,SumW2Error(kFALSE),Extended(kTRUE),InitialHesse(kTRUE),Hesse(kTRUE),Minos(kTRUE),Save(kTRUE));
    } 

    // Draw
    RooArgSet obs(x);
    RooArgSet* flparams = (RooArgSet*) model.getParameters(obs)->selectByAttrib("Constant",kFALSE) ; 

    RooPlot* xframe = x.frame() ; 
    //xframe->SetTitle("In signal region");
    xframe->SetTitle("");
    float ymax_ = 1.5*h_TemplatesExtractedFromData[3]->GetMaximum();
    data.plotOn(xframe);
    model.plotOn(xframe,LineColor(kBlue-2));
    model.plotOn(xframe,Components(Cate0),LineStyle(kDashed),LineColor(kRed));
    //model.plotOn(xframe,Components(Cate3),LineStyle(kDashed),LineColor(kBlue));
    model.plotOn(xframe,Components(Cate4),LineStyle(kDashed),LineColor(kYellow-2));
    model.plotOn(xframe,Components(Signal),LineStyle(kDashed),LineColor(432+2));
    xframe->GetYaxis()->SetRangeUser(0.,ymax_);
    if(UsingLogyScaleInMass)
        xframe->GetYaxis()->SetRangeUser(0.001,pow(10,2.0*log10(ymax_)));
    xframe->GetYaxis()->SetTitle("Events / (150 GeV/c^{2})");
    xframe->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
    xframe->Draw() ;
    drawCMSLumi();

    // polish the plot here
    float offsetLegend = 0.0;
    if(UsingLogyScaleInMass)
        offsetLegend=0.06;
    TLegend *legend_nmFit = new TLegend(0.5077,0.605+offsetLegend,0.9,0.87+offsetLegend);
    TH1F *htmp[6];
    TrimChar(SAMPLE[Tgamma500GeV+sidx_].tag, TrimChar_);
    sprintf(buffer,"t* (m_{t*}= %1.0d GeV/c^{2})",atoi(TrimChar_));
    string nametemp[6] = {"data","Fit model","Cat M (data)","Cat 3 (data)","Cat 4 (data)",buffer};
    if(usingCate1insteadofCate3) nametemp[3] = "Cat 1(data)";

    if(unc_on_icate_mode_tag[1]!= -1){
        if(unc_on_icate_mode_tag[1]==1)
            nametemp[unc_on_icate_mode_tag[0]+2] += " - unc.";
        if(unc_on_icate_mode_tag[1]==0)
            nametemp[unc_on_icate_mode_tag[0]+2] += " + unc.";
    }

    int Lcolors[6] = {1,kBlue-2,kRed,kBlue,kYellow-2,432+2};
    int Lstyle[6] = {1,1,kDashed,kDashed,kDashed,kDashed};
    for(int i=0;i<6;i++){
        if(i==3) continue;
        sprintf(buffer,"%s",nametemp[i].c_str());
        htmp[i] = new TH1F(buffer,"",10,0,10);
        htmp[i]->SetLineColor(Lcolors[i]);
        htmp[i]->SetLineStyle(Lstyle[i]);
        htmp[i]->SetLineWidth(3);

        if(i==0){
            htmp[i]->SetMarkerStyle(20);
            htmp[i]->SetLineWidth(1);
            legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"lpe");
        }else{
            legend_nmFit->AddEntry(htmp[i],nametemp[i].c_str(),"l");
        }   
    }   

    legend_nmFit->SetBorderSize(0);
    legend_nmFit->SetFillColor(0);
    legend_nmFit->SetFillStyle(0);
    legend_nmFit->SetNColumns(1);
    legend_nmFit->SetTextSize(0.04);
    legend_nmFit->SetTextSizePixels(25);
    legend_nmFit->Draw();

    //old : model = fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)( fbkg4 * Cate4 + (1- fbkg4) * Signal))
    // model = fbkg ( fbkg0 * Cate0 + (1-fbkg0) (fbkg3* Cate3 + (1-fbkg3)* Cate4) ) + (1- fbkg) * Signal
    //         Nbkg                                                                   NSign
    float TextOffsetX = 0.01;
    float TextOffsetY = -0.01;
    TLatex *latexComp;
    if(UsingLogyScaleInMass)
        latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.51+0.42-0.04+TextOffsetY,"Fit model = ");
    else
        latexComp = new TLatex(0.52,0.51,"Fit model = ");
    latexComp->SetTextSize(0.03);
    latexComp->SetTextColor(kBlue-2);
    latexComp->SetNDC();
    latexComp->Draw();

    float fbkg = Nbkgall.getVal()/(Nbkgall.getVal()+NSign.getVal());

    sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 4 +}",
            fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal()));
    if(usingCate1insteadofCate3)
        sprintf(buffer,"#splitline{%1.2fxCat M +}{%1.2fxCat 1 +}",
                fbkg*fbkg0.getVal(), fbkg*(1-fbkg0.getVal())*fbkg3.getVal());

    if(UsingLogyScaleInMass)
        latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.50+0.42-0.04-0.005+TextOffsetY,buffer);
    latexComp->SetTextSize(0.03);
    latexComp->SetTextColor(kBlue-2);
    latexComp->SetNDC();
    latexComp->Draw();

    sprintf(buffer,"%1.2fxSig",
            //        fbkg*(1 - fbkg0.getVal())*(1 - fbkg3.getVal()),
            1.-fbkg
           );
    latexComp = new TLatex(0.63,0.46,buffer);
    if(UsingLogyScaleInMass)
        latexComp = new TLatex(0.63-0.49+0.12+0.04+TextOffsetX,0.46+0.42-0.04-0.007-0.005+TextOffsetY,buffer);
    latexComp->SetTextSize(0.03);
    latexComp->SetTextColor(kBlue-2);
    latexComp->SetNDC();
    latexComp->Draw();

    sprintf(buffer,"#splitline{N(Fitted-Bkg) : %1.1f#pm%1.1f}{N(Fitted-Sig) : %1.1f#pm%1.1f }",
            Nbkgall.getVal(),
            Nbkgall.getError(),
            NSign.getVal(),
            NSign.getError()//, hdata_toy->Integral()
           );

    latexComp = new TLatex(0.52,0.36,buffer);
    if(UsingLogyScaleInMass)
        latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.36+0.43-0.04+0.028+TextOffsetY,buffer);
    latexComp->SetTextSize(0.03);
    latexComp->SetTextColor(kBlue-2);
    latexComp->SetNDC();
    latexComp->Draw();

    sprintf(buffer,"#chi^{2} / ndf = %1.2f",xframe->chiSquare()/flparams->getSize());
    latexComp = new TLatex(0.52,0.56,buffer);
    if(UsingLogyScaleInMass)
        latexComp = new TLatex(0.52-0.49+0.12+TextOffsetX,0.56+0.43-0.04-0.03+TextOffsetY,buffer);
    latexComp->SetTextSize(0.03);
    latexComp->SetTextColor(kBlack);
    latexComp->SetNDC();
    latexComp->Draw();

    //fracs[0] = fbkg0.getVal();
    //fracs[1] = fbkg3.getVal();
    //fracs[2] = fbkg4.getVal();
    fracs[0] = (1./total_data)* Nbkgall.getVal()*fbkg0.getVal();
    //fracs[1] = (1./total_data)* Nbkgall.getVal()*(1. - fbkg0.getVal())* fbkg3.getVal();
    //fracs[2] = (1./total_data)* Nbkgall.getVal()*(1. - fbkg0.getVal())* (1. - fbkg3.getVal());
    fracs[1] = 0.;
    fracs[2] = (1./total_data)* Nbkgall.getVal()*(1. - fbkg0.getVal());
}

void DrawRATIO(TH1F *numeritor, TH1F *denominator){

    TH1F *ratio = (TH1F*)numeritor->Clone();
    ratio->Sumw2();
    float kst = ratio->KolmogorovTest(denominator);
    float chi2test = ratio->Chi2Test(denominator,"CHI2/NDF");
    ratio->Divide(denominator);
    ratio->GetYaxis()->SetLabelSize(0.055);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetTitleSize(0.055);
    ratio->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
    ratio->GetXaxis()->SetLabelSize(0.055);

    ratio->GetYaxis()->SetTitle("ratio");
    ratio->GetYaxis()->SetTitleOffset(0.23);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetRangeUser(0.0,3);
    ratio->GetXaxis()->SetTitleOffset(0.8);
    ratio->GetXaxis()->SetTitleSize(0.21);
    ratio->GetXaxis()->SetLabelSize(0.20);
    ratio->GetYaxis()->SetLabelSize(0.14);
    ratio->GetXaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->SetLineColor(kBlack);
    ratio->Draw();

    TH1F *unc = new TH1F("unc_tmp","",denominator->GetNbinsX(),denominator->GetXaxis()->GetXmin(),
            denominator->GetXaxis()->GetXmax());
    //TGraphAsymmErrors *unc = new TGraphAsymmErrors(denominator->GetNbinsX());
    for(int i=1;i<=denominator->GetNbinsX();i++){ 
        if(denominator->GetBinContent(i)==0) continue;
        unc->SetBinError(i, denominator->GetBinError(i)/denominator->GetBinContent(i));
        unc->SetBinContent(i, 1.0);
        //unc->SetPointEYlow(i,  
        //        1. - (denominator->GetBinContent(i)-denominator->GetBinError(i))/denominator->GetBinContent(i)
        //        ); 
        //unc->SetPoint(i,denominator->GetBinLowEdge(i),1.);
        //unc->SetPointEYhigh(i,  
        //        (denominator->GetBinContent(i)+denominator->GetBinError(i))/denominator->GetBinContent(i) - 1.
        //        );
        printf("[ERROR DEBUG] %i th bin : (x, y+/- e) = (%f, %f +/- %f)\n",
                i,
                denominator->GetBinLowEdge(i),
                denominator->GetBinContent(i),
                denominator->GetBinError(i));
    }
    //unc->SetLineColor(kWhite);
    unc->SetFillColor(kYellow-9);
    unc->SetMarkerSize(0.001);
    unc->Draw("same,E2");

    ratio->Draw("same");

    TLine *line_ = new TLine(ratio->GetXaxis()->GetXmax(),1,ratio->GetXaxis()->GetXmin(),1);
    line_->SetLineStyle(2);
    line_->SetLineWidth(2);
    line_->SetLineColor(kMagenta);
    line_->Draw();

    sprintf(buffer,"Kolmogorov-Smirnov test : %1.3f",kst);
    TLatex *latex_ = new TLatex(0.187,0.86,buffer);
    latex_->SetTextSize(0.12);
    latex_->SetTextColor(kBlack);
    latex_->SetNDC();
    latex_->Draw();

    sprintf(buffer,"#chi^{2} test : %1.3f",chi2test);
    latex_ = new TLatex(0.187,0.74,buffer);
    latex_->SetTextSize(0.12);
    latex_->SetTextColor(kBlack);
    latex_->SetNDC();
    latex_->Draw();
}
void DrawPULL(TH1F *numeritor, TH1F *denominator, string pullname){

    // (numeritor - denominator)/sigma(denominator) +/-  sqrt(sigma(denominator)^2 + sigma(numeritor)^2) /sigma(denominator)

    //TH1F *pull = (TH1F*)numeritor->Clone();
    //sprintf(buffer,"pull_%s",numeritor->GetTitle());
    sprintf(buffer,"pull_%f",rnd->Uniform(0,1));
    TH1F *pull = new TH1F(buffer,"",numeritor->GetNbinsX(),numeritor->GetXaxis()->GetXmin(),numeritor->GetXaxis()->GetXmax());
    for(int i=1;i<=numeritor->GetNbinsX();i++){
        if(denominator->GetBinError(i)==0) continue;
        if(numeritor->GetBinContent(i)==0) continue;
        float pull_ = (numeritor->GetBinContent(i) - denominator->GetBinContent(i)) / denominator->GetBinError(i);
        float error_ = sqrt( denominator->GetBinError(i)*denominator->GetBinError(i) + 
                numeritor->GetBinError(i)*numeritor->GetBinError(i)) / denominator->GetBinError(i);
        printf("[PULL] %i th bin, %f +/- %f \n", i, pull_, error_);
        pull->SetBinContent(i,pull_);
        //pull->SetBinError(i,error_);
        pull->SetBinError(i,0);
    }

    //pull->Sumw2();
    //pull->Divide(denominator);
    //pull->KolmogorovTest(denominator);
    pull->GetYaxis()->SetLabelSize(0.055);
    pull->GetYaxis()->SetTitleOffset(1.2);
    pull->GetYaxis()->SetTitleSize(0.055);
    pull->GetXaxis()->SetTitle("M_{t#gamma} [GeV/c^{2}]");
    pull->GetXaxis()->SetLabelSize(0.055);

    pull->GetYaxis()->SetTitle(pullname.c_str());
    pull->GetYaxis()->SetTitleOffset(0.23);
    pull->GetYaxis()->SetTitleSize(0.2);
    //pull->GetYaxis()->SetRangeUser(0.0,3);
    pull->GetYaxis()->SetRangeUser(-3,3);
    pull->GetXaxis()->SetTitleOffset(0.8);
    pull->GetXaxis()->SetTitleSize(0.21);
    pull->GetXaxis()->SetLabelSize(0.20);
    pull->GetYaxis()->SetLabelSize(0.14);
    pull->GetXaxis()->SetNdivisions(505);
    pull->GetYaxis()->SetNdivisions(505);
    pull->SetLineColor(kBlack);
    pull->SetLineWidth(2);
    pull->SetFillColor(kYellow-9);
    pull->Draw();

    TLine *line_ = new TLine(pull->GetXaxis()->GetXmax(),0,pull->GetXaxis()->GetXmin(),0);
    line_->SetLineStyle(2);
    line_->SetLineWidth(2);
    line_->SetLineColor(kMagenta);
    line_->Draw();
}
double MaxUnc(double NBgks[]){
    double max_ = fabs(NBgks[1] - NBgks[0])/NBgks[1];
    if( max_< fabs(NBgks[1] - NBgks[2])/NBgks[1] )
        max_ = fabs(NBgks[1] - NBgks[2])/NBgks[1];
    return 100.*max_;   //[%]
}

void TrimChar(char *OrignalChar_, char *TrimChar_){
    memset(TrimChar_, 0, sizeof(TrimChar_));
    int matchIndex_ = -1; 
    for(int i=Tgamma500GeV; i<=TGluonTgamma1500GeV; i++){
	if(!(strcmp(SAMPLE[i].tag,OrignalChar_)))
	    matchIndex_=i;
    }   
    if(matchIndex_>=Tgamma500GeV && matchIndex_<=TGluon1500GeV){
	int Nchar = strlen(OrignalChar_) -6 -3;  // remove Tgamma and GeV
	for(int i=0;i<Nchar;i++)
	    TrimChar_[i] = OrignalChar_[6+i];
    }else if(matchIndex_>=TGluonTgamma500GeV && matchIndex_<=TGluonTgamma1500GeV){
	int Nchar = strlen(OrignalChar_) -12 -3;  // remove TGluonTgamma and GeV
	for(int i=0;i<Nchar;i++)
	    TrimChar_[i] = OrignalChar_[12+i];
    }   
}

void PrintUncTable(){
    char TrimChar_[128];
    string UncNames[4+2+2] = {
        // Reduction
        "Reduction (Cat M-0)",
        "Reduction (Cat M-1)",
        "Reduction (Cat M-2)",
        "Reduction (Cat M-3)",
        // SR-CR 
        "SR-CR (Cat 0)",
            //"SR-CR (Cate.3)",
        "SR-CR (Cat 4)",
        // DATA-MC 
        "DATA-MC (Cat 0)",
            //"DATA-MC (Cate.3)",
        "DATA-MC (Cat 4)"
    };
    double NBgkWithUncTemp[Tgamma1500GeV-Tgamma500GeV+1][4+2+2][3];  //[sample][4Merge, 2SR-CR, 2Data-MC][down, normal, up]
    for(int isample=0;isample<Tgamma1500GeV-Tgamma500GeV+1;isample++){
        for(int imode=0;imode<3;imode++){
            for(int iunc=0;iunc<4+2+2;iunc++){
                if(iunc<4){
                    NBgkWithUncTemp[isample][iunc][imode] = NBgkWithUncOnMerge[isample][iunc][imode]; 
                }else if(iunc<4+2 && iunc>=4){
                    NBgkWithUncTemp[isample][iunc][imode] = NBgkWithUncOnSRandCR[isample][2*(iunc-4)][imode]; 
                }else if(iunc<4+2+2 && iunc>= 4+2){
                    NBgkWithUncTemp[isample][iunc][imode] = NBgkWithUncOnDATAandMC[isample][2*(iunc-4-2)][imode]; 
                }
            }
        }
    }

    printf("\n\n\n Latex for uncertainty table \n\n");
    printf("\\begin{tabular}{|l|| "); 
    for(int isample=Tgamma700GeV-Tgamma500GeV;isample<Tgamma1200GeV-Tgamma500GeV+1;isample++)
        printf("c");
    printf("|}\n");
    printf("\\hline\n");

    printf("\\multirow{2}{*}{Source of shape unc.} & \\multicolumn{%i}{c|}{ $\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) [\\GeVcc] =$}\\\\ \n",
            Tgamma1200GeV-Tgamma700GeV+1);

    for(int isample=Tgamma700GeV-Tgamma500GeV;isample<Tgamma1200GeV-Tgamma500GeV+1;isample++){
        TrimChar(SAMPLE[isample + Tgamma500GeV].tag, TrimChar_);
        printf("& %5s ",TrimChar_); 
    }
    printf("\\\\ \n");


    printf("\\hline\n");
    printf("\\hline\n");
    for(int iunc=0;iunc<4+2+2;iunc++){
        printf("%s ",UncNames[iunc].c_str());
        for(int isample=Tgamma700GeV-Tgamma500GeV;isample<Tgamma1200GeV-Tgamma500GeV+1;isample++){
            printf("& %1.2f ",MaxUnc(NBgkWithUncTemp[isample][iunc]) );
        }
        printf("\\\\ \n");
    }
    printf("\\hline\n");
    printf("\\end{tabular}\n");
}

TH1F *hshape(int M_, int iUncName){
    sprintf(buffer,"tree_Tgamma%iGeV%s",M_,RunStatusNames[iUncName].c_str());
    TTree *tree_TgammaGeV_tmp = (TTree *) file->Get(buffer);

    sprintf(buffer,"h_Tgamma%iGeV%s",M_,RunStatusNames[iUncName].c_str());
    TH1F *hist = new TH1F(buffer,"",numberbins,Xregions[0],Xregions[1]);
    hist->Sumw2();

    float weight_;
    float Mass_;
    float category_;
    tree_TgammaGeV_tmp->SetBranchAddress("weight",&weight_);
    tree_TgammaGeV_tmp->SetBranchAddress("Mass",&Mass_);
    tree_TgammaGeV_tmp->SetBranchAddress("category",&category_);

    for(int ientry = 0; ientry < tree_TgammaGeV_tmp->GetEntries(); ientry++){
	tree_TgammaGeV_tmp->GetEntry(ientry);
	if(category_==10)
	    hist->Fill(Mass_,weight_);
    }   

    hist->GetXaxis()->SetTitle("Mass(top+#gamma)");
    delete tree_TgammaGeV_tmp;
    return hist;
}

void drawCMSLumi(){

    TLatex * latexCMS;
    // Lumi
    latexCMS = new TLatex(0.94,0.96,"19.7 fb^{-1} (8 TeV)");
    latexCMS->SetTextFont(42);
    latexCMS->SetTextAlign(31);
    latexCMS->SetTextSize(0.04);
    latexCMS->SetNDC();
    latexCMS->Draw();

    // cms
    latexCMS = new TLatex(0.13,0.96,"CMS");
    latexCMS->SetTextFont(61);
    latexCMS->SetTextAlign(11);
    latexCMS->SetTextSize(0.05);
    latexCMS->SetNDC();
    latexCMS->Draw();

    // preliminary
    latexCMS = new TLatex(0.13+0.01 + 0.12*(1-0.13),0.96,"Preliminary");
    latexCMS->SetTextSize(0.05*0.76);
    latexCMS->SetTextFont(52);
    latexCMS->SetTextAlign(11);
    latexCMS->SetNDC();
    latexCMS->Draw();
}
