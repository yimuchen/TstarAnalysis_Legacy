#ifndef DRAWSTACK_H
#define DRAWSTACK_H

#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLatex.h"
#include "ConstantNumbers.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"

int skipData=0;
int skipQCD=1;
bool IsRatioOrPull = false; // [false : true] = [pull, ratio]
// hmerge[0,1,2,3,4] = [boson(s), TTbar, singleTop, ttX , QCD]
void DRAWSTACK(TCanvas *c1, THStack *Stacks1_1,
        TH1F *hmerge[], TH1F *signal,TH1F *data, TLegend *legend_nm,string Ytitle, string Xtitle, string Unit, int rebin_,
        TH1F *hMainSysUnc){

    gStyle->SetOptStat(0);
    char buffer[128];

    float padtop_ = 0.15;
    float padbottom_ = padtop_ + 0.043;
    TPad *tp2 = new TPad("tp2","",0,0,1,padbottom_);
    TPad *tp1 = new TPad("tp1","",0,padtop_,1,1);

    c1->cd();
    tp1->Draw();
    tp1->cd();
    tp1->SetBottomMargin(0.05);
    tp1->SetTopMargin(0.06);


    /*  
        enum EColor { kWhite =0,   kBlack =1,   kGray=920,
        kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
        kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
     */
    const int N_file = 7;
    int hColor[N_file]={
        800-7,
        600-6,
        kAzure+2,//432+2,
        616+2,
        900+2,
        1,
        1};
    int hFillstyle[N_file]={
        3013,
        1001,
        3005,
        3007,
        3006,
        0,
        0};

    for(int i=0;i<5;i++) {
        hmerge[i]->SetFillColor(hColor[i]);
        hmerge[i]->SetLineColor(hColor[i]);
        hmerge[i]->SetFillStyle(hFillstyle[i]);
        hmerge[i]->Rebin(rebin_);
    }
    signal->SetFillColor(hColor[5]);
    signal->SetLineColor(hColor[5]);
    signal->SetFillStyle(hFillstyle[5]);
    signal->Rebin(rebin_);
    data->SetFillColor(hColor[6]);
    data->SetLineColor(hColor[6]);
    data->SetFillStyle(hFillstyle[6]);
    data->Rebin(rebin_);
    data->SetMarkerStyle(21);
    data->SetMarkerSize(1.0);
    data->Sumw2();

    legend_nm = new TLegend(0.577,0.605,0.9,0.92);
    legend_nm->AddEntry(signal,"t*#bar{t*} (950 GeV/c^{2})","f");
    legend_nm->AddEntry(hmerge[0],"boson(s)","f");
    legend_nm->AddEntry(hmerge[1],"t#bar{t}+j","f");
    legend_nm->AddEntry(hmerge[2],"singleTop","f");
    //legend_nm->AddEntry(hmerge[2],"boson(s)","f");
    legend_nm->AddEntry(hmerge[3],"t#bar{t}X","f");
    if(!skipQCD) legend_nm->AddEntry(hmerge[4],"QCD","f");
    if(!skipData) legend_nm->AddEntry(data,"data","lep");
    legend_nm->SetBorderSize(0);
    legend_nm->SetFillColor(0);
    legend_nm->SetFillStyle(0);
    legend_nm->SetNColumns(1);
    legend_nm->SetTextSize(0.04);
    legend_nm->SetTextSizePixels(25);

    TH1F *unc; 
    bool given=0;

    for(int i=0;i<5;i++) {
        if(skipQCD&&i==0) continue;  // skip QCD
        Stacks1_1->Add(hmerge[4-i]);
        if(!given){
            unc = (TH1F*)hmerge[4-i]->Clone();
            given=1;
        }else{
            unc->Add(hmerge[4-i]);
        }   

    }
    Stacks1_1->Add(signal);

    Stacks1_1->Draw();

    float max_ = data->GetMaximum();
    if(Stacks1_1->GetMaximum()>max_) max_=Stacks1_1->GetMaximum();
    int secondConcern_ = data->GetXaxis()->GetNbins()*0.6;
    if(data->GetBinContent(secondConcern_)>0.605*max_) max_ = 1.1*data->GetBinContent(secondConcern_)/0.605;
    //secondConcern_ = Stacks1_1->GetXaxis()->GetNbins()*0.6;
    max_ = pow(max_, 1.2);
    //if(Stacks1_1->GetBinContent(secondConcern_)>0.605*max_) max_ = Stacks1_1->GetBinContent(secondConcern_)/0.605;
    if(!skipData) Stacks1_1->SetMaximum(1.1*max_);

    float min_ = 10001;
    for(int ibin=1;ibin<=data->GetNbinsX();ibin++)
        if(min_> data->GetBinContent(ibin) && data->GetBinContent(ibin)!=0)
            min_ = data->GetBinContent(ibin);
    Stacks1_1->SetMinimum(pow(min_,0.9));

    sprintf(buffer,"%s / (%1.1f %s)",Ytitle.c_str(),signal->GetBinWidth(1),Unit.c_str());
    if(!strcmp(Unit.c_str(),"-"))
        sprintf(buffer,"%s / (%1.1f)",Ytitle.c_str(),signal->GetBinWidth(1));
    Stacks1_1->GetYaxis()->SetTitle(buffer);
    Stacks1_1->GetYaxis()->SetLabelSize(0.055);
    Stacks1_1->GetYaxis()->SetTitleOffset(0.9);
    Stacks1_1->GetYaxis()->SetTitleSize(0.055);
    sprintf(buffer,"%s [%s]",Xtitle.c_str(),Unit.c_str());
    if(!strcmp(Unit.c_str(),"-"))
        sprintf(buffer,"%s",Xtitle.c_str());
    Stacks1_1->GetXaxis()->SetTitle(buffer);
    Stacks1_1->GetXaxis()->SetLabelSize(0.055);
    Stacks1_1->GetXaxis()->SetTitleOffset(0.95);
    Stacks1_1->GetXaxis()->SetNdivisions(505);
    Stacks1_1->Draw("hist");
    tp1->SetLogy(1);
    if(!skipData) data->Draw("P,same");

    if(RunStatus_ == UncQsquare){
        printf("\n");
        float unc_sys= 0.;
        float unc_stat= 0.;
        for(int ibin = 1;ibin <= unc->GetNbinsX(); ibin++){
            unc_sys = hMainSysUnc->GetBinContent(ibin) * unc->GetBinContent(ibin); // hMainSysUnc :RelUnc; unc : absolut
            unc_stat= unc->GetBinError(ibin);
            unc->SetBinError(ibin, sqrt( unc_sys*unc_sys +
                        unc->GetBinError(ibin) * unc->GetBinError(ibin) ));
            printf("[%i th bin] %f (stat, syst) = (%f, %f)\n", ibin,unc->GetBinError(ibin), unc_stat, unc_sys);
        }
    }
    unc->SetMarkerSize(0.001);
    unc->SetLineColor(kWhite);
    //unc->SetFillColor(kGray+2);
    unc->SetFillColor(kBlack);
    //unc->SetFillStyle(3002);
    unc->SetFillStyle(3005);
    unc->Draw("E2,same");

    if(RunStatus_ == UncQsquare){
        legend_nm->AddEntry(unc,"unc.(syst.+stat.)","fe");
    }else{
        legend_nm->AddEntry(unc,"unc.(stat.)","fe");
    }


    legend_nm->Draw();

    char bufferLumi[128];
    sprintf(bufferLumi,"2012 CMS %1.1ffb^{-1} #sqrt{s}=8TeV",LUMINOSITY/1000.);
    TLatex *latexCMS = new TLatex(0.48,0.96,bufferLumi);
    latexCMS->SetTextSize(0.04);
    latexCMS->SetNDC();
    latexCMS->Draw();

    c1->cd();
    tp2->Draw();
    tp2->cd();
    tp2->SetTopMargin(0.0);
    tp2->SetBottomMargin(0.4);

    TH1F *ratio;
    int baseline = 1;
    double yranges[2] = {0.,2.};
    if(IsRatioOrPull){
        ratio = (TH1F*)data->Clone();
        ratio->Sumw2();
        ratio->Divide(unc);
        ratio->GetYaxis()->SetTitle("Data/MC");
    }else{
        yranges[0] = 999.;   // min
        yranges[1] = -999.;    // max
        baseline = 0;
        ratio = new TH1F("ratio_tmp","", data->GetNbinsX(),data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax());
        for(int ibin=1;ibin <= data->GetNbinsX();ibin++){
            if(data->GetBinContent(ibin)==0) continue;
            if(unc->GetBinError(ibin)!=0){
                ratio->SetBinContent(ibin, (data->GetBinContent(ibin) - unc->GetBinContent(ibin))/unc->GetBinError(ibin));
                //ratio->SetBinError(ibin, 
                //        sqrt(data->GetBinError(ibin)*data->GetBinError(ibin) + 
                //            unc->GetBinError(ibin)*unc->GetBinError(ibin)) /unc->GetBinError(ibin));
                ratio->SetBinError(ibin, 0);
                if (ratio->GetBinContent(ibin) > yranges[1] )
                    yranges[1] = ratio->GetBinContent(ibin);
                if (ratio->GetBinContent(ibin) < yranges[0] )
                    yranges[0] = ratio->GetBinContent(ibin);
            }
        }
        if(fabs(yranges[0]) > fabs(yranges[1])){
            yranges[1] = fabs(yranges[0]);
        }else{
            yranges[0] = -1.* fabs(yranges[1]);
        }
        if(fabs(yranges[0])>3.){
            yranges[0] = -3.;
            yranges[1] = 3.;
        }
        //ratio->GetYaxis()->SetTitle("(Data-MC)/#sigma");
        ratio->GetYaxis()->SetTitle("pull");
        ratio->SetLineColor(kBlack);
        ratio->SetLineWidth(2);
        ratio->SetFillColor(kYellow-9);

    }

    ratio->GetYaxis()->SetLabelSize(0.055);
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetTitleSize(0.055);
    sprintf(buffer,"%s [%s]",Xtitle.c_str(),Unit.c_str());
    if(!strcmp(Unit.c_str(),"-"))
        sprintf(buffer,"%s",Xtitle.c_str());
    ratio->GetXaxis()->SetTitle(buffer);
    ratio->GetXaxis()->SetLabelSize(0.055);

    ratio->GetYaxis()->SetTitleOffset(0.23);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetRangeUser(yranges[0],yranges[1]);
    ratio->GetXaxis()->SetTitleOffset(0.8);
    ratio->GetXaxis()->SetTitleSize(0.21);
    ratio->GetXaxis()->SetLabelSize(0.20);
    ratio->GetYaxis()->SetLabelSize(0.14);
    ratio->GetXaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->Draw();

    TLine *line_ = new TLine(ratio->GetXaxis()->GetXmax(),baseline,ratio->GetXaxis()->GetXmin(),baseline);
    line_->SetLineStyle(2);
    line_->SetLineWidth(2);
    line_->SetLineColor(kMagenta);
    line_->Draw();
}

// hmerge[0,1,2,3,4] = [boson(s), TTbar, singleTop, ttX , QCD]
void DRAWSTACK(THStack *Stacks1_1,
        TH2D *hmerge[], TH2D *signal,TH2D *data, TLegend *legend_nm,string Ytitle, string Xtitle, 
        string XUnit, int Xrebin_,
        string YUnit, int Yrebin_
        ){

    gStyle->SetPalette(1);
    char buffer[128];

    /*  
        enum EColor { kWhite =0,   kBlack =1,   kGray=920,
        kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
        kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
     */
    const int N_file = 7;
    int hColor[N_file]={
        800-7,
        600-6,
        432+2,
        616+2,
        900+2,
        1,
        1};
    int hFillstyle[N_file]={
        3013,
        1001,
        3005,
        3007,
        3006,
        0,
        0};

    for(int i=0;i<5;i++) {
        hmerge[i]->SetFillColor(hColor[i]);
        hmerge[i]->SetLineColor(hColor[i]);
        hmerge[i]->SetFillStyle(hFillstyle[i]);
        hmerge[i]->Rebin2D(Xrebin_,Yrebin_);
    }
    signal->SetFillColor(hColor[5]);
    signal->SetLineColor(hColor[5]);
    signal->SetFillStyle(hFillstyle[5]);
    signal->Rebin2D(Xrebin_,Yrebin_);
    data->SetFillColor(hColor[6]);
    data->SetLineColor(hColor[6]);
    data->SetFillStyle(hFillstyle[6]);
    data->Rebin2D(Xrebin_,Yrebin_);
    data->SetMarkerStyle(21);
    data->SetMarkerSize(1.0);
    data->Sumw2();

    legend_nm = new TLegend(0.577,0.605,0.9,0.92);
    //legend_nm->AddEntry(data,"data","lep");
    //legend_nm->AddEntry(signal,"T#bar{T} (500 GeV/c^{2})","f");
    legend_nm->AddEntry(hmerge[0],"boson(s)","f");
    legend_nm->AddEntry(hmerge[1],"t#bar{t}+j","f");
    legend_nm->AddEntry(hmerge[2],"singleTop","f");
    legend_nm->AddEntry(hmerge[3],"t#bar{t}X","f");
    legend_nm->AddEntry(hmerge[4],"QCD","f");
    legend_nm->SetBorderSize(0);
    legend_nm->SetFillColor(0);
    legend_nm->SetFillStyle(0);
    legend_nm->SetNColumns(1);
    legend_nm->SetTextSize(0.04);
    legend_nm->SetTextSizePixels(25);

    for(int i=0;i<5;i++) {
        Stacks1_1->Add(hmerge[4-i]);
    }
    //Stacks1_1->Add(signal);

    Stacks1_1->Draw();

    float max_ = data->GetMaximum();
    if(Stacks1_1->GetMaximum()>max_) max_=Stacks1_1->GetMaximum();

    int secondConcern_ = data->GetXaxis()->GetNbins()*(0.55);
    for(int n_=secondConcern_;n_<=data->GetXaxis()->GetNbins();n_++){
        if(data->GetBinContent(n_)>0.605*max_) max_ = 1.1*data->GetBinContent(n_)/0.605;
        if(Stacks1_1->GetHistogram()->GetBinContent(n_)>0.605*max_) 
            max_ = Stacks1_1->GetHistogram()->GetBinContent(n_)/0.605;
    }

    //Stacks1_1->SetMaximum(20*max_);
    Stacks1_1->SetMaximum(1.2*max_);

    sprintf(buffer,"%s [%s]",Ytitle.c_str(),YUnit.c_str());
    //sprintf(buffer,"%s / (%1.1f %s)",Ytitle.c_str(),signal->GetYaxis()->GetBinWidth(1),YUnit.c_str());
    if(!strcmp(XUnit.c_str(),"-"))
        sprintf(buffer,"%s",Ytitle.c_str());
        //sprintf(buffer,"%s / (%1.1f)",Ytitle.c_str(),signal->GetYaxis()->GetBinWidth(1));
    Stacks1_1->GetYaxis()->SetTitle(buffer);
    Stacks1_1->GetYaxis()->SetLabelSize(0.055);
    Stacks1_1->GetYaxis()->SetTitleOffset(1.2);
    Stacks1_1->GetYaxis()->SetTitleSize(0.055);
    //sprintf(buffer,"%s / (%1.1f %s)",Xtitle.c_str(),signal->GetXaxis()->GetBinWidth(1),XUnit.c_str());
    sprintf(buffer,"%s [%s]",Xtitle.c_str(),XUnit.c_str());
    if(!strcmp(XUnit.c_str(),"-"))
        sprintf(buffer,"%s",Xtitle.c_str());
        //sprintf(buffer,"%s / (%1.1f)",Xtitle.c_str(),signal->GetXaxis()->GetBinWidth(1));
    Stacks1_1->GetXaxis()->SetTitle(buffer);
    Stacks1_1->GetXaxis()->SetLabelSize(0.055);
    Stacks1_1->GetXaxis()->SetTitleOffset(0.95);
    Stacks1_1->GetXaxis()->SetNdivisions(505);
    Stacks1_1->Draw();
    //gPad->SetLogy(1);
    //Stacks1_1->Draw("lego1");
    Stacks1_1->Draw("colz");
    //data->Draw("P,same");

    //legend_nm->Draw();

    char bufferLumi[128];
    sprintf(bufferLumi,"2011 CMS %1.1ffb^{-1} #sqrt{s}=7TeV",LUMINOSITY/1000.);
    TLatex *latexCMS = new TLatex(0.48,0.96,bufferLumi);
    latexCMS->SetTextSize(0.04);
    latexCMS->SetNDC();
    latexCMS->Draw();
}

void DRAWSTACK(THStack *Stacks1_1,
        TH1F *hmerge[], TH1F *signal,TH1F *data, TLegend *legend_nm,string Ytitle, string Xtitle, string Unit, int rebin_){

    gStyle->SetOptStat(0);
    char buffer[128];

    /*  
        enum EColor { kWhite =0,   kBlack =1,   kGray=920,
        kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
        kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
     */
    const int N_file = 7;
    int hColor[N_file]={
        800-7,
        600-6,
        kAzure+2,//432+2,
        616+2,
        900+2,
        1,
        1};
    int hFillstyle[N_file]={
        3013,
        1001,
        3005,
        3007,
        3006,
        0,
        0};

    for(int i=0;i<5;i++) {
        hmerge[i]->SetFillColor(hColor[i]);
        hmerge[i]->SetLineColor(hColor[i]);
        hmerge[i]->SetFillStyle(hFillstyle[i]);
        hmerge[i]->Rebin(rebin_);
    }
    signal->SetFillColor(hColor[5]);
    signal->SetLineColor(hColor[5]);
    signal->SetFillStyle(hFillstyle[5]);
    signal->Rebin(rebin_);
    data->SetFillColor(hColor[6]);
    data->SetLineColor(hColor[6]);
    data->SetFillStyle(hFillstyle[6]);
    data->Rebin(rebin_);
    data->SetMarkerStyle(21);
    data->SetMarkerSize(1.0);
    data->Sumw2();

    legend_nm = new TLegend(0.577,0.605,0.9,0.92);
    legend_nm->AddEntry(signal,"t*#bar{t*} (950 GeV/c^{2})","f");
    legend_nm->AddEntry(hmerge[0],"boson(s)","f");
    legend_nm->AddEntry(hmerge[1],"t#bar{t}+j","f");
    legend_nm->AddEntry(hmerge[2],"singleTop","f");
    //legend_nm->AddEntry(hmerge[2],"boson(s)","f");
    legend_nm->AddEntry(hmerge[3],"t#bar{t}X","f");
    if(!skipQCD) legend_nm->AddEntry(hmerge[4],"QCD","f");
    if(!skipData) legend_nm->AddEntry(data,"data","lep");
    legend_nm->SetBorderSize(0);
    legend_nm->SetFillColor(0);
    legend_nm->SetFillStyle(0);
    legend_nm->SetNColumns(1);
    legend_nm->SetTextSize(0.04);
    legend_nm->SetTextSizePixels(25);

    TH1F *unc; 
    bool given=0;

    for(int i=0;i<5;i++) {
        if(skipQCD&&i==0) continue;  // skip QCD
        Stacks1_1->Add(hmerge[4-i]);
        if(!given){
            unc = (TH1F*)hmerge[4-i]->Clone();
            given=1;
        }else{
            unc->Add(hmerge[4-i]);
        }   

    }
    Stacks1_1->Add(signal);

    Stacks1_1->Draw();

    float max_ = data->GetMaximum();
    if(Stacks1_1->GetMaximum()>max_) max_=Stacks1_1->GetMaximum();
    int secondConcern_ = data->GetXaxis()->GetNbins()*0.6;
    if(data->GetBinContent(secondConcern_)>0.605*max_) max_ = 1.1*data->GetBinContent(secondConcern_)/0.605;
    //secondConcern_ = Stacks1_1->GetXaxis()->GetNbins()*0.6;
    //if(Stacks1_1->GetBinContent(secondConcern_)>0.605*max_) max_ = Stacks1_1->GetBinContent(secondConcern_)/0.605;
    if(!skipData) Stacks1_1->SetMaximum(1.1*max_);

    sprintf(buffer,"%s / (%1.1f %s)",Ytitle.c_str(),signal->GetBinWidth(1),Unit.c_str());
    if(!strcmp(Unit.c_str(),"-"))
        sprintf(buffer,"%s / (%1.1f)",Ytitle.c_str(),signal->GetBinWidth(1));
    Stacks1_1->GetYaxis()->SetTitle(buffer);
    Stacks1_1->GetYaxis()->SetLabelSize(0.055);
    Stacks1_1->GetYaxis()->SetTitleOffset(0.9);
    Stacks1_1->GetYaxis()->SetTitleSize(0.055);
    sprintf(buffer,"%s [%s]",Xtitle.c_str(),Unit.c_str());
    if(!strcmp(Unit.c_str(),"-"))
        sprintf(buffer,"%s",Xtitle.c_str());
    Stacks1_1->GetXaxis()->SetTitle(buffer);
    Stacks1_1->GetXaxis()->SetLabelSize(0.055);
    Stacks1_1->GetXaxis()->SetTitleOffset(0.95);
    Stacks1_1->GetXaxis()->SetNdivisions(505);
    Stacks1_1->Draw("hist");
    if(!skipData) data->Draw("P,same");
    unc->SetMarkerSize(0.001);
    unc->SetLineColor(kWhite);
    //unc->SetFillColor(kGray+2);
    unc->SetFillColor(kBlack);
    unc->SetFillStyle(3002);
    unc->Draw("E2,same");
    legend_nm->AddEntry(unc,"uncertainty","fe");


    legend_nm->Draw();

    char bufferLumi[128];
    sprintf(bufferLumi,"2012 CMS %1.1ffb^{-1} #sqrt{s}=8TeV",LUMINOSITY/1000.);
    TLatex *latexCMS = new TLatex(0.48,0.96,bufferLumi);
    latexCMS->SetTextSize(0.04);
    latexCMS->SetNDC();
    latexCMS->Draw();

}

#endif
