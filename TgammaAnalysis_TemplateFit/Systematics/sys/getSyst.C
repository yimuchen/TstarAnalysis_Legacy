#include "interface/ConstantNumbers.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <string>

const int Ncate = 5;
float MaxAllowedUnc = 0.4; // 0.1 means 10% --> because Qscale or Matching has larger uncertainty due to limited statistics
TH1F *htemplates[Ncate][RunStatusSize];
TH1F *htemplatesUnc[Ncate];
char buffer_[128];
void loadTemplates(string RegionName);
TH1F *hdiff(string hname, TH1F *hNorm, TH1F *hplus, TH1F *hminus);
void getallsys(TH1F *hall, int N, TH1F *hunc[], int mode/*0 : including all; 1 : without Q^2 and match*/);

const int Nregion = 2; // CR, SR
string RegionsName[Nregion] = {
    "", // CR
    " in SR"    // SR
};

void getSyst(){

    TFile *output = new TFile("sysforshape.root","recreate");
    const int Nsys = 10;
    TH1F *hSumSys[Ncate][Nregion];
    for(int irg = 0; irg < Nregion ; irg ++){

        loadTemplates(RegionsName[irg]);
        for(int icate=0;icate<Ncate;icate++){
            sprintf(buffer_,"SumSys_cate%i%s",icate,RegionsName[irg].c_str());
            hSumSys[icate][irg] = new TH1F(buffer_,"", htemplates[icate][Normal]->GetXaxis()->GetNbins(), 
                    htemplates[icate][Normal]->GetXaxis()->GetXmin(), 
                    htemplates[icate][Normal]->GetXaxis()->GetXmax());
            TH1F *hUncs[Nsys];

            printf("[  Cate %i loop]\n", icate);
            hUncs[0] = hdiff("XsecUnc", htemplates[icate][Normal], htemplates[icate][UncXsecPlus], 
                    htemplates[icate][UncXsecMinus]);
            hUncs[1] = hdiff("PUUnc", htemplates[icate][Normal], htemplates[icate][UncPUPlus], 
                    htemplates[icate][UncPUMinus]);
            hUncs[2] = hdiff("JESUnc", htemplates[icate][Normal], htemplates[icate][UncJESPlus], 
                    htemplates[icate][UncJESMinus]);
            hUncs[3] = hdiff("JERUnc", htemplates[icate][Normal], htemplates[icate][UncJERPlus], 
                    htemplates[icate][UncJERMinus]);
            hUncs[4] = hdiff("TrigUnc", htemplates[icate][Normal], htemplates[icate][UncTrigPlus], 
                    htemplates[icate][UncTrigMinus]);
            hUncs[5] = hdiff("PhoIDUnc", htemplates[icate][Normal], htemplates[icate][UncPhoIDPlus], 
                    htemplates[icate][UncPhoIDMinus]);
            hUncs[6] = hdiff("LepIDUnc", htemplates[icate][Normal], htemplates[icate][UncLepIDPlus], 
                    htemplates[icate][UncLepIDMinus]);
            hUncs[7] = hdiff("TopPtUnc", htemplates[icate][Normal], htemplates[icate][UncTopPtPlus], 
                    htemplates[icate][UncTopPtMinus]);
            hUncs[8] = hdiff("QsquareUnc", htemplates[icate][Normal], htemplates[icate][UncQsquarePlus], 
                    htemplates[icate][UncQsquareMinus]);
            hUncs[9] = hdiff("MatchingUnc", htemplates[icate][Normal], htemplates[icate][UncMatchingPlus], 
                    htemplates[icate][UncMatchingMinus]);

            getallsys(hSumSys[icate][irg], Nsys, hUncs, irg);
            TCanvas *c1 = new TCanvas("c1","",640,640);
            c1->Divide(5,2);

            for(int iu = 0;iu < Nsys;iu++){
                c1->cd(iu+1);
                hUncs[iu]->Draw();
            }
            sprintf(buffer_,"Merge%i%s_sys.pdf",icate,RegionsName[irg].c_str());
            c1->SaveAs(buffer_);

            delete c1;
        }

        TCanvas *c1 = new TCanvas("c1","",640,640);
        c1->Divide(3,2);
        for(int icate=0;icate < Ncate;icate++){
            c1->cd(icate + 1);
            hSumSys[icate][irg]->Draw();
            output->cd();
            hSumSys[icate][irg]->Write();

        }
        sprintf(buffer_,"SumSys%s.pdf", RegionsName[irg].c_str());
        c1->SaveAs(buffer_);
    }

    output->Close();
}

void loadTemplates(string RegionName){
    for(int irun=0;irun<RunStatusSize;irun++){
        if(irun==UncQsquare) continue;
        if(irun!=Normal)
            sprintf(buffer_,"../%s/FittingOutput/output_theta_1.0.root",RunStatusNames[irun].c_str());
        else
            sprintf(buffer_,"../Normal/FittingOutput/output_theta_1.0.root");

        TFile *file_ = new TFile(buffer_);
        for(int icate=0;icate<5;icate++){
            sprintf(buffer_,"Merge_cate-%i%s",icate,RegionName.c_str());
            printf("[Loading] %s in %s\n",buffer_, RunStatusNames[irun].c_str());
            htemplates[icate][irun] = (TH1F*) file_->Get(buffer_);
        }
    }
}

TH1F *hdiff(string hname, TH1F *hNorm, TH1F *hplus, TH1F *hminus){
    TH1F *hdifftmp = new TH1F(hname.c_str(),"", hNorm->GetNbinsX(), hNorm->GetXaxis()->GetXmin(), hNorm->GetXaxis()->GetXmax());
    float unc_ = 0.;
    float stat_normal = 0.;
    float stat_up = 0.;
    float stat_down = 0.;
    float sys_up = 0.;
    float sys_down = 0.;
    for(int ibin=1;ibin<=hNorm->GetNbinsX();ibin++){
        unc_ = 0;
        stat_normal = 0;
        stat_up = 0;
        stat_down = 0;
        sys_up = 0;
        sys_down = 0;
        if(hNorm->GetBinContent(ibin)==0){
	    hdifftmp->SetBinContent(ibin,unc_);
	    continue;
	}

        stat_normal = hNorm->GetBinError(ibin)/hNorm->GetBinContent(ibin);
        if(hplus->GetBinContent(ibin)!=0)
            stat_up = hplus->GetBinError(ibin)/hplus->GetBinContent(ibin);

        if(hminus->GetBinContent(ibin)!=0)
            stat_down = hminus->GetBinError(ibin)/hminus->GetBinContent(ibin);

	sys_up = fabs(hNorm->GetBinContent(ibin) - hplus->GetBinContent(ibin))/ hNorm->GetBinContent(ibin);
	sys_down = fabs(hNorm->GetBinContent(ibin) - hminus->GetBinContent(ibin))/ hNorm->GetBinContent(ibin);

        if(stat_up <  MaxAllowedUnc)
            unc_ = sys_up;

        if(( unc_< sys_down && stat_down <  MaxAllowedUnc ) )   
            unc_ = sys_down;

	//if(unc_>1.00){
	//    unc_ = sys_up;
	//    if(sys_down<sys_up)
	//	unc_ = sys_down;
	//}

        if(hminus->GetBinContent(ibin) == hplus->GetBinContent(ibin))
            unc_ = 0.; 

        printf("    [Sys %s (bin %i) ] %f +/- %f (+ %f +/- %f, - %f +/- %f) with max error of %f \n",
                hname.c_str(), 
                ibin,
                hNorm->GetBinContent(ibin),
                hNorm->GetBinError(ibin),
                hplus->GetBinContent(ibin),
                hplus->GetBinError(ibin),
                hminus->GetBinContent(ibin),
                hminus->GetBinError(ibin),
                unc_
              );
        hdifftmp->SetBinContent(ibin,unc_);
    }
    return hdifftmp;
}

void getallsys(TH1F *hall, int N, TH1F *hunc[], int mode/*0 : including all; 1 : without Q^2 and match*/){

    for(int ibin = 1;ibin <= hall->GetXaxis()->GetNbins();ibin++){
        float sum=0.;
        for(int isys = 0; isys < N; isys++){
            //if((isys == 2 ))	// no JES due to limited statistics causing larger unc
            //    continue;
            //if((isys == 3 ))	// no JER due to limited statistics causing larger unc
            //    continue;
            if(mode==1&&(isys == N-1 || isys == N-2))
                continue;

            sum += pow(hunc[isys]->GetBinContent(ibin),2);
        }
        printf("        [Total sys (bin %i)]  max error of %f \n",ibin, sqrt(sum) );
        hall->SetBinContent(ibin, sqrt(sum) );
    }

}
