#ifndef SAMPLES_tr_H
#define SAMPLES_tr_H

#include "ConstantNumbers.h"
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV?skin=drupal

//Top-14-016 (LHC dilep combination 241.5 +/- 5.869) Top-006 (semi-lep 235.4+/-14.314),  Top-007 (di-lep 239), Top-026 ( tau 257)
double ttbarXsec = 235.4; 
double ttbarXsecUnc = 14.314; 

enum samples_tr_order{
    Data,
    //Tgamma450GeV,
    Tgamma500GeV,
    Tgamma550GeV,
    Tgamma600GeV,
    Tgamma650GeV,
    Tgamma700GeV,
    Tgamma750GeV,
    Tgamma800GeV,
    Tgamma850GeV,
    Tgamma900GeV,
    Tgamma950GeV,
    Tgamma1000GeV,
    Tgamma1100GeV,
    Tgamma1200GeV,
    Tgamma1300GeV,
    Tgamma1400GeV,
    Tgamma1500GeV,
    // Q^2 --start--
    TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_powheg_v2_scaledown,
    TTJets_powheg_v2_scaleup,
    DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,

    // Q^2 --end--
    // more statistics --start--
    ///*
    TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2,
    TTJets_HadronicMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_SemiLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_ext_v1,

    TTJets_MSDecays,
    TTJets_powheg_v2,
    TT_700_1000_powheg,
    TT_1000_inf_powheg,

    DYJetsToLL_PtZ_100_TuneZ2star_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2,
    DYJetsToLL_PtZ_50To70_TuneZ2star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_PtZ_70To100_TuneZ2star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2,
    DYToEE_M_20_CT10_TuneZ2star_v2_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYToMuMu_M_20_CT10_TuneZ2star_v2_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,

    ZZJetsTo2L2Q_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo2e2mu_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo2e2tau_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo2mu2tau_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo4e_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo4mu_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZTo4tau_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,

    GluGluToHToZZTo4L_M_126_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    GluGluToHToZZTo2L2Q_M_125_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    VBF_HToZZTo2L2Q_M_125_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTbarH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1,

    VBF_HToZZTo4L_M_126_8TeV_powheg_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1,
    ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1,
    //*/
    // more statistics --end--
    QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2,
    T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_10To50filter_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    //GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    GJets_HT_40To100_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V19_v1,
    GJets_HT_100To200_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V19_v1,
    GJets_HT_200To400_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    GJets_HT_400ToInf_8TeV_madgraph_v3_Summer12_DR53X_PU_S10_START53_V7C_v1,
    DiPhotonJets_8TeV_madgraph_tarball_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WGToLNuG_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZGToLLG_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2,
    WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WWGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    samples_order_size
};

struct sample_t{
    char filename[512];
    char latex[128];
    char tag[64];
    double ngen;
    double xsec;
    double unc; // 1 : 100%, 0.11 : 11%, 0 : no uncertainty
};


struct sample_t SAMPLE[samples_order_size] = {
    {"./REDUCE_DATA/PhotonData/1lepton4jets_*_Run2012*.root",
        "Data",
        "data", LUMINOSITY , 1 , 0},
    //{"./REDUCE_DATA/1lepton4jets_Tgamma450GeV.root",
    //    "Tgamma450GeV", 28546 , 11.984291 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 500\\GeVcc$",
        "Tgamma500GeV", 97662 , 5.087496 , 0},  // LO from excited quark  // okay
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-550_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 550\\GeVcc$",
        "Tgamma550GeV", 98650 , 2.315665 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-600_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 600\\GeVcc$",
        "Tgamma600GeV", 99408 , 1.115097 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-650_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 650\\GeVcc$",
        "Tgamma650GeV", 101131 , 0.562451 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-700_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 700\\GeVcc$",
        "Tgamma700GeV", 98736 , 0.294853 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-750_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 750\\GeVcc$",
        "Tgamma750GeV", 99847 , 0.159653 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-800_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 800\\GeVcc$",
        "Tgamma800GeV", 97373 , 0.088853 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-850_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 850\\GeVcc$",
        "Tgamma850GeV", 98312 , 5.06221252E-02 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-900_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 900\\GeVcc$",
        "Tgamma900GeV", 96254 , 2.94248249E-02 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-950_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 950\\GeVcc$",
        "Tgamma950GeV", 97617 , 1.74004431E-02 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1000\\GeVcc$",
        "Tgamma1000GeV", 95614 , 1.04431023E-02 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1100\\GeVcc$",
        "Tgamma1100GeV", 94642 , 3.90115915E-03 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1200_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1200\\GeVcc$",
        "Tgamma1200GeV", 93570 , 1.51328287E-03 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1300_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1300\\GeVcc$",
        "Tgamma1300GeV", 95859 , 6.03060062E-04 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1400_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1400\\GeVcc$",
        "Tgamma1400GeV", 95396 , 2.44765697E-04 , 0},  // LO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1500\\GeVcc$",
        "Tgamma1500GeV", 94563 , 1.00460294E-04 , 0},  // LO from excited quark

    // Q^2 --start--
    {"./REDUCE_DATA/1lepton4jets_TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets(Matchingdown)",
        "TTJetsMatchingdown", 5303083 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"./REDUCE_DATA/1lepton4jets_TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets(Matchingup)",
        "TTJetsMatchingup", 5306766 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"./REDUCE_DATA/1lepton4jets_TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets(Scaledown)",
        "TTJetsScaledown", 5317172 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"./REDUCE_DATA/1lepton4jets_TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets(Scaleup)",
        "TTJetsScaleup", 4963854 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay

    {"./REDUCE_DATA/1lepton4jets_TTJets_powheg_v2_scaledown.root",
        "$\\ttbar+$jets(POWHEG:scale down)",
        "TTJets_powhegScaledown", 14998602 , ttbarXsec , ttbarXsecUnc/ttbarXsec},  
    {"./REDUCE_DATA/1lepton4jets_TTJets_powheg_v2_scaleup.root",
        "$\\ttbar+$jets(POWHEG:scale up)",
        "TTJets_powhegScaleup", 14998717 , ttbarXsec , ttbarXsecUnc/ttbarXsec},  

    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_matchingdown_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}(Matchingdown)$",
        "DYJetsToLL_TuneZ2_M-50Matchingdown", 2028774  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_matchingup_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}(Matchingup)$",
        "DYJetsToLL_TuneZ2_M-50Matchingup",   1981995  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_scaledown_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}(Scaledown)$",
        "DYJetsToLL_TuneZ2_M-50Scaledown",    1931192  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_scaleup_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}(Scaleup)$",
        "DYJetsToLL_TuneZ2_M-50Scaleup",      2146505  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    // Q^2 --end--
    // more statistics --start--
    ///*
    {"./REDUCE_DATA/1lepton4jets_TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "$\\ttbar+$jets(FullLeptMGDecays)",
        "TTJets_FullLeptMGDecays", 11964351 , 234.2*(13.43/(53.4+53.2+13.43)) , 0},  
    {"./REDUCE_DATA/1lepton4jets_TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets(HadronicMGDecays)",
        "TTJets_HadronicMGDecays", 10087972 , 234.2*(53.4/(53.4+53.2+13.43)) , 0},  
    {"./REDUCE_DATA/1lepton4jets_TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1.root",
        "$\\ttbar+$jets(SemiLeptMGDecays)",
        "TTJets_SemiLeptMGDecays", 22656256 , 234.2*(53.2/(53.4+53.2+13.43)) , 0},  

    {"./REDUCE_DATA/1lepton4jets_TTJets_MSDecays*.root",
        "$\\ttbar+$jets(MSDecays)",
        "TTJets_MSDecays", 53308506 , ttbarXsec , ttbarXsecUnc/ttbarXsec},  
    {"./REDUCE_DATA/1lepton4jets_TTJets_powheg_v2_part*.root",
        "$\\ttbar+$jets(POWHEG)",
        "TTJets_powheg", 21461120 , ttbarXsec , ttbarXsecUnc/ttbarXsec},  
    {"./REDUCE_DATA/1lepton4jets_TT_700_1000_powheg*.root",
        "$\\ttbar+$jets(POWHEG:TT_700_1000)",
        "TT_700_1000_powheg", 3058249 , ttbarXsec*0.074 , ttbarXsecUnc/ttbarXsec},  
    {"./REDUCE_DATA/1lepton4jets_TT_1000_inf_powheg*.root",
        "$\\ttbar+$jets(POWHEG:TT_1000_inf)",
        "TT_1000_inf_powheg", 1159788 , ttbarXsec*0.014 , ttbarXsecUnc/ttbarXsec},  

    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "DYJetsToLL_PtZ-100",
        "DYJetsToLL_PtZ-100", 2626905 , 34.1 , 0},  
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DYJetsToLL_PtZ-50To70",
        "DYJetsToLL_PtZ-50To70", 4876674 , 93.8 , 0},  
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "DYJetsToLL_PtZ-70To100",
        "DYJetsToLL_PtZ-70To100", 776772 , 52.31, 0},  

    {"./REDUCE_DATA/1lepton4jets_DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DYToEE_M-20",
        "DYToEE_M-20", 34265939 , 1871. , 0},  
    {"./REDUCE_DATA/1lepton4jets_DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DYToMuMu_M-20",
        "DYToMuMu_M-20", 47006919 , 1871. , 0},  


    {"./REDUCE_DATA/1lepton4jets_ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZJetsTo2L2Q",
        "ZZJetsTo2L2Q", 1790107 , 0.91 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo2e2mu_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo2e2mu",
        "ZZTo2e2mu", 1390178 , 0.1767 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo2e2tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo2e2tau",
        "ZZTo2e2tau", 723132 , 0.1767 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo2mu2tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo2mu2tau",
        "ZZTo2mu2tau", 721675 , 0.1767 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo4e_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo4e",
        "ZZTo4e", 1237568 , 0.07691 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo4mu_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo4mu",
        "ZZTo4mu", 1469235 , 0.07691 , 0},  
    {"./REDUCE_DATA/1lepton4jets_ZZTo4tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZTo4tau",
        "ZZTo4tau", 823403 , 0.07691 , 0},  


    {"./REDUCE_DATA/1lepton4jets_GluGluToHToZZTo4L_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ggToH2ZZTo4l",
        "ggToH2ZZTo4l", 290510 , 18.97*3.02*pow(10.,-4), 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_GluGluToHToZZTo2L2Q_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ggToH2ZZTo2l2q"  ,
        "ggToH2ZZTo2l2q"  , 298986, 18.97*3.02*pow(10.,-4)*2 , 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_VBF_HToZZTo2L2Q_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "VBFToH2ZZTo2l2q" ,
        "VBFToH2ZZTo2l2q" , 49674 , 1.568*3.02*pow(10.,-4)*2 , 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_TTbarH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "TTbarH_HToZZTo4L",
        "TTbarH_HToZZTo4L", 97650 , 0.13 , 0.50},   // okay (HIG-13-015)

    {"./REDUCE_DATA/1lepton4jets_VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "VBFToH2ZZTo4l",
        "VBFToH2ZZTo4l", 49614 , 1.568*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_WH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "WHToH2ZZTo4l",
        "WHToH2ZZTo4l", 151309 , 0.6860*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_ZH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "ZHToH2ZZTo4l",
        "ZHToH2ZZTo4l", 140098 , 0.4050*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    //*/
    // more statistics --end--

    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-15to30",
        "QCD_Pt-15to30", 361396 , 988287420 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-30to50",
        "QCD_Pt-30to50", 392963 , 66285328 , 1},  // LO // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-120to170",
        "QCD_Pt-120to170", 823243 , 156293.3 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-170to300",
        "QCD_Pt-170to300", 235702 , 34138.15 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-300to470",
        "QCD_Pt-300to470", 314037 , 1759.549 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-470to600",
        "QCD_Pt-470to600", 263895 , 113.8791 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-50to80",
        "QCD_Pt-50to80", 385008 , 8148778.0 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-600to800",
        "QCD_Pt-600to800", 263403 , 26.9921 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-80to120",
        "QCD_Pt-80to120", 843700 , 1033680.0 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-800to1000",
        "QCD_Pt-800to1000", 309388 , 3.550036 , 1},  // LO  // okay
    {"./REDUCE_DATA/1lepton4jets_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\ttbar+$jets",
        //"TTJets", 6866660 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
        "TTJets", 6866660 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-006   // okay
    {"./REDUCE_DATA/1lepton4jets_T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\cPqt$ (s-channel)",
        "T_TuneZ2_s-channel", 259180 , 3.79 , 0.07/3.79}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\cPqt$ (t-channel)",
        "T_TuneZ2_t-channel", 3728236 , 56.4 , 2.1/56.4},  // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\cPqt$ ($\\cPqt\\PW$-channel)",
        "T_TuneZ2_tW-channel", 495571 , 11.1 , 0.3/11.1},    // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\overline{\\cPqt}$ (s-channel)",
        "Tbar_TuneZ2_s-channel", 139604 , 1.76 , 0.01/1.76},  // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\overline{\\cPqt}$ (t-channel)",
        "Tbar_TuneZ2_t-channel", 1910249 , 30.7 , 0.7/30.7}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Single $\\overline{\\cPqt}$ ($\\cPqt\\PW$-channel)",
        "Tbar_TuneZ2_tW-channel", 491471 , 11.1 , 0.3/11.1}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}$ (M-50)",
        "DYJetsToLL_TuneZ2_M-50", 19028293  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay

    // need to update unc on xsec
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-10To50filter_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ+{\\rm jets}$ (M-10to50)",
        "DYJetsToLL_TuneZ2_M-10To50", 6855679/0.069  , 11050.0 , 0.0701/3.5121},  // SMP-12-011 // okay
    //{"./REDUCE_DATA/1lepton4jets_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
    //    "GJet_Pt40", 4781024/0.05387 , 8884.0 ,0.},   // CTEQ // okay

    // need to update xsec
    {"./REDUCE_DATA/1lepton4jets_GJets_HT-40To100_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V19-v1.root",
        "$\\gamma+{\\rm jets}$ (HT-40to100)",
        "GJets_HT-40To100", 19800506 , 20930.0 ,0.01},   // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_GJets_HT-100To200_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V19-v1.root",
        "$\\gamma+{\\rm jets}$ (HT-100to200)",
        "GJets_HT-100To200", 9246629 , 5331 , 36.11/5331}, // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1/1/1/1/1/1/1.html
    {"./REDUCE_DATA/1lepton4jets_GJets_HT-200To400_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\gamma+{\\rm jets}$ (HT-200to400)",
        "GJets_HT-200To400", 53029964 , 960.5 ,7./960.5},   // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1.html
    {"./REDUCE_DATA/1lepton4jets_GJets_HT-400ToInf_8TeV-madgraph_v3_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "$\\gamma+{\\rm jets}$ (HT-400toInf)",
        "GJets_HT-400ToInf", 38972947 , 107.5 ,0.9/102.9},  // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1.html 

    {"./REDUCE_DATA/1lepton4jets_DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\gamma\\gamma+{\\rm jets}$",
        "DiPhotonJets", 1152389 , 75.39 ,0.01},   // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\PW\\gamma$",
        "WGToLNuG", 4736074 , 461.6 , 0.11},   // LO // okay    // unc from EWK-11-009
    {"./REDUCE_DATA/1lepton4jets_ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\cPZ\\gamma$",
        "ZGToLLG", 5277997 , 132.6 , 0.049},   // LO // okay // NLO : 159.12    // unc from EWK-11-009
    {"./REDUCE_DATA/1lepton4jets_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "$\\PW+{\\rm jets}$",
        "WJetsToLNu", 54487115 , 35640 , 666.1/35640 },   // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\PW\\PW",
        "WW", 8345262 , 56 , 2.3/56},  // MCFM 6.6

    // need to update xsec
    {"./REDUCE_DATA/1lepton4jets_WWGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "$\\PW\\PW\\gamma+{\\rm jets}$",
        "WWGJets", 214319 , 0.528 , 6.26/69.9},  // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_WWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\PW\\PW\\PW + jets",
        "WWWJets", 159718 , 0.08058 , 0.047},  // MCFM 6.6

    {"./REDUCE_DATA/1lepton4jets_WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\PW\\cPZ",
        "WZ", 9681965 , 33.21 ,0.51/33.21},  // CTEQ  // okay
    {"./REDUCE_DATA/1lepton4jets_ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\cPZ\\cPZ",
        "ZZ", 9560761 , 8.4 ,1.22/8.4},   // CTEQ // okay

    {"./REDUCE_DATA/1lepton4jets_TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\ttbar\\PW",
        "ttW", 194231 , 0.232 , 0.067/0.232}, // okay
    // need to update xsec
    {"./REDUCE_DATA/1lepton4jets_TTWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\ttbar\\PW\\PW",
        "ttWW", 214805 ,0.002037  , 0.067/0.232}, // okay
    //{"./REDUCE_DATA/1lepton4jets_TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
    //    "\\ttbar$\\gamma$",
    //    "ttG", 71049 , ttbarXsec*1.07*pow(10.,-2.), 0.27/1.07},   // TOP-13-011
    {"./REDUCE_DATA/1lepton4jets_TTGamma_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_RD1_START53_V7N-v1.root",
        "\\ttbar$\\gamma$",
        //"ttG", 872053 , 2.073, fabs(1.-ttbarXsec*1.07*pow(10.,-2.) - 2.073)},   // TOP-13-011
        "ttG", 872053 ,  ttbarXsec*1.07*pow(10.,-2.), 0.27/1.07},   // TOP-13-011
    {"./REDUCE_DATA/1lepton4jets_TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\ttbar H",
        "ttH", 925373 , 0.13 , 0.50},   // okay (HIG-13-015)


    {"./REDUCE_DATA/1lepton4jets_TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "\\ttbar\\cPZ",
        "ttZ", 208118 , 0.2057 ,0.024/0.2057}    // okay
};

struct sample_t SAMPLE_PDFnoSelection[samples_order_size] = {
    {"./REDUCE_DATA/PhotonData/1lepton4jets_ElectronHad_Run2012A-22Jan2013-v1_190645-193621.root",
        "Data",
        "data", LUMINOSITY , 1 , 0},
    //{"./REDUCE_DATA/1lepton4jets_Tgamma450GeV.root",
    //    "Tgamma450GeV", 28546 , 11.984291 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 500\\GeVcc$",
        "Tgamma500GeV", 97662 , 5.087496 , 0},  // LO from excited quark  // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-550_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 550\\GeVcc$",
        "Tgamma550GeV", 98650 , 2.315665 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-600_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 600\\GeVcc$",
        "Tgamma600GeV", 99408 , 1.115097 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-650_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 650\\GeVcc$",
        "Tgamma650GeV", 101131 , 0.562451 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-700_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 700\\GeVcc$",
        "Tgamma700GeV", 98736 , 0.294853 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-750_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 750\\GeVcc$",
        "Tgamma750GeV", 99847 , 0.159653 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-800_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 800\\GeVcc$",
        "Tgamma800GeV", 97373 , 0.088853 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-850_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 850\\GeVcc$",
        "Tgamma850GeV", 98312 , 5.06221252E-02 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-900_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 900\\GeVcc$",
        "Tgamma900GeV", 96254 , 2.94248249E-02 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-950_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 950\\GeVcc$",
        "Tgamma950GeV", 97617 , 1.74004431E-02 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1000\\GeVcc$",
        "Tgamma1000GeV", 95614 , 1.04431023E-02 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1100\\GeVcc$",
        "Tgamma1100GeV", 94642 , 3.90115915E-03 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1200_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1200\\GeVcc$",
        "Tgamma1200GeV", 93570 , 1.51328287E-03 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1300_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1300\\GeVcc$",
        "Tgamma1300GeV", 95859 , 6.03060062E-04 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1400_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1400\\GeVcc$",
        "Tgamma1400GeV", 95396 , 2.44765697E-04 , 0},  // LO from excited quark
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TprimeTprimeToTgammaTgammainc_M-1500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) = 1500\\GeVcc$",
        "Tgamma1500GeV", 94563 , 1.00460294E-04 , 0},  // LO from excited quark

    // Q^2 --start--
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_matchingdown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets(Matchingdown)",
        "TTJetsMatchingdown", 5303083 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_matchingup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets(Matchingup)",
        "TTJetsMatchingup", 5306766 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_scaledown_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets(Scaledown)",
        "TTJetsScaledown", 5317172 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_scaleup_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets(Scaleup)",
        "TTJetsScaleup", 4963854 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-007   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-50_matchingdown_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}(Matchingdown)$",
        "DYJetsToLL_TuneZ2_M-50Matchingdown", 2028774  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-50_matchingup_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}(Matchingup)$",
        "DYJetsToLL_TuneZ2_M-50Matchingup",   1981995  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-50_scaledown_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}(Scaledown)$",
        "DYJetsToLL_TuneZ2_M-50Scaledown",    1931192  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-50_scaleup_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}(Scaleup)$",
        "DYJetsToLL_TuneZ2_M-50Scaleup",      2146505  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    // Q^2 --end--
    // more statistics --start--
    ///*
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "$\\ttbar+$jets(FullLeptMGDecays)",
        "TTJets_FullLeptMGDecays", 11964351 , 234.2*(13.43/(53.4+53.2+13.43)) , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets(HadronicMGDecays)",
        "TTJets_HadronicMGDecays", 10087972 , 234.2*(53.4/(53.4+53.2+13.43)) , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/*root",
        "$\\ttbar+$jets(SemiLeptMGDecays)",
        "TTJets_SemiLeptMGDecays", 22656256 , 234.2*(53.2/(53.4+53.2+13.43)) , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "DYJetsToLL_PtZ-100",
        "DYJetsToLL_PtZ-100", 2626905 , 34.1 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "DYJetsToLL_PtZ-50To70",
        "DYJetsToLL_PtZ-50To70", 4876674 , 93.8 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "DYJetsToLL_PtZ-70To100",
        "DYJetsToLL_PtZ-70To100", 776772 , 52.31, 0},  

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "DYToEE_M-20",
        "DYToEE_M-20", 34265939 , 1871. , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "DYToMuMu_M-20",
        "DYToMuMu_M-20", 47006919 , 1871. , 0},  


    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZJetsTo2L2Q",
        "ZZJetsTo2L2Q", 1790107 , 0.91 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo2e2mu_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo2e2mu",
        "ZZTo2e2mu", 1390178 , 0.1767 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo2e2tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo2e2tau",
        "ZZTo2e2tau", 723132 , 0.1767 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo2mu2tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo2mu2tau",
        "ZZTo2mu2tau", 721675 , 0.1767 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo4e_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo4e",
        "ZZTo4e", 1237568 , 0.07691 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo4mu_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo4mu",
        "ZZTo4mu", 1469235 , 0.07691 , 0},  
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZTo4tau_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ZZTo4tau",
        "ZZTo4tau", 823403 , 0.07691 , 0},  


    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GluGluToHToZZTo4L_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ggToH2ZZTo4l",
        "ggToH2ZZTo4l", 290510 , 18.97*3.02*pow(10.,-4), 0.50},   // okay (HIG-13-015)
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GluGluToHToZZTo2L2Q_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "ggToH2ZZTo2l2q"  ,
        "ggToH2ZZTo2l2q"  , 298986, 18.97*3.02*pow(10.,-4)*2 , 0.50},   // okay (HIG-13-015)
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/VBF_HToZZTo2L2Q_M-125_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "VBFToH2ZZTo2l2q" ,
        "VBFToH2ZZTo2l2q" , 49674 , 1.568*3.02*pow(10.,-4)*2 , 0.50},   // okay (HIG-13-015)
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTbarH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "TTbarH_HToZZTo4L",
        "TTbarH_HToZZTo4L", 97650 , 0.13 , 0.50},   // okay (HIG-13-015)

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/VBF_HToZZTo4L_M-126_8TeV-powheg-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "VBFToH2ZZTo4l",
        "VBFToH2ZZTo4l", 49614 , 1.568*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "WHToH2ZZTo4l",
        "WHToH2ZZTo4l", 151309 , 0.6860*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZH_HToZZTo4L_M-126_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "ZHToH2ZZTo4l",
        "ZHToH2ZZTo4l", 140098 , 0.4050*3.02*pow(10.,-4) , 0.50},   // okay (HIG-13-015)
    //*/
    // more statistics --end--

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "QCD_Pt-15to30",
        "QCD_Pt-15to30", 361396 , 988287420 , 1},  // LO    // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-30to50",
        "QCD_Pt-30to50", 392963 , 66285328 , 1},  // LO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "QCD_Pt-120to170",
        "QCD_Pt-120to170", 823243 , 156293.3 , 1},  // LO   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-170to300",
        "QCD_Pt-170to300", 235702 , 34138.15 , 1},  // LO   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-300to470",
        "QCD_Pt-300to470", 314037 , 1759.549 , 1},  // LO   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-470to600",
        "QCD_Pt-470to600", 263895 , 113.8791 , 1},  // LO   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-50to80",
        "QCD_Pt-50to80", 385008 , 8148778.0 , 1},  // LO    // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-600to800",
        "QCD_Pt-600to800", 263403 , 26.9921 , 1},  // LO    // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "QCD_Pt-80to120",
        "QCD_Pt-80to120", 843700 , 1033680.0 , 1},  // LO   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "QCD_Pt-800to1000",
        "QCD_Pt-800to1000", 309388 , 3.550036 , 1},  // LO  // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\ttbar+$jets",
        //"TTJets", 6866660 , ttbarXsec , 0.049},    // TOP-12-007   // okay
        "TTJets", 6866660 , ttbarXsec , ttbarXsecUnc/ttbarXsec},    // TOP-12-006   // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\cPqt$ (s-channel)",
        "T_TuneZ2_s-channel", 259180 , 3.79 , 0.07/3.79}, // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\cPqt$ (t-channel)",
        "T_TuneZ2_t-channel", 3728236 , 56.4 , 2.1/56.4},  // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\cPqt$ ($\\cPqt\\PW$-channel)",
        "T_TuneZ2_tW-channel", 495571 , 11.1 , 0.3/11.1},    // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\overline{\\cPqt}$ (s-channel)",
        "Tbar_TuneZ2_s-channel", 139604 , 1.76 , 0.01/1.76},  // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\overline{\\cPqt}$ (t-channel)",
        "Tbar_TuneZ2_t-channel", 1910249 , 30.7 , 0.7/30.7}, // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "Single $\\overline{\\cPqt}$ ($\\cPqt\\PW$-channel)",
        "Tbar_TuneZ2_tW-channel", 491471 , 11.1 , 0.3/11.1}, // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}$ (M-50)",
        "DYJetsToLL_TuneZ2_M-50", 19028293  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay

    // need to update unc on xsec
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DYJetsToLL_M-10To50filter_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ+{\\rm jets}$ (M-10to50)",
        "DYJetsToLL_TuneZ2_M-10To50", 6855679/0.069  , 11050.0 , 0.0701/3.5121},  // SMP-12-011 // okay
    //{"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
    //    "GJet_Pt40", 4781024/0.05387 , 8884.0 ,0.},   // CTEQ // okay

    // need to update xsec
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GJets_HT-40To100_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V19-v1/*root",
        "$\\gamma+{\\rm jets}$ (HT-40to100)",
        "GJets_HT-40To100", 19800506 , 20930.0 ,0.01},   // CTEQ // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GJets_HT-100To200_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V19-v1/*root",
        "$\\gamma+{\\rm jets}$ (HT-100to200)",
        "GJets_HT-100To200", 9246629 , 5331 , 36.11/5331}, // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1/1/1/1/1/1/1.html
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GJets_HT-200To400_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\gamma+{\\rm jets}$ (HT-200to400)",
        "GJets_HT-200To400", 53029964 , 960.5 ,7./960.5},   // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1.html
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/GJets_HT-400ToInf_8TeV-madgraph_v3_Summer12_DR53X-PU_S10_START53_V7C-v1/*root",
        "$\\gamma+{\\rm jets}$ (HT-400toInf)",
        "GJets_HT-400ToInf", 38972947 , 107.5 ,0.9/102.9},  // https://hypernews.cern.ch/HyperNews/CMS/get/generators/2135/1.html 

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\gamma\\gamma+{\\rm jets}$",
        "DiPhotonJets", 1152389 , 75.39 ,0.01},   // CTEQ // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\PW\\gamma$",
        "WGToLNuG", 4736074 , 461.6 , 0.11},   // LO // okay    // unc from EWK-11-009
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\cPZ\\gamma$",
        "ZGToLLG", 5277997 , 132.6 , 0.049},   // LO // okay // NLO : 159.12    // unc from EWK-11-009
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2/*root",
        "$\\PW+{\\rm jets}$",
        "WJetsToLNu", 54487115 , 35640 , 666.1/35640 },   // NNLO // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\PW\\PW",
        "WW", 8345262 , 56 , 2.3/56},  // MCFM 6.6

    // need to update xsec
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WWGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "$\\PW\\PW\\gamma+{\\rm jets}$",
        "WWGJets", 214319 , 0.528 , 6.26/69.9},  // CTEQ // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\PW\\PW\\PW + jets",
        "WWWJets", 159718 , 0.08058 , 0.047},  // MCFM 6.6

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\PW\\cPZ",
        "WZ", 9681965 , 33.21 ,0.51/33.21},  // CTEQ  // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\cPZ\\cPZ",
        "ZZ", 9560761 , 8.4 ,1.22/8.4},   // CTEQ // okay

    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\ttbar\\PW",
        "ttW", 194231 , 0.232 , 0.067/0.232}, // okay
    // need to update xsec
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTWWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\ttbar\\PW\\PW",
        "ttWW", 214805 ,0.002037  , 0.067/0.232}, // okay
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\ttbar$\\gamma$",
        //"ttG", 71049 , ttbarXsec*1.07*pow(10.,-2.), 0.27/1.07},   // TOP-13-011
        "ttG", 872053 ,  ttbarXsec*1.07*pow(10.,-2.), 0.27/1.07},   // TOP-13-011
    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\ttbar H",
        "ttH", 925373 , 0.13 , 0.50},   // okay (HIG-13-015)


    {"/wk3/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1/*root",
        "\\ttbar\\cPZ",
        "ttZ", 208118 , 0.2057 ,0.024/0.2057}    // okay
};
#endif
