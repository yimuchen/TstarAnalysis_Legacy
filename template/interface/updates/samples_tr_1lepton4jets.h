#ifndef SAMPLES_tr_H
#define SAMPLES_tr_H

#include "ConstantNumbers.h"

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
    Tgamma950GeV,
    Tgamma1100GeV,
    Tgamma1300GeV,
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
    //GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DiPhotonJets_8TeV_madgraph_tarball_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WGToLNuG_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZGToLLG_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2,
    WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    samples_order_size
};

struct sample_t{
    char filename[256];
    char tag[64];
    double ngen;
    double xsec;
    double unc; // 1 : 100%, 0.11 : 11%, 0 : no uncertainty
};


struct sample_t SAMPLE[samples_order_size] = {
    {"./REDUCE_DATA/PhotonData/1lepton4jets_*_Run2012*.root",
        "data", LUMINOSITY , 1 , 0},
    //{"./REDUCE_DATA/1lepton4jets_Tgamma450GeV.root",
    //    "Tgamma450GeV", 28546 , 11.984291 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma500GeV", 97662 , 5.087496 , 0},  // NNLO from excited quark  // okay
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-550_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma550GeV", 98650 , 2.315665 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-600_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma600GeV", 99408 , 1.115097 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-650_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma650GeV", 101131 , 0.562451 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-700_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma700GeV", 98736 , 0.294853 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-750_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma750GeV", 92102 , 0.159653 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-800_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tgamma800GeV", 97373 , 0.088853 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-850_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "Tgamma850GeV", 98306 , 5.06221252E-02 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-950_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "Tgamma950GeV", 97604 , 1.74004431E-02 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "Tgamma1100GeV", 94618 , 3.90115915E-03 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_TprimeTprimeToTgammaTgammainc_M-1300_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1.root",
        "Tgamma1300GeV", 95833 , 6.03060062E-04 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-15to30", 361396 , 988287420 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-30to50", 404954 , 66285328 , 1},  // LO // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-120to170", 883075 , 156293.3 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-170to300", 287436 , 34138.15 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-300to470", 314037 , 1759.549 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-470to600", 322824 , 113.8791 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-50to80", 385008 , 8148778.0 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-600to800", 263403 , 26.9921 , 1},  // LO    // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-80to120", 827723 , 1033680.0 , 1},  // LO   // okay
    {"./REDUCE_DATA/1lepton4jets_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-800to1000", 301619 , 3.550036 , 1},  // LO  // okay
    {"./REDUCE_DATA/1lepton4jets_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "TTJets", 6775135 , 227 , 0.11},    // TOP-12-007   // okay
    {"./REDUCE_DATA/1lepton4jets_T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "T_TuneZ2_s-channel", 259180 , 3.79 , 0.07/3.79}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "T_TuneZ2_t-channel", 23712 , 56.4 , 2.1/56.4},  // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "T_TuneZ2_tW-channel", 479630 , 11.1 , 0.3/11.1},    // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_s-channel", 123657 , 1.76 , 0.01/1.76},  // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_t-channel", 1882319 , 30.7 , 0.7/30.7}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_tW-channel", 491471 , 11.1 , 0.3/11.1}, // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DYJetsToLL_TuneZ2_M-50", 23315575  , 1120.*3*3503.71/3351.97 , 0.0701/3.5121},  // SMP-12-011 // okay
    //{"./REDUCE_DATA/1lepton4jets_GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
    //    "GJet_Pt40", 4781024/0.05387 , 8884.0 ,0.},   // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DiPhotonJets", 1152389 , 75.39 ,0.},   // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "WGToLNuG", 4732092 , 461.6 , 0.11},   // LO // okay    // unc from EWK-11-009
    {"./REDUCE_DATA/1lepton4jets_ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZGToLLG", 6575596 , 132.6 , 0.049},   // LO // okay // NLO : 159.12    // unc from EWK-11-009
    {"./REDUCE_DATA/1lepton4jets_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "WJetsToLNu", 50341708 , 35640 , 666.1/35640 },   // NNLO // okay
    {"./REDUCE_DATA/1lepton4jets_WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "WW", 3019619 , 69.9 , 6.26/69.9},  // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "WZ", 8376430 , 33.21 ,0.51/33.21},  // CTEQ  // okay
    {"./REDUCE_DATA/1lepton4jets_ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZ", 9042572 , 8.4 ,1.22/8.4},   // CTEQ // okay
    {"./REDUCE_DATA/1lepton4jets_TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttW", 194231 , 0.232 , 0.067/0.232}, // okay
    {"./REDUCE_DATA/1lepton4jets_TTGJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttG", 71049 , 1.444 , 0.50},   // okay // assume Unc of 50%
    {"./REDUCE_DATA/1lepton4jets_TTH_Inclusive_M-125_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttH", 984751 , 0.13 , 0.50},   // okay (HIG-13-015)
    {"./REDUCE_DATA/1lepton4jets_TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttZ", 208118 , 0.2057 ,0.024/0.2057}    // okay
};

#endif
