#ifndef CorrectionOnTstar_H
#define CorrectionOnTstar_H
float CorrectionOnTstar(int index, int MCTruthChannel_){
float MCTruthChannelsDB[48][4] ={ // [sample][index, hardronic, semilepton, dilepton]
//TprimeTprimeToTgammaTgammainc_M-500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(50617,39396,7649)
{1, 0.85752, 1.10177, 1.41866},
//TprimeTprimeToTgammaTgammainc_M-550_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(51048,39856,7746)
{2, 0.85889, 1.10007, 1.41507},
//TprimeTprimeToTgammaTgammainc_M-600_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(51165,40328,7915)
{3, 0.86351, 1.09555, 1.39549},
//TprimeTprimeToTgammaTgammainc_M-650_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(52509,40666,7956)
{4, 0.85599, 1.10527, 1.41237},
//TprimeTprimeToTgammaTgammainc_M-700_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(51264,39722,7750)
{5, 0.85601, 1.10474, 1.41557},
//TprimeTprimeToTgammaTgammainc_M-750_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(51767,40210,7870)
{6, 0.85723, 1.10362, 1.40967},
//TprimeTprimeToTgammaTgammainc_M-800_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(50274,39482,7617)
{7, 0.86082, 1.09612, 1.42040},
//TprimeTprimeToTgammaTgammainc_M-850_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(50967,39555,7790)
{8, 0.85730, 1.10464, 1.40225},
//TprimeTprimeToTgammaTgammainc_M-900_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(49608,38892,7754)
{9, 0.86235, 1.09996, 1.37927},
//TprimeTprimeToTgammaTgammainc_M-950_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(50321,39602,7694)
{10, 0.86217, 1.09553, 1.40971},
//TprimeTprimeToTgammaTgammainc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(49281,38742,7591)
{11, 0.86230, 1.09687, 1.39952},
//TprimeTprimeToTgammaTgammainc_M-1100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(48697,38344,7601)
{12, 0.86377, 1.09699, 1.38347},
//TprimeTprimeToTgammaTgammainc_M-1200_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(48409,37672,7489)
{13, 0.85907, 1.10391, 1.38826},
//TprimeTprimeToTgammaTgammainc_M-1300_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(49502,38592,7765)
{14, 0.86065, 1.10396, 1.37167},
//TprimeTprimeToTgammaTgammainc_M-1400_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(49326,38449,7621)
{15, 0.85955, 1.10271, 1.39084},
//TprimeTprimeToTgammaTgammainc_M-1500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(48889,38272,7402)
{16, 0.85966, 1.09814, 1.41948},
//TprimeTprimeToTGluonTGluoninc_M-500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(40859,32046,6220)
{17, 0.86068, 1.09738, 1.41345},
//TprimeTprimeToTGluonTGluoninc_M-550_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(51340,40182,7858)
{18, 0.86032, 1.09922, 1.40522},
//TprimeTprimeToTGluonTGluoninc_M-600_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(50239,39080,7592)
{19, 0.85733, 1.10214, 1.41832},
//TprimeTprimeToTGluonTGluoninc_M-650_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(50168,38944,7625)
{20, 0.85700, 1.10400, 1.40965},
//TprimeTprimeToTGluonTGluoninc_M-700_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(49525,38626,7410)
{21, 0.85758, 1.09956, 1.43291},
//TprimeTprimeToTGluonTGluoninc_M-750_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(42046,32562,6370)
{22, 0.85597, 1.10528, 1.41249},
//TprimeTprimeToTGluonTGluoninc_M-800_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1(48912,37966,7400)
{23, 0.85667, 1.10365, 1.41559},
//TprimeTprimeToTGluonTGluoninc_M-850_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(51180,39594,7864)
{24, 0.85657, 1.10722, 1.39366},
//TprimeTprimeToTGluonTGluoninc_M-900_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(51255,40460,7982)
{25, 0.86450, 1.09515, 1.38780},
//TprimeTprimeToTGluonTGluoninc_M-950_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(40211,31931,6183)
{26, 0.86571, 1.09020, 1.40753},
//TprimeTprimeToTGluonTGluoninc_M-1000_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(19671,15465,3018)
{27, 0.86205, 1.09650, 1.40468},
//TprimeTprimeToTGluonTGluoninc_M-1100_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(47422,37155,7213)
{28, 0.86027, 1.09798, 1.41396},
//TprimeTprimeToTGluonTGluoninc_M-1200_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(48005,37697,7311)
{29, 0.86114, 1.09662, 1.41359},
//TprimeTprimeToTGluonTGluoninc_M-1300_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(18706,14550,2837)
{30, 0.85755, 1.10250, 1.41358},
//TprimeTprimeToTGluonTGluoninc_M-1400_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(19658,15524,3076)
{31, 0.86497, 1.09531, 1.38195},
//TprimeTprimeToTGluonTGluoninc_M-1500_TuneZ2star_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7C-v1(38853,30043,6021)
{32, 0.85699, 1.10829, 1.38251},
//TprimeTprimeToTGluonTgammainc_M_500_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(65420,50814,10060)
{33, 0.85800, 1.10463, 1.39490},
//TprimeTprimeToTGluonTgammainc_M-550_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(66330,51851,10164)
{34, 0.85998, 1.10012, 1.40305},
//TprimeTprimeToTGluonTgammainc_M_600_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(66672,51560,9942)
{35, 0.85443, 1.10485, 1.43246},
//TprimeTprimeToTGluonTgammainc_M-650_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(66836,51864,10261)
{36, 0.85756, 1.10512, 1.39645},
//TprimeTprimeToTGluonTgammainc_M_700_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(61464,47943,9278)
{37, 0.85821, 1.10024, 1.42134},
//TprimeTprimeToTGluonTgammainc_M-750_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(70038,54201,10699)
{38, 0.85628, 1.10648, 1.40136},
//TprimeTprimeToTGluonTgammainc_M_800_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(64762,50600,9975)
{39, 0.86015, 1.10090, 1.39612},
//TprimeTprimeToTGluonTgammainc_M-850_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(71411,55807,10912)
{40, 0.85969, 1.10006, 1.40650},
//TprimeTprimeToTGluonTgammainc_M_900_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(126035,98322,19148)
{41, 0.85869, 1.10071, 1.41300},
//TprimeTprimeToTGluonTgammainc_M-950_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(73416,57411,10838)
{42, 0.85761, 1.09669, 1.45235},
//TprimeTprimeToTGluonTgammainc_M-1000_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(73520,57219,11379)
{43, 0.85913, 1.10389, 1.38772},
//TprimeTprimeToTGluonTgammainc_M-1100_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(73627,57509,11331)
{44, 0.85999, 1.10102, 1.39702},
//TprimeTprimeToTGluonTgammainc_M-1200_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(77324,60525,11821)
{45, 0.86028, 1.09905, 1.40682},
//TprimeTprimeToTGluonTgammainc_M-1300_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(79012,61906,11935)
{46, 0.85980, 1.09738, 1.42301},
//TprimeTprimeToTGluonTgammainc_M-1400_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(71316,56027,10818)
{47, 0.86103, 1.09599, 1.41904},
//TprimeTprimeToTGluonTgammainc_M-1500_TuneZ2star_8TeV-madgraph_Summer12DR53X-PU_S10_START53_V19-v1(82607,64570,12615)
{48, 0.85972, 1.09987, 1.40743}
};   // MCTruthChannelsDB
for(int iscan=0;iscan<48;iscan++){
if(index == (int)MCTruthChannelsDB[iscan][0]){
return MCTruthChannelsDB[iscan][MCTruthChannel_+1];
}// if
}// for
return -1;
}// CorrectionOnTstar
#endif
