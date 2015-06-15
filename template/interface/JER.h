#ifndef JER_H
#define JER_H

// https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
// pT->max[0.,pTgen+c*(pT–pTgen)]
// mode : [-1, 0, 1] = [minus, normal, plus]
float JER(float eta,int mode){

    float scalingfactor = -1.;
    /*
    eta bin  Data/MC Ratio (factor +-stat. +syst.- syst.)
        0.0–0.5  1.052+-0.012+0.062-0.061
        0.5–1.1  1.057+-0.012+0.056-0.055
        1.1–1.7  1.096+-0.017+0.063-0.062
        1.7–2.3  1.134+-0.035+0.087-0.085
        2.3–5.0  1.288+-0.127+0.155-0.153
        */
    const int Nbins = 5;
    float EtaRegins[Nbins][2] = {
        {0.0,0.5},
        {0.5,1.1},
        {1.1,1.7},
        {1.7,2.3},
        {2.3,5.0}
    };

    float scale[Nbins] = {
        1.052,
        1.057,
        1.096,
        1.134,
        1.288
    };

    // +-stat. +syst.- syst.
    float unc[Nbins][3] = {
        {0.012,0.062,0.061},
        {0.012,0.056,0.055},
        {0.017,0.063,0.062},
        {0.035,0.087,0.085},
        {0.127,0.155,0.153}
    };

    for(int i=0;i<Nbins;i++){
        if( fabs(eta)>=EtaRegins[i][0] && fabs(eta)<EtaRegins[i][1] ){
            if(mode==0){
                scalingfactor = scale[i];
            }else if(mode ==-1){
                scalingfactor = scale[i] - sqrt(unc[i][0]*unc[i][0] + unc[i][2]*unc[i][2]);
            }else if(mode == 1){
                scalingfactor = scale[i] + sqrt(unc[i][0]*unc[i][0] + unc[i][1]*unc[i][1]);
            }
            break;
        }
    }

    return scalingfactor;
}

#endif
