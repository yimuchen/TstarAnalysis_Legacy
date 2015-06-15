#ifndef EffectiveAreaPhoton_H
#define EffectiveAreaPhoton_H



float EffectiveAreaPhoton(float eta, int cate){
	// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
	// EA charged hadrons	EA neutral hadrons	EA photons
	float EA[7][3] = {
		{0.012,0.030,0.148},
		{0.010,0.057,0.130},
		{0.014,0.039,0.112},
		{0.012,0.015,0.216},
		{0.016,0.024,0.262},
		{0.020,0.039,0.260},
		{0.012,0.072,0.266}
	};



	if (fabs(eta)<1.0)
		return EA[0][cate];
	if (1.0<=fabs(eta)&&fabs(eta)<1.479)
		return EA[1][cate];
	if (1.479<=fabs(eta)&&fabs(eta)<2.0)
		return EA[2][cate];
	if (2.0<=fabs(eta)&&fabs(eta)<2.2)
		return EA[3][cate];
	if (2.2<=fabs(eta)&&fabs(eta)<2.3)
		return EA[4][cate];
	if (2.3<=fabs(eta)&&fabs(eta)<2.4)
		return EA[5][cate];
	if (fabs(eta)>2.4)
		return EA[6][cate];

}


#endif
