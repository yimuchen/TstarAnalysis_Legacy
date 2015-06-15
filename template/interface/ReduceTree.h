void ReduceTree(TTree *root)
{
    root->SetBranchStatus("*"        ,0);

    //root->SetBranchStatus("EvtInfo.*"        ,1);
    //root->SetBranchStatus("GenInfo.*"			   ,1);
    //root->SetBranchStatus("PhotonInfo*"			   ,1);
    //root->SetBranchStatus("PFLepInfo.*"			   ,1);
    //root->SetBranchStatus("PFJetInfo.*"			 ,1);
    //root->SetBranchStatus("VertexInfo.*"        ,1);


    ///*
    root->SetBranchStatus("EvtInfo.PDFid1",1);
    root->SetBranchStatus("EvtInfo.PDFid2",1);
    root->SetBranchStatus("EvtInfo.PDFx1",1);
    root->SetBranchStatus("EvtInfo.PDFx2",1);
    root->SetBranchStatus("EvtInfo.PDFscale",1);
    root->SetBranchStatus("EvtInfo.PDFv1",1);
    root->SetBranchStatus("EvtInfo.PDFv2",1);

	root->SetBranchStatus("EvtInfo.RhoPU",1);
	root->SetBranchStatus("EvtInfo.RunNo",1);
	root->SetBranchStatus("EvtInfo.LumiNo",1);
	root->SetBranchStatus("EvtInfo.TrueIT",1);
	root->SetBranchStatus("EvtInfo.PFMET",1);
	root->SetBranchStatus("EvtInfo.PFMETPhi",1);
	root->SetBranchStatus("EvtInfo.McFlag",1);
	root->SetBranchStatus("EvtInfo.EvtNo",1);

	root->SetBranchStatus("GenInfo.Size",1);
	root->SetBranchStatus("GenInfo.Phi",1);
	root->SetBranchStatus("GenInfo.Eta",1);
	root->SetBranchStatus("GenInfo.Pt",1);
	root->SetBranchStatus("GenInfo.Status",1);
	root->SetBranchStatus("GenInfo.PdgID",1);
	root->SetBranchStatus("GenInfo.PhotonFlag",1);
	root->SetBranchStatus("GenInfo.Da*",1);

	root->SetBranchStatus("VertexInfo.isValid",1);
	root->SetBranchStatus("VertexInfo.Size",1);
	root->SetBranchStatus("VertexInfo.Type",1);
	root->SetBranchStatus("VertexInfo.Ndof",1);
	root->SetBranchStatus("VertexInfo.Rho",1);
	root->SetBranchStatus("VertexInfo.z",1);
	root->SetBranchStatus("VertexInfo.y",1);
	root->SetBranchStatus("VertexInfo.x",1);

	root->SetBranchStatus("PFLepInfo.Size",1);
	root->SetBranchStatus("PFLepInfo.Charge",1);
	root->SetBranchStatus("PFLepInfo.LeptonType",1);
	root->SetBranchStatus("PFLepInfo.Pt",1);
	root->SetBranchStatus("PFLepInfo.Eta",1);
	root->SetBranchStatus("PFLepInfo.Phi",1);
	root->SetBranchStatus("PFLepInfo.Px",1);
	root->SetBranchStatus("PFLepInfo.Py",1);
	root->SetBranchStatus("PFLepInfo.Pz",1);
	root->SetBranchStatus("PFLepInfo.GenPt",1);
	root->SetBranchStatus("PFLepInfo.GenEta",1);
	root->SetBranchStatus("PFLepInfo.GenPhi",1);
    root->SetBranchStatus("PFLepInfo.EgammaMVANonTrig*",1);
    root->SetBranchStatus("PFLepInfo.EgammaMVATrig*",1);
    root->SetBranchStatus("PFLepInfo.isGoodMuonTMOneStationTight*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdTRIGGERTIGHT*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdTRIGGERWP70*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdVETO*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdLOOSE*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdMEDIUM*",1);
    root->SetBranchStatus("PFLepInfo.EgammaCutBasedEleIdTIGHT*",1);


	root->SetBranchStatus("PFLepInfo.MuInnerTrackDxy_PV",1);
	root->SetBranchStatus("PFLepInfo.MuInnerTrackDz",1);
	root->SetBranchStatus("PFLepInfo.MuInnerTrackNHits",1);
	root->SetBranchStatus("PFLepInfo.MuGlobalNormalizedChi2",1);
	root->SetBranchStatus("PFLepInfo.MuNPixelLayers",1);
	root->SetBranchStatus("PFLepInfo.MuNChambersMatchesSegment",1);
	root->SetBranchStatus("PFLepInfo.MuNMuonhits",1);
	root->SetBranchStatus("PFLepInfo.MuNTrackLayersWMeasurement",1);
	root->SetBranchStatus("PFLepInfo.ChargedHadronIsoR04",1);
	root->SetBranchStatus("PFLepInfo.NeutralHadronIsoR04",1);
	root->SetBranchStatus("PFLepInfo.PhotonIsoR04",1);
	root->SetBranchStatus("PFLepInfo.sumPUPtR04",1);
	root->SetBranchStatus("PFLepInfo.MuType",1);
	root->SetBranchStatus("PFLepInfo.Energy",1);
	root->SetBranchStatus("PFLepInfo.NeutralHadronIsoR03",1);
	root->SetBranchStatus("PFLepInfo.PhotonIsoR03",1);
	root->SetBranchStatus("PFLepInfo.ChargedHadronIsoR03",1);
	root->SetBranchStatus("PFLepInfo.EldeltaEta",1);
	root->SetBranchStatus("PFLepInfo.EldeltaPhi",1);
	root->SetBranchStatus("PFLepInfo.ElsigmaIetaIeta",1);
	root->SetBranchStatus("PFLepInfo.ElHadoverEm",1);
	root->SetBranchStatus("PFLepInfo.ElEcalE",1);
	root->SetBranchStatus("PFLepInfo.ElEoverP",1);
	root->SetBranchStatus("PFLepInfo.ElEcalE",1);
	root->SetBranchStatus("PFLepInfo.ElTrackDxy_PV",1);
	root->SetBranchStatus("PFLepInfo.ElTrackDz",1);
	root->SetBranchStatus("PFLepInfo.ElhasConv",1);
	root->SetBranchStatus("PFLepInfo.NumberOfExpectedInnerHits",1);
	root->SetBranchStatus("PFLepInfo.GenMCTag",1);

	root->SetBranchStatus("PFJetInfo.Pt",1);
	root->SetBranchStatus("PFJetInfo.Eta",1);
	root->SetBranchStatus("PFJetInfo.Phi",1);
	root->SetBranchStatus("PFJetInfo.GenPt",1);
	root->SetBranchStatus("PFJetInfo.GenEta",1);
	root->SetBranchStatus("PFJetInfo.GenPhi",1);
	root->SetBranchStatus("PFJetInfo.GenJet*",1);
	root->SetBranchStatus("PFJetInfo.TrackCountHiPurBJetTags",1);
	root->SetBranchStatus("PFJetInfo.CombinedSVMVABJetTags",1);
	root->SetBranchStatus("PFJetInfo.CEF",1);
	root->SetBranchStatus("PFJetInfo.CHF",1);
	root->SetBranchStatus("PFJetInfo.NEF",1);
	root->SetBranchStatus("PFJetInfo.NHF",1);
	root->SetBranchStatus("PFJetInfo.Py",1);
	root->SetBranchStatus("PFJetInfo.Px",1);
	root->SetBranchStatus("PFJetInfo.Pz",1);
	root->SetBranchStatus("PFJetInfo.Energy",1);
	root->SetBranchStatus("PFJetInfo.Size",1);
	root->SetBranchStatus("PFJetInfo.JetIDLOOSE",1);
	root->SetBranchStatus("PFJetInfo.GenMCTag",1);
	root->SetBranchStatus("PFJetInfo.Unc",1);

	root->SetBranchStatus("PhotonInfo*Size",1);
	root->SetBranchStatus("PhotonInfo*Pt",1);
	root->SetBranchStatus("PhotonInfo*Phi",1);
	root->SetBranchStatus("PhotonInfo*Eta",1);
	root->SetBranchStatus("PhotonInfo*HoverE",1);
	root->SetBranchStatus("PhotonInfo*SigmaIetaIeta",1);
	root->SetBranchStatus("PhotonInfo*hadTowOverEm",1);
	root->SetBranchStatus("PhotonInfo*phoPFChIsoDR03",1);
	root->SetBranchStatus("PhotonInfo*phoPFPhoIsoDR03",1);
	root->SetBranchStatus("PhotonInfo*phoPFNeuIsoDR03",1);
	root->SetBranchStatus("PhotonInfo*r9",1);
	root->SetBranchStatus("PhotonInfo*passelectronveto",1);
	root->SetBranchStatus("PhotonInfo*GenPt",1);
	root->SetBranchStatus("PhotonInfo*GenPhi",1);
	root->SetBranchStatus("PhotonInfo*GenEta",1);
//*/


}
