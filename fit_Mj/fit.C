void mj_fit_Landau_Gaus(TString file="ofile_WJets.root")
{
	gSystem->Load("libRooFit");
	using namespace RooFit;

	TFile * fp = new TFile(file);
	TTree * tree = (TTree *)fp->Get("otree");

	RooRealVar x("x","x",30,200);
	RooRealVar mean_landau("mean_landau","mean_landau",80,50,100);
	//RooRealVar mean_landau("mean_landau","mean_landau",40,0,100);
	RooRealVar c_landau("c_landau","c_landau",25,0,100);
	RooLandau landau("landau","landau",x,mean_landau,c_landau) ;

    RooRealVar mean_gaus("mean_gaus","mean_gaus",60,30,200);
    RooRealVar sigma_gaus("sigma_gaus","sigma_gaus",10,0.1,100);
    RooGaussian gaussian("gaussian","gaussian",x,mean_gaus,sigma_gaus);

    RooRealVar f1("f1","f1",0.5,0,1);

    RooAddPdf model("model","model",RooArgList(landau,gaussian),f1);

	float jet_mass_pr;
	float jet_massdrop_pr;
	float jet_pt_pr;

	tree->SetBranchAddress("jet_mass_pr",&jet_mass_pr);
	tree->SetBranchAddress("jet_massdrop_pr",&jet_massdrop_pr);
	tree->SetBranchAddress("jet_pt_pr",&jet_pt_pr);

	RooDataSet data("data","Dataset With X",x);
	for(Int_t j=0; j<tree->GetEntries(); j++){
		tree->GetEntry(j); 
		if(j%1000==0)
		cout<<"Processing the events: "<<j<<endl;
	//	cout<<jet_massdrop_pr<<endl;
	//	cout<<jet_mass_pr<<endl;
		if(jet_massdrop_pr<0.25 && jet_mass_pr <200 && jet_mass_pr >30 && jet_pt_pr > 200.)
		{
			x=jet_mass_pr;
			data.add(x);
		}
	}

	//landau.fitTo(data);
	model.fitTo(data);
	model.fitTo(data);
	//landau.fitTo(data,Range(30,200));
	RooPlot* frame = x.frame();
	data.plotOn(frame,Binning(100));
	model.plotOn(frame);
    model.plotOn(frame, Components(gaussian),LineStyle(kDashed));
    model.plotOn(frame, Components(landau),LineStyle(kDashed),LineColor(kRed));

    TString rlt_file=file;
    rlt_file.ReplaceAll(".root","_LG_rlt.png");
    TCanvas *cfit=new TCanvas("cfit","cfit");
	frame->Draw();
    cfit->Print(rlt_file);
    rlt_file.ReplaceAll(".png",".eps"); cfit->Print(rlt_file);
}


void mj_fit_Arg_Gaus(TString file="ofile_WJets.root")
{
	gSystem->Load("libRooFit");
	using namespace RooFit;

	TFile * fp = new TFile(file);
	TTree * tree = (TTree *)fp->Get("otree");

	RooRealVar x("x","x",30,200);
	RooRealVar m0_arg("m0_arg","m0_arg",80,30,150);
	RooRealVar c_arg("c_arg","c_arg",-2,-50,0);
	RooArgusBG arg("arg","arg",x,m0_arg,c_arg) ;

    RooRealVar mean_gaus("mean_gaus","mean_gaus",80,70,90);
    RooRealVar sigma_gaus("sigma_gaus","sigma_gaus",7,5,10);
    RooGaussian gaussian("gaussian","gaussian",x,mean_gaus,sigma_gaus);

    RooRealVar f1("f1","f1",0.1,0,1);

    RooAddPdf model("model","model",RooArgList(arg,gaussian),f1);

	float jet_mass_pr;
	float jet_massdrop_pr;
	float jet_pt_pr;

	tree->SetBranchAddress("jet_mass_pr",&jet_mass_pr);
	tree->SetBranchAddress("jet_massdrop_pr",&jet_massdrop_pr);
	tree->SetBranchAddress("jet_pt_pr",&jet_pt_pr);

	RooDataSet data("data","Dataset With X",x);
	for(Int_t j=0; j<tree->GetEntries(); j++){
		tree->GetEntry(j); 
		if(j%1000==0)
		cout<<"Processing the events: "<<j<<endl;
	//	cout<<jet_massdrop_pr<<endl;
	//	cout<<jet_mass_pr<<endl;
		if(jet_massdrop_pr<0.25 && jet_mass_pr <200 && jet_mass_pr >30 && jet_pt_pr > 200.)
		{
			x=jet_mass_pr;
			data.add(x);
		}
	}

	//landau.fitTo(data);
	model.fitTo(data);
	model.fitTo(data);
	RooPlot* frame = x.frame();
	data.plotOn(frame,Binning(100));
	model.plotOn(frame);
    model.plotOn(frame, Components(gaussian),LineStyle(kDashed));
    model.plotOn(frame, Components(arg),LineStyle(kDashed),LineColor(kRed));

    TString rlt_file=file;
    rlt_file.ReplaceAll(".root","_ARG_rlt.png");
    TCanvas *cfit=new TCanvas("cfit","cfit");
	frame->Draw();
    cfit->Print(rlt_file);
    rlt_file.ReplaceAll(".png",".eps"); cfit->Print(rlt_file);
}


void mj_fit_ErfExp(TString file="ofile_WJets.root")
{
	gSystem->Load("libRooFit");
	using namespace RooFit;

	TFile * fp = new TFile(file);
	TTree * tree = (TTree *)fp->Get("otree");

	RooRealVar x("x","x",30,200);
	RooRealVar c_ErfExp("c_ErfExp","c_ErfExp",-0.1,-10.,0.);
	RooRealVar off_ErfExp("off_ErfExp","off_ErfExp",50.,10.,100.);
	RooRealVar sigma_ErfExp("sigma_ErfExp","sigma_ErfExp",10.,0,100.);
    RooErfExpPdf erfExp("erfExp","erfExp",x,c_ErfExp,off_ErfExp,sigma_ErfExp);

    RooRealVar mean_gaus("mean_gaus","mean_gaus",85,75,95);
    RooRealVar sigma_gaus("sigma_gaus","sigma_gaus",5,0.1,100);
    RooGaussian gaussian("gaussian","gaussian",x,mean_gaus,sigma_gaus);

    RooRealVar f1("f1","f1",0.5,0,1);

    RooAddPdf model("model","model",RooArgList(erfExp,gaussian),f1);

	float jet_mass_pr;
	float jet_massdrop_pr;
	float jet_pt_pr;

	tree->SetBranchAddress("jet_mass_pr",&jet_mass_pr);
	tree->SetBranchAddress("jet_massdrop_pr",&jet_massdrop_pr);
	tree->SetBranchAddress("jet_pt_pr",&jet_pt_pr);

	RooDataSet data("data","Dataset With X",x);
	for(Int_t j=0; j<tree->GetEntries(); j++){
		tree->GetEntry(j); 
		if(j%1000==0)
		cout<<"Processing the events: "<<j<<endl;
	//	cout<<jet_massdrop_pr<<endl;
	//	cout<<jet_mass_pr<<endl;
		if(jet_massdrop_pr<0.25 && jet_mass_pr <200 && jet_mass_pr >30 && jet_pt_pr > 200.)
		{
			x=jet_mass_pr;
			data.add(x);
		}
	}

	//landau.fitTo(data);
	model.fitTo(data);
	model.fitTo(data);
	//landau.fitTo(data,Range(30,200));
	RooPlot* frame = x.frame();
	data.plotOn(frame,Binning(100));
	model.plotOn(frame);
    model.plotOn(frame, Components(gaussian),LineStyle(kDashed));
    model.plotOn(frame, Components(erfExp),LineStyle(kDashed),LineColor(kRed));

    TString rlt_file=file;
    rlt_file.ReplaceAll(".root","_ErfExp_rlt.png");
    TCanvas *cfit=new TCanvas("cfit","cfit");
	frame->Draw();
    cfit->Print(rlt_file);
    rlt_file.ReplaceAll(".png",".eps"); cfit->Print(rlt_file);
}

void mj_fit_CB(TString file="ofile_WJets.root")
{
	gSystem->Load("libRooFit");
	using namespace RooFit;

	TFile * fp = new TFile(file);
	TTree * tree = (TTree *)fp->Get("otree");

	RooRealVar x("x","x",30,200);

    RooRealVar mean_gaus("mean_gaus","mean_gaus",80,70,95);
    RooRealVar sigma_gaus("sigma_gaus","sigma_gaus",5,4,6);
    RooRealVar alpha_CB("alpha_CB","alpha_CB",2.,1,10);
    RooRealVar n_CB("n_CB","n_CB",3,0,10);
    RooCBShape model("cbshape","cbshape",x,mean_gaus,sigma_gaus,alpha_CB,n_CB);
    //RooGaussian gaussian("gaussian","gaussian",x,mean_gaus,sigma_gaus);


    //RooAddPdf model("model","model",RooArgList(erfExp,gaussian),f1);

	float jet_mass_pr;
	float jet_massdrop_pr;
	float jet_pt_pr;

	tree->SetBranchAddress("jet_mass_pr",&jet_mass_pr);
	tree->SetBranchAddress("jet_massdrop_pr",&jet_massdrop_pr);
	tree->SetBranchAddress("jet_pt_pr",&jet_pt_pr);

	RooDataSet data("data","Dataset With X",x);
	for(Int_t j=0; j<tree->GetEntries(); j++){
		tree->GetEntry(j); 
		if(j%1000==0)
		cout<<"Processing the events: "<<j<<endl;
	//	cout<<jet_massdrop_pr<<endl;
	//	cout<<jet_mass_pr<<endl;
		if(jet_massdrop_pr<0.25 && jet_mass_pr <200 && jet_mass_pr >30 && jet_pt_pr > 200.)
		{
			x=jet_mass_pr;
			data.add(x);
		}
	}

	//landau.fitTo(data);
	model.fitTo(data);
	model.fitTo(data);
	//landau.fitTo(data,Range(30,200));
	RooPlot* frame = x.frame();
	data.plotOn(frame,Binning(100));
	model.plotOn(frame);
    //model.plotOn(frame, Components(gaussian),LineStyle(kDashed));
    //model.plotOn(frame, Components(erfExp),LineStyle(kDashed),LineColor(kRed));

    TString rlt_file=file;
    rlt_file.ReplaceAll(".root","_CB_rlt.png");
    TCanvas *cfit=new TCanvas("cfit","cfit");
	frame->Draw();
    cfit->Print(rlt_file);
    rlt_file.ReplaceAll(".png",".eps"); cfit->Print(rlt_file);
}




void fit(){
	using namespace RooFit;
    gROOT->ProcessLine(".L ./PDFs/RooErfExpPdf.cxx+");

    //mj_fit_ErfExp("ofile_TTbar.root");
    //mj_fit_ErfExp("ofile_WW.root");
    //mj_fit_Landau_Gaus("ofile_WW.root");
   mj_fit_Arg_Gaus("ofile_WW.root");
    //mj_fit_ErfExp("ofile_WJets.root");
   // mj_fit_Landau_Gaus("ofile_WJets.root");
}
