/************
 * PRex Carbon  Crystal Ball Fit Functions
 ************/
#include <string>
#include <TChain.h>
#include <iostream>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TF1NormSum.h>

void CarbonCrystalBallFit(std::string HRSArm="L"){
	TChain *t=new TChain("T");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/prexLHRS_1676_-1_0_optics_check_1.root");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/prexLHRS_1676_-1_0_optics_check_2.root");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/prexLHRS_1676_-1_0_optics_check_3.root");

	std::string basicCuts =
			Form(   "%s.vdc.u1.nclust==1 && %s.vdc.v1.nclust==1 && %s.vdc.u2.nclust==1 && %s.vdc.v2.nclust==1 &&  %s.gold.dp<1 && %s.gold.dp > -0.1 && fEvtHdr.fEvtType==1",
					HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(),HRSArm.c_str());

//	std::string fiducials = Form("abs(%s.gold.th)<0.004 && abs(%s.gold.ph)<0.002",
//			HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

	// select the large hole on the left side
//	std::string fiducials = Form("(%s.gold.th)<0.022 &&(%s.gold.th)>0.014 && (%s.gold.ph)<-0.007 && (%s.gold.ph)>-0.011",
//				HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

	// select the righ bottom hole
	std::string fiducials = Form("(%s.gold.th)<-0.02 &&(%s.gold.th)>-0.033 && (%s.gold.ph)<0.008 && (%s.gold.ph)>0.004",
				HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

	std::cout<<"basic cuts : "<<basicCuts.c_str()<<std::endl;
	std::cout<<"fiducial cuts : "<<fiducials.c_str()<<std::endl;

	TCanvas *c0=new TCanvas("c0","c0",1600,1600);
	c0->Divide(3,2);
	c0->cd(1);

	TH2F *htgthph = new TH2F("htgthph","Target theta-phi (theta vertical)",500,-0.03,0.035,500,-0.07,0.07);
	t->Project(htgthph->GetName(),Form("%s.gold.th:%s.gold.ph",HRSArm.c_str(),HRSArm.c_str()),Form("%s && %s",basicCuts.c_str(),fiducials.c_str()));
	htgthph->DrawCopy("colz");

	c0->cd(3);
	TH1F *goldYhist=new TH1F("gold.y","gold.y",200,-0.03,0.03);
    t->Project(goldYhist->GetName(),Form("%s.gold.ph",HRSArm.c_str()),Form("%s && %s",basicCuts.c_str(),fiducials.c_str()));
	goldYhist->DrawCopy();

	c0->cd(5);
	TH1F *trXhist= new TH1F("tr.x","tr.x",1000,-0.15,0.0);
	t->Project(trXhist->GetName(),Form("%s.tr.x",HRSArm.c_str()),Form("%s && %s",basicCuts.c_str(),fiducials.c_str()));
	trXhist->DrawCopy();

	// fit the two peak with the gaus functions
	double_t fHtrxGausPar[3];
	TF1 *fHtrxGaus=new TF1("trxHGaus","gaus",-0.105,-0.085);
	trXhist->Fit("trxHGaus","R","ep",-0.105,-0.085);
	fHtrxGaus->GetParameters(fHtrxGausPar);

	double_t fHtrxCrystalPar[5];
	TF1 *fHtrxCrystal=new TF1("trxHCrystal","crystalball",-0.13,fHtrxGaus->GetXmax());
	fHtrxCrystal->SetParameters(fHtrxGausPar[0],fHtrxGausPar[1],fHtrxGausPar[2],1,1);
	trXhist->Fit("trxHCrystal","R","ep",-0.13,fHtrxGaus->GetXmax());
	fHtrxCrystal->GetParameters(fHtrxCrystalPar);

	// fit the O peak with gaus
	double_t  fOtrxGausPar[3];
	TF1 *fOtrxGaus=new TF1("trxOGaus","gaus",-0.055,-0.03);
	trXhist->Fit("trxOGaus","R","ep",-0.055,-0.03);
	fOtrxGaus->GetParameters(fOtrxGausPar);

	double_t fOtrxCrystalPar[5];
	TF1 *fOtrxCrystal=new TF1("fOtrxCrystal","crystalball",-0.07,-0.03);
	fOtrxCrystal->SetParameters(fOtrxGausPar[0],fOtrxGausPar[1],fOtrxGausPar[2],1,1);
	trXhist->Fit("fOtrxCrystal","R","ep",-0.07,-0.03);
	fOtrxCrystal->GetParameters(fOtrxCrystalPar);


	double_t fH2OtrxCrystalPar[10];
	std::copy(fHtrxCrystalPar,fHtrxCrystalPar+5,fH2OtrxCrystalPar);
	std::copy(fOtrxCrystalPar,fOtrxCrystalPar+5,fH2OtrxCrystalPar+5);

	TF1 *fH2OtrxCrystal=new TF1("fH2OtrxCrystal","crystalball(0)+crystalball(5)",-0.14,-0.0);
	fH2OtrxCrystal->SetParameters(fH2OtrxCrystalPar);
	fH2OtrxCrystal->SetParNames("H_Constant","H_Mean","H_sigma","H_Alpha","H_N","O_Constant","O_Mean","O_sigma","O_Alpha","O_N");
	std::cout<<"\n tr.x double crystal ball fit"<<std::endl;
	trXhist->Fit("fH2OtrxCrystal","R","ep",-0.14,-0.0);

	gPad->SetLogy(1);

	c0->cd(4);
	TH2F *htgtThPhNoCut = new TH2F("htgthphNoCut","Target theta-phi (theta vertical)",
					 500,-0.03,0.035,500,-0.07,0.07);
	t->Project(htgtThPhNoCut->GetName(),Form("%s.gold.th:%s.gold.ph",HRSArm.c_str(),HRSArm.c_str()),basicCuts.c_str());
	htgtThPhNoCut->DrawCopy("colz");

	gPad->SetGridx(1);
	gPad->SetGridy(1);

	c0->cd(2);

	TH1F *golddPhist = new TH1F("gold.dp", "gold.dp", 1000, -0.01, -0.00); //1676
	t->Project(golddPhist->GetName(), Form("%s.gold.dp", HRSArm.c_str()),
			Form("%s && %s", basicCuts.c_str(), fiducials.c_str()));
	golddPhist->DrawCopy();

	//Get the fit result for gaus, and use the parameter for the crystal ball fit
	// will need to add the code for peak search


	TF1 *fOGaus=new TF1("fOGaus","gaus(0)",-0.005,-0.003);  // get the gaus fit for each individual peak, and use the paramater for the crystal ball
	golddPhist->Fit("fOGaus","R","ep",-0.005,-0.003);
	double_t a[3];
	fOGaus->GetParameters(a);
	std::cout<<a[0]<<"  "<<a[1]<<"   "<<a[2]<<std::endl;

	// add the single crystal ball function fit, used for extract the fit paramaters;
	TF1 *fOCrystal=new TF1("fOCrystal","crystalball",-0.0065,-0.003);  // get the gaus fit for each individual peak, and use the paramater for the crystal ball
	fOCrystal->SetParameters(a[0],a[1],a[2],1.64,1.1615);
	golddPhist->Fit("fOCrystal","R","ep",-0.0065,-0.003);
	double_t fOCrystalPar[5];
	fOCrystal->GetParameters(fOCrystalPar);

	//  get the fit result for the H peak
	TF1 *fHGaus=new TF1("fHGaus","gaus",-0.008,-0.0068);
	golddPhist->Fit("fHGaus","R","ep",-0.008,-0.0068);
	double_t fHGausPar[3];
	fHGaus->GetParameters(fHGausPar);

	TF1 *fHCrystal=new TF1("fHCrystal","crystalball",-0.01,-0.0068);
	fHCrystal->SetParameters(fHGausPar[0],fHGausPar[1],fHGausPar[2],1.64,1.1615);
	golddPhist->Fit("fHCrystal","R","ep",-0.01,-0.0068);
	double_t fHCrystalPar[5];
	fHCrystal->GetParameters(fHCrystalPar);

	// fit for the H2O crystal ball
	TF1  *fH2OSpecCrystal=new  TF1("fH2OSpecCrystal","crystalball(0)+crystalball(5)",-0.01,-0.0);
	double_t fH2OSpecCrystalPar[10];
	std::copy(fHCrystalPar,fHCrystalPar+5,fH2OSpecCrystalPar);
	std::copy(fOCrystalPar,fOCrystalPar+5,fH2OSpecCrystalPar+5);
	fH2OSpecCrystal->SetParameters(fH2OSpecCrystalPar);
	fH2OSpecCrystal->SetParNames("H_Constant","H_Mean","H_sigma","H_Alpha","H_N","O_Constant","O_Mean","O_sigma","O_Alpha","O_N");
	std::cout<<"\n gold.dp double crystal ball fit"<<std::endl;
	golddPhist->Fit("fH2OSpecCrystal","R","ep",-0.01,0);

	gPad->SetLogy(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
}
