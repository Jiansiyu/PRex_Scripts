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
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/Dp/prexLHRS_1696_-1_0_optics_check.root");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/Dp/prexLHRS_1696_-1_0_optics_check_1.root");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/Dp/prexLHRS_1696_-1_0_optics_check_2.root");
	t->Add("/home/newdriver/PRex/PRex_Data/PRex_Data/rootfile/Dp/prexLHRS_1696_-1_0_optics_check_3.root");

	std::string basicCuts =
			Form(   "%s.vdc.u1.nclust==1 && %s.vdc.v1.nclust==1 && %s.vdc.u2.nclust==1 && %s.vdc.v2.nclust==1 &&  %s.gold.dp<1 && %s.gold.dp > -0.1 && fEvtHdr.fEvtType==1",
					HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(), HRSArm.c_str(),HRSArm.c_str());

	// centrol hole cut
//	std::string fiducials = Form("abs(%s.gold.th)<0.004 && abs(%s.gold.ph)<0.002",
//				HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

	// select the large hole on the left side
	std::string fiducials = Form("(%s.gold.th)<0.022 &&(%s.gold.th)>0.014 && (%s.gold.ph)<-0.007 && (%s.gold.ph)>-0.011",
				HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

//	// select the righ bottom hole
//	std::string fiducials = Form("(%s.gold.th)<-0.02 &&(%s.gold.th)>-0.033 && (%s.gold.ph)<0.008 && (%s.gold.ph)>0.004",
//				HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str(),HRSArm.c_str());

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
	TH1F *trXhist= new TH1F("tr.x","tr.x",1000,-0.25,0.15);
	t->Project(trXhist->GetName(),Form("%s.tr.x",HRSArm.c_str()),Form("%s && %s",basicCuts.c_str(),fiducials.c_str()));
	trXhist->DrawCopy();


	// fit fit fit fit
	double_t fCtrxGausPar[5];
	TF1 *fCtrxGaus=new TF1("fCtrxGaus","gaus",-0.02,-0.008);
	trXhist->Fit("fCtrxGaus","RQ","ep",fCtrxGaus->GetXmin(),fCtrxGaus->GetXmax());
	fCtrxGaus->GetParameters(fCtrxGausPar);

	double_t fCtrxCrystalPar[5];
	fCtrxGausPar[3]=1.64;
	fCtrxGausPar[4]=1.16;
	TF1 *fCtrxCrystal=new TF1("fCtrxCrystal","crystalball",-0.20,fCtrxGaus->GetXmax());
	fCtrxCrystal->SetParameters(fCtrxGausPar);
	std::cout<<"\n tr.x single crystal ball fit"<<std::endl;
	trXhist->Fit("fCtrxCrystal","R","ep",fCtrxCrystal->GetXmin(),fCtrxCrystal->GetXmax());
	fCtrxCrystal->Draw("same");


	gPad->SetLogy(1);

	c0->cd(4);
	TH2F *htgtThPhNoCut = new TH2F("htgthphNoCut","Target theta-phi (theta vertical)",
					 500,-0.03,0.035,500,-0.07,0.07);
	t->Project(htgtThPhNoCut->GetName(),Form("%s.gold.th:%s.gold.ph",HRSArm.c_str(),HRSArm.c_str()),basicCuts.c_str());
	htgtThPhNoCut->DrawCopy("colz");

	gPad->SetGridx(1);
	gPad->SetGridy(1);

	c0->cd(2);

	TH1F *golddPhist = new TH1F("gold.dp", "gold.dp", 1000,-0.01, 0.0); //1676
	t->Project(golddPhist->GetName(), Form("%s.gold.dp", HRSArm.c_str()),
			Form("%s && %s", basicCuts.c_str(), fiducials.c_str()));
	golddPhist->DrawCopy();

	//Get the fit result for gaus, and use the parameter for the crystal ball fit
	// will need to add the code for peak search

	double_t fCgoldDpGausPar[5];

    TF1 *fCgoldDpGaus=new TF1("fCgoldDpGaus","gaus",-0.0032,-0.002);
    golddPhist->Fit("fCgoldDpGaus","RQ","",fCgoldDpGaus->GetXmin(),fCgoldDpGaus->GetXmax());
    fCgoldDpGaus->GetParameters(fCgoldDpGausPar);

    // fit the peak with crystal ball
    fCgoldDpGausPar[3]=1.64;
    fCgoldDpGausPar[4]=1.16;
    TF1 *fCgoldDpCrystal=new TF1("fCgoldDpCrystal","crystalball",-0.007,fCgoldDpGaus->GetXmax());  // get the gaus fit for each individual peak, and use the paramater for the crystal ball
	fCgoldDpCrystal->SetParameters(fCgoldDpGausPar);
	std::cout<<"\n gold.dp single crystal ball fit"<<std::endl;
	golddPhist->Fit("fCgoldDpCrystal","R","",fCgoldDpCrystal->GetXmin(),fCgoldDpCrystal->GetXmax());
	fCgoldDpCrystal->Draw("same");



	gPad->SetLogy(1);
	gPad->SetGridx(1);
	gPad->SetGridy(1);
}
