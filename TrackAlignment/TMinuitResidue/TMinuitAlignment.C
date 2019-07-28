/*
 * TMinuitAlignment.C
 *
 *  Created on: Jul 26, 2019
 *      Author: newdriver
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string.h>
#include <vector>

#include "TMinuit.h"
#include "TF1.h"
#include <TH2F.h>
#include <TH1F.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include <TLine.h>

int TargetGEMID=6;

enum DetectorID{
	vdc,
	GEM1,
	GEM2,
	GEM3,
	GEM4,
	GEM5,
	GEM6
};

// CREATE global paramaters that used for buffer the Hit for VDC and the GEMs
struct HitStruct{
public:
	HitStruct(int8_t detID, double_t hitX, double_t hitY, double_t hitZ, double_t hitTh=0.0, double_t hitPh=0.0){
		this->detectorID=detID;
		this->x=hitX;
		this->y=hitY;
		this->z=hitZ;
		this->theta=hitTh;
		this->phi=hitPh;
	};

	~HitStruct(){};

	inline int8_t GetDetectorID(){return detectorID;};
	inline double_t GetX(){return x;};
	inline double_t GetY(){return y;};
	inline double_t GetZ(){return z;};
	inline double_t GetTheta(){return theta;};
	inline double_t GetPhi(){return phi;};



	void Print(){
		std::cout<<"==> ("<<(int)detectorID<<")"<<std::endl;
		std::cout<<"	x    :"<<x<<std::endl
				 <<"	y    :"<<y<<std::endl
				 <<"	z    :"<<z<<std::endl
				 <<"	theta:"<<theta<<std::endl
				 <<"	phi  :"<<phi<<std::endl;};

	inline bool operator == (const HitStruct &y){
		return this->detectorID==y.detectorID;}
private:
	int8_t detectorID;
	double_t x;
	double_t y;
	double_t z;
	double_t theta;   // x'
	double_t phi;     // y'

};

std::vector<std::vector<HitStruct>> DetHitBuff;    // buffers all the hit on the detector VDC & GEM

int LoadDetectorHit(std::string fname="trackxyz.txt") {

	DetHitBuff.clear();

	// VDC result ID, theta, phi, x, y ,z
	// GEM       : ID, x, y, z
	std::cout<<"Input filename: "<<fname.c_str()<<std::endl;
	std::cout<<"==> Data structure requirement"<<std::endl;
	std::ifstream infile(fname.c_str());

	while(infile){
		std::vector<double_t> line_elements;
		std::string line;
		if(!getline(infile,line)) break;
		std::istringstream ss(line);
		while(ss){
			std::string s;
			if(!getline(ss,s,','))break;
			line_elements.push_back(atof(s.c_str()));
		}

		std::vector<HitStruct> DetHit;
	    DetHit.clear();

		//get VDC values
		int8_t vdcID      = line_elements[0];
		double_t vdctheta = line_elements[1];
		double_t vdcphi   = line_elements[2];
		double_t vdcx     = line_elements[3];
		double_t vdcy     = line_elements[4];
		double_t vdcz     = line_elements[5];

		HitStruct vdcHit(vdcID,vdcx,vdcy,vdcz,vdctheta,vdcphi);

		DetHit.push_back(vdcHit);
		// loop on the line elements
		for(int8_t GEMCount=0; GEMCount < (line_elements.size()-6)/4; GEMCount++)
		{
			int8_t gemID=(int8_t)line_elements[6+GEMCount*4];
			double_t gemX=line_elements[7+GEMCount*4];
			double_t gemY=line_elements[8+GEMCount*4];
			double_t gemZ=line_elements[9+GEMCount*4];

			HitStruct gemHit(gemID,gemX,gemY,gemZ);
			DetHit.push_back(gemHit);
		}
		assert(DetHit.size()==7);
		DetHitBuff.push_back(DetHit);
	}
	return 0;
}


// make an invert fit would be better
double_t linearFit(std::vector<HitStruct> DetHit, double_t *x,std::string dimension="xz"){

	//get the average
	assert(DetHit.size()=6);

	std::vector<HitStruct> subGEM(DetHit.begin()++,DetHit.end());  // get the second GEM to the end . the first one will be take as reference, will not going to the fit

	double_t Xaverage, Yaverage, Zaverage;
	double_t Xsum=0.0,Ysum=0.0, Zsum=0.0;
	double_t HitNum=(double_t)subGEM.size();

	for(auto Hit : subGEM){
		Xsum+=Hit.GetX();
		Ysum+=Hit.GetY();
		Zsum+=Hit.GetZ();
	}

	Xaverage=Xsum/HitNum;
	Yaverage=Ysum/HitNum;
	Zaverage=Zsum/HitNum;

	// Calculate the slop
	double_t fitSlop=0.0;

	// loop on the hit
	double_t SlopXZSum=0.0;
	double_t SlopXSum=0.0;
	for(auto Hit : subGEM){
		// this is a X-Z fit
		SlopXZSum+=(Hit.GetX()-Xaverage)*(Hit.GetZ()-Zaverage);
		SlopXSum+=(Hit.GetX()-Xaverage)*(Hit.GetX()-Xaverage);
	}
    fitSlop=SlopXZSum/SlopXSum;

    // generate the functions
    // y-y0=m(x-x0)=> y=m(x-x0)+y0
    return fitSlop*(x[0]-DetHit[0].GetX())+DetHit[0].GetZ();
}

double_t linearInvertFit(std::vector<HitStruct> DetHit, double_t *x,double_t &fitpar, std::string dimension="xz"){

	//get the average
	assert(DetHit.size()==6);

	std::vector<HitStruct> subGEM(DetHit.begin()+1,DetHit.end());  // get the second GEM to the end . the first one will be take as reference, will not going to the fit

	double_t Xaverage, Yaverage, Zaverage;
	double_t Xsum=0.0,Ysum=0.0, Zsum=0.0;
	double_t HitNum=(double_t)subGEM.size();

	for(auto Hit : subGEM){
		Xsum+=Hit.GetX();
		Ysum+=Hit.GetY();
		Zsum+=Hit.GetZ();
	}

	Xaverage=Xsum/HitNum;
	Yaverage=Ysum/HitNum;
	Zaverage=Zsum/HitNum;

	// Calculate the slop
	double_t fitSlop=0.0;

	// loop on the hit
	double_t SlopXZSum=0.0;
	double_t SlopXSum=0.0;
	double_t SlopZSum=0.0;

	for(auto Hit : subGEM){
		// this is a X-Z fit
		SlopXZSum+=(Hit.GetX()-Xaverage)*(Hit.GetZ()-Zaverage);
		SlopXSum+=(Hit.GetX()-Xaverage)*(Hit.GetX()-Xaverage);
		SlopZSum+=(Hit.GetZ()-Zaverage)*(Hit.GetZ()-Zaverage);
	}
    fitSlop=SlopXZSum/SlopZSum;
    fitpar=fitSlop;

    // generate the functions
    // x-x0=m(z-z0)
    return fitSlop*(x[0]-DetHit[0].GetZ())+DetHit[0].GetX();
}


// all the fit applied on the Z-X dimension, in order to solve the INFINIT issues
// create custom residue functions
double_t GetResidual(std::vector<HitStruct> DetHit,std::string dimension="xz"){
	// calculate the distance from the fit functions
	//get the average
	assert(DetHit.size()==6);

	double_t residualF=0.0;
	double_t residualParSq=0.0;

	double_t residue=0.0;
	for(auto Hit : DetHit){
		double_t x[1]={Hit.GetZ()};
		double_t slop=0.0;
		residualF=linearInvertFit(DetHit,x,slop)-Hit.GetX();

		double_t a=slop;
		double_t b=-1.0;

		double_t deltaD=residualF/(std::sqrt(a*a+b*b));
		residue+=deltaD*deltaD;
//		residue+=std::abs(deltaD);

//		std::cout<<"Distance "<<std::abs(deltaD)<<std::endl;
	}
return residue;
}




double_t FitTest(){

	TCanvas *a=new TCanvas("Fit","Fit",1000,1000);
	a->cd();
	a->Draw();

	TH2F *DetHist;

	TH2F *FitHist;
	for (auto Event : DetHitBuff){
		DetHist=new TH2F("Initial Hit","Initial Hit",1000, -0, 3.0, 2000, -0.4, 0.4);
		DetHist->SetMarkerSize(1);
		DetHist->SetMarkerStyle(20);

		FitHist=new TH2F("Fit Hit","Fit Hit",1000, -0, 3.0, 2000, -0.4, 0.4);
		FitHist->SetMarkerSize(1);
		FitHist->SetMarkerStyle(20);
		FitHist->SetMarkerColor(3);

		for(auto Hit:Event){
			DetHist->Fill(Hit.GetZ(),Hit.GetX());
		}
		DetHist->Draw();

		std::vector<HitStruct> gemHit(Event.begin()+1,Event.end());
		for(auto Hit:gemHit){
			std::cout<<gemHit.size()<<std::endl;
			double_t a[1]={Hit.GetZ()};
			double_t slop=0.0;
			FitHist->Fill(Hit.GetZ(),linearInvertFit(gemHit,a,slop));
			GetResidual(gemHit);
			std::cout<<"slop:"<<slop<<std::endl;
		}
		FitHist->Draw("same");
		a->Update();
		getchar();

	}
	return 0.0;
}

void fcn_residual(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// correction matrix     [0] [1] [2]
	// firt ste, only take the tranlation into consideration
	//
	//   6 GEM detetors  5 need to be aligned
	//   GEM 2  X par[0]  Y par[1]
	//   GEM 3  X par[2]  Y par[3]
	//   GEM 4  X par[4]  Y par[5]
	//   GEM 5  X par[6]  Y par[7]
	//   GEM 6  X par[8]  Y par[9]
	Double_t chisq=0.0;

	TH2F *FitFunc;
	// loop on the event and apply the corretion
	int counter=0;
	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		assert(gemEvent.size()==6);
		// appli the correction matrix on the GEM detectors
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());
		assert(gemcorrEvent.size()==5);

		// create vector of GEMs that contain the correction parameters
		std::vector<HitStruct> gemCorrected;

		gemCorrected.push_back(gemEvent[0]);
		for(auto Hit : gemcorrEvent){
			// start from GEM 2
			HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
			gemCorrected.push_back(a);
		}

		// test the funct
//		std::cout<<gemCorrected.size()<<std::endl;
/*		for (auto Hit : gemCorrected){
			FitFunc->Fill(Hit.GetZ(),Hit.GetX());
		}
		FitFunc=new TH2F("zx","zx",1000,-0,3.0,2000,-0.4,0.4);
		FitFunc->Fit("pol1","Q");
		chisq+=FitFunc->GetFunction("pol1")->GetChisquare();

		std::cout<<counter++<<std::endl;

		FitFunc->Delete();*/
		chisq+=GetResidual(gemCorrected);
	}


	f=chisq;

}

void fcn_residual_tranlationRotation(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// take the translation and rotation into consideration

}

double_t GetDeltD(std::vector<HitStruct> DetHit,std::string dimension="xz",int detectorID=0){
	//
	// calculate the distance from the fit functions
	//get the average
	assert(DetHit.size()==6);

	double_t residualF=0.0;
	double_t residualParSq=0.0;

	double_t residue=0.0;
	for(auto Hit : DetHit){
		double_t x[1]={Hit.GetZ()};
		double_t slop=0.0;
		residualF=linearInvertFit(DetHit,x,slop)-Hit.GetX();

		double_t a=slop;
		double_t b=-1.0;

		double_t deltaD=residualF/(std::sqrt(a*a+b*b));
		residue+=(deltaD);
		if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltaD);
//		std::cout<<"Distance "<<std::abs(deltaD)<<std::endl;
	}
return residue;
}

double_t MinimizerCheck(double_t *par, TCanvas *a){
	a->Divide(2,1);
	a->cd(1);

	for(int i =0 ; i < 10 ; i ++) std::cout<<par[i]<<std::endl;
	TH1F *residualBefore=new TH1F("ResdualBefore","ResdualBefore",100,-0.01,0.01);
	residualBefore->GetYaxis()->SetRangeUser(0,3000);
	TH1F *residualAfter=new TH1F("ResdualAfter","ResdualAfter",100,residualBefore->GetXaxis()->GetXmin(),residualBefore->GetXaxis()->GetXmax());
	residualAfter->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
	    	gemCorrected.push_back(a);
		}
	    residualBefore->Fill(GetResidual(gemEvent));
	    residualAfter->Fill(GetResidual(gemCorrected));

	}
	residualBefore->Draw();
	residualAfter->Draw("same");

	//distance check
	a->cd(2)->Divide(1,2);
	a->cd(2)->cd(1);
	TH1F *deltDBeforAlignSum=new TH1F("DeltaDBefore","DeltaDBefore",100,-0.1,0.1);
	TH1F *deltDAfterAlignSum=new TH1F("DeltaDAfter","DeltadAfter",100,-0.1,0.1);
	deltDAfterAlignSum->SetLineColor(3);

	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
	    	gemCorrected.push_back(a);
		}
	    deltDBeforAlignSum->Fill(GetDeltD(gemEvent));
	    deltDAfterAlignSum->Fill(GetDeltD(gemCorrected));
	}
	deltDBeforAlignSum->Draw();
	deltDAfterAlignSum->Draw("same");

	a->cd(2)->cd(2)->Divide(3,2);
	TH1F *deltDBeforAlign[6];
	TH1F *deltDAfterAlign[6];
	for (int i =0 ; i < 6 ; i ++){
		deltDBeforAlign[i]=new TH1F(Form("%d_before_deltaD",i),Form("%d_before_deltaD",i),100,-0.1,0.1);
		deltDAfterAlign[i]=new TH1F(Form("%d_After_deltaD",i),Form("%d_After_deltaD",i),100,-0.1,0.1);
	}
	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
	    	gemCorrected.push_back(a);
		}
	    for (auto Hit : gemEvent){
	    	deltDBeforAlign[Hit.GetDetectorID()-1]->Fill(GetDeltD(gemEvent,"",Hit.GetDetectorID()));
	    	deltDAfterAlign[Hit.GetDetectorID()-1]->Fill(GetDeltD(gemCorrected,"",Hit.GetDetectorID()));
	    }
	    deltDBeforAlignSum->Fill(GetDeltD(gemEvent));
	    deltDAfterAlignSum->Fill(GetDeltD(gemCorrected));
	}

	for (int i =0 ; i < 6 ; i ++){
		a->cd(2)->cd(2)->cd(i+1);
		deltDBeforAlign[i]->Draw();
		deltDAfterAlign[i]->SetLineColor(3);
		deltDAfterAlign[i]->Draw("same");
	}

	a->Update();

	// calculate the prediced position vs the real pos

	return 0.0;
}

void TMinimer(){
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(10);  // initialize the minuit for 10 parameters
	gMinuit->SetFCN(fcn_residual);

	//   6 GEM detetors  5 need to be aligned
	//   GEM 2  X par[0]  Y par[1]
	//   GEM 3  X par[2]  Y par[3]
	//   GEM 4  X par[4]  Y par[5]
	//   GEM 5  X par[6]  Y par[7]
	//   GEM 6  X par[8]  Y par[9]

	double_t vstart[10]={0.0, 0.0, 0.0, 0.0, 0.0,
			             0.0, 0.0, 0.0, 0.0, 0.0};

	double_t step[10]  ={1e-06, 1e-06, 1e-06, 1e-06, 1e-06,
			             1e-06, 1e-06, 1e-06, 1e-06, 1e-06};

	                 // about the 10 degree,  for translation about 20cm

//	double_t bmin[10]={vstart[0]-0.1,vstart[1]-0.1,vstart[2]-0.1,vstart[3]-0.1,vstart[4]-0.1,
//					   vstart[5]-0.1,vstart[6]-0.1,vstart[7]-0.1,vstart[8]-0.1,vstart[9]-0.1};
//
//	double_t bmax[10]={vstart[0]+0.1,vstart[1]+0.1,vstart[2]+0.1,vstart[3]+0.1,vstart[4]+0.1,
//					   vstart[5]+0.1,vstart[6]+0.1,vstart[7]+0.1,vstart[8]+0.1,vstart[9]+0.1};

	double_t Binrange=0.01;
	double_t bmin[10]={vstart[0]-Binrange,vstart[1]-Binrange,vstart[2]-Binrange,vstart[3]-Binrange,vstart[4]-Binrange,
					   vstart[5]-Binrange,vstart[6]-Binrange,vstart[7]-Binrange,vstart[8]-Binrange,vstart[9]-Binrange};

	double_t bmax[10]={vstart[0]+Binrange,vstart[1]+Binrange,vstart[2]+Binrange,vstart[3]+Binrange,vstart[4]+Binrange,
					   vstart[5]+Binrange,vstart[6]+Binrange,vstart[7]+Binrange,vstart[8]+Binrange,vstart[9]+Binrange};

	double_t arglist[10];
	int ierflg=0;

	gMinuit->mnparm( 0, "a", vstart[0],step[0],bmin[0],bmax[0],ierflg);
	gMinuit->mnparm( 1, "b", vstart[1],step[1],bmin[1],bmax[1],ierflg);
	gMinuit->mnparm( 2, "c", vstart[2],step[2],bmin[2],bmax[2],ierflg);
	gMinuit->mnparm( 3, "d", vstart[3],step[3],bmin[3],bmax[3],ierflg);
	gMinuit->mnparm( 4, "e", vstart[4],step[4],bmin[4],bmax[4],ierflg);
	gMinuit->mnparm( 5, "f", vstart[5],step[5],bmin[5],bmax[5],ierflg);
	gMinuit->mnparm( 6, "g", vstart[6],step[6],bmin[6],bmax[6],ierflg);
	gMinuit->mnparm( 7, "h", vstart[7],step[7],bmin[7],bmax[7],ierflg);
	gMinuit->mnparm( 8, "l", vstart[8],step[8],bmin[8],bmax[8],ierflg);
	gMinuit->mnparm( 9, "m", vstart[9],step[9],bmin[9],bmax[9],ierflg);


	// Set the output
	// set the print level
	// -1 no output
	// 1 standdard output
	gMinuit->SetPrintLevel(1);

	//minimization strategy
	// 1 standard
	// 2 try to improve minimum (slower)
	arglist[0]=2;
	gMinuit->mnexcm("SET STR",arglist, 1, ierflg);

	// Call the minimizer
	arglist[0]=5e6;

	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
//	gMinuit->mnsimp();

	// read out the parameters.
	double_t MiniPars[10];
	double_t MiniParsErr[10];
	for(int i =0 ; i < 10 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
		vstart[i]=MiniPars[i];
	}

	// changed the coordination and refit
//	Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]
	std::vector<std::vector<HitStruct>> newDetHitBuff;
	for (auto Event : DetHitBuff){
		std::vector<HitStruct> EventBuff;
		EventBuff.push_back(Event[0]);
		EventBuff.push_back(Event[1]);

		std::vector<HitStruct> GEMCorr(Event.begin()+2,Event.end());

		for (auto Hit : GEMCorr){

			if(Hit.GetDetectorID()>=2){
				HitStruct a(Hit.GetDetectorID(),Hit.GetX()-MiniPars[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-MiniPars[Hit.GetDetectorID()-1]);
				EventBuff.push_back(a);
			}
		}
		assert(EventBuff.size()==7);
		newDetHitBuff.push_back(EventBuff);
	}
	DetHitBuff.clear();
	DetHitBuff.insert(DetHitBuff.begin(),newDetHitBuff.begin(),newDetHitBuff.end());

	gMinuit->mnparm( 0, "a", vstart[0],step[0],bmin[0],bmax[0],ierflg);
	gMinuit->mnparm( 1, "b", vstart[1],step[1],bmin[1],bmax[1],ierflg);
	gMinuit->mnparm( 2, "c", vstart[2],step[2],bmin[2],bmax[2],ierflg);
	gMinuit->mnparm( 3, "d", vstart[3],step[3],bmin[3],bmax[3],ierflg);
	gMinuit->mnparm( 4, "e", vstart[4],step[4],bmin[4],bmax[4],ierflg);
	gMinuit->mnparm( 5, "f", vstart[5],step[5],bmin[5],bmax[5],ierflg);
	gMinuit->mnparm( 6, "g", vstart[6],step[6],bmin[6],bmax[6],ierflg);
	gMinuit->mnparm( 7, "h", vstart[7],step[7],bmin[7],bmax[7],ierflg);
	gMinuit->mnparm( 8, "l", vstart[8],step[8],bmin[8],bmax[8],ierflg);
	gMinuit->mnparm( 9, "m", vstart[9],step[9],bmin[9],bmax[9],ierflg);
	// Set the output
	// set the print level
	// -1 no output
	// 1 standdard output
	gMinuit->SetPrintLevel(1);

	//minimization strategy
	// 1 standard
	// 2 try to improve minimum (slower)
	arglist[0]=2;
	gMinuit->mnexcm("SET STR",arglist, 1, ierflg);

	// Call the minimizer
	arglist[0]=5e8;

	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
//	gMinuit->mnsimp();

	// read out the parameters.
//	double_t MiniPars[10];
//	double_t MiniParsErr[10];
	for(int i =0 ; i < 10 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
		vstart[i]=MiniPars[i];
	}


	TCanvas *a=new TCanvas("a","a",1000,1000);
	MinimizerCheck(MiniPars,a);
	/*
	// load the parameter to draw
	if(true){
		TCanvas *a = new TCanvas("Minimize Check", "Minimize Check", 1000, 1000);
		a->Divide(2,1);
		a->Draw();

		a->cd(1);
		// plot the initial result before the correction
		TH1F *DeltaX=new TH1F("DeltaX","DeltaX",100,-0.01,0.01);
		DeltaX->GetYaxis()->SetRangeUser(0,500);
		TH1F *DeltaXCorr=new TH1F("DeltaXCorr","DeltaXCorr",100,-0.01,0.01);
		DeltaXCorr->SetLineColor(3);

		for(auto Event : DetHitBuff){
			auto vdc=Event[0];
			auto gem=Event[1];

			double_t residue=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
			if(residue<-0.004 || residue>0.01) continue;

			// correct the position for the X and Y
			double_t Xcorr=gem.GetX() * MiniPars[0] + gem.GetZ() * MiniPars[1] + MiniPars[2];
			double_t Zcorr=gem.GetX() * MiniPars[3] + gem.GetZ() * MiniPars[4] + MiniPars[5];

			double_t delta=vdc.GetX()+Zcorr*vdc.GetTheta()-Xcorr;
			DeltaXCorr->Fill(delta);
			DeltaX->Fill(vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX());
		}
	DeltaX->Draw();
	DeltaXCorr->Draw("same");
	a->Modified();
	a->Update();
	}
*/

}

/*
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// correction matrix     [0] [1] [2]
	//						 [3] [4] [5]
	// Correct Matrix Range
	//   Xcorr = [0]*X + [1] * Z + [2]
	//   Zcorr = [3]*X + [4] * Z + [5]
	//   [0]    : 1
	//   [1]    : 0
	//   [2]    : 0
	//
	//	 [3]    : 0
	//   [4]    : 1
	//   [5]    : 0

	Double_t chisq=0.0;
	// loop on the Event, and get the hit
	for( auto Event : DetHitBuff){
		// loop on the event
		// just get the first Detector to have a try
		HitStruct vdc=Event[0];
		HitStruct gem=Event[TargetGEMID];


		double_t residue=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
		if((TargetGEMID==1)&&(residue<-0.004 || residue>0.01)) continue;
		if((TargetGEMID==4)&&(residue<-0.02 || residue>0.02))continue;

		assert(gem.GetDetectorID()==TargetGEMID);
		// Get the corrected
		double_t Xcorr=gem.GetX() * par[0] + gem.GetZ() * par[1] + par[2];
		double_t Zcorr=gem.GetX() * par[3] + gem.GetZ() * par[4] + par[5];

		// Project the VDC hit to the GEM plane
		double_t vdcProjectX=vdc.GetX() + Zcorr * vdc.GetTheta();    //Z should use the corrected Z to do the projection
		double_t vdcProjectY;
		double_t vdcProjectZ;

		//Since here we use the corrected Z, in the residue, just take the x into consideration
		double_t delta= vdcProjectX-Xcorr;
		chisq=delta*delta;
	}
	f=chisq;
}*/


int ProjectVDC(int id=0){
//	LoadDetectorHit();
	TCanvas *a = new TCanvas("Project", "Project", 1000, 1000);
	a->Divide(2,1);

	a->Draw();
	a->cd(1);
	TH1F *DeltaX=new TH1F("DeltaX","DeltaX",100,-0.05,0.05);
	TH1F *ResidueX=new TH1F("ResidueX","ResidueX",100,-0.05,0.05);
	ResidueX->SetLineColor(3);
	for (auto Event : DetHitBuff)
	{
		HitStruct vdc=Event[0];
		HitStruct gem=Event[TargetGEMID];
		double_t residue=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
		if(residue<-0.02 || residue>0.02) continue;
		double_t delta=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
		DeltaX->Fill(delta);
		ResidueX->Fill(delta*delta);

	}
	ResidueX->Draw();
	DeltaX->Draw("same");

	a->cd(2);
	for (auto Event : DetHitBuff)
	{
		HitStruct vdc=Event[0];
		HitStruct gem=Event[TargetGEMID];
		TH2F *vdcHist = new TH2F("Projectvdc", "Projectvdc", 2000, -0.4, 0.4,1000, -0, 3.0);
		TH2F *vdcProjectHist = new TH2F("ProjectvdcGEM", "ProjectvdcGEM", 2000,-0.4, 0.4, 1000, -0, 3.0);
		TH2F *gemHist = new TH2F("Projectgem", "Projectgem", 2000, -0.4, 0.4,1000, -0, 3.0);

		vdcHist->SetMarkerSize(1);
		vdcHist->SetMarkerStyle(20);
		vdcProjectHist->SetMarkerSize(1);
		vdcProjectHist->SetMarkerStyle(20);

		gemHist->SetMarkerSize(1);
		gemHist->SetMarkerStyle(20);
		gemHist->SetMarkerColor(3);

//		double_t residue=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
//		if(residue<-0.004 || residue>0.01) continue;
        vdcHist->Fill(vdc.GetX(), vdc.GetZ());
		vdcProjectHist->Fill(vdc.GetX() + gem.GetZ() * vdc.GetTheta(),gem.GetZ());
		gemHist->Fill(gem.GetX(), gem.GetZ());


		vdcHist->Draw();
		vdcProjectHist->Draw("same");
		gemHist->Draw("same");

		a->Modified();
		a->Update();
		getchar();
	}

	return 0;
}

/*
void TMinimer(){
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(6);  // initialize the minuit for 6 parameters
	gMinuit->SetFCN(fcn);

	// correction matrix     [0] [1] [2]
	//						 [3] [4] [5]
	// Correct Matrix Range
	//   Xcorr = [0]*X + [1] * Z + [2]
	//   Zcorr = [3]*X + [4] * Z + [5]
	//   [0]    : 1    rotation
	//   [1]    : 0    rotation
	//   [2]    : 0    translation
	//
	//	 [3]    : 0    rotation
	//   [4]    : 1    rotation
	//   [5]    : 0    translation
	// setup for the initial values
	double_t vstart[6]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
	double_t step[6]  ={1e-04, 1e-04, 1e-05, 1e-04, 1e-04, 1e-04};
	                 // about the 10 degree,  for translation about 20cm
	double_t bmin[6]={vstart[0]-0.1,vstart[1]-0.1,vstart[2]-0.1,vstart[3]-0.1,vstart[4]-0.1,vstart[5]-0.1};
	double_t bmax[6]={vstart[0]+0.1,vstart[1]+0.1,vstart[2]+0.1,vstart[3]+0.1,vstart[4]+0.1,vstart[5]+0.1};

	double_t arglist[10];
	int ierflg=0;

	gMinuit->mnparm( 0, "a", vstart[0],step[0],bmin[0],bmax[0],ierflg);
	gMinuit->mnparm( 1, "b", vstart[1],step[1],bmin[1],bmax[1],ierflg);
	gMinuit->mnparm( 2, "c", vstart[2],step[2],bmin[2],bmax[2],ierflg);
	gMinuit->mnparm( 3, "d", vstart[3],step[3],bmin[3],bmax[3],ierflg);
	gMinuit->mnparm( 4, "e", vstart[4],step[4],bmin[4],bmax[4],ierflg);
	gMinuit->mnparm( 5, "f", vstart[5],step[5],bmin[5],bmax[5],ierflg);

	// Set the output
	// set the print level
	// -1 no output
	// 1 standdard output
	gMinuit->SetPrintLevel(1);

	//minimization strategy
	// 1 standard
	// 2 try to improve minimum (slower)
	arglist[0]=2;
	gMinuit->mnexcm("SET STR",arglist, 1, ierflg);

	// Call the minimizer
	arglist[0]=50000;

	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
//	gMinuit->mnsimp();

	// read out the parameters.
	double_t MiniPars[6];
	double_t MiniParsErr[6];
	for(int i =0 ; i < 6 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
	}
	for(int i =0 ; i < 6 ; i ++){
		std::cout<<"Par "<< i <<"   "<< MiniPars[i]<<"   "<<MiniParsErr[i]<<std::endl;
		vstart[i]= MiniPars[i];
	}

	// load the parameter to draw
	if(true){
		TCanvas *a = new TCanvas("Minimize Check", "Minimize Check", 1000, 1000);
		a->Divide(2,1);
		a->Draw();

		a->cd(1);
		// plot the initial result before the correction
		TH1F *DeltaX=new TH1F("DeltaX","DeltaX",100,-0.01,0.01);
		DeltaX->GetYaxis()->SetRangeUser(0,500);
		TH1F *DeltaXCorr=new TH1F("DeltaXCorr","DeltaXCorr",100,-0.01,0.01);
		DeltaXCorr->SetLineColor(3);

		for(auto Event : DetHitBuff){
			auto vdc=Event[0];
			auto gem=Event[1];

			double_t residue=vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX();
			if(residue<-0.004 || residue>0.01) continue;

			// correct the position for the X and Y
			double_t Xcorr=gem.GetX() * MiniPars[0] + gem.GetZ() * MiniPars[1] + MiniPars[2];
			double_t Zcorr=gem.GetX() * MiniPars[3] + gem.GetZ() * MiniPars[4] + MiniPars[5];

			double_t delta=vdc.GetX()+Zcorr*vdc.GetTheta()-Xcorr;
			DeltaXCorr->Fill(delta);
			DeltaX->Fill(vdc.GetX() + gem.GetZ() * vdc.GetTheta()-gem.GetX());
		}
	DeltaX->Draw();
	DeltaXCorr->Draw("same");
	a->Modified();
	a->Update();

	}


}

// the shift parameter should close to 0
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

	Double_t chisq=0.0;
	Double_t delta=0.0;
	for(auto dhit : dd){
		chisq=(dhit.GetX()-par[0])*(dhit.GetX()-par[0])+
				(dhit.GetY()-par[1])*(dhit.GetY()-par[1])+
				(dhit.GetZ()-par[2])*(dhit.GetZ()-par[2])/0.001;
	}
	f=chisq;
}


void TMinimizer(int8_t detectorID=0){

	dd.clear();
	GetAlignmentMatrix(detectorID);

	gStyle->SetOptFile(1111111);
//	TGraphErrors *gr_f1 = new TGraphErrors();
//	gr_f1->SetTitle("Linear y=a*x+b;x;y");
//	gr_f1->SetMarkerSize(1.0);
//	gr_f1->SetMarkerColor(kBlue);
	TMinuit *gMinuit=new TMinuit(3);  // initialize the minuit for 3 parameters
	gMinuit->SetFCN(fcn);

    Double_t arglist[10];//declare flags
    Int_t ierflg = 0;
    arglist[0] = 2;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);// just for information
    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);// just for information
    // set the start value and the step size for parameters
    static Double_t vstart[3]={-0.2,-0.2,-0.5};
    static Double_t step[3]={0.00001,0.00001,0.000001};

    gMinuit->mnparm(0, "shiftX", vstart[0], step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "shiftY", vstart[1], step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "shiftZ", vstart[1], step[1], 0,0,ierflg);

    // now ready for minimization step
    arglist[0]=5000000; // 5000 step
    arglist[1]=0.00001;    // stop when it reach a condition. If you want to do more change 1 to 0.1
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);//call MIGRAD to do optimization



    // The following part of the code used for generate the plot
    double_t shiftX,shiftY,shiftZ;
    double_t shiftXErr,shiftYErr,shiftZErr;

    gMinuit->GetParameter(0,shiftX,shiftXErr);
    gMinuit->GetParameter(1,shiftY,shiftYErr);
    gMinuit->GetParameter(2,shiftZ,shiftZErr);

    TCanvas *a=new TCanvas("test","test",1000,1000);
    a->Divide(1,3);

    TH1F  *residueBeforeX=new TH1F("alignX0","alignX0",100,0,0.001);
    TH1F  *residueBeforeY=new TH1F("alignY0","alignY0",100,0,0.001);
    TH1F  *residueBeforeZ=new TH1F("alignZ0","alignZ0",100,0,0.001);
    residueBeforeX->SetLineColor(2);
	residueBeforeY->SetLineColor(2);
	residueBeforeZ->SetLineColor(2);

    TH1F  *residueAfterX=new TH1F("alignX","alignX",100,0,0.001);
    TH1F  *residueAfterY=new TH1F("alignY","alignY",100,0,0.001);
    TH1F  *residueAfterZ=new TH1F("alignZ","alignZ",100,0,0.001);
    residueAfterX->SetLineColor(3);
	residueAfterY->SetLineColor(3);
	residueAfterZ->SetLineColor(3);

    for (auto dhit : dd){
    	residueBeforeX->Fill((dhit.GetX()*dhit.GetX()));
    	residueBeforeY->Fill((dhit.GetY()*dhit.GetY()));
    	residueBeforeZ->Fill((dhit.GetZ()*dhit.GetZ()));

    	residueAfterX->Fill((dhit.GetX()-shiftX)*(dhit.GetX()-shiftX));
    	residueAfterY->Fill((dhit.GetY()-shiftY)*(dhit.GetY()-shiftY));
    	residueAfterZ->Fill((dhit.GetZ()-shiftZ)*(dhit.GetZ()-shiftZ));
    }
    a->cd(1);
    residueBeforeX->GetYaxis()->SetRangeUser(0,3000);
    residueBeforeX->Draw();
    residueAfterX->Draw("same");
    a->cd(2);
    residueBeforeY->GetYaxis()->SetRangeUser(0,3000);
    residueBeforeY->Draw();
    residueAfterY->Draw("same");
    a->cd(3);
    residueBeforeZ->GetYaxis()->SetRangeUser(0,3000);
    residueBeforeZ->Draw();
    residueAfterZ->Draw("same");


//    TH1F *residueBefore=new TH1F("alignment1","alignment1",100,0,0.001);
//    TH1F *residueAfter=new TH1F("alignment","alignment",100,0,0.001);
//    for(auto dhit : dd){
//    	residueBefore->Fill((dhit.GetX())*(dhit.GetX())+
//				(dhit.GetY())*(dhit.GetY())+
//				(dhit.GetZ())*(dhit.GetZ()));
//    	residueAfter->Fill((dhit.GetX()-shiftX)*(dhit.GetX()-shiftX)+
//				(dhit.GetY()-shiftY)*(dhit.GetY()-shiftY)+
//				(dhit.GetZ()-shiftZ)*(dhit.GetZ()-shiftZ));
//    }
//    residueAfter->SetLineColor(3);
//    residueAfter->Draw();
//    residueBefore->Draw("same");

}

void AligmentCheck(std::string fname="trackxyz.txt"){
	// VDC result ID, theta, phi, x, y ,z
	// GEM       : ID, x, y, z
	//
	double correction_x[]={0,5.47824e-03,1.86478e-03,2.27380e-03,-7.04919e-03,-9.74938e-03,-1.29634e-02 };
	double correction_y[]={0,6.34000e-03,9.61117e-03,9.75368e-03,3.72736e-03, 2.19724e-03,2.63852e-04 };
	double correction_z[]={0,1.64880e-06,1.64880e-06,1.64880e-06,1.64880e-06,1.64880e-06,1.64880e-06};
	std::cout<<"Input filename: "<<fname.c_str()<<std::endl;
	std::cout<<"==> Data structure requirement"<<std::endl;
	std::ifstream infile(fname.c_str());

	TCanvas *a=new TCanvas("a","a",1000,1000);
	a->cd();



	while(infile){
		TH2F *beforexz=new TH2F("axz","axz",2000,-0.4,0.4,1000,-0,3.0);
		TH2F *beforeyz=new TH2F("ayz","ayz",2000,-0.4,0.4,1000,-0,3.0);
		beforexz->SetMarkerSize(1);
		beforexz->SetMarkerColor(2);
		beforexz->SetMarkerStyle(20);

		beforeyz->SetMarkerColor(2);
		beforeyz->SetMarkerSize(1);
		beforeyz->SetMarkerStyle(20);

		TH2F *afterxz=new TH2F("axz","axz",2000,-0.4,0.4,1000,-0,3.0);
		TH2F *afteryz=new TH2F("ayz","ayz",2000,-0.4,0.4,1000,-0,3.0);

		afterxz->SetMarkerSize(1);
		afterxz->SetMarkerColor(3);
		afterxz->SetMarkerStyle(20);

		afteryz->SetMarkerColor(3);
		afteryz->SetMarkerSize(1);
		afteryz->SetMarkerStyle(20);

		std::vector<double_t> line_elements;
		std::string line;
		if(!getline(infile,line)) break;
		std::istringstream ss(line);
		while(ss){
			std::string s;
			if(!getline(ss,s,','))break;
			line_elements.push_back(atof(s.c_str()));
		}

		//get VDC values
		int8_t vdcID      = line_elements[0];
		double_t vdctheta = line_elements[1];
		double_t vdcphi   = line_elements[2];
		double_t vdcx     = line_elements[3];
		double_t vdcy     = line_elements[4];
		double_t vdcz     = line_elements[5];

		HitStruct vdcHit(vdcID,vdcx,vdcy,vdcz,vdctheta,vdcphi);
		beforexz->Fill(vdcx,vdcz);
		beforeyz->Fill(vdcy,vdcz);
		afterxz->Fill(vdcx,vdcz);
		afteryz->Fill(vdcy,vdcz);

		std::vector<HitStruct> GEMHit;
		// loop on the line elements
		for(int8_t GEMCount=0; GEMCount < (line_elements.size()-6)/4; GEMCount++)
		{
			int8_t gemID=(int8_t)line_elements[6+GEMCount*4];
			double_t gemX=line_elements[7+GEMCount*4];
			double_t gemY=line_elements[8+GEMCount*4];
			double_t gemZ=line_elements[9+GEMCount*4];
			beforexz->Fill(gemX,gemZ);
			beforeyz->Fill(gemY,gemZ);
			afterxz->Fill(gemX+correction_x[gemID],gemZ+correction_z[gemID]);
			afteryz->Fill(gemY+correction_y[gemID],gemZ+correction_z[gemID]);
		}
		TLine *vdctrackXZ=new TLine(vdcx,vdcz,vdcx+vdctheta*2.60,2.60);

		vdctrackXZ->SetLineWidth(1);
		vdctrackXZ->SetLineColor(6);
		beforexz->Draw();
		vdctrackXZ->Draw("same");
		afterxz->Draw("same");
		a->Update();

		getchar();
		a->Clear();
		beforexz->Clear();
		beforeyz->Clear();
		afterxz->Clear();
		afteryz->Clear();

	}
}*/
