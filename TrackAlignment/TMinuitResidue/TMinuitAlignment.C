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
#include <assert.h>
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
	assert(DetHit.size()==6);

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
//	double_t residualParSq=0.0;

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

double GetChisq(std::vector<HitStruct> DetHit,std::string dimension="xz"){
	assert(DetHit.size()==6);

	double_t residualF=0.0;
	double_t residualPar=0.0;

	double_t residue=0.0;
	double_t Chisq=0.0;
	for(auto Hit : DetHit){
		double_t x[1]={Hit.GetZ()};
		double_t slop=0.0;
		double f=linearInvertFit(DetHit,x,slop);
		Chisq+=(f-Hit.GetX())*(f-Hit.GetX())/Hit.GetX();
//		residualF=linearInvertFit(DetHit,x,slop)-Hit.GetX();

	}
return Chisq;
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
			assert((Hit.GetDetectorID()==2)||(Hit.GetDetectorID()==3)||
					(Hit.GetDetectorID()==4)||(Hit.GetDetectorID()==5)||
					(Hit.GetDetectorID()==6));
			double_t  Xcorr=Hit.GetX()-par[(Hit.GetDetectorID()-2)*2];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=Hit.GetZ()-par[(Hit.GetDetectorID()-2)*2+1];
			HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);
			//HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
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

		FitFunc->Delete();
		*/

		chisq+=GetResidual(gemCorrected);
	}
	f=chisq;
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

	}
return residue;
}

double_t MinimizerCheck(double_t *par, TCanvas *a){
	a->Divide(2,1);
	a->cd(1)->Divide(1,2);
	a->cd(1)->cd(1);
	for(int i =0 ; i < 10 ; i ++) std::cout<<par[i]<<std::endl;
	TH1F *residualBefore=new TH1F("ResdualBefore","ResdualBefore",100,-0.001,0.001);
	residualBefore->GetYaxis()->SetRangeUser(0,3000);
	TH1F *residualAfter=new TH1F("ResdualAfter","ResdualAfter",100,residualBefore->GetXaxis()->GetXmin(),residualBefore->GetXaxis()->GetXmax());
	residualAfter->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
	    	double_t  Xcorr=Hit.GetX()-par[(Hit.GetDetectorID()-2)*2];
	    	double_t  Ycorr=Hit.GetY();
	    	double_t  Zcorr=Hit.GetZ()-par[(Hit.GetDetectorID()-2)*2+1];
	    	HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

	    	gemCorrected.push_back(a);
		}
	    residualBefore->Fill(GetResidual(gemEvent));
	    residualAfter->Fill(GetResidual(gemCorrected));

	}
	residualBefore->Draw();
	residualAfter->Draw("same");

	// check Theta compare with VDC
	a->cd(1)->cd(2)->Divide(2,1);
	a->cd(1)->cd(2)->cd(1);
	TH2F *thetaCheckAfter=new TH2F("thetaCheckAfter","thetaCheckAfter",1000,-0.03,0.03,1000,-0.03,0.03);
	TH2F *thetaCheckBefore=new TH2F("thetaCheckBefore","thetaCheckBefore",1000,-0.03,0.03,1000,-0.03,0.03);
	thetaCheckAfter->SetMarkerColor(3);
	TH1F *thetaVDCHisto=new TH1F("vdc","vdc",100,-0.03,0.03);
	thetaVDCHisto->SetLineColor(1);
	TH1F *thetaGEMBeforeHisto=new TH1F("GEMBefore","GEMBefore",100,-0.03,0.03);
	thetaGEMBeforeHisto->SetLineColor(2);
	TH1F *thetaGEMAfterHisto=new TH1F("thetaGEMAfterHisto","thetaGEMAfterHisto",100,-0.03,0.03);
	thetaGEMAfterHisto->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
	    	double_t  Xcorr=Hit.GetX()-par[(Hit.GetDetectorID()-2)*2];
	    	double_t  Ycorr=Hit.GetY();
	    	double_t  Zcorr=Hit.GetZ()-par[(Hit.GetDetectorID()-2)*2+1];
	    	HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

	    	gemCorrected.push_back(a);
		}
	    double x[1]={0.0};
	    double slop=0.0;
	    linearInvertFit(gemCorrected,x,slop);
	    thetaCheckAfter->Fill(Event[0].GetTheta(),slop);
	    thetaGEMAfterHisto->Fill(slop);

//	    thetaCheckBefore
	    linearInvertFit(gemcorrEvent,x,slop);
	    thetaCheckBefore->Fill(Event[0].GetTheta(),slop);
	    thetaVDCHisto->Fill(Event[0].GetTheta());
	    thetaGEMBeforeHisto->Fill(slop);

	}
	thetaCheckAfter->Draw();
	thetaCheckBefore->Draw("same");
	a->cd(1)->cd(2)->cd(2);
	thetaVDCHisto->Draw();
	thetaGEMBeforeHisto->Draw("same");
	thetaGEMAfterHisto->Draw("same");

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
	    	double_t  Xcorr=Hit.GetX()-par[(Hit.GetDetectorID()-2)*2];
	    				double_t  Ycorr=Hit.GetY();
	    				double_t  Zcorr=Hit.GetZ()-par[(Hit.GetDetectorID()-2)*2+1];
	    				HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
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
		deltDBeforAlign[i]=new TH1F(Form("%d_before_deltaD",i),Form("%d_before_deltaD",i),100,-0.02,0.02);
		deltDAfterAlign[i]=new TH1F(Form("%d_After_deltaD",i),Form("%d_After_deltaD",i),100,-0.02,0.02);
	}
	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
	    	double_t  Xcorr=Hit.GetX()-par[(Hit.GetDetectorID()-2)*2];
	    				double_t  Ycorr=Hit.GetY();
	    				double_t  Zcorr=Hit.GetZ()-par[(Hit.GetDetectorID()-2)*2+1];
	    				HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
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

	double_t Binrange=0.1;
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

	TCanvas *a=new TCanvas("a","a",1000,1000);
	MinimizerCheck(MiniPars,a);
}




void fcn_residual_tranlationRotation(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// take the translation and rotation into consideration
	// correction matrix     [0] [1] [2]
	//
	//   6 GEM detetors  5 need to be aligned
	//   GEM 2  X par[0]  Y par[1]
	//    X=[0]*x+[1]*z+[2]   1 0 0
	//    z=[3]*x+[4]*z+[5]   0 1 0
	//
	// correction matrix
	// [0] -sqrt[1-[0][0]] [1]
	// sqrt[1-[0][0]] [0]  [2]
	// 3 parameter is enough

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
			assert((Hit.GetDetectorID()==2)||(Hit.GetDetectorID()==3)||
					(Hit.GetDetectorID()==4)||(Hit.GetDetectorID()==5)||
					(Hit.GetDetectorID()==6));

			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
										   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
										   +par[(Hit.GetDetectorID()-2)*3+1];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];

/*
			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
													   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
													   +par[(Hit.GetDetectorID()-2)*3+1];
						double_t  Ycorr=Hit.GetY();
						double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];
*/

			HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);
			//HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
			gemCorrected.push_back(a);
		}
		chisq+=GetResidual(gemCorrected);
	}
	f=chisq;
}

double_t MinimizerCheck2DRT(double_t *par, TCanvas *a){
	a->Divide(2,1);
	a->cd(1)->Divide(1,2);
	a->cd(1)->cd(1);
	TH1F *residualBefore=new TH1F("ResdualBefore","ResdualBefore",100,-0.001,0.001);
	residualBefore->GetYaxis()->SetRangeUser(0,3000);
	TH1F *residualAfter=new TH1F("ResdualAfter","ResdualAfter",100,residualBefore->GetXaxis()->GetXmin(),residualBefore->GetXaxis()->GetXmax());
	residualAfter->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
										   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
										   +par[(Hit.GetDetectorID()-2)*3+1];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];

	    	HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

	    	gemCorrected.push_back(a);
		}
	    residualBefore->Fill(GetResidual(gemEvent));
	    residualAfter->Fill(GetResidual(gemCorrected));

	}
	residualBefore->Draw();
	residualAfter->Draw("same");

	// check Theta compare with VDC
	a->cd(1)->cd(2)->Divide(2,1);
	a->cd(1)->cd(2)->cd(1);
	TH2F *thetaCheckAfter=new TH2F("thetaCheckAfter","thetaCheckAfter",1000,-0.03,0.03,1000,-0.03,0.03);
	TH2F *thetaCheckBefore=new TH2F("thetaCheckBefore","thetaCheckBefore",1000,-0.03,0.03,1000,-0.03,0.03);
	thetaCheckAfter->SetMarkerColor(3);
	TH1F *thetaVDCHisto=new TH1F("vdc","vdc",100,-0.03,0.03);
	thetaVDCHisto->SetLineColor(1);
	TH1F *thetaGEMBeforeHisto=new TH1F("GEMBefore","GEMBefore",100,-0.03,0.03);
	thetaGEMBeforeHisto->SetLineColor(2);
	TH1F *thetaGEMAfterHisto=new TH1F("thetaGEMAfterHisto","thetaGEMAfterHisto",100,-0.03,0.03);
	thetaGEMAfterHisto->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
										   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
										   +par[(Hit.GetDetectorID()-2)*3+1];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];

	    	HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

	    	gemCorrected.push_back(a);
		}
	    double x[1]={0.0};
	    double slop=0.0;
	    linearInvertFit(gemCorrected,x,slop);
	    thetaCheckAfter->Fill(Event[0].GetTheta(),slop);
	    thetaGEMAfterHisto->Fill(slop);

//	    thetaCheckBefore
	    linearInvertFit(gemcorrEvent,x,slop);
	    thetaCheckBefore->Fill(Event[0].GetTheta(),slop);
	    thetaVDCHisto->Fill(Event[0].GetTheta());
	    thetaGEMBeforeHisto->Fill(slop);

	}
	thetaCheckAfter->Draw();
	thetaCheckBefore->Draw("same");
	a->cd(1)->cd(2)->cd(2);
	thetaVDCHisto->Draw();
	thetaGEMBeforeHisto->Draw("same");
	thetaGEMAfterHisto->Draw("same");

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
			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
										   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
										   +par[(Hit.GetDetectorID()-2)*3+1];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];

	    	HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

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
		deltDBeforAlign[i]=new TH1F(Form("%d_before_deltaD",i),Form("%d_before_deltaD",i),100,-0.02,0.02);
		deltDAfterAlign[i]=new TH1F(Form("%d_After_deltaD",i),Form("%d_After_deltaD",i),100,-0.02,0.02);
	}
	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected;
		gemCorrected.push_back(gemEvent[0]);
	    for(auto Hit : gemcorrEvent){
			double_t  Xcorr=Hit.GetX()*par[(Hit.GetDetectorID()-2)*3]
										   +Hit.GetZ()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])
										   +par[(Hit.GetDetectorID()-2)*3+1];
			double_t  Ycorr=Hit.GetY();
			double_t  Zcorr=-Hit.GetX()*std::sqrt(1-par[(Hit.GetDetectorID()-2)*3]*par[(Hit.GetDetectorID()-2)*3])+Hit.GetZ()*par[(Hit.GetDetectorID()-2)*3]+par[(Hit.GetDetectorID()-2)*3+2];

			HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);

//	    	HitStruct a(Hit.GetDetectorID(),Hit.GetX()-par[Hit.GetDetectorID()-2],Hit.GetY(),Hit.GetZ()-par[Hit.GetDetectorID()-1]);
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




	return 0.0;
}

// rotation tranlation 2D
void TMinimizer2DRT(){
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(15);  // initialize the minuit for 10 parameters
	gMinuit->SetFCN(fcn_residual_tranlationRotation);
	double_t vstart[15]={1,0,0,
			1,0,0,
			1,0,0,
			1,0,0,
			1,0,0};
	double_t step[15]  ={
			1e-06, 1e-06, 1e-06,
			1e-06, 1e-06,1e-06,
			1e-06, 1e-06, 1e-06,
			1e-06, 1e-06,1e-06,
			1e-06, 1e-06, 1e-06};
	double_t Binrange=0.1;
	double_t bmin[15];
	double_t bmax[15];
	for(int i = 0 ; i < 15 ; i ++){
		bmin[i]=vstart[i]-Binrange;
		bmax[i]=vstart[i]+Binrange;
	}
	double_t arglist[10];
	int ierflg=0;
	for(int i =0 ; i < 15 ; i ++){
		gMinuit->mnparm(i,Form("Par%d",i),vstart[i],step[i],bmin[i],bmax[i],ierflg);
	}
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
//	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
	gMinuit->mnsimp();
	// read out the parameters.
	double_t MiniPars[15];
	double_t MiniParsErr[15];
	for(int i =0 ; i < 15 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
	}
	TCanvas *a=new TCanvas("a","a",1000,1000);
	MinimizerCheck2DRT(MiniPars,a);
}


// a more general code used for fit
//
//
// DetCut: exclude the detector with ID number
double_t LineFit2D(std::vector<HitStruct> DetHit, double_t *x,double_t &fitpar, std::string dimension="xz",int8_t DetCut=-1){

	assert(DetHit.size()==6);

	std::vector<HitStruct> subGEM(DetHit.begin(),DetHit.end());  // get the second GEM to the end . the first one will be take as reference, will not going to the fit
	double_t Xaverage, Yaverage, Zaverage;
	double_t Xsum=0.0,Ysum=0.0, Zsum=0.0;

	if(subGEM.size()!=6){exit(-1);};

	double_t HitNum=0.0;

	for(auto Hit : subGEM){
		if(Hit.GetDetectorID()==DetCut) continue;
		Xsum+=Hit.GetX();
		Ysum+=Hit.GetY();
		Zsum+=Hit.GetZ();
		HitNum+=1.0;
	}

	Xaverage=Xsum/HitNum;
	Yaverage=Ysum/HitNum;
	Zaverage=Zsum/HitNum;

	// start for individual fit
	double_t fitSlop=0.0;

	double_t fitXZSum=0.0;
	double_t fitYZSum=0.0;
	double_t fitXYSum=0.0;
	double_t SlopXSum=0.0;
	double_t SlopYSum=0.0;
	double_t SlopZSum=0.0;

	for(auto Hit : subGEM){
		if(Hit.GetDetectorID()==DetCut) continue;
		fitXYSum+=(Hit.GetX()-Xaverage)*(Hit.GetY()-Yaverage);
		fitXZSum+=(Hit.GetX()-Xaverage)*(Hit.GetZ()-Zaverage);
		fitYZSum+=(Hit.GetY()-Yaverage)*(Hit.GetZ()-Zaverage);

		SlopXSum+=(Hit.GetX()-Xaverage)*(Hit.GetX()-Xaverage);
		SlopYSum+=(Hit.GetY()-Yaverage)*(Hit.GetY()-Yaverage);
		SlopZSum+=(Hit.GetZ()-Zaverage)*(Hit.GetZ()-Zaverage);
	}


	if (dimension == "xz"){
		//xz:  z-z0=m()(x-x0)
		fitSlop=fitXZSum/SlopXSum;
		fitpar=fitSlop;
		double_t Fita=Zaverage-Xaverage*fitSlop;
		return fitSlop*(x[0])+Fita;

	}else if (dimension == "zx") {
		//zx: x-x0=m(z-z0)
		fitSlop=fitXZSum/SlopZSum;
		fitpar=fitSlop;
		double_t Fita=Xaverage-Zaverage*fitSlop;
		return fitSlop * (x[0])+Fita;

	}else if ( dimension == "yz") {
		// YZ : Z-Z0=M(Y-Y0)
		fitSlop=fitYZSum/SlopYSum;
		fitpar=fitSlop;
		double_t Fita=Zaverage - Yaverage*fitSlop;
		return fitSlop * (x[0])+Fita;

	}else if (dimension == "zy"){
		//ZY : y-y0=m(z-z0)
		fitSlop=fitYZSum/SlopZSum;
		fitpar=fitSlop;
		double_t Fita=Yaverage-Zaverage*fitSlop;
		return fitSlop * (x[0])+Fita;
	}else{
		std::cout<<"Undefined fit"<<std::endl;
		exit(-1);
		return -1.0;
	}
}

std::vector<HitStruct> HitPosCorrection3DRT(std::vector<HitStruct> gemHit, Double_t *par){
	std::vector<HitStruct> gemHitCorrected;
	// here only use the parameter for the GEM position correction
	if(gemHit.size()==6){
		// The first GEM detector is taken to be the reference
		std::vector<HitStruct> gemcorrEvent(gemHit.begin()+1,gemHit.end());
		if(gemcorrEvent.size()!=5){
			std::cout<<__FUNCTION__<<"("<<__LINE__<<") Paramater need to be 6"<<std::endl;
			exit(-1);
		}
		gemHitCorrected.push_back(gemHit[0]);
		// start correct the matrix
		// for 3D RT correction need 30 parameters
		for(auto Hit : gemcorrEvent){
			if(Hit.GetDetectorID()<2){
				std::cout<<__FUNCTION__<<"("<<__LINE__<<") detectorID does not match"<<std::endl;
				exit(-1);
			}
			int16_t parBaseID=(Hit.GetDetectorID()-2)*6; //six paramaters on each plane
			// generate the Eular angle correction for the paramaters
			// here chose the sin(Enlar Angle) for the minimize parameters
			double_t SinAlpha=par[parBaseID+0] ; // sin(alpha)
			double_t SinBeta=par[parBaseID+1] ; // sin(beta)
			double_t SinGamma=par[parBaseID+2] ; // sin gamma

			double_t CosAlpha = std::sqrt(1- SinAlpha * SinAlpha);
			double_t CosBeta  = std::sqrt(1-SinBeta * SinBeta);
			double_t CosGamma = std::sqrt(1- SinGamma * SinGamma);

			// generate the Eular Matrix for the alignment

			double_t RotationX0 =  CosAlpha * CosGamma - CosBeta * SinAlpha * SinGamma;
			double_t RotationX1 = -CosBeta* CosGamma * SinAlpha - CosAlpha * SinGamma;
			double_t RotationX2 =  SinAlpha * SinBeta;
			double_t TranslationX=par[parBaseID+3];

			double_t RotationY0 = CosGamma * SinAlpha + CosAlpha * CosBeta * SinGamma;
			double_t RotationY1 = CosAlpha * CosBeta * CosGamma - SinAlpha * SinGamma;
			double_t RotationY2 = - CosAlpha * SinBeta;
			double_t TranslationY=par[parBaseID+4];

			double_t RotationZ0 = SinBeta * SinGamma;
			double_t RotationZ1 = CosGamma * SinBeta;
			double_t RotationZ2 = CosBeta;
			double_t TranslationZ=par[parBaseID+5];

			double_t Xcorr=Hit.GetX() * RotationX0 + Hit.GetY() * RotationX1 + Hit.GetZ() * RotationX2 + TranslationX;
			double_t Ycorr=Hit.GetX() * RotationY0 + Hit.GetY() * RotationY1 + Hit.GetZ() * RotationY2 + TranslationY;
			double_t Zcorr=Hit.GetX() * RotationZ0 + Hit.GetY() * RotationZ1 + Hit.GetZ() * RotationZ2 + TranslationZ;

			HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);
			gemHitCorrected.push_back(a);
		}
	}else{
		std::cout<<__FUNCTION__<<"("<<__LINE__<<") Paramater need to be 6"<<std::endl;
		exit(-1);
	}

	if(gemHit.size()!=gemHitCorrected.size()){
		std::cout<<__FUNCTION__<<"("<<__LINE__<<") vector size does not match !!!"<<std::endl;
		exit(-1);
	}
	return gemHitCorrected;
}



// start the 3d line fit
void fcn_3D_residual(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// fit 3D
	// with Euler angles we just need 6 parameters to discribe the rotation and translation
	// alpha [0]   0
	// beta  [1]   0
	// gamma [2]   0
	// translation X [3]  0
	// translation Y [4]   0
	// translation Z [5]   0

	double_t chisq=0.0;

	for(auto Event : DetHitBuff){  // loop on all the event and get the chisq
		std::vector<HitStruct> gemEvent(Event.begin()+1,Event.end()); // all the GEM detectors
	    HitStruct gemReference=Event[1]; // take the first GEM as reference
		std::vector<HitStruct> gemFitEvent(Event.begin()+2,Event.end());

		auto gemCorrected = HitPosCorrection3DRT(gemEvent, par);

		// fit this event
		// get the X parameter
		double_t fitM=0.0;
		double_t fitN=0.0;
		double_t fitX[1]={0.0};

		// get the zx fit
		LineFit2D(gemCorrected,fitX,fitM,"zx");
		LineFit2D(gemCorrected,fitX,fitN,"zy");

		// generate the residue
		// (m,n, 1) is the vector along the fit line
		// (Hit.x - x0, hit.y - y0, hit.z -z0) project to the line vector
		for(auto Hit : gemCorrected){ // loop on the event, and get the residual
			if(gemReference.GetDetectorID()!=1){
				exit(-1);
			}
			double delta=((Hit.GetX()-gemReference.GetX())*fitM + (Hit.GetY()-gemReference.GetY()) * fitN + (Hit.GetZ()-gemReference.GetZ()) * 1.0 )/std::sqrt( fitM * fitM + fitN * fitN +1.0);
			chisq+=delta*delta;
		}
	}
	f=chisq;
}

double_t MinimizerCheck3DRT(double_t *par, TCanvas *a){
	a->Divide(2,1);
	a->cd(1)->Divide(1,2);
	a->cd(1)->cd(1);
	TH1F *residualBefore=new TH1F("ResdualBefore","ResdualBefore",100,-10,100);
	residualBefore->GetYaxis()->SetRangeUser(0,3000);
	TH1F *residualAfter=new TH1F("ResdualAfter","ResdualAfter",100,residualBefore->GetXaxis()->GetXmin(),residualBefore->GetXaxis()->GetXmax());
	residualAfter->SetLineColor(3);

	TH2F *DetHitHistXZ;
	TH2F *DetCorrHitHistXZ;

	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());
		HitStruct gemReference=Event[1];
		std::vector<HitStruct> gemCorrected=HitPosCorrection3DRT(gemEvent,par);

	    // fit this event
		// get the X parameter
		double_t fitM=0.0;
		double_t fitN=0.0;
		double_t fitX[1]={0.0};

		// get the zx fit
		LineFit2D(gemCorrected,fitX,fitM,"zx");
		LineFit2D(gemCorrected,fitX,fitN,"zy");
		double_t chisq=0.0;
		for(auto Hit : gemCorrected){ // loop on the event, and get the residual
			 double delta=((Hit.GetX()-gemReference.GetX())*fitM + (Hit.GetY()-gemReference.GetY()) * fitN + (Hit.GetZ()-gemReference.GetZ()) * 1.0 )/std::sqrt( fitM * fitM + fitN * fitN +1.0);
			chisq+=delta*delta;
		}
	    residualAfter->Fill(chisq);
	    // get the zx fit
	    LineFit2D(gemEvent,fitX,fitM,"zx");
	    LineFit2D(gemEvent,fitX,fitN,"zy");
	    chisq=0.0;


	    //
	    std::cout<<"Fit Param(m.n) : (" << fitM <<" ,"<<fitN<<")  ";
	    for(auto Hit : gemEvent){ // loop on the event, and get the residual
			    double delta=((Hit.GetX()-gemReference.GetX())*fitM + (Hit.GetY()-gemReference.GetY()) * fitN + (Hit.GetZ()-gemReference.GetZ()) * 1.0 )/std::sqrt( fitM * fitM + fitN * fitN +1.0);
	    		chisq+=delta*delta;
	    		std::cout<<"Detector: "<<(int)Hit.GetDetectorID() <<"  ("<<Hit.GetX()<<", "<<Hit.GetY()<<", "<<Hit.GetZ()<<")   ";
	    	}
	    std::cout<<chisq<<std::endl;

	    residualBefore->Fill(chisq);
	}
	residualBefore->Draw();
	residualAfter->Draw("same");

	a->cd(2);

	for(auto Event : DetHitBuff){

			std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
			std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());
			HitStruct gemReference=Event[1];
			std::vector<HitStruct> gemCorrected=HitPosCorrection3DRT(gemEvent,par);
			DetHitHistXZ=new TH2F("XZ","XZ", 2000, -0.4, 0.4,1000, -0, 3.0);
			DetCorrHitHistXZ=new TH2F("XZcorr","XZcorr", 2000, -0.4, 0.4,1000, -0, 3.0);
		    for (auto Hit : gemEvent){
		    	DetHitHistXZ->Fill(Hit.GetX(), Hit.GetZ());
		    }
		    for(auto Hit : gemCorrected){
		    	DetCorrHitHistXZ->Fill(Hit.GetX(),Hit.GetZ());
		    }
			DetHitHistXZ->SetMarkerStyle(20);
			DetHitHistXZ->SetMarkerSize(1);
			DetHitHistXZ->Draw();
			DetCorrHitHistXZ->SetMarkerStyle(20);
			DetCorrHitHistXZ->SetMarkerSize(1);
			DetCorrHitHistXZ->SetMarkerColor(3);
			DetCorrHitHistXZ->Draw("same");
			a->Update();
			getchar();
	}
/*
	// check Theta compare with VDC
	a->cd(1)->cd(2)->Divide(2,1);
	a->cd(1)->cd(2)->cd(1);
	TH2F *thetaCheckAfter=new TH2F("thetaCheckAfter","thetaCheckAfter",1000,-0.03,0.03,1000,-0.03,0.03);
	TH2F *thetaCheckBefore=new TH2F("thetaCheckBefore","thetaCheckBefore",1000,-0.03,0.03,1000,-0.03,0.03);
	thetaCheckAfter->SetMarkerColor(3);
	TH1F *thetaVDCHisto=new TH1F("vdc","vdc",100,-0.03,0.03);
	thetaVDCHisto->SetLineColor(1);
	TH1F *thetaGEMBeforeHisto=new TH1F("GEMBefore","GEMBefore",100,-0.03,0.03);
	thetaGEMBeforeHisto->SetLineColor(2);
	TH1F *thetaGEMAfterHisto=new TH1F("thetaGEMAfterHisto","thetaGEMAfterHisto",100,-0.03,0.03);
	thetaGEMAfterHisto->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected=HitPosCorrection3DRT(gemEvent,par);

	    double x[1]={0.0};
	    double slop=0.0;
	    linearInvertFit(gemCorrected,x,slop);
	    thetaCheckAfter->Fill(Event[0].GetTheta(),slop);
	    thetaGEMAfterHisto->Fill(slop);

//	    thetaCheckBefore
	    linearInvertFit(gemcorrEvent,x,slop);
	    thetaCheckBefore->Fill(Event[0].GetTheta(),slop);
	    thetaVDCHisto->Fill(Event[0].GetTheta());
	    thetaGEMBeforeHisto->Fill(slop);

	}
	thetaCheckAfter->Draw();
	thetaCheckBefore->Draw("same");
	a->cd(1)->cd(2)->cd(2);
	thetaVDCHisto->Draw();
	thetaGEMBeforeHisto->Draw("same");
	thetaGEMAfterHisto->Draw("same");

	//distance check
	a->cd(2)->Divide(1,2);
	a->cd(2)->cd(1);
	TH1F *deltDBeforAlignSum=new TH1F("DeltaDBefore","DeltaDBefore",100,-0.1,0.1);
	TH1F *deltDAfterAlignSum=new TH1F("DeltaDAfter","DeltadAfter",100,-0.1,0.1);
	deltDAfterAlignSum->SetLineColor(3);

	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected=HitPosCorrection3DRT(gemEvent,par);

	    deltDBeforAlignSum->Fill(GetDeltD(gemEvent));
	    deltDAfterAlignSum->Fill(GetDeltD(gemCorrected));
	}
	deltDBeforAlignSum->Draw();
	deltDAfterAlignSum->Draw("same");

	a->cd(2)->cd(2)->Divide(3,2);
	TH1F *deltDBeforAlign[6];
	TH1F *deltDAfterAlign[6];
	for (int i =0 ; i < 6 ; i ++){
		deltDBeforAlign[i]=new TH1F(Form("%d_before_deltaD",i),Form("%d_before_deltaD",i),100,-0.02,0.02);
		deltDAfterAlign[i]=new TH1F(Form("%d_After_deltaD",i),Form("%d_After_deltaD",i),100,-0.02,0.02);
	}
	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());

		std::vector<HitStruct> gemCorrected=HitPosCorrection3DRT(gemEvent,par);

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
	}*/

	a->Update();

	return 0.0;
}


// Minimize for 3drt
void TMinimizer3DRT(){
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(30);  // for 3d fit. 3 Eular angle + 3 translation = 6 for each detector, 6* 5= 30
	gMinuit->SetFCN(fcn_3D_residual);  //

	double_t  vstart[30]={
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0};
	double_t step[30]={
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,};
	double_t Binrange=0.1;
	double_t bmin[30];
	double_t bmax[30];
	for(int i =0 ; i < 30; i ++){
		bmin[i]=vstart[i]-Binrange;
		bmax[i]=vstart[i]+Binrange;
	}
	double_t arglist[10];
	int ierflg=0;
	for(int i =0 ; i < 30 ; i ++){
		gMinuit->mnparm(i,Form("Par%d",i),vstart[i],step[i],bmin[i],bmax[i],ierflg);
	}
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
	//arglist[1] = 1e-5;// stop when it reach a condition. If you want to do more change 1 to 0.1
	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
	// read out the parameters.
	double_t MiniPars[30];
	double_t MiniParsErr[30];
	for(int i =0 ; i < 30 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
	}
	TCanvas *a=new TCanvas("a","a",1000,1000);
	MinimizerCheck3DRT(MiniPars,a);
}



// start fit for 3DT
// function for 3D translation only correction

std::vector<HitStruct> HitPosCorrection3DT(std::vector<HitStruct> gemHit, Double_t *par){
	// only for the 3D tranlation correction
	// x [0]
	// y [1]
	// z [3]
//	for (int i = 0 ; i < 15 ; i ++)std::cout<<par[i]<<std::endl;
	std::vector<HitStruct> gemHitCorrected;
	// here only use the parameter for the GEM position correction
	if(gemHit.size()==6){
		// The first GEM detector is taken to be the reference
		std::vector<HitStruct> gemcorrEvent(gemHit.begin()+1,gemHit.end());
		if((gemcorrEvent.size()!=5 )||gemHit[0].GetDetectorID()!=1){
			std::cout<<__FUNCTION__<<"("<<__LINE__<<") Paramater need to be 6"<<std::endl;
			exit(-1);
		}

		gemHitCorrected.push_back(gemHit[0]);
		// start correct the matrix
		// for 3D RT correction need 30 parameters
		for(auto Hit : gemcorrEvent){
			if(Hit.GetDetectorID()<2){
				std::cout<<__FUNCTION__<<"("<<__LINE__<<") detectorID does not match"<<std::endl;
				exit(-1);
			}

			int16_t parBaseID=(Hit.GetDetectorID()-2)*3; //three parameters on each plane

			double_t Xcorr=Hit.GetX() + par[parBaseID + 0];
			double_t Ycorr=Hit.GetY() + par[parBaseID + 1];
			double_t Zcorr=Hit.GetZ() + par[parBaseID + 2];
			HitStruct a(Hit.GetDetectorID(),Xcorr,Ycorr,Zcorr);
			gemHitCorrected.push_back(a);
		}
	}else{
		std::cout<<__FUNCTION__<<"("<<__LINE__<<") Paramater need to be 6"<<std::endl;
		exit(-1);
	}

	if(gemHit.size()!=gemHitCorrected.size()){
		std::cout<<__FUNCTION__<<"("<<__LINE__<<") vector size does not match !!!"<<std::endl;
		exit(-1);
	}
	return gemHitCorrected;

}

void TestFit3DT(){
	LoadDetectorHit();
	TCanvas *a=new TCanvas("a","a",1000,1000);
	a->Divide(2,2);
	TH2F *xzFit;
	TH2F *xzHist;
	TH2F *yzFit;
	TH2F *yzHist;

	TH2F *zxFit;
	TH2F *zxHist;
	TH2F *zyFit;
	TH2F *zyHist;
	for (auto Event : DetHitBuff){
		xzFit = new TH2F("xzFit", "xzFit", 2000, -0.4, 0.4, 1000, -0, 3.0);
		xzFit->SetMarkerSize(1);
		xzFit->SetMarkerStyle(20);
		xzFit->SetMarkerColor(3);
		xzHist = new TH2F("xzHist", "xzHist", 2000, -0.4, 0.4, 1000, -0, 3.0);
		xzHist->SetMarkerStyle(20);
		xzHist->SetMarkerSize(1);
		yzFit = new TH2F("yzFit", "yzFit", 2000, -0.4, 0.4, 1000, -0, 3.0);
		yzFit->SetMarkerSize(1);
		yzFit->SetMarkerStyle(20);
		yzFit->SetMarkerColor(3);
		yzHist = new TH2F("yzHist", "yzHist", 2000, -0.4, 0.4, 1000, -0, 3.0);
		yzHist->SetMarkerStyle(20);
		yzHist->SetMarkerSize(1);


		zxFit = new TH2F("xzFit", "xzFit", 1000, -0, 3.0, 2000, -0.4, 0.4);
		zxFit->SetMarkerSize(1);
		zxFit->SetMarkerStyle(20);
		zxFit->SetMarkerColor(3);
		zxHist = new TH2F("xzHist", "xzHist", 1000, -0, 3.0, 2000, -0.4, 0.4);
		zxHist->SetMarkerStyle(20);
		zxHist->SetMarkerSize(1);
		zyFit = new TH2F("yzFit", "yzFit", 1000, -0, 3.0, 2000, -0.4, 0.4);
		zyFit->SetMarkerSize(1);
		zyFit->SetMarkerStyle(20);
		zyFit->SetMarkerColor(3);
		zyHist = new TH2F("yzHist", "yzHist", 1000, -0, 3.0, 2000, -0.4, 0.4);
		zyHist->SetMarkerStyle(20);
		zyHist->SetMarkerSize(1);

		std::vector<HitStruct> gemEvent(Event.begin()+1,Event.end());
		for(auto Hit : Event){
			xzHist->Fill(Hit.GetX(),Hit.GetZ());
			double_t slop;
			double_t x[1]={Hit.GetX()};

			xzFit->Fill(Hit.GetX(),LineFit2D(gemEvent,x,slop,"xz"));

			double_t y[1]={Hit.GetY()};
			yzHist->Fill(Hit.GetY(),Hit.GetZ());
			yzFit->Fill(Hit.GetY(),LineFit2D(gemEvent,y,slop,"yz"));

			double_t z[1]={Hit.GetZ()};
			zxHist->Fill(Hit.GetZ(),Hit.GetX());
			zxFit->Fill(Hit.GetZ(),LineFit2D(gemEvent,z,slop,"zx"));

			zyHist->Fill(Hit.GetZ(),Hit.GetY());
			zyFit->Fill(Hit.GetZ(),LineFit2D(gemEvent,z,slop,"zy"));
		}
		a->cd(1);
		xzHist->Draw();
		xzFit->Draw("same");
		a->cd(2);
		yzHist->Draw();
		yzFit->Draw("same");
		a->cd(3);
		zxHist->Draw();
		zxFit->Draw("same");
		a->cd(4);
		zyFit->Draw();
		zyHist->Draw("same");

		a->Update();
		getchar();
	}
}

double GetDeltaD3D(std::vector<HitStruct> DetHit,std::string dimension="xz", int detectorID=0){
	double_t residualF=0.0;
	double_t residue=0.0;

	if(dimension == "zx"){
		for (auto Hit : DetHit){
			double_t fitM=0.0;
			double_t fitX[1]={Hit.GetZ()};
			residualF=LineFit2D(DetHit,fitX,fitM,"zx")-Hit.GetX();
			double_t deltaD=residualF/(std::sqrt(fitM*fitM+1));
			residue+=deltaD;
			if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltaD);
		}

	}else if (dimension == "zy") {
		for (auto Hit : DetHit) {
			double_t fitM = 0.0;
			double_t fitX[1] = { Hit.GetZ() };
			residualF = LineFit2D(DetHit, fitX, fitM, "zy") - Hit.GetY();
			double_t deltaD = residualF / (std::sqrt(fitM * fitM + 1));
			residue += deltaD;
			if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltaD);
		}
	}
return residue;
}

double_t GetResidual3D(std::vector<HitStruct> DetHit,std::string dimension="xz", int detectorID=0){
	// calculate the distance from the fit functions
	//get the average
	double_t residualF=0.0;
	double_t residue=0.0;

	// when culclate the residue. need to excluse the target detector when doing the fit

	if(dimension == "zx"){
		for (auto Hit : DetHit){
			double_t fitM=0.0;
			double_t fitX[1]={Hit.GetZ()};
			residualF=LineFit2D(DetHit,fitX,fitM,"zx",Hit.GetDetectorID())-Hit.GetX();
			double_t deltaD=residualF/(std::sqrt(fitM*fitM+1));
			residue+=deltaD*deltaD;
			if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltaD*deltaD);
		}

	}else if (dimension == "zy") {
		for (auto Hit : DetHit) {
			double_t fitM = 0.0;
			double_t fitX[1] = { Hit.GetZ() };
			residualF = LineFit2D(DetHit, fitX, fitM, "zy",Hit.GetDetectorID()) - Hit.GetY();
			double_t deltaD = residualF / (std::sqrt(fitM * fitM + 1));
			residue += deltaD * deltaD;
			if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltaD*deltaD);
		}
	}else if (dimension == "zx-zy") {

	}else if (dimension == "xyz") {
		double_t deltazx=0.0;
		double_t deltazy=0.0;
		for (auto Hit : DetHit){
			double_t fitM=0.0;
			double_t fitX[1]={Hit.GetZ()};
			residualF=LineFit2D(DetHit,fitX,fitM,"zx",Hit.GetDetectorID())-Hit.GetX();
			double_t deltaD=residualF/(std::sqrt(fitM*fitM+1));
			residue+=deltaD*deltaD;
			deltazx=deltaD*deltaD;

			residualF = LineFit2D(DetHit, fitX, fitM, "zy",Hit.GetDetectorID()) - Hit.GetY();
			deltaD = residualF / (std::sqrt(fitM * fitM + 1));
			residue += deltaD * deltaD;
			deltazy=deltaD * deltaD;
			if((detectorID!=0) && (Hit.GetDetectorID()==detectorID)) return (deltazx+deltazy);
		}

	}else if (dimension == "xyz-sin") {
		for (auto Hit : DetHit){
			double_t fitM=0.0;
			double_t fitX[1]={Hit.GetZ()};
			residualF=LineFit2D(DetHit,fitX,fitM,"zx",Hit.GetDetectorID())-Hit.GetX();
			double_t deltaD=residualF/(std::sqrt(fitM*fitM+1));
			residue+=deltaD;
		}
		for (auto Hit : DetHit) {
			double_t fitM = 0.0;
			double_t fitX[1] = { Hit.GetZ() };
			residualF = LineFit2D(DetHit, fitX, fitM, "zy",Hit.GetDetectorID()) - Hit.GetY();
			double_t deltaD = residualF / (std::sqrt(fitM * fitM + 1));
			residue += deltaD ;
		}
	}
return residue;
}

// start the 3d line fit
// 3D translation optimization
void fcn_3DT_residual(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// fit 3D
	// with Euler angles we just need 6 parameters to discribe the rotation and translation
	// alpha [0]   0
	// beta  [1]   0
	// gamma [2]   0
	// translation X [3]  0
	// translation Y [4]   0
	// translation Z [5]   0

	double_t chisq=0.0;

	for(auto Event : DetHitBuff){  // loop on all the event and get the chisq
		std::vector<HitStruct> gemEvent(Event.begin()+1,Event.end()); // all the GEM detectors
		std::vector<HitStruct> gemCorrected = HitPosCorrection3DT(gemEvent, par);

/*		// fit this event
		// get the X parameter
		double_t fitM=0.0;
		double_t fitN=0.0;
		double_t fitX[1]={1.0};

		// get the zx fit
		LineFit2D(gemCorrected,fitX,fitM,"zx");
		LineFit2D(gemCorrected,fitX,fitN,"zy");

		// generate the residue
		// (m,n, 1) is the vector along the fit line
		// (Hit.x - x0, hit.y - y0, hit.z -z0) project to the line vector
		for(auto Hit : gemCorrected){ // loop on the event, and get the residual
			double delta=((Hit.GetX()-gemEvent[0].GetX())*fitM + (Hit.GetY()-gemEvent[0].GetY()) * fitN + (Hit.GetZ()-gemEvent[0].GetZ()) * 1.0 )/std::sqrt( fitM * fitM + fitN * fitN +1.0);
			chisq+=delta*delta;

		}*/
		chisq+=GetResidual3D(gemCorrected,"xyz");
	}
//	std::cout<<chisq<<std::endl;
	f=chisq;
}

void fcn_3DRT_residual(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
	// fit 3D
	// with Euler angles we just need 6 parameters to discribe the rotation and translation
	// alpha [0]   0
	// beta  [1]   0
	// gamma [2]   0
	// translation X [3]  0
	// translation Y [4]   0
	// translation Z [5]   0


	double_t chisq=0.0;

	for(auto Event : DetHitBuff){  // loop on all the event and get the chisq
		std::vector<HitStruct> gemEvent(Event.begin()+1,Event.end()); // all the GEM detectors
		std::vector<HitStruct> gemCorrected = HitPosCorrection3DRT(gemEvent, par);


		// fit this event
		// get the X parameter
		double_t fitM=0.0;
		double_t fitN=0.0;
		double_t fitX[1]={1.0};

		// get the zx fit
		LineFit2D(gemCorrected,fitX,fitM,"zx");
		LineFit2D(gemCorrected,fitX,fitN,"zy");
		chisq+=GetResidual3D(gemCorrected,"xyz");
	}
	f=chisq;
}

double_t MinimizerCheck3DT(double_t *par, TCanvas *a){
	a->Divide(2,1);
	a->cd(1)->Divide(1,2);
	a->cd(1)->cd(1);
	TH1F *residualBefore=new TH1F("ResdualBefore","ResdualBefore",100,0,0.0003);
	residualBefore->GetYaxis()->SetRangeUser(0,5000);
	TH1F *residualAfter=new TH1F("ResdualAfter","ResdualAfter",100,residualBefore->GetXaxis()->GetXmin(),residualBefore->GetXaxis()->GetXmax());
	residualAfter->SetLineColor(3);



	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemCorrected=HitPosCorrection3DT(gemEvent,par);

	    residualAfter->Fill(GetResidual3D(gemCorrected,"xyz"));
	    residualBefore->Fill(GetResidual3D(gemEvent,"xyz"));
	}
	residualBefore->Draw();
	residualAfter->Draw("same");

	// check Theta compare with VDC
	a->cd(1)->cd(2)->Divide(2,1);
	a->cd(1)->cd(2)->cd(1);
	TH2F *thetaCheckAfter=new TH2F("thetaCheckAfter","thetaCheckAfter",1000,-0.03,0.03,1000,-0.03,0.03);
	TH2F *thetaCheckBefore=new TH2F("thetaCheckBefore","thetaCheckBefore",1000,-0.03,0.03,1000,-0.03,0.03);
	thetaCheckAfter->SetMarkerColor(3);
	TH1F *thetaVDCHisto=new TH1F("vdc","vdc",100,-0.03,0.03);
	thetaVDCHisto->SetLineColor(1);
	TH1F *thetaGEMBeforeHisto=new TH1F("GEMBefore","GEMBefore",100,-0.03,0.03);
	thetaGEMBeforeHisto->SetLineColor(2);
	TH1F *thetaGEMAfterHisto=new TH1F("thetaGEMAfterHisto","thetaGEMAfterHisto",100,-0.03,0.03);
	thetaGEMAfterHisto->SetLineColor(3);
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemCorrected=HitPosCorrection3DT(gemEvent,par);

	    double x[1]={0.0};
	    double slop=0.0;
	    LineFit2D(gemCorrected,x,slop,"zx");
	    thetaCheckAfter->Fill(Event[0].GetTheta(),slop);
	    thetaGEMAfterHisto->Fill(slop);

//	    thetaCheckBefore
	    LineFit2D(gemEvent,x,slop,"zx");
	    thetaCheckBefore->Fill(Event[0].GetTheta(),slop);
	    thetaVDCHisto->Fill(Event[0].GetTheta());
	    thetaGEMBeforeHisto->Fill(slop);

	}

	thetaCheckAfter->Draw();

	thetaCheckBefore->Draw("same");
	a->cd(1)->cd(2)->cd(2);
	thetaVDCHisto->GetYaxis()->SetRangeUser(0,150);
	thetaVDCHisto->Draw();
	thetaGEMBeforeHisto->Draw("same");
	thetaGEMAfterHisto->Draw("same");

	//distance check
	a->cd(2)->Divide(1,2);
	a->cd(2)->cd(1);

	//check the chisq

	TH1F *deltDBeforAlignSum=new TH1F("DeltaDBefore","DeltaDBefore",100,-0.1,100);
	TH1F *deltDAfterAlignSum=new TH1F("DeltaDAfter","DeltadAfter",100,-0.1,100);
	deltDAfterAlignSum->SetLineColor(3);

	deltDBeforAlignSum->Draw();
	deltDAfterAlignSum->Draw("same");

	a->cd(2)->cd(2)->Divide(3,2);
	TH1F *deltDBeforAlign[6];
	TH1F *deltDAfterAlign[6];
	for (int i =0 ; i < 6 ; i ++){
		deltDBeforAlign[i]=new TH1F(Form("GEM%d_before_deltaD",i),Form("GEM%d_before_deltaD",i),100,-0.005,0.005);
		deltDAfterAlign[i]=new TH1F(Form("GEM%d_After_deltaD",i),Form("GEM%d_After_deltaD",i),100,-0.005,0.005);
	}

	for(auto Event : DetHitBuff){

		std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
		std::vector<HitStruct> gemCorrected=HitPosCorrection3DT(gemEvent,par);

	    for (auto Hit : gemEvent){
	    	deltDBeforAlign[Hit.GetDetectorID()-1]->Fill(GetDeltaD3D(gemEvent,"zx",Hit.GetDetectorID()));
	    	deltDAfterAlign[Hit.GetDetectorID()-1]->Fill(GetDeltaD3D(gemCorrected,"zx",Hit.GetDetectorID()));
	    }
	    deltDBeforAlignSum->Fill(GetDeltD(gemEvent));
	    deltDAfterAlignSum->Fill(GetDeltD(gemCorrected));

	}

	for (int i =0 ; i < 6 ; i ++){
		a->cd(2)->cd(2)->cd(i+1);
		deltDAfterAlign[i]->SetLineColor(3);
		deltDAfterAlign[i]->Draw();
		deltDBeforAlign[i]->Draw("same");
	}

/*	a->cd(2);

	for(auto Event : DetHitBuff){

			std::vector<HitStruct> gemEvent(Event.begin()+1, Event.end());
			std::vector<HitStruct> gemcorrEvent(Event.begin()+2, Event.end());
			HitStruct gemReference=Event[1];
			std::vector<HitStruct> gemCorrected=HitPosCorrection3DT(gemEvent,par);
			DetHitHistXZ=new TH2F("XZ","XZ", 2000, -0.4, 0.4,1000, -0, 3.0);
			DetCorrHitHistXZ=new TH2F("XZcorr","XZcorr", 2000, -0.4, 0.4,1000, -0, 3.0);
		    for (auto Hit : gemEvent){
		    	DetHitHistXZ->Fill(Hit.GetX(), Hit.GetZ());
		    }
		    for(auto Hit : gemCorrected){
		    	DetCorrHitHistXZ->Fill(Hit.GetX(),Hit.GetZ());
		    }
			DetHitHistXZ->SetMarkerStyle(20);
			DetHitHistXZ->SetMarkerSize(1);
			DetHitHistXZ->Draw();
			DetCorrHitHistXZ->SetMarkerStyle(20);
			DetCorrHitHistXZ->SetMarkerSize(1);
			DetCorrHitHistXZ->SetMarkerColor(3);
			DetCorrHitHistXZ->Draw("same");
			a->Update();
			getchar();
	}*/
	a->Update();

	return 0.0;
}

// Minimize for 3drt
void TMinimizer3DRT_NEW(){
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(30);  // for 3d fit. 3 Eular angle + 3 translation = 6 for each detector, 6* 5= 30
	gMinuit->SetFCN(fcn_3DRT_residual);  //

	double_t  vstart[30]={
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0};
	double_t step[30]={
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6,};
	double_t Binrange=0.1;
	double_t bmin[30];
	double_t bmax[30];
	for(int i =0 ; i < 30; i ++){
		bmin[i]=vstart[i]-Binrange;
		bmax[i]=vstart[i]+Binrange;
	}
	double_t arglist[10];
	int ierflg=0;
	for(int i =0 ; i < 30 ; i ++){
		gMinuit->mnparm(i,Form("Par%d",i),vstart[i],step[i],bmin[i],bmax[i],ierflg);
	}
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
	//arglist[1] = 1e-5;// stop when it reach a condition. If you want to do more change 1 to 0.1
	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
	// read out the parameters.
	double_t MiniPars[30];
	double_t MiniParsErr[30];
	for(int i =0 ; i < 30 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
	}
	TCanvas *a=new TCanvas("a","a",3000,3000);
	MinimizerCheck3DT(MiniPars,a);
}





int fitstep=0;
// Minimize for 3drt
void TMinimizer3DT(){
	if(DetHitBuff.size()==0)
	LoadDetectorHit();      // load the raw hit from the txt file
	gStyle->SetOptFile(1111111);

	TMinuit *gMinuit=new TMinuit(15);  // for 3d fit. 3 Eular angle + 3 translation = 6 for each detector, 6* 5= 30
	gMinuit->SetFCN(fcn_3DT_residual);  //

	double_t  vstart[15]={
			0, 0, 0,
			0, 0, 0,
			0, 0, 0,
			0, 0, 0,
			0, 0, 0};
	double_t step[15]={
			1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6,
			1e-6, 1e-6, 1e-6};
	double_t Binrange=0.1;
	double_t bmin[15];
	double_t bmax[15];
	for(int i =0 ; i < 15; i ++){
		bmin[i]=vstart[i]-Binrange;
		bmax[i]=vstart[i]+Binrange;
	}
	double_t arglist[10];
	int ierflg=0;
	for(int i =0 ; i < 15 ; i ++){
		gMinuit->mnparm(i,Form("Par%d",i),vstart[i],step[i],bmin[i],bmax[i],ierflg);
	}
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
	//arglist[1] = 1e-5;// stop when it reach a condition. If you want to do more change 1 to 0.1
	gMinuit->mnexcm("MIGRAD",arglist,1,ierflg);
//	gMinuit->mnsimp();
	// read out the parameters.
	double_t MiniPars[15];
	double_t MiniParsErr[15];
	for(int i =0 ; i < 15 ; i ++){
		gMinuit->GetParameter(i,MiniPars[i],MiniParsErr[i]);
	}
	for(int i =0 ; i < 15 ; i ++){
		std::cout<<MiniPars[i]<<std::endl;
	}


	TCanvas *a=new TCanvas("a","a",3000,3000);
	MinimizerCheck3DT(MiniPars,a);
	a->Update();
	a->SaveAs(Form("Fitstep%d.jpg",fitstep++));

	std::vector<std::vector<HitStruct>> newHitBuff;
	for(auto Event : DetHitBuff){
		std::vector<HitStruct> evenHit;
		std::vector<HitStruct> gemHit(Event.begin()+1,Event.end());

		evenHit.push_back(Event[0]);
		auto corrected=HitPosCorrection3DRT(gemHit,MiniPars);
		evenHit.insert(evenHit.end(),corrected.begin(),corrected.end());
		newHitBuff.push_back(evenHit);
	}
	DetHitBuff.clear();
	DetHitBuff.insert(DetHitBuff.begin(),newHitBuff.begin(),newHitBuff.end());
//	TMinimizer3DT();


}


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

