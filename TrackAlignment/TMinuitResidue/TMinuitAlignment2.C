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
				 <<"	phi  :"<<phi<<std::endl;
	};


//	inline bool operator == (const HitStruct &x, const HitStruct &y){
//		return x.detectorID==y.detectorID;}
private:
	int8_t detectorID;
	double_t x;
	double_t y;
	double_t z;
	double_t theta;   // x'
	double_t phi;     // y'

};

// CREATE global paramaters that used for buffer the Hit for VDC and the GEMs
std::vector<HitStruct> dd; // used for buffer the distance difference between measured value and the predicted value

//
int GetAlignmentMatrix(int8_t detectorID=0,std::string fname="trackxyz.txt") {

	// VDC result ID, theta, phi, x, y ,z
	// GEM       : ID, x, y, z
	//
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
//		std::cout<<"Number of Detectors::"<<(line_elements.size()-6)/4<<std::endl;

		//get VDC values
		int8_t vdcID      = line_elements[0];
		double_t vdctheta = line_elements[1];
		double_t vdcphi   = line_elements[2];
		double_t vdcx     = line_elements[3];
		double_t vdcy     = line_elements[4];
		double_t vdcz     = line_elements[5];

		HitStruct vdcHit(vdcID,vdcx,vdcy,vdcz,vdctheta,vdcphi);
//		vdcHit.Print();
		std::vector<HitStruct> GEMHit;
		// loop on the line elements
		int8_t GEMCount=detectorID;
//		for(int8_t GEMCount=0; GEMCount < (line_elements.size()-6)/4; GEMCount++)
		{
			int8_t gemID=(int8_t)line_elements[6+GEMCount*4];
			double_t gemX=line_elements[7+GEMCount*4];
			double_t gemY=line_elements[8+GEMCount*4];
			double_t gemZ=line_elements[9+GEMCount*4];

			double_t dgemX=vdcHit.GetX()+gemZ*vdcHit.GetTheta()-gemX;
			double_t dgemY=vdcHit.GetY()+gemZ*vdcHit.GetPhi()-gemY;


			HitStruct gemHit(gemID,dgemX,dgemY,0.0);
			dd.push_back(gemHit);
//			gemHit.Print();
//			GEMHit.push_back(gemHit);
		}

		// loop on each GEM detector, Get the optimize result
	}

	return 0;
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

/*
    TH1F *residueBefore=new TH1F("alignment1","alignment1",100,0,0.001);
    TH1F *residueAfter=new TH1F("alignment","alignment",100,0,0.001);
    for(auto dhit : dd){
    	residueBefore->Fill((dhit.GetX())*(dhit.GetX())+
				(dhit.GetY())*(dhit.GetY())+
				(dhit.GetZ())*(dhit.GetZ()));
    	residueAfter->Fill((dhit.GetX()-shiftX)*(dhit.GetX()-shiftX)+
				(dhit.GetY()-shiftY)*(dhit.GetY()-shiftY)+
				(dhit.GetZ()-shiftZ)*(dhit.GetZ()-shiftZ));
    }
    residueAfter->SetLineColor(3);
    residueAfter->Draw();
    residueBefore->Draw("same");
*/
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
}
