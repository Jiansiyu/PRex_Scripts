/*
 * MillepedeAlignment.C
 *
 *  Created on: Jul 28, 2019
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


