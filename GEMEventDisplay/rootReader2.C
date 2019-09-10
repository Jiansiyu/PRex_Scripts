/*
 * rootReader2.C
 *
 *  Created on: Aug 4, 2019
 *      Author: newdriver
 */

#include <TString.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrix.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TLine.h>
#include <string>
#include <stdlib.h>
#include <cctype>
#include <numeric>
#include <TPaveText.h>
#include <TStyle.h>
//#pragma cling load("libTreeSearch-GEM.so");
//#pragma cling load("/adaqfs/home/a-onl/siyu/PRex_MPD/libprexCounting.so"); // load

#define _ROOTREADER_MAX_GEM_CHAMBER 7
#define _ROOTREADER_MAX_TSAMPLE 3


#define _ROOTREADER_MAX_GEM_CHAMBER 7
#define _ROOTREADER_MAX_TSAMPLE 3

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

std::vector<std::vector<Int_t>> splitCluster(std::vector<Int_t> StripsVec,Int_t clusterSizeCut=2){
	std::sort(StripsVec.begin(),StripsVec.end());
	std::vector<std::vector<Int_t>> result;

	if(StripsVec.size()!=0){
		std::vector<Int_t> buff;
		for(auto value=StripsVec.begin();value!=StripsVec.end();value++){
			if((value+1)!=StripsVec.end()){
				//std::cout<<*value<<"  next:"<<*value+1<<"  next in buff:"<< *(value+1);

				if((*value+1)==*(value+1)){
					//std::cout<<"  Match"<<std::endl;
					buff.push_back(*value);
				}else{
					buff.push_back(*value);
					if(buff.size()>=clusterSizeCut)
					result.push_back(buff);
					buff.clear();
				}
			}else{
				buff.push_back(*value);
				if(buff.size()>=clusterSizeCut)
				result.push_back(buff);
				buff.clear();
			}
		}

	}
	return result;
}

/// add the code for the                                     // RGEM or left GEM
void rootReader(TString fname="test_20532.root", std::string HRSarm="RGEM.rgems"){
	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);
	eventCanvas->Divide(1,2);
    eventCanvas->cd(1)->Divide(3,1);
	eventCanvas->cd(2)->Divide(3,1);
	TLine *detPlaneline[6];

	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

	// read the left HRS or the right HRS
	TFile *fileio=TFile::Open(fname.Data());
	assert(fileio);
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);
	if(PRex_GEM_tree->IsZombie()){
		std::cout<<"[Error]: can not find tree in the file !!!"<<std::endl;
	}else{
		std::cout<<"Total Entries in the file:"<< (PRex_GEM_tree->GetEntries())<<std::endl;
	}

	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::map<int16_t,std::vector<int16_t>> TsampleLit;

	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc1",HRSarm.c_str(),chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// initialize the buffers
	for (auto chamberID : chamberList){

		std::cout<<"Reading chamber :"<<chamberID<<"\n	sample:";
		//std::vector<int16_t> TsampleLit;
		for (int adc_sample =0; adc_sample <_ROOTREADER_MAX_TSAMPLE; adc_sample++){
			std::string getstring;
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc%d",HRSarm.c_str(),chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;
	}

	std::cout<<"******************  Reading Done  *********************"<<std::endl;

	// loop on the chamber and read out the data for the position

	std::string HRSarmTag="R";
	if(HRSarm=="LGEM.lgems"){
		HRSarmTag="L";
	}

	//----------------------------------------------------------------------------
	// search for the vdcs
	Int_t fEvtNum=0;
	Int_t fRun=0;
	Int_t fvdcXNum=0;
	double_t fvdcX[100];
	Int_t fvdcYNum=0;
	double_t fvdcY[100];
	double_t fvdc_th[10];
	double_t fvdc_ph[10];

	std::string fEvtNumForm("fEvtHdr.fEvtNum");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fEvtNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fEvtNumForm.c_str(),&fEvtNum);
	}else{
		std::cout<<"[Warning]:: fEvtNum data did not find in the replay resuly, skip it"<<std::endl;
	}
	std::string fRunForm("Event_Branch/fEvtHdr/fEvtHdr.fRun");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fRunForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fRunForm.c_str(),&fRun);
	}else{
		std::cout<<"[Warning]:: fRun data did not find in the replay resuly, skip it"<<std::endl;
	}

	std::string fvdcXNumForm(Form("Ndata.R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXNumForm.c_str(),&fvdcXNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcXForm(Form("R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXForm.c_str(),fvdcX);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_th_Form(Form("R.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_th_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_th_Form.c_str(),fvdc_th);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_ph_Form(Form("R.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_ph_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_ph_Form.c_str(),fvdc_ph);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdcYNumForm(Form("Ndata.R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYNumForm.c_str(),&fvdcYNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcYForm(Form("R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYForm.c_str(),fvdcY);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	// get the detector position from the root file

	// load the data from  root file

	// load the Z position
	double_t DetectorZpos[]={0.0,  0.9, 1.5982485,  1.8558464, 2.3839658, 2.5141378, 2.6395974};

	for (int chamberID=0 ; chamberID<=6 ; chamberID++){
		if(chamberID <=3){
			detPlaneline[chamberID]= new TLine(-0.1,DetectorZpos[chamberID],0.1,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}else{
			detPlaneline[chamberID]= new TLine(-0.3,DetectorZpos[chamberID],0.3,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}
	}

	//
	std::map<Int_t,Int_t> GEMNDetX;
	std::map<Int_t,Int_t> GEMNDetY;
	std::map<Int_t,double_t  *> GEMDetX;
	std::map<Int_t,double_t  *> GEMDetY;


	// read out the coord parameters
	//  coord.pos parameters
	std::map<Int_t, Int_t> GEM_NCoordPosX;
	std::map<Int_t, Int_t> GEM_NCoordPosY;
	std::map<Int_t, double *>	GEM_CoordPosX;
	std::map<Int_t, double *>	GEM_CoordPosY;

	//trkPos
	std::map<Int_t, Int_t> GEM_NCoordTrackPosX;
	std::map<Int_t, Int_t> GEM_NCoordTrackPosY;
	std::map<Int_t, double *>	GEM_CoordTrackPosX;
	std::map<Int_t, double *>	GEM_CoordTrackPosY;

	// load the GEM th-ph result

	Int_t NGEMTracktheta=0;
	Int_t NGEMTrackphi=0;
	double	GEMTracktheta[50];
	double  GEMTrackphi[50];

	// load the data for the theta-phi for GEM detectors
	std::string NGEMTracktheta_str(Form("Ndata.RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTracktheta_str.c_str(), &NGEMTracktheta);
	}
	std::string GEMTracktheta_str(Form("RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTracktheta_str.c_str(), GEMTracktheta);
	}

	std::string NGEMTrackphi_str(Form("Ndata.RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTrackphi_str.c_str(),&NGEMTrackphi);
	}
	std::string GEMTrackphi_str(Form("RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTrackphi_str.c_str(), GEMTrackphi);
	}

	// loop on the chamber list
	for (auto chamberID : chamberList){
		std::cout<<chamberID<<std::endl;

		double_t detX[20];
		double_t detY[20];
		Int_t NdetX=0;
		Int_t NdetY=0;
		double_t detZ=DetectorZpos[chamberID];

		// load the branch data
		// load the number of the Size X
		std::string NdetX_str(Form("Ndata.RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetX_str.c_str(),&GEMNDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the data
		GEMDetX[chamberID]=new double_t [100];
		std::string detX_str(Form("RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detX_str.c_str(),GEMDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the Y dimension
		std::string NdetY_str(Form("Ndata.RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetY_str.c_str(),&GEMNDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}
		GEMDetY[chamberID]=new double_t [100];
		std::string detY_str(Form("RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detY_str.c_str(),GEMDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		//finish load the data

		// start load the coordi parameter in the root file

		std::string NCoordPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordPosX_str.c_str(), &GEM_NCoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}
		GEM_CoordPosX[chamberID]=new double_t [100];
		std::string CoordPosX_str(Form("RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosX_str.c_str(), GEM_CoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordPosY_str.c_str(), &GEM_NCoordPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordPosY[chamberID]=new double_t [100];
		std::string CoordPosY_str(Form("RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosY_str.c_str(), GEM_CoordPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}


		// start load the coordi parameter in the root file
		std::string NCoordTrackPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordTrackPosX_str.c_str(), &GEM_NCoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it "<<__LINE__<<std::endl;
		}

		GEM_CoordTrackPosX[chamberID]=new double_t [100];
		std::string CoordTrackPosX_str(Form("RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosX_str.c_str(), GEM_CoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordTrackPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordTrackPosY_str.c_str(), &GEM_NCoordTrackPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordTrackPosY[chamberID]=new double_t [100];
		std::string CoordTrackPosY_str(Form("RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosY_str.c_str(), GEM_CoordTrackPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}

		// load the theta-phi result
	}

	// loop on the data
    std::vector<std::vector<HitStruct>> DetEventBuff;
	std::vector<HitStruct> DetHitArr;
	std::vector<HitStruct> DetHitArrX;
	std::vector<HitStruct> DetHitArrY;
	std::cout<<"Total Entries:"<<PRex_GEM_tree->GetEntries()<<std::endl;

	TH2F *DetHist2DXZ;//= new TH2F("XZ","XZ",2000,-0.4,0.4,1000,-0,3.0);
	TH2F *DetHist2DYZ;
	TH2F *DetHist2DXZCorr;//= new TH2F("XZ","XZ",2000,-0.4,0.4,1000,-0,3.0);
	TH2F *DetHist2DYZCorr;

	//used for plot the coord positions
	TH2F *DetCoordPosXZ;
	TH2F *DetCoordPosYZ;
	TH2F *DetCoordTrackPosXZ;
	TH2F *DetCoordTrackPosYZ;


	for(auto entry=1;entry<(PRex_GEM_tree->GetEntries()) && entry<2000;entry++)
	{

		PRex_GEM_tree->GetEntry(entry);
		//PRex_GEM_tree->GetEntry(1);
		DetHitArr.clear();
		DetHitArrX.clear();
		DetHitArrY.clear();

		// load the data to data buff
		for (auto chamberID : chamberList){
			if (GEMNDetX.find(chamberID) != GEMNDetX.end()) {
				for (int i = 0; i < GEMNDetX[chamberID]; i++) {
					HitStruct hittemp(chamberID, GEMDetX[chamberID][i], 0.0,
							DetectorZpos[chamberID]);
					DetHitArrX.push_back(hittemp);
				}
			}

			// load the Y dimension
			if (GEMNDetY.find(chamberID) != GEMNDetY.end()) {
				for (int i = 0; i < GEMNDetY[chamberID]; i++) {
					double_t y = GEMDetY[chamberID][i];
//					if(chamberID ==6) y+=0.05;
					HitStruct hittemp(chamberID, 0.0, y,
							DetectorZpos[chamberID]);
					DetHitArrY.push_back(hittemp);
				}
			}

		}
		// write all the good event to the buff
		// check the VDC
		Bool_t goodHitFlag = kTRUE;
		if((fvdcYNum==1)&&(fvdcXNum==1)){
			HitStruct hit(0,fvdcX[0],fvdcY[0],0.0,fvdc_th[0],fvdc_ph[0]);
			DetHitArr.push_back(hit);
		}else{
			goodHitFlag = kFALSE;
		}


		for (auto chamberID : chamberList) {
			if ((GEMNDetX.find(chamberID) != GEMNDetX.end())
					&& (GEMNDetY.find(chamberID) != GEMNDetY.end())) {
				if ((GEMNDetX[chamberID] == 1)
						&& (GEMNDetY[chamberID] == 1)) {

				} else {
					goodHitFlag = kFALSE;
				}
			} else {
				goodHitFlag = kFALSE;
			}
		}
		if(goodHitFlag){
			for (auto chamberID : chamberList) {
				HitStruct hit(chamberID,GEMDetX[chamberID][0],GEMDetY[chamberID][0],DetectorZpos[chamberID]);
				DetHitArr.push_back(hit);
			}
		}
		if(goodHitFlag){
			DetEventBuff.push_back(DetHitArr);
		}



		for (auto hit : DetHitArrX) hit.Print();
		for (auto hit : DetHitArrY) hit.Print();


		// finish load the data
		// load the data for display
		double_t CorrectionMatrix[]={
		0.0,
		0.0,
		0.0,
		0.00484111,
		0.00446876,
		-0.0612685,

		0.00674022,
		0.00656245,
		-0.0266333,

		0.00587216,
		0.00606263,
		0.0715993,

		0.00696359,
		0.0073248,
		0.0997759,

		0.00765592,
		0.0078796,
		0.0970067};


		DetHist2DXZ= new TH2F(Form("XZ_Evt%d",entry),Form("XZ_Evt%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetHist2DYZ= new TH2F(Form("YZ_Evt%d",entry),Form("YZ_Evt%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetHist2DXZCorr= new TH2F(Form("XZ_CEvt%d",entry),Form("XZ_CEvt%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetHist2DYZCorr= new TH2F(Form("YZ_CEvt%d",entry),Form("YZ_CEvt%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetHist2DXZ->SetMarkerSize(1);
		DetHist2DXZ->SetMarkerColor(2);
		DetHist2DXZ->SetMarkerStyle(20);

		DetHist2DYZ->SetMarkerSize(1);
		DetHist2DYZ->SetMarkerColor(2);
		DetHist2DYZ->SetMarkerStyle(20);


		DetHist2DXZCorr->SetMarkerSize(1);
		DetHist2DXZCorr->SetMarkerColor(3);
		DetHist2DXZCorr->SetMarkerStyle(20);

		DetHist2DYZCorr->SetMarkerSize(1);
		DetHist2DYZCorr->SetMarkerColor(3);
		DetHist2DYZCorr->SetMarkerStyle(20);


		std::cout<<"--------------Hit Pos----------------"<<std::endl;
		std::cout<<"====> Hit Position "<<std::endl;
		for(auto Hit : DetHitArrX){
			DetHist2DXZ->Fill(Hit.GetX(),Hit.GetZ());
			DetHist2DXZCorr->Fill(Hit.GetX()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
			Hit.Print();
		}

		for (auto Hit : DetHitArrY){
			DetHist2DYZ->Fill(Hit.GetY(),Hit.GetZ());
			DetHist2DYZCorr->Fill(Hit.GetY()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+1],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
		}

		// plot the data
		DetCoordPosXZ= new TH2F(Form("XZ_Pos%d",entry),Form("XZ_Pos%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetCoordPosYZ= new TH2F(Form("YZ_Pos%d",entry),Form("YZ_Pos%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetCoordTrackPosXZ= new TH2F(Form("XZ_CTRackPos%d",entry),Form("XZ_CTrackPos%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetCoordTrackPosYZ= new TH2F(Form("YZ_CTrackPos%d",entry),Form("YZ_CTrackPos%d",entry),2000,-0.4,0.4,1000,-0,3.0);


		DetCoordPosXZ->SetMarkerSize(1);
		DetCoordPosXZ->SetMarkerColor(4);
		DetCoordPosXZ->SetMarkerStyle(20);

		DetCoordPosYZ->SetMarkerSize(1);
		DetCoordPosYZ->SetMarkerColor(4);
		DetCoordPosYZ->SetMarkerStyle(20);


		DetCoordTrackPosXZ->SetMarkerSize(1);
		DetCoordTrackPosXZ->SetMarkerColor(4);
		DetCoordTrackPosXZ->SetMarkerStyle(20);

        DetCoordTrackPosYZ->SetMarkerSize(1);
		DetCoordTrackPosYZ->SetMarkerColor(4);
		DetCoordTrackPosYZ->SetMarkerStyle(20);

		std::cout<<"--------------Coord Pos----------------"<<std::endl;
		for (auto chamberID : chamberList){

			std::cout<<"====>"<<chamberID<<std::endl;
			if((GEM_NCoordPosX.find(chamberID)!=GEM_NCoordPosX.end())){
				for(int i = 0; i < GEM_NCoordPosX[chamberID]; i++){

				    std::cout<< i<<std::endl;
					DetCoordPosXZ->Fill(GEM_CoordPosX[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  xzPos: ("<<GEM_CoordPosX[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
			if((GEM_NCoordPosY.find(chamberID)!=GEM_NCoordPosY.end())){
				for(int i = 0; i < GEM_NCoordPosY[chamberID]; i++){
					DetCoordPosYZ->Fill(GEM_CoordPosY[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  yzPos: ("<<GEM_CoordPosY[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}

			if((GEM_NCoordTrackPosX.find(chamberID)!=GEM_NCoordTrackPosX.end())){
				for(int i = 0; i < GEM_NCoordTrackPosX[chamberID]; i++){
					DetCoordTrackPosXZ->Fill(GEM_CoordTrackPosX[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  xzTrackPos: ("<<GEM_CoordTrackPosX[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
			if((GEM_NCoordTrackPosY.find(chamberID)!=GEM_NCoordTrackPosY.end())){
				for(int i = 0; i < GEM_NCoordTrackPosY[chamberID]; i++){
					DetCoordTrackPosYZ->Fill(GEM_CoordTrackPosY[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  yzTrackPos: ("<<GEM_CoordTrackPosY[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
		std::cout<<std::endl;
		}



		eventCanvas->cd(1);


		eventCanvas->cd(1)->cd(1);
		DetHist2DXZ->Draw();
//		{
//			TLine *yztrack=new //TLine(DetHitArrX.front().GetX(),DetHitArrX.front().GetZ(),DetHitArrX.back().GetX(),DetHitArrX.back().GetZ());
//			yztrack->Draw("same");
//		}

		// draw the plane pannel
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}


		eventCanvas->cd(1)->cd(2);
		DetHist2DXZ->Draw();
		DetCoordPosXZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

        eventCanvas->cd(1)->cd(3);
		DetHist2DXZ->Draw();
		DetCoordTrackPosXZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		{
			TLine *xztrack=new TLine(fvdcX[0],0.0,fvdcX[0]+fvdc_th[0]*2.7,2.7);
			xztrack->Draw("same");
		}

		eventCanvas->cd(2)->cd(1);
		DetHist2DYZ->Draw();
//		{
//			TLine *yztrack=new TLine(fvdcY[0],0.0,fvdcY[0]+fvdc_ph[0]*2.7,2.7);
//			yztrack->Draw("same");
//		}
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		eventCanvas->cd(2)->cd(2);
		DetHist2DYZ->Draw();
		DetCoordPosYZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		eventCanvas->cd(2)->cd(3);
		DetHist2DYZ->Draw();
        DetCoordTrackPosYZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		{
			TLine *yztrack=new TLine(fvdcY[0],0.0,fvdcY[0]+fvdc_ph[0]*2.7,2.7);
			yztrack->Draw("same");
		}

		eventCanvas->Update();
		//if(goodHitFlag)

		//-------------------------------------
		std::cout<<"------------ Located the reconstructed angles -------------------"<<std::endl;
		std::cout<<"===> "<<std::endl;
		std::cout<<"vdc    ::  ";
		for(int i =0 ; i < fvdcXNum; i ++){
			std::cout<<" ("<<fvdc_th[i]<<",  "<<fvdc_ph[i]<<"),   ";
		}
		std::cout<<std::endl;

		// load the GEM detectors
        std::cout<<"===> "<<std::endl;
		std::cout<<"GEM    ::  theta-> ";
		for(int i = 0 ; i < NGEMTracktheta; i ++){
			std::cout<<"  "<< GEMTracktheta[i]<<",  ";
		}
		std::cout<<std::endl;

		std::cout<<"       ::  phi-> ";

		for (int i = 0 ; i < NGEMTrackphi; i++){
			std::cout<<"  "<<GEMTrackphi[i]<<",   ";
		}
		std::cout<<std::endl;
		//if(DetHitArrX.size()<=4)
		eventCanvas->SaveAs(Form("result/HitEvt%d.jpg",entry));
		getchar();

		DetHist2DYZ->Delete();
		DetHist2DXZ->Delete();
		DetCoordPosXZ->Delete();
		DetCoordPosYZ->Delete();
	}
	FILE *trackXYZ=fopen("trackxyz.txt","w");
	// write the data to file
	for (auto Event : DetEventBuff){
		std::string line;
		for(auto Hit : Event){
			if(Hit.GetDetectorID()==0){
				line.append(Form("%d, %f, %f, %f, %f, %f, ",Hit.GetDetectorID(),Hit.GetTheta(),Hit.GetPhi(), Hit.GetX(),Hit.GetY(),Hit.GetZ() ));
			}else{
				line.append(Form("%d, %f, %f, %f,",Hit.GetDetectorID(), Hit.GetX(),Hit.GetY(),Hit.GetZ() ));
			}
		}
		line.append("\n");
		fprintf(trackXYZ,"%s",line.c_str());
		line.clear();
	}
	fclose(trackXYZ);
}




/// add the code for the                                     // RGEM or left GEM
void hitMatchCheck(TString fname="test_20532.root", std::string HRSarm="RGEM.rgems"){
	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);
	eventCanvas->Divide(1,2);
    eventCanvas->cd(1)->Divide(3,1);
	eventCanvas->cd(2)->Divide(3,1);
	eventCanvas->Draw();

	TCanvas *hitMatchCanvas=new TCanvas("ClusterMatch","ClusterMatch",1000,1000);
	hitMatchCanvas->Divide(6,1);
	hitMatchCanvas->Draw();
	TLine *detPlaneline[6];

	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

	// read the left HRS or the right HRS
	TFile *fileio=TFile::Open(fname.Data());
	assert(fileio);
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);
	if(PRex_GEM_tree->IsZombie()){
		std::cout<<"[Error]: can not find tree in the file !!!"<<std::endl;
	}else{
		std::cout<<"Total Entries in the file:"<< (PRex_GEM_tree->GetEntries())<<std::endl;
	}

	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::map<int16_t,std::vector<int16_t>> TsampleLit;

	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc1",HRSarm.c_str(),chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// initialize the buffers
	for (auto chamberID : chamberList){

		std::cout<<"Reading chamber :"<<chamberID<<"\n	sample:";
		//std::vector<int16_t> TsampleLit;
		for (int adc_sample =0; adc_sample <_ROOTREADER_MAX_TSAMPLE; adc_sample++){
			std::string getstring;
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc%d",HRSarm.c_str(),chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;
	}

	std::cout<<"******************  Reading Done  *********************"<<std::endl;

	// loop on the chamber and read out the data for the position

	std::string HRSarmTag="R";
	if(HRSarm=="LGEM.lgems"){
		HRSarmTag="L";
	}

	//----------------------------------------------------------------------------
	// search for the vdcs
	Int_t fEvtNum=0;
	Int_t fRun=0;
	Int_t fvdcXNum=0;
	double_t fvdcX[100];
	Int_t fvdcYNum=0;
	double_t fvdcY[100];
	double_t fvdc_th[10];
	double_t fvdc_ph[10];

	std::string fEvtNumForm("fEvtHdr.fEvtNum");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fEvtNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fEvtNumForm.c_str(),&fEvtNum);
	}else{
		std::cout<<"[Warning]:: fEvtNum data did not find in the replay resuly, skip it"<<std::endl;
	}
	std::string fRunForm("Event_Branch/fEvtHdr/fEvtHdr.fRun");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fRunForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fRunForm.c_str(),&fRun);
	}else{
		std::cout<<"[Warning]:: fRun data did not find in the replay resuly, skip it"<<std::endl;
	}

	std::string fvdcXNumForm(Form("Ndata.R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXNumForm.c_str(),&fvdcXNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcXForm(Form("R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXForm.c_str(),fvdcX);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_th_Form(Form("R.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_th_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_th_Form.c_str(),fvdc_th);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_ph_Form(Form("R.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_ph_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_ph_Form.c_str(),fvdc_ph);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdcYNumForm(Form("Ndata.R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYNumForm.c_str(),&fvdcYNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcYForm(Form("R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYForm.c_str(),fvdcY);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	// get the detector position from the root file

	// load the data from  root file

	// load the Z position
	double_t DetectorZpos[]={0.0,  0.9, 1.5982485,  1.8558464, 2.3839658, 2.5141378, 2.6395974};

	for (int chamberID=0 ; chamberID<=6 ; chamberID++){
		if(chamberID <=3){
			detPlaneline[chamberID]= new TLine(-0.1,DetectorZpos[chamberID],0.1,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}else{
			detPlaneline[chamberID]= new TLine(-0.3,DetectorZpos[chamberID],0.3,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}
	}

	//
	std::map<Int_t,Int_t> GEMNDetX;
	std::map<Int_t,Int_t> GEMNDetY;
	std::map<Int_t,double_t  *> GEMDetX;
	std::map<Int_t,double_t  *> GEMDetY;


	// read out the coord parameters
	//  coord.pos parameters
	std::map<Int_t, Int_t> GEM_NCoordPosX;
	std::map<Int_t, Int_t> GEM_NCoordPosY;
	std::map<Int_t, double *>	GEM_CoordPosX;
	std::map<Int_t, double *>	GEM_CoordPosY;

	//trkPos
	std::map<Int_t, Int_t> GEM_NCoordTrackPosX;
	std::map<Int_t, Int_t> GEM_NCoordTrackPosY;
	std::map<Int_t, double *>	GEM_CoordTrackPosX;
	std::map<Int_t, double *>	GEM_CoordTrackPosY;

	// load the GEM th-ph result

	Int_t NGEMTracktheta=0;
	Int_t NGEMTrackphi=0;
	double	GEMTracktheta[50];
	double  GEMTrackphi[50];

	// load the data for the theta-phi for GEM detectors
	std::string NGEMTracktheta_str(Form("Ndata.RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTracktheta_str.c_str(), &NGEMTracktheta);
	}
	std::string GEMTracktheta_str(Form("RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTracktheta_str.c_str(), GEMTracktheta);
	}

	std::string NGEMTrackphi_str(Form("Ndata.RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTrackphi_str.c_str(),&NGEMTrackphi);
	}
	std::string GEMTrackphi_str(Form("RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTrackphi_str.c_str(), GEMTrackphi);
	}

	// loop on the chamber list
	for (auto chamberID : chamberList){
		std::cout<<chamberID<<std::endl;

		double_t detX[20];
		double_t detY[20];
		Int_t NdetX=0;
		Int_t NdetY=0;
		double_t detZ=DetectorZpos[chamberID];

		// load the branch data
		// load the number of the Size X
		std::string NdetX_str(Form("Ndata.RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetX_str.c_str(),&GEMNDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the data
		GEMDetX[chamberID]=new double_t [100];
		std::string detX_str(Form("RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detX_str.c_str(),GEMDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the Y dimension
		std::string NdetY_str(Form("Ndata.RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetY_str.c_str(),&GEMNDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}
		GEMDetY[chamberID]=new double_t [100];
		std::string detY_str(Form("RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detY_str.c_str(),GEMDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		//finish load the data

		// start load the coordi parameter in the root file

		std::string NCoordPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordPosX_str.c_str(), &GEM_NCoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}
		GEM_CoordPosX[chamberID]=new double_t [100];
		std::string CoordPosX_str(Form("RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosX_str.c_str(), GEM_CoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordPosY_str.c_str(), &GEM_NCoordPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordPosY[chamberID]=new double_t [100];
		std::string CoordPosY_str(Form("RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosY_str.c_str(), GEM_CoordPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}


		// start load the coordi parameter in the root file
		std::string NCoordTrackPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordTrackPosX_str.c_str(), &GEM_NCoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it "<<__LINE__<<std::endl;
		}

		GEM_CoordTrackPosX[chamberID]=new double_t [100];
		std::string CoordTrackPosX_str(Form("RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosX_str.c_str(), GEM_CoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordTrackPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordTrackPosY_str.c_str(), &GEM_NCoordTrackPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordTrackPosY[chamberID]=new double_t [100];
		std::string CoordTrackPosY_str(Form("RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosY_str.c_str(), GEM_CoordTrackPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}

		// load the theta-phi result
	}

	// loop on the data
    std::vector<std::vector<HitStruct>> DetEventBuff;
	std::vector<HitStruct> DetHitArr;
	std::vector<HitStruct> DetHitArrX;
	std::vector<HitStruct> DetHitArrY;
	std::cout<<"Total Entries:"<<PRex_GEM_tree->GetEntries()<<std::endl;

	TH2F *DetHist2DXZ;//= new TH2F("XZ","XZ",2000,-0.4,0.4,1000,-0,3.0);
	TH2F *DetHist2DYZ;
	TH2F *DetHist2DXZCorr;//= new TH2F("XZ","XZ",2000,-0.4,0.4,1000,-0,3.0);
	TH2F *DetHist2DYZCorr;

	//used for plot the coord positions
	TH2F *DetCoordPosXZ;
	TH2F *DetCoordPosYZ;
	TH2F *DetCoordTrackPosXZ;
	TH2F *DetCoordTrackPosYZ;


	for(auto entry=1;entry<(PRex_GEM_tree->GetEntries()) && entry<2000;entry++)
	{

		PRex_GEM_tree->GetEntry(entry);
		//PRex_GEM_tree->GetEntry(1);
		DetHitArr.clear();
		DetHitArrX.clear();
		DetHitArrY.clear();

		// load the data to data buff
		for (auto chamberID : chamberList){
			if (GEMNDetX.find(chamberID) != GEMNDetX.end()) {
				for (int i = 0; i < GEMNDetX[chamberID]; i++) {
					HitStruct hittemp(chamberID, GEMDetX[chamberID][i], 0.0,
							DetectorZpos[chamberID]);
					DetHitArrX.push_back(hittemp);
				}
			}

			// load the Y dimension
			if (GEMNDetY.find(chamberID) != GEMNDetY.end()) {
				for (int i = 0; i < GEMNDetY[chamberID]; i++) {
					double_t y = GEMDetY[chamberID][i];
//					if(chamberID ==6) y+=0.05;
					HitStruct hittemp(chamberID, 0.0, y,
							DetectorZpos[chamberID]);
					DetHitArrY.push_back(hittemp);
				}
			}

		}
		// write all the good event to the buff
		// check the VDC
		Bool_t goodHitFlag = kTRUE;
		if((fvdcYNum==1)&&(fvdcXNum==1)){
			HitStruct hit(0,fvdcX[0],fvdcY[0],0.0,fvdc_th[0],fvdc_ph[0]);
			DetHitArr.push_back(hit);
		}else{
			goodHitFlag = kFALSE;
		}


		for (auto chamberID : chamberList) {
			if ((GEMNDetX.find(chamberID) != GEMNDetX.end())
					&& (GEMNDetY.find(chamberID) != GEMNDetY.end())) {
				if ((GEMNDetX[chamberID] == 1)
						&& (GEMNDetY[chamberID] == 1)) {

				} else {
					goodHitFlag = kFALSE;
				}
			} else {
				goodHitFlag = kFALSE;
			}
		}
		if(goodHitFlag){
			for (auto chamberID : chamberList) {
				HitStruct hit(chamberID,GEMDetX[chamberID][0],GEMDetY[chamberID][0],DetectorZpos[chamberID]);
				DetHitArr.push_back(hit);
			}
		}
		if(goodHitFlag){
			DetEventBuff.push_back(DetHitArr);
		}



		for (auto hit : DetHitArrX) hit.Print();
		for (auto hit : DetHitArrY) hit.Print();


		// finish load the data
		// load the data for display
		double_t CorrectionMatrix[]={
		0.0,
		0.0,
		0.0,
		0.00484111,
		0.00446876,
		-0.0612685,

		0.00674022,
		0.00656245,
		-0.0266333,

		0.00587216,
		0.00606263,
		0.0715993,

		0.00696359,
		0.0073248,
		0.0997759,

		0.00765592,
		0.0078796,
		0.0970067};


		DetHist2DXZ= new TH2F(Form("XZ_Evt%d",entry),Form("XZ_Evt%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetHist2DYZ= new TH2F(Form("YZ_Evt%d",entry),Form("YZ_Evt%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetHist2DXZCorr= new TH2F(Form("XZ_CEvt%d",entry),Form("XZ_CEvt%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetHist2DYZCorr= new TH2F(Form("YZ_CEvt%d",entry),Form("YZ_CEvt%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetHist2DXZ->SetMarkerSize(1);
		DetHist2DXZ->SetMarkerColor(2);
		DetHist2DXZ->SetMarkerStyle(20);

		DetHist2DYZ->SetMarkerSize(1);
		DetHist2DYZ->SetMarkerColor(2);
		DetHist2DYZ->SetMarkerStyle(20);


		DetHist2DXZCorr->SetMarkerSize(1);
		DetHist2DXZCorr->SetMarkerColor(3);
		DetHist2DXZCorr->SetMarkerStyle(20);

		DetHist2DYZCorr->SetMarkerSize(1);
		DetHist2DYZCorr->SetMarkerColor(3);
		DetHist2DYZCorr->SetMarkerStyle(20);


		std::cout<<"--------------Hit Pos----------------"<<std::endl;
		std::cout<<"====> Hit Position "<<std::endl;
		for(auto Hit : DetHitArrX){
			DetHist2DXZ->Fill(Hit.GetX(),Hit.GetZ());
			DetHist2DXZCorr->Fill(Hit.GetX()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
			Hit.Print();
		}

		for (auto Hit : DetHitArrY){
			DetHist2DYZ->Fill(Hit.GetY(),Hit.GetZ());
			DetHist2DYZCorr->Fill(Hit.GetY()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+1],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
		}

		// plot the data
		DetCoordPosXZ= new TH2F(Form("XZ_Pos%d",entry),Form("XZ_Pos%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetCoordPosYZ= new TH2F(Form("YZ_Pos%d",entry),Form("YZ_Pos%d",entry),2000,-0.4,0.4,1000,-0,3.0);

		DetCoordTrackPosXZ= new TH2F(Form("XZ_CTRackPos%d",entry),Form("XZ_CTrackPos%d",entry),2000,-0.4,0.4,1000,-0,3.0);
		DetCoordTrackPosYZ= new TH2F(Form("YZ_CTrackPos%d",entry),Form("YZ_CTrackPos%d",entry),2000,-0.4,0.4,1000,-0,3.0);


		DetCoordPosXZ->SetMarkerSize(1);
		DetCoordPosXZ->SetMarkerColor(4);
		DetCoordPosXZ->SetMarkerStyle(20);

		DetCoordPosYZ->SetMarkerSize(1);
		DetCoordPosYZ->SetMarkerColor(4);
		DetCoordPosYZ->SetMarkerStyle(20);


		DetCoordTrackPosXZ->SetMarkerSize(1);
		DetCoordTrackPosXZ->SetMarkerColor(4);
		DetCoordTrackPosXZ->SetMarkerStyle(20);

        DetCoordTrackPosYZ->SetMarkerSize(1);
		DetCoordTrackPosYZ->SetMarkerColor(4);
		DetCoordTrackPosYZ->SetMarkerStyle(20);

		std::cout<<"--------------Coord Pos----------------"<<std::endl;
		for (auto chamberID : chamberList){

			std::cout<<"====>"<<chamberID<<std::endl;
			if((GEM_NCoordPosX.find(chamberID)!=GEM_NCoordPosX.end())){
				for(int i = 0; i < GEM_NCoordPosX[chamberID]; i++){

				    std::cout<< i<<std::endl;
					DetCoordPosXZ->Fill(GEM_CoordPosX[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  xzPos: ("<<GEM_CoordPosX[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
			if((GEM_NCoordPosY.find(chamberID)!=GEM_NCoordPosY.end())){
				for(int i = 0; i < GEM_NCoordPosY[chamberID]; i++){
					DetCoordPosYZ->Fill(GEM_CoordPosY[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  yzPos: ("<<GEM_CoordPosY[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}

			if((GEM_NCoordTrackPosX.find(chamberID)!=GEM_NCoordTrackPosX.end())){
				for(int i = 0; i < GEM_NCoordTrackPosX[chamberID]; i++){
					DetCoordTrackPosXZ->Fill(GEM_CoordTrackPosX[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  xzTrackPos: ("<<GEM_CoordTrackPosX[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
			if((GEM_NCoordTrackPosY.find(chamberID)!=GEM_NCoordTrackPosY.end())){
				for(int i = 0; i < GEM_NCoordTrackPosY[chamberID]; i++){
					DetCoordTrackPosYZ->Fill(GEM_CoordTrackPosY[chamberID][i],DetectorZpos[chamberID]);
					std::cout<<"  yzTrackPos: ("<<GEM_CoordTrackPosY[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
		std::cout<<std::endl;
		}



		eventCanvas->cd(1);


		eventCanvas->cd(1)->cd(1);
		DetHist2DXZ->Draw();


		// draw the plane pannel
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}


		eventCanvas->cd(1)->cd(2);
		DetHist2DXZ->Draw();
		DetCoordPosXZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

        eventCanvas->cd(1)->cd(3);
		DetHist2DXZ->Draw();
		DetCoordTrackPosXZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		{
			TLine *xztrack=new TLine(fvdcX[0],0.0,fvdcX[0]+fvdc_th[0]*2.7,2.7);
			xztrack->Draw("same");
		}

		eventCanvas->cd(2)->cd(1);
		DetHist2DYZ->Draw();

		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		eventCanvas->cd(2)->cd(2);
		DetHist2DYZ->Draw();
		DetCoordPosYZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		eventCanvas->cd(2)->cd(3);
		DetHist2DYZ->Draw();
        DetCoordTrackPosYZ->Draw("same");
		for(int i = 0 ; i <=6 ; i++){
			detPlaneline[i]->Draw("same");
		}

		{
			TLine *yztrack=new TLine(fvdcY[0],0.0,fvdcY[0]+fvdc_ph[0]*2.7,2.7);
			yztrack->Draw("same");
		}

		eventCanvas->Update();



		std::cout<<"------------ Located the reconstructed angles -------------------"<<std::endl;
		std::cout<<"===> "<<std::endl;
		// Plot the Hit and the matched Hit
		std::map<Int_t, TH2F *>gemChamberCluster2D;
		for(auto chamberID : chamberList){
			if(gemChamberCluster2D.find(chamberID)==gemChamberCluster2D.end()){
				gemChamberCluster2D[chamberID]=new TH2F(Form("gem%d_xy_matchedCluster",chamberID),Form("gem%d_xy_matchedCluster",chamberID),80,-0.4,0.4,80,-0.4,0.4);
			}

			if((GEM_NCoordPosY.find(chamberID)!=GEM_NCoordPosY.end()) && (GEM_NCoordPosX.find(chamberID)!=GEM_NCoordPosX.end())){


				for(int i = 0; (i < GEM_NCoordPosY[chamberID]) &&(i < GEM_NCoordPosX[chamberID]); i++){

					gemChamberCluster2D[chamberID]->Fill(GEM_CoordPosX[chamberID][i],GEM_CoordPosY[chamberID][i]);
					//DetCoordPosYZ->Fill(GEM_CoordPosY[chamberID][i],DetectorZpos[chamberID]);
					//std::cout<<"  yzPos: ("<<GEM_CoordPosY[chamberID][i]<<",  "<<DetectorZpos[chamberID]<<");   ";
				}
			}
			hitMatchCanvas->cd(chamberID);
			gemChamberCluster2D[chamberID]->SetMarkerSize(1);
			gemChamberCluster2D[chamberID]->SetMarkerColor(3);
			gemChamberCluster2D[chamberID]->SetMarkerStyle(20);
			gemChamberCluster2D[chamberID]->Draw();
		}


		std::map<Int_t,std::map<Int_t, TLine *>> XLine;
		std::map<Int_t,std::map<Int_t, TLine *>> YLine;
		Int_t counter_temp=0;
		for(auto Hit : DetHitArrX){
			XLine[Hit.GetDetectorID()][counter_temp]=new TLine(Hit.GetX(),-0.4,Hit.GetX(),0.4);
			hitMatchCanvas->cd(Hit.GetDetectorID());

			XLine[Hit.GetDetectorID()][counter_temp]->Draw("same");
			counter_temp++;
		}
		counter_temp=0;
		for(auto Hit : DetHitArrY){
			YLine[Hit.GetDetectorID()][counter_temp]=new TLine(-0.4,Hit.GetY(),0.4,Hit.GetY());
			hitMatchCanvas->cd(Hit.GetDetectorID());
			YLine[Hit.GetDetectorID()][counter_temp]->Draw("same");
			counter_temp++;
		}
		hitMatchCanvas->Modified();
		hitMatchCanvas->Update();



		//-------------------------------------
		std::cout<<"------------ Located the reconstructed angles -------------------"<<std::endl;
		std::cout<<"===> "<<std::endl;
		std::cout<<"vdc    ::  ";
		for(int i =0 ; i < fvdcXNum; i ++){
			std::cout<<" ("<<fvdc_th[i]<<",  "<<fvdc_ph[i]<<"),   ";
		}
		std::cout<<std::endl;

		// load the GEM detectors
        std::cout<<"===> "<<std::endl;
		std::cout<<"GEM    ::  theta-> ";
		for(int i = 0 ; i < NGEMTracktheta; i ++){
			std::cout<<"  "<< GEMTracktheta[i]<<",  ";
		}
		std::cout<<std::endl;

		std::cout<<"       ::  phi-> ";

		for (int i = 0 ; i < NGEMTrackphi; i++){
			std::cout<<"  "<<GEMTrackphi[i]<<",   ";
		}
		std::cout<<std::endl;
		//if(DetHitArrX.size()<=4)
		//eventCanvas->SaveAs(Form("result/HitEvt%d.jpg",entry));
		getchar();

		DetHist2DYZ->Delete();
		DetHist2DXZ->Delete();
		DetCoordPosXZ->Delete();
		DetCoordPosYZ->Delete();
	}


}





///---------------------------------------------------------------------------------------------------------------
///---------------------------------------------------------------------------------------------------------------
/// add the code for the                                     // RGEM or left GEM
///
///
///
void gemEfficiency(TString fname="/home/newdriver/PRex/PRex_Data/GEMRootFile/prexRHRS_20862_00_test.root", std::string HRSarm="RGEM.rgems"){
	/*TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);
	eventCanvas->Divide(1,2);
    eventCanvas->cd(1)->Divide(3,1);
	eventCanvas->cd(2)->Divide(3,1);*/
	TLine *detPlaneline[6];


	std::map<Int_t, TCanvas *> gemEffCanvas;
	std::map<Int_t, TH2F *> gemProjectedHist2D;
	std::map<Int_t, TH2F *> gemRealHist2D;
	std::map<Int_t, TH2F *> gemEfficiencyHist2D;



	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

	// read the left HRS or the right HRS
	TFile *fileio=TFile::Open(fname.Data());
	assert(fileio);
	TTree *PRex_GEM_tree;
	fileio->GetObject("T",PRex_GEM_tree);
	if(PRex_GEM_tree->IsZombie()){
		std::cout<<"[Error]: can not find tree in the file !!!"<<std::endl;
	}else{
		std::cout<<"Total Entries in the file:"<< (PRex_GEM_tree->GetEntries())<<std::endl;
	}

	// check the detector are listed in the tree
	std::vector<int16_t> chamberList;
	std::map<int16_t,std::vector<int16_t>> TsampleLit;

	std::cout<<"List of Chambers:"<<std::endl;
	for(int16_t chambercount=0; chambercount<=_ROOTREADER_MAX_GEM_CHAMBER;chambercount++){
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc1",HRSarm.c_str(),chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// for run 20862, the Z-distance is {1.161, 1.7979800, 2.0902131, 2.7165651, 2.8749137, 2.9976041}
	const double_t zpos[]={0,1.161, 1.7979800, 2.0902131, 2.7165651, 2.8749137, 2.9976041};

	const double_t xMinCut[]={0.0,  -0.0594110320,  -0.0549617020,  -0.0534276410,   -0.46155132,  -0.46307054,  -0.46354599  };
	const double_t xMaxCut[]={0.0,   0.040588968,    0.045038298,    0.046572359,     0.13844868,   0.13692946,   0.13645401  };
	const double_t yMinCut[]={0.0,  -0.096555290,   -0.091156590,   -0.089864033,    -0.27825711,  -0.27738908,  -0.22982879 };
	const double_t yMaxCut[]={0.0,   0.1034447100,   0.10884341,     0.11013597,      0.22174289,   0.22261092,   0.27017121 };

	// searching area (square)
	const double_t xCut[]={0.02,  0.02,   0.02,  0.02,  0.10,  0.10,  0.10};
	const double_t yCut[]={0.02,  0.02,   0.02,  0.02,  0.10,  0.10,  0.10};


	// innitialize the canvas
	for(auto chamberID: chamberList){
		gemEffCanvas[chamberID] = new TCanvas(Form("chamber%d_eff",chamberID),Form("chamber%d_eff",chamberID),600,600);
		gemEffCanvas[chamberID]->Divide(3,1);
		gemEffCanvas[chamberID]->Draw();

		if(chamberID<=3){
			gemProjectedHist2D[chamberID]  =  new TH2F(Form("vdcPredicted_ch%d",chamberID),Form("vdcPredicted_ch%d",chamberID),20,xMinCut[chamberID],xMaxCut[chamberID],40,yMinCut[chamberID],yMaxCut[chamberID]);
			gemRealHist2D[chamberID]       =  new TH2F(Form("gemDetected_ch%d",chamberID),Form("gemDetected_ch%d",chamberID),20,xMinCut[chamberID],xMaxCut[chamberID],40,yMinCut[chamberID],yMaxCut[chamberID]);
		}else{
			gemProjectedHist2D[chamberID]  =  new TH2F(Form("vdcPredicted_ch%d",chamberID),Form("vdcPredicted_ch%d",chamberID),120,xMinCut[chamberID],xMaxCut[chamberID],100,yMinCut[chamberID],yMaxCut[chamberID]);
			gemRealHist2D[chamberID]       =  new TH2F(Form("gemDetected_ch%d",chamberID),Form("gemDetected_ch%d",chamberID),120,xMinCut[chamberID],xMaxCut[chamberID],100,yMinCut[chamberID],yMaxCut[chamberID]);
		}
	}



	// initialize the buffers
	for (auto chamberID : chamberList){

		std::cout<<"Reading chamber :"<<chamberID<<"\n	sample:";
		//std::vector<int16_t> TsampleLit;
		for (int adc_sample =0; adc_sample <_ROOTREADER_MAX_TSAMPLE; adc_sample++){
			std::string getstring;
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("%s.x%d.adc%d",HRSarm.c_str(),chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;
	}

	std::cout<<"******************  Reading Done  *********************"<<std::endl;

	// loop on the chamber and read out the data for the position

	std::string HRSarmTag="R";
	if(HRSarm=="LGEM.lgems"){
		HRSarmTag="L";
	}

	//----------------------------------------------------------------------------
	// search for the vdcs
	Int_t fEvtNum=0;
	Int_t fRun=0;
	Int_t fvdcXNum=0;
	double_t fvdcX[100];
	Int_t fvdcYNum=0;
	double_t fvdcY[100];
	double_t fvdc_th[10];
	double_t fvdc_ph[10];

	std::string fEvtNumForm("fEvtHdr.fEvtNum");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fEvtNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fEvtNumForm.c_str(),&fEvtNum);
	}else{
		std::cout<<"[Warning]:: fEvtNum data did not find in the replay resuly, skip it"<<std::endl;
	}
	std::string fRunForm("Event_Branch/fEvtHdr/fEvtHdr.fRun");
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fRunForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fRunForm.c_str(),&fRun);
	}else{
		std::cout<<"[Warning]:: fRun data did not find in the replay resuly, skip it"<<std::endl;
	}

	std::string fvdcXNumForm(Form("Ndata.R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXNumForm.c_str(),&fvdcXNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcXForm(Form("R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXForm.c_str(),fvdcX);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_th_Form(Form("R.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_th_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_th_Form.c_str(),fvdc_th);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdc_ph_Form(Form("R.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_ph_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_ph_Form.c_str(),fvdc_ph);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::string fvdcYNumForm(Form("Ndata.R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYNumForm.c_str(),&fvdcYNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::string fvdcYForm(Form("R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYForm.c_str(),fvdcY);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	// get the detector position from the root file

	// load the data from  root file

	// load the Z position
	double_t DetectorZpos[]={0.0,  0.9, 1.5982485,  1.8558464, 2.3839658, 2.5141378, 2.6395974};

	for (int chamberID=0 ; chamberID<=6 ; chamberID++){
		if(chamberID <=3){
			detPlaneline[chamberID]= new TLine(-0.1,DetectorZpos[chamberID],0.1,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}else{
			detPlaneline[chamberID]= new TLine(-0.3,DetectorZpos[chamberID],0.3,DetectorZpos[chamberID]);
			//detPlaneline[chamberID]->Draw("same");
		}
	}

	//
	std::map<Int_t,Int_t> GEMNDetX;
	std::map<Int_t,Int_t> GEMNDetY;
	std::map<Int_t,double_t  *> GEMDetX;
	std::map<Int_t,double_t  *> GEMDetY;


	// read out the coord parameters
	//  coord.pos parameters
	std::map<Int_t, Int_t> GEM_NCoordPosX;
	std::map<Int_t, Int_t> GEM_NCoordPosY;
	std::map<Int_t, double *>	GEM_CoordPosX;
	std::map<Int_t, double *>	GEM_CoordPosY;

	//trkPos
	std::map<Int_t, Int_t> GEM_NCoordTrackPosX;
	std::map<Int_t, Int_t> GEM_NCoordTrackPosY;
	std::map<Int_t, double *>	GEM_CoordTrackPosX;
	std::map<Int_t, double *>	GEM_CoordTrackPosY;

	// load the GEM th-ph result

	Int_t NGEMTracktheta=0;
	Int_t NGEMTrackphi=0;
	double	GEMTracktheta[50];
	double  GEMTrackphi[50];

	// load the data for the theta-phi for GEM detectors
	std::string NGEMTracktheta_str(Form("Ndata.RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTracktheta_str.c_str(), &NGEMTracktheta);
	}
	std::string GEMTracktheta_str(Form("RGEM.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTracktheta_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTracktheta_str.c_str(), GEMTracktheta);
	}

	std::string NGEMTrackphi_str(Form("Ndata.RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(NGEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(NGEMTrackphi_str.c_str(),&NGEMTrackphi);
	}
	std::string GEMTrackphi_str(Form("RGEM.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(GEMTrackphi_str.c_str())){
		PRex_GEM_tree->SetBranchAddress(GEMTrackphi_str.c_str(), GEMTrackphi);
	}

	// loop on the chamber list
	for (auto chamberID : chamberList){
		std::cout<<chamberID<<std::endl;

		double_t detX[20];
		double_t detY[20];
		Int_t NdetX=0;
		Int_t NdetY=0;
		double_t detZ=DetectorZpos[chamberID];

		// load the branch data
		// load the number of the Size X
		std::string NdetX_str(Form("Ndata.RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetX_str.c_str(),&GEMNDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the data
		GEMDetX[chamberID]=new double_t [100];
		std::string detX_str(Form("RGEM.rgems.x%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detX_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detX_str.c_str(),GEMDetX[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		// load the Y dimension
		std::string NdetY_str(Form("Ndata.RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NdetY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NdetY_str.c_str(),&GEMNDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}
		GEMDetY[chamberID]=new double_t [100];
		std::string detY_str(Form("RGEM.rgems.y%d.hit.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(detY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(detY_str.c_str(),GEMDetY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
			}

		//finish load the data

		// start load the coordi parameter in the root file

		std::string NCoordPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordPosX_str.c_str(), &GEM_NCoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}
		GEM_CoordPosX[chamberID]=new double_t [100];
		std::string CoordPosX_str(Form("RGEM.rgems.x%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosX_str.c_str(), GEM_CoordPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordPosY_str.c_str(), &GEM_NCoordPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordPosY[chamberID]=new double_t [100];
		std::string CoordPosY_str(Form("RGEM.rgems.y%d.coord.pos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordPosY_str.c_str(), GEM_CoordPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}


		// start load the coordi parameter in the root file
		std::string NCoordTrackPosX_str(Form("Ndata.RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(NCoordTrackPosX_str.c_str(), &GEM_NCoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it "<<__LINE__<<std::endl;
		}

		GEM_CoordTrackPosX[chamberID]=new double_t [100];
		std::string CoordTrackPosX_str(Form("RGEM.rgems.x%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosX_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosX_str.c_str(), GEM_CoordTrackPosX[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosX data did not find in the replay resuly, skip it"<<std::endl;
		}

		std::string NCoordTrackPosY_str(Form("Ndata.RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(NCoordTrackPosY_str.c_str())){
				PRex_GEM_tree->SetBranchAddress(NCoordTrackPosY_str.c_str(), &GEM_NCoordTrackPosY[chamberID]);
		}else{
				std::cout<<"[Warning]:: NdetX data did not find in the replay resuly, skip it"<<std::endl;
		}

		GEM_CoordTrackPosY[chamberID]=new double_t [100];
		std::string CoordTrackPosY_str(Form("RGEM.rgems.y%d.coord.trkpos",chamberID));
		if(PRex_GEM_tree->GetListOfBranches()->Contains(CoordTrackPosY_str.c_str())){
			PRex_GEM_tree->SetBranchAddress(CoordTrackPosY_str.c_str(), GEM_CoordTrackPosY[chamberID]);
		}else{
			std::cout<<"[Warning]:: GEM_CoordPosY data did not find in the replay resuly, skip it"<<std::endl;
		}

		// load the theta-phi result
	}

	// loop on the data
    std::vector<std::vector<HitStruct>> DetEventBuff;
	std::vector<HitStruct> DetHitArr;
	std::vector<HitStruct> DetHitArrX;
	std::vector<HitStruct> DetHitArrY;
	std::cout<<"Total Entries:"<<PRex_GEM_tree->GetEntries()<<std::endl;

	for(auto entry=1;entry<(PRex_GEM_tree->GetEntries());entry++)
	{

		PRex_GEM_tree->GetEntry(entry);
		//PRex_GEM_tree->GetEntry(1);
		DetHitArr.clear();
		DetHitArrX.clear();
		DetHitArrY.clear();

		// load the data to data buff
		for (auto chamberID : chamberList){
			if (GEMNDetX.find(chamberID) != GEMNDetX.end()) {
				for (int i = 0; i < GEMNDetX[chamberID]; i++) {
					HitStruct hittemp(chamberID, GEMDetX[chamberID][i], 0.0,
							DetectorZpos[chamberID]);
					DetHitArrX.push_back(hittemp);
				}
			}

			// load the Y dimension
			if (GEMNDetY.find(chamberID) != GEMNDetY.end()) {
				for (int i = 0; i < GEMNDetY[chamberID]; i++) {
					double_t y = GEMDetY[chamberID][i];
//					if(chamberID ==6) y+=0.05;
					HitStruct hittemp(chamberID, 0.0, y,
							DetectorZpos[chamberID]);
					DetHitArrY.push_back(hittemp);
				}
			}

		}
		// write all the good event to the buff
		// check the VDC
		Bool_t goodHitFlag = kTRUE;
		if((fvdcYNum==1)&&(fvdcXNum==1)){
			HitStruct hit(0,fvdcX[0],fvdcY[0],0.0,fvdc_th[0],fvdc_ph[0]);
			DetHitArr.push_back(hit);
		}else{
			goodHitFlag = kFALSE;
		}


		for (auto chamberID : chamberList) {
			if ((GEMNDetX.find(chamberID) != GEMNDetX.end())
					&& (GEMNDetY.find(chamberID) != GEMNDetY.end())) {
				if ((GEMNDetX[chamberID] == 1)
						&& (GEMNDetY[chamberID] == 1)) {

				} else {
					goodHitFlag = kFALSE;
				}
			} else {
				goodHitFlag = kFALSE;
			}
		}
		if(goodHitFlag){
			for (auto chamberID : chamberList) {
				HitStruct hit(chamberID,GEMDetX[chamberID][0],GEMDetY[chamberID][0],DetectorZpos[chamberID]);
				DetHitArr.push_back(hit);
			}
		}
		if(goodHitFlag){
			DetEventBuff.push_back(DetHitArr);
		}


		// check the data within the range of GEM detector
		// check wether there is any hit candidate with in cerntain range

		//

		// loop on the GEM detector
		for(auto chamberID : chamberList){
			//check the projection within the range

			double_t vdcProjectX=fvdcX[0]+fvdc_th[0]*zpos[chamberID];
			double_t vdcProjectY=fvdcY[0]+fvdc_ph[0]*zpos[chamberID];

			if(((fvdcX[0]+fvdc_th[0]*zpos[chamberID])>=xMinCut[chamberID])&&((fvdcX[0]+fvdc_th[0]*zpos[chamberID])<=xMaxCut[chamberID])&&
					(fvdcY[0]+fvdc_ph[0]*zpos[chamberID]>=yMinCut[chamberID])&&(fvdcY[0]+fvdc_ph[0]*zpos[chamberID]<=yMaxCut[chamberID]))
			{
				// TODO

				gemProjectedHist2D[chamberID]->Fill(vdcProjectX,vdcProjectY);


				//check the GEM hit with certain range
				//loop on the hit
				Bool_t match_Flag=false;
				for(auto hitX :DetHitArrX ){
					for(auto hitY : DetHitArrY){
						if((hitY.GetDetectorID()==chamberID)&&(hitX.GetDetectorID()==chamberID)){
							if((std::abs(fvdcX[0]+fvdc_th[0]*zpos[chamberID]-hitX.GetX())<=xCut[chamberID]) &&(std::abs(fvdcY[0]+fvdc_ph[0]*zpos[chamberID]-hitY.GetY())<=yCut[chamberID])){
								// within range cut
								// TODO
								if(match_Flag==false){
									gemRealHist2D[chamberID]->Fill(vdcProjectX,vdcProjectY);
									match_Flag=true;
								}

							}
						}

					}

				}
				match_Flag=false;
			}
		}

	}

	TPaveText *pt;
	// plot the result
	for (auto chamberID : chamberList){
		gemEffCanvas[chamberID]->cd(1);

//		gemEffCanvas[chamberID]->cd(1)->SetBorderSize(20);
		gemRealHist2D[chamberID]->GetXaxis()->SetTitle("X");
		gemRealHist2D[chamberID]->GetYaxis()->SetTitle("Y");
		gemRealHist2D[chamberID]->GetXaxis()->SetLabelSize(0.03);
		gemRealHist2D[chamberID]->GetYaxis()->SetLabelSize(0.03);
		gemRealHist2D[chamberID]->GetZaxis()->SetLabelSize(0.03);
		gemRealHist2D[chamberID]->SetStats(0);
		gStyle->SetOptStat("e");
		gemRealHist2D[chamberID]->Draw("zcol");

		gemEffCanvas[chamberID]->cd(2);

		gemProjectedHist2D[chamberID]->GetXaxis()->SetTitle("X");
		gemProjectedHist2D[chamberID]->GetYaxis()->SetTitle("Y");
		gemProjectedHist2D[chamberID]->GetXaxis()->SetLabelSize(0.03);
		gemProjectedHist2D[chamberID]->GetYaxis()->SetLabelSize(0.03);
		gemProjectedHist2D[chamberID]->GetZaxis()->SetLabelSize(0.03);
		gStyle->SetOptStat("e");
		gemProjectedHist2D[chamberID]->Draw("zcol");

		gemEffCanvas[chamberID]->cd(3);
		gemEfficiencyHist2D[chamberID]=(TH2F *)gemRealHist2D[chamberID]->Clone(Form("efficiency_ch%d",chamberID));
		gemEfficiencyHist2D[chamberID]->SetTitle(Form("efficiency_ch%d",chamberID));

		gemEfficiencyHist2D[chamberID]->Divide(gemProjectedHist2D[chamberID]);
		gemEfficiencyHist2D[chamberID]->SetMinimum(0.5);

		gemEfficiencyHist2D[chamberID]->GetXaxis()->SetTitle("X");
		gemEfficiencyHist2D[chamberID]->GetYaxis()->SetTitle("Y");
		gemEfficiencyHist2D[chamberID]->GetXaxis()->SetLabelSize(0.03);
		gemEfficiencyHist2D[chamberID]->GetYaxis()->SetLabelSize(0.03);
		gemEfficiencyHist2D[chamberID]->GetZaxis()->SetLabelSize(0.03);

		for (auto binx =0 ;(binx< gemRealHist2D[chamberID]->GetXaxis()->GetNbins()); binx++){
			for (auto biny=0; biny<(gemRealHist2D[chamberID]->GetYaxis()->GetNbins());biny++){
				if(gemRealHist2D[chamberID]->GetBinContent(binx,biny)<=30){
					gemEfficiencyHist2D[chamberID]->SetBinContent(binx,biny,0.0);
				}
			}
		}

		gemEfficiencyHist2D[chamberID]->Draw("zcol");
		gStyle->SetOptStat("e");
		pt = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
		pt->AddText(Form("Efficiency=%f",(double_t)(gemRealHist2D[chamberID]->GetEntries()/(double_t)(gemProjectedHist2D[chamberID]->GetEntries()))));
		pt->SetLineColor(0);
		pt->SetFillColor(0);
		pt->SetShadowColor(0);
		pt->SetTextSize(0.04);
		pt->Draw();
	}

}















