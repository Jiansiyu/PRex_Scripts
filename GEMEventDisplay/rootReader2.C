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


#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexAnalyzer/libprexCounting.so"); // load

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


	//
	std::map<Int_t,Int_t> GEMNDetX;
	std::map<Int_t,Int_t> GEMNDetY;
	std::map<Int_t,double_t  *> GEMDetX;
	std::map<Int_t,double_t  *> GEMDetY;

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



		for(auto Hit : DetHitArrX){
			DetHist2DXZ->Fill(Hit.GetX(),Hit.GetZ());
			DetHist2DXZCorr->Fill(Hit.GetX()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
		}

		for (auto Hit : DetHitArrY){
			DetHist2DYZ->Fill(Hit.GetY(),Hit.GetZ());
			DetHist2DYZCorr->Fill(Hit.GetY()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+1],Hit.GetZ()+CorrectionMatrix[3*(Hit.GetDetectorID()-1)+2]);
		}
		eventCanvas->cd(1);
		DetHist2DXZ->Draw();
		DetHist2DXZCorr->Draw("same");

		/*{
			TLine *xztrack=new TLine(fvdcX[0],0.0,fvdcX[0]+fvdc_th[0]*3.0,3.0);
			xztrack->Draw("same");
		}*/
		{
					TLine *yztrack=new TLine(DetHitArrX.front().GetX(),DetHitArrX.front().GetZ(),DetHitArrX.back().GetX(),DetHitArrX.back().GetZ());
					yztrack->Draw("same");
		}


		eventCanvas->cd(2);
		DetHist2DYZ->Draw();
		DetHist2DYZCorr->Draw("same");
		{
			TLine *yztrack=new TLine(fvdcY[0],0.0,fvdcY[0]+fvdc_ph[0]*3.0,3.0);
			yztrack->Draw("same");
		}
		eventCanvas->Update();
		if(goodHitFlag)
		getchar();

		DetHist2DYZ->Delete();
		DetHist2DXZ->Delete();
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
