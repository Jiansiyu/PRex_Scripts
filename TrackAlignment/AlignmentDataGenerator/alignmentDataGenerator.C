/*
 *scripts used for generate the txt data for the alignment algorithm
 *scripts
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
#include <TMatrixD.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1.h>
#include <stdio.h>
#include <fstream>
#pragma cling load("libTreeSearch-GEM.so");
#pragma cling load("/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexAnalyzer/libprexCounting.so");

#define _ROOTREADER_MAX_GEM_CHAMBER 7
#define _ROOTREADER_MAX_TSAMPLE 3

//	#define __DEBUG

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


void rootReader(TString fname="test_20532.root"){



	FILE *trackXZ=fopen("trackxz.txt","w");
	FILE *trackYZ=fopen("trackyz.txt","w");
	FILE *trackXYZ=fopen("trackxyz.txt","w");
#ifdef __DEBUG
	TCanvas *eventCanvas=new TCanvas("CanvasDisplay","CanvasDisplay",1000,1000);
	eventCanvas->Divide(1,2);

	TCanvas *eventCanvas2D=new TCanvas("CanvasDisplay2D","CanvasDisplay2D",1000,1000);
	eventCanvas2D->Divide(2,2);
#endif
	if(fname.IsNull()){
		std::cout<<"Please input the file name"<<std::endl;
	}

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
		if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("RGEM.rgems.x%d.adc1",chambercount))){
			std::cout<<"	->"<< chambercount<<std::endl;
			chamberList.push_back(chambercount);
		}
	}

	// load the data
	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	// initialize the buffers
	for (auto chamberID : chamberList){

		std::cout<<"Reading chamber :"<<chamberID<<"\n	sample:";
		//std::vector<int16_t> TsampleLit;
		for (int adc_sample =0; adc_sample <_ROOTREADER_MAX_TSAMPLE; adc_sample++){
			if(PRex_GEM_tree->GetListOfBranches()->Contains(Form("RGEM.rgems.x%d.adc%d",chamberID,adc_sample))){
				std::cout<<adc_sample<<"  ";
				//TsampleLit.push_back(adc_sample);
				TsampleLit[chamberID].push_back(adc_sample);
			}
		}
		std::cout<<std::endl;
	}
	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	// allocate the memory to the tree

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

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdcXNumForm(Form("Ndata.R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXNumForm.c_str(),&fvdcXNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdcXForm(Form("R.tr.x"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcXForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcXForm.c_str(),fvdcX);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdc_th_Form(Form("R.tr.th"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_th_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_th_Form.c_str(),fvdc_th);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdc_ph_Form(Form("R.tr.ph"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdc_ph_Form.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdc_ph_Form.c_str(),fvdc_ph);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdcYNumForm(Form("Ndata.R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYNumForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYNumForm.c_str(),&fvdcYNum);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay resuly"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;

	std::string fvdcYForm(Form("R.tr.y"));
	if(PRex_GEM_tree->GetListOfBranches()->Contains(fvdcYForm.c_str())){
		PRex_GEM_tree->SetBranchAddress(fvdcYForm.c_str(),fvdcY);
	}else{
		std::cout<<"[Warning]:: VDC data did not find in the replay result"<<std::endl;
	}

	std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;
    //  chamber / value
	std::map<int16_t,double_t *>fstrip;
	std::map<int16_t,Int_t> fstripNum;
	// chamber / adcID / value
	std::map<Int_t,std::map<Int_t,Int_t>> fadcNum;
	std::map<Int_t,std::map<Int_t,double_t *>> fadc;
	std::string gem_root_header("RGEM.rgems");

	for (auto chamberID : chamberList){
		for(auto adc_sample : TsampleLit[chamberID]){
			//std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<" Working on : Chamber"<< chamberID<<"	sample"<<adc_sample<<std::endl;

			std::string fstripNumFormat(Form("Ndata.%s.x%d.strip.number",gem_root_header.c_str(),chamberID));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fstripNumFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripNumFormat.c_str(),&fstripNum[chamberID]);
			}
			//std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;
			// load the strip number informations
			fstrip[chamberID]=new double_t [1000];//[(fstripNum[chamberID])];
			std::string fstripFormat(Form("%s.x%d.strip.number",gem_root_header.c_str(),chamberID));
			if (PRex_GEM_tree->GetListOfBranches()->Contains(fstripFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripFormat.c_str(),fstrip[chamberID]);
			}

			std::string fadcNumformat(Form("Ndata.%s.x%d.adc%d",gem_root_header.c_str(),chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcNumformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcNumformat.c_str(),&fadcNum[chamberID][adc_sample]);
			}

			fadc[chamberID][adc_sample]=new double_t [5000];//[(fadcNum[chamberID][adc_sample])];
			std::string fadcformat(Form("%s.x%d.adc%d",gem_root_header.c_str(),chamberID,adc_sample));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fadcformat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fadcformat.c_str(),fadc[chamberID][adc_sample]);
			}
		}
	}

	//std::cout<<"[Debug]::"<<__FUNCTION__<<"  LINE::"<<__LINE__<<std::endl;
	// read the data for Y dimension
    //  chamber / value
	std::map<int16_t,double_t *>fstrip_y;
	std::map<int16_t,Int_t> fstripNum_y;
	// chamber / adcID / value
	std::map<Int_t,std::map<Int_t,Int_t>> fadcNum_y;
	std::map<Int_t,std::map<Int_t,double_t *>> fadc_y;
	for (auto chamberID : chamberList){
		for(auto adc_sample : TsampleLit[chamberID]){
			std::string fstripNumFormat(Form("Ndata.%s.y%d.strip.number",gem_root_header.c_str(),chamberID));
			if(PRex_GEM_tree->GetListOfBranches()->Contains(fstripNumFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripNumFormat.c_str(),&fstripNum_y[chamberID]);
			}

			fstrip_y[chamberID]=new double_t [1000];//[(fstripNum[chamberID])];
			std::string fstripFormat(Form("%s.y%d.strip.number",gem_root_header.c_str(),chamberID));
			if (PRex_GEM_tree->GetListOfBranches()->Contains(fstripFormat.c_str())){
				PRex_GEM_tree->SetBranchAddress(fstripFormat.c_str(),fstrip_y[chamberID]);
			}
		}
	}

	// read the vdc values
	std::cout<<"Total Entries:"<<PRex_GEM_tree->GetEntries()<<std::endl;
	PRex_GEM_tree->Show(1);

	for(auto entry=1;entry<(PRex_GEM_tree->GetEntries());entry++){

	// load the data to the buff
	PRex_GEM_tree->GetEntry(entry);

	// load the data to plot
	const double_t positionshift[]={0,256.0,256.0,256.0,768.0,768.0,768.0};   //  the center of the detector
	double_t positionZpos[]={0.0,  0.9, 1.5858, 1.789, 2.34135, 2.44305, 2.54};   // updated version of the z-position. take the vdc to be 0,

	const double_t positionshift_x[]={0,256.0,256.0,256.0,768.0,768.0,768.0};   //  the center of the detector
	const double_t positionshift_y[]={0,128.0,128.0,128.0,640.0,640.0,640.0};

	//---------------------------------------------------------------
	// select hit for alignment algorithm
	// requirement: all the chamber have and only have one hit
	// fit cut, fit the hit, apply cut on the X
	// require hit on x and on Y
	//---------------------------------------------------------------


	// get hit for fit

	TH2F * trackFitXZ=new TH2F("X-z","X-z",2000,-0.4,0.4,1000,-0,3.0);
	TH2F * trackFitYZ=new TH2F("Y-z","Y-z",2000,-0.4,0.4,1000,-0,3.0);

	std::map<Int_t,std::vector<TVector2>> HitArrayX_Z;
	std::map<Int_t,std::vector<TVector2>> HitArrayY_Z;
	std::vector<TVector2> HitVectX_Z;
	std::vector<TVector2> HitVectY_Z;
	HitArrayX_Z.clear();
	HitArrayY_Z.clear();
	HitVectX_Z.clear();
	HitVectY_Z.clear();

	// write the data for y-z dimension
	std::map<int16_t,TH1F *> GEMHisto_xz;
	std::map<int16_t,TH1F *> GEMHisto_yz;
	GEMHisto_xz.clear();
	GEMHisto_yz.clear();

	for(auto chamberID : chamberList){

		// add the cluster cut for the GEM detectors
		std::vector<Int_t> strips_buff_x;
		std::vector<Int_t> strips_buff_y;
		if(fstripNum[chamberID]>0){ // if there is hit on this dimension
			for(auto strips_iter=0;strips_iter<fstripNum[chamberID];strips_iter++){
				strips_buff_x.push_back(fstrip[chamberID][strips_iter]);
			}
		}

		if(fstripNum_y[chamberID]>0){
			for(auto strips_iter=0;strips_iter<fstripNum_y[chamberID];strips_iter++){
				strips_buff_y.push_back(fstrip_y[chamberID][strips_iter]);
			}
		}

		auto clusterVec_x=splitCluster(strips_buff_x);
		auto clusterVec_y=splitCluster(strips_buff_y);
		strips_buff_x.clear();
		strips_buff_y.clear();

		if(clusterVec_x.size()!=0){
			if(!(GEMHisto_xz.find(chamberID)!=GEMHisto_xz.end())){
			GEMHisto_xz[chamberID]=new TH1F(Form("chamber%d_xz",chamberID),Form("chamber%d_xz",chamberID),6/0.00004,-3,3);
			GEMHisto_xz[chamberID]->SetMarkerStyle(20);
			GEMHisto_xz[chamberID]->SetMarkerSize(1);
			}
			for(auto cluster : clusterVec_x){
				auto strip=std::accumulate(cluster.begin(),cluster.end(),0.0)/cluster.size();
				{
					double_t x=(double_t)(strip-positionshift[chamberID])*0.0004;
					double_t z=positionZpos[chamberID];

					// add shift to the detectors
					if(chamberID==1){
						x-=0.01;
					}
					if(chamberID==2){
						x-=0.005;
					}
					if(chamberID==3){
						x-=0.005;
					}
					if(chamberID==4){
						x-=0.155;
					}
					if(chamberID==5){
						x-=0.155;
					}
					if(chamberID==6){
						x-=0.155;
					}


					GEMHisto_xz[chamberID]->Fill(x,z);
					// save the hit into buffer
					// only select the event that only have one hit
//					if(clusterVec_x.size()==1)
					{
						TVector2 a(x,z);
						if(clusterVec_x.size()==1)
							HitArrayX_Z[chamberID].push_back(a);
						HitVectX_Z.push_back(a);
					}
				}
			}
		}

		if(clusterVec_y.size()!=0){
			if(!(GEMHisto_yz.find(chamberID)!=GEMHisto_yz.end())){
				GEMHisto_yz[chamberID]=new TH1F(Form("chamber%d_yz",chamberID),Form("chamber%d_yz",chamberID),6/0.00004,-3,3);
				GEMHisto_yz[chamberID]->SetMarkerStyle(20);
				GEMHisto_yz[chamberID]->SetMarkerSize(1);
			}
			for(auto cluster: clusterVec_y){
				auto strip=std::accumulate(cluster.begin(),cluster.end(),0.0)/cluster.size();
				{
					double_t y=(double_t)(strip-positionshift_y[chamberID])*0.0004;
					double_t z=positionZpos[chamberID];

					if(chamberID==4){
						y-=0.03;
					}
					if(chamberID==5){
						y-=0.028;
					}
					if(chamberID==6){
						y-=0.03;
					}
					GEMHisto_yz[chamberID]->Fill(y,z);
//					if(clusterVec_y.size()==1)
					{  // select the event only have one hit on each plane
						TVector2 a(y,z);
						if(clusterVec_y.size()==1)
							HitArrayY_Z[chamberID].push_back(a);
						HitVectY_Z.push_back(a);
					}

				}
			}
		}
	}

	HitVectX_Z.shrink_to_fit();
	HitVectY_Z.shrink_to_fit();
#ifdef __DEBUG
	TLine *beamcenter=new TLine(0,positionZpos[6],0,0);
	beamcenter->SetLineWidth(1);
	beamcenter->SetLineColor(45);

	eventCanvas->SetTitle(Form("Tracking_Detector_Run%d_Evt%d",fEvtNum,fRun));
	eventCanvas->cd(1);
#endif
	TH1F *trackingHut_xz=new TH1F("Tracking X-Z ","Tracking X-Z",2/0.00004,-1.1,1.1);
	trackingHut_xz->GetYaxis()->SetRangeUser(0,2.8);
#ifdef __DEBUG
	trackingHut_xz->Draw("histp");
//	trackingHut_xz->Delete();
	// draw beam center
	beamcenter->Draw("same");
#endif
	// draw VDC
	TLine *vdcplane_xz=new TLine(-1.059,0,1.059,0);
	vdcplane_xz->SetLineWidth(2);
#ifdef __DEBUG
	vdcplane_xz->Draw("same");
#endif

	// draw GEM planestd::map<int16_t,TLine *>
	std::map<int16_t,TLine *>GEMPlane_xz;
	for (int i =1; i <=6;i++){
		GEMPlane_xz[i]=new TLine(-positionshift_x[i]*0.0004,positionZpos[i],positionshift_x[i]*0.0004,positionZpos[i]);
		GEMPlane_xz[i]->SetLineWidth(2);
#ifdef __DEBUG

		GEMPlane_xz[i]->Draw("same");
#endif
	}

#ifdef __DEBUG
	for(auto histo=GEMHisto_xz.begin();histo!=GEMHisto_xz.end();histo++){
		histo->second->Draw("HISTPsame");
//		histo->second->Delete();
	}
#endif

	// draw VDC
	if(fvdcXNum>0){
		TH1F *vdchisto=new TH1F("vdc_xz","vdc_xz",3.8/0.00004,-2.5,1.3);
		vdchisto->Clear();
		vdchisto->Fill(fvdcX[0],0.00001);     // change to the first track
		vdchisto->SetMarkerStyle(20);
		vdchisto->SetMarkerColor(2);
		vdchisto->SetMarkerSize(1);
#ifdef __DEBUG
		vdchisto->Draw("HISTPsame");
#endif
		//		vdchisto->Delete();
		TLine *vdcTrack=new TLine(fvdcX[0],0.001,fvdcX[0]+positionZpos[6]*fvdc_th[0],positionZpos[6]);
		vdcTrack->SetLineWidth(1);
		vdcTrack->SetLineColor(6);
#ifdef __DEBUG
		vdcTrack->Draw("HISTPsame");
#endif
	}

	// draw the y-z plane
#ifdef __DEBUG
	eventCanvas->cd(2);
#endif
	TH1F *trackingHut_yz=new TH1F("Tracking Y-Z ","Tracking Y-Z",2/0.00004,-1.0,1.0);
	trackingHut_yz->GetYaxis()->SetRangeUser(0,2.8);
#ifdef __DEBUG
	trackingHut_yz->Draw("histp");

//	trackingHut_yz->Delete();
	// draw beam center
	beamcenter->Draw("same");
#endif
	// draw VDC
	TLine *vdcplane_yz=new TLine(-0.2,0,0.2,0);
	vdcplane_yz->SetLineWidth(2);
#ifdef __DEBUG
	vdcplane_yz->Draw("same");
#endif
	std::map<int16_t,TLine *>GEMPlane_yz;
	for (int i =1; i <=6;i++){
		GEMPlane_yz[i]=new TLine(-positionshift_y[i]*0.0004,positionZpos[i],positionshift_y[i]*0.0004,positionZpos[i]);
		GEMPlane_yz[i]->SetLineWidth(2);
#ifdef __DEBUG
		GEMPlane_yz[i]->Draw("same");
#endif
	}

#ifdef __DEBUG
	//draw GEM
	for(auto histo=GEMHisto_yz.begin();histo!=GEMHisto_yz.end();histo++){
		histo->second->Draw("HISTPsame");
//		histo->second->Delete();
	}

	// draw VDC
	if(fvdcYNum>0){
		TH1F *vdchisto=new TH1F("vdc_yz","vdc_yz",3.8/0.00004,-2.5,1.3);
		vdchisto->Clear();
		vdchisto->Fill(fvdcY[0],0.00001);     // change to the first track
		vdchisto->SetMarkerStyle(20);
		vdchisto->SetMarkerColor(2);
		vdchisto->SetMarkerSize(1);
		vdchisto->Draw("HISTPsame");
//		vdchisto->Delete();
		TLine *vdcTrack=new TLine(fvdcY[0],0.001,fvdcY[0]+positionZpos[6]*fvdc_ph[0],positionZpos[6]);
		vdcTrack->SetLineWidth(1);
		vdcTrack->SetLineColor(6);
		vdcTrack->Draw("Psame");
	}
	eventCanvas->Update();
#endif
	// release the mem
	for(auto hit = GEMHisto_xz.begin();hit != GEMHisto_xz.end();hit++){
		hit->second->Delete();
	}
	for(auto hit = GEMHisto_yz.begin();hit != GEMHisto_yz.end();hit++){
		hit->second->Delete();
	}
	trackingHut_xz->Delete();
	trackingHut_yz->Delete();

	// plot all the 2D staff in anothe canvas
	// all the data should have been saved in the vector
	//if(true)



	Bool_t GoodHitFlag_XZ=kFALSE;
	Bool_t GoodHitFlag_YZ=kFALSE;

	{
#ifdef __DEBUG
		eventCanvas2D->cd(1);
		TH2F *tracking2D_XZ=new TH2F("Tracking2D","Tracking2D",2000,-0.4,0.4,1000,-0,3.0);
		tracking2D_XZ->SetMarkerSize(1);
		tracking2D_XZ->SetMarkerColor(2);
		tracking2D_XZ->SetMarkerStyle(20);
		for(auto hit : HitVectX_Z){
			tracking2D_XZ->Fill(hit.X(),hit.Y());
		}
		if(fvdcXNum!=0)
		tracking2D_XZ->Fill(fvdcX[0],0.00001);
		tracking2D_XZ->Draw();
#endif

		if((HitArrayX_Z.size()==6)&& (fvdcXNum==1)){
			for (auto chamberIter=HitArrayX_Z.begin();chamberIter!=HitArrayX_Z.end();chamberIter++){
				trackFitXZ->Fill(chamberIter->second[0].X(),chamberIter->second[0].Y());
			}
			trackFitXZ->Fill(fvdcX[0],0.0);
			trackFitXZ->SetMarkerSize(1);
			trackFitXZ->SetMarkerColor(1);
			trackFitXZ->SetMarkerStyle(20);
			trackFitXZ->Fit("pol1");
#ifdef __DEBUG
			trackFitXZ->Draw("same");

			eventCanvas2D->cd(2);
#endif
			TH2F *trackFitXZ_invert=new TH2F("zx","zx",1000,-0,3.0,2000,-0.4,0.4);
			trackFitXZ_invert->Clear();
			for (auto chamberIter=HitArrayX_Z.begin();chamberIter!=HitArrayX_Z.end();chamberIter++){
				trackFitXZ_invert->Fill(chamberIter->second[0].Y(),chamberIter->second[0].X());
			}
			trackFitXZ_invert->Fill(0.0,fvdcX[0]);
			trackFitXZ_invert->SetMarkerSize(1);
			trackFitXZ_invert->SetMarkerStyle(20);
			trackFitXZ_invert->Fit("pol1");
			TF1 *lFitFunctionXZ=trackFitXZ_invert->GetFunction("pol1");
			TF1 *lFitFunctionZX=trackFitXZ->GetFunction("pol1");

			Double_t a=lFitFunctionXZ->GetChisquare() > lFitFunctionZX->GetChisquare() ? lFitFunctionZX->GetChisquare():lFitFunctionXZ->GetChisquare();

			// write candidate
			if(a<1000)
			{
				GoodHitFlag_XZ=kTRUE;
//				std::string a(Form("%d,",entry));
//				for (auto chamberIter=HitArrayX_Z.begin();chamberIter!=HitArrayX_Z.end();chamberIter++){
//					a.append(Form(" %d, %0.8f, %0.8f,",chamberIter->first,chamberIter->second[0].X(),chamberIter->second[0].Y()));
//				}
//				a.append("\n");
//				fprintf(trackXZ,"%s",a.c_str());
//				a.clear();
			}else{
				GoodHitFlag_XZ=kFALSE;
			}
			std::cout<< a<<std::endl;
#ifdef __DEBUG
			trackFitXZ_invert->Draw();
#endif
			//eventCanvas2D->Update();

			lFitFunctionXZ->Delete();
			lFitFunctionZX->Delete();

		}
#ifdef __DEBUG
		tracking2D_XZ->Delete();
#endif
	}


	// check the Y-Z dimension
	{
#ifdef __DEBUG
		eventCanvas2D->cd(3);

		TH2F *tracking2D_YZ=new TH2F("Tracking2D_YZ","Tracking2D_YZ",2000,-0.4,0.4,1000,-0,3.0);
		tracking2D_YZ->SetMarkerSize(1);
		tracking2D_YZ->SetMarkerColor(2);
		tracking2D_YZ->SetMarkerStyle(20);

		for(auto hit : HitVectY_Z){   // save the esult into the yz histogram
			tracking2D_YZ->Fill(hit.X(),hit.Y());
		}
		// if there is a hit on the vdc
		if(fvdcYNum!=0)
			tracking2D_YZ->Fill(fvdcY[0],0.0);

	    tracking2D_YZ->Draw();
#endif
		if((HitArrayY_Z.size()==6)&& (fvdcYNum==1))
		{
			for (auto chamberIter=HitArrayY_Z.begin();chamberIter!=HitArrayY_Z.end();chamberIter++){
				trackFitYZ->Fill(chamberIter->second[0].X(),chamberIter->second[0].Y());
			}
			trackFitYZ->Fill(fvdcY[0],0.0);
			trackFitYZ->SetMarkerSize(1);
			trackFitYZ->SetMarkerStyle(20);
			trackFitYZ->Fit("pol1");
#ifdef __DEBUG
			trackFitYZ->Draw("same");

			eventCanvas2D->cd(4);
#endif
			TH2F *trackFitYZ_invert=new TH2F("zY","zY",1000,-0,3.0,2000,-0.4,0.4);
			for (auto chamberIter=HitArrayY_Z.begin();chamberIter!=HitArrayY_Z.end();chamberIter++){
				trackFitYZ_invert->Fill(chamberIter->second[0].Y(),chamberIter->second[0].X());
			}
			trackFitYZ_invert->Fill(0.0,fvdcY[0]);
			trackFitYZ_invert->SetMarkerSize(1);
			trackFitYZ_invert->SetMarkerStyle(20);
			trackFitYZ_invert->Fit("pol1");
			TF1 *lFitFunctionYZ=trackFitYZ_invert->GetFunction("pol1");
			TF1 *lFitFunctionZY=trackFitYZ->GetFunction("pol1");

			Double_t a=lFitFunctionYZ->GetChisquare() > lFitFunctionZY->GetChisquare() ? lFitFunctionZY->GetChisquare():lFitFunctionYZ->GetChisquare();

			// write candidate
			if(a<1500){
				GoodHitFlag_YZ=kTRUE;
//				std::string a(Form("%d,	",entry));
//				for (auto chamberIter=HitArrayY_Z.begin();chamberIter!=HitArrayY_Z.end();chamberIter++){
//					a.append(Form("	%d,	%f,	%f,",chamberIter->first,chamberIter->second[0].X(),chamberIter->second[0].Y()));
//				}
//				a.append("\n");
//				fprintf(trackYZ,"%s",a.c_str());
//				a.clear();
			}else{
				GoodHitFlag_YZ=kFALSE;
			}

			std::cout<< a<<std::endl;
#ifdef __DEBUG
			trackFitYZ_invert->Draw();
#endif
			lFitFunctionYZ->Delete();
			lFitFunctionZY->Delete();
			trackFitYZ_invert->Clear();

		}
#ifdef __DEBUG
		tracking2D_YZ->Delete();
#endif
	}

#ifdef __DEBUG
	eventCanvas->Modified();
	eventCanvas->Update();
#endif
	//	eventCanvas2D->Modified();
//	eventCanvas2D->Update();
		if ((fvdcXNum > 0 || fvdcYNum > 0) && GoodHitFlag_XZ
				&& GoodHitFlag_YZ) {

			//std::cout << "Entry:" << entry << std::endl;
			//eventCanvas->SaveAs(Form("result/PRex_Evt%d.jpg", entry));

			// save the result in files
			// need to write the data for vdc position, angle, and the position for all the GEM detectors
			// [VDC[code 0] th ph x y z] [GEMID x y z] etc

			// write the VDC result
			std::string line;
			assert(fvdcXNum==1);
			assert(fvdcYNum==1);

			line.append(Form("0, %f, %f, %f, %f, %f, ",fvdc_th[0],fvdc_ph[0],fvdcX[0],fvdcY[0],(double_t)0.0));

			// form the GEM data
			assert(HitVectX_Z.size()==6);
			assert(HitVectY_Z.size()==6);
			std::vector<TVector3> Hit3D;
			std::map<int8_t, TVector3> Hit3DMap;
			for(auto hitXZ : HitVectX_Z){
				for(auto hitYZ : HitVectY_Z){
					if(hitXZ.Y()==hitYZ.Y()){
						TVector3 a(hitXZ.X(),hitYZ.X(),hitXZ.Y());
						Hit3D.push_back(a);
						for(int i =0; i <7;i++){
							if(positionZpos[i]==hitXZ.Y()){
							  Hit3DMap[i]=a;
							}
						}

					}
				}
			}

			// check all the six detector have hit
			assert(Hit3D.size()==6);
			assert(Hit3DMap.size()==6);
			//write the data to the string
			for(auto hitIter=Hit3DMap.begin(); hitIter != Hit3DMap.end() ; hitIter++){
				int16_t detectorID = hitIter->first;
				double_t x=hitIter->second.X();
				double_t y=hitIter->second.Y();
				double_t z=hitIter->second.Z();
				line.append(Form(" %d, %f, %f, %f, ",detectorID,x,y,z)); // detectorID, x, y, z
			}
			line.append("\n");
			fprintf(trackXYZ,"%s",line.c_str());
			line.clear();
			Hit3D.clear();
			Hit3DMap.clear();
			HitVectX_Z.clear();
			HitVectY_Z.clear();
		}
		// end of write txt file
#ifdef __DEBUG
		eventCanvas2D->Clear();
		eventCanvas->Clear();
		eventCanvas2D->Divide(2, 2);
		eventCanvas->Divide(1, 2);
#endif
		// collect memorys
		trackFitXZ->Delete();
		trackFitYZ->Delete();
#ifdef __DEBUG
		beamcenter->Delete();
#endif
//		for(auto i = GEMHisto_xz.begin(); i != GEMHisto_xz.end();i++) i->second->Delete();
//		for(auto i = GEMHisto_yz.begin(); i != GEMHisto_yz.end();i++) i->second->Delete();
		vdcplane_xz->Delete();
		HitArrayX_Z.clear();
		HitArrayY_Z.clear();
		HitVectX_Z.clear();
		HitVectY_Z.clear();

	}


	fclose(trackXZ);
	fclose(trackYZ);
	fclose(trackXYZ);
}
