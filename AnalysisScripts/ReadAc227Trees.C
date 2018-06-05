#include <fstream>

#include "PROSPECT_Style.cc"
#include "TROOT.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TPave.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TColor.h"
#include "TExec.h"
#include <sstream>
#include "TLatex.h"
#include "TMath.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include <TError.h>
#include "TChain.h"
#include <tuple>

using namespace std;

tuple<int, double, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>>ReadAc227Trees(const int NUMCELLS, int IDX, int NUMTREES, int nLoop, int PLOTFLAG, double PSDCUTLOW, double PSDCUTHIGH){

//	double ZFIDCUT = 448;
	double ZFIDCUT = 1000;

	const int NUMEXCLUDECELLS = 31;
	int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};
	bool exclude;
	gStyle->SetPalette(kTemperatureMap);


//===============================================================
//Set variables for Ac-227 tree
	
	printf("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+ \n");
	printf("Number of cells: %i \n",NUMCELLS);

	TChain AcChain("TAc");

	AcChain.Add(Form("%s/WetCommissioning/Series015_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series016_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series017_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series018_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series019_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series020_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series021_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series022_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set0.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set1.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set2.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180316_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series001_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180417_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180420_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180501_ReactorOn/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series001_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series002_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series004_AcTrees.root",getenv("P2X_ANALYZED")));


	Double_t Ac_tstamp ;
	Float_t  Ac_promptECut[2];
	Float_t  Ac_delayECut[2];
	Float_t  Ac_zCut;
	Float_t  Ac_promptPSDCut[2];
	Float_t  Ac_delayPSDCut[2];
	Long64_t Ac_evt[3];
	Int_t    Ac_seg[3];
	Float_t  Ac_E[3];
	Double_t Ac_t[3];
	Float_t  Ac_z[3];
	Float_t  Ac_PSD[3];

	AcChain.SetBranchAddress("Ac_tstamp", 		&Ac_tstamp);	//[s] in epoch time
	AcChain.SetBranchAddress("Ac_promptECut", 	Ac_promptECut);
	AcChain.SetBranchAddress("Ac_delayECut", 	Ac_delayECut);
	AcChain.SetBranchAddress("Ac_zCut", 		&Ac_zCut);
	AcChain.SetBranchAddress("Ac_promptPSDCut", Ac_promptPSDCut);
	AcChain.SetBranchAddress("Ac_delayPSDCut", 	Ac_delayPSDCut);
	AcChain.SetBranchAddress("Ac_evt",    		Ac_evt);		//evt number
	AcChain.SetBranchAddress("Ac_seg",    		Ac_seg);		//seg number
	AcChain.SetBranchAddress("Ac_E",      		Ac_E);		//[MeVee]
	AcChain.SetBranchAddress("Ac_t",      		Ac_t);		//[ns]
	AcChain.SetBranchAddress("Ac_z",      		Ac_z);		//[mm]
	AcChain.SetBranchAddress("Ac_PSD",    		Ac_PSD);		//[arb]

	Int_t nAcEvents = AcChain.GetEntries();
	printf("Number of candidate Ac-227 events: %d \n",nAcEvents);

	if(nAcEvents==0){ 
		printf("NO Ac227 EVENTS! :( STOPPING ANALYSIS \n");
		IDX = -1;
		double ts = 0.0;
		vector<double> vNull;
		vNull.push_back(0.0);
		vector<vector<double>> vvNull;
		vvNull.push_back(vNull); 
		return make_tuple(IDX, ts, 
			   vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, 
			   vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull);
	}
	else if((IDX+1)>=nAcEvents){ 
		printf("FINISHED ANALYSIS! :D \n");
		IDX = -1;
		double ts = 0.0;
		vector<double> vNull;
		vNull.push_back(0.0);
		vector<vector<double>> vvNull;
		vvNull.push_back(vNull); 
		return make_tuple(IDX, ts, 
			   vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, 
			   vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull);
	}

//===============================================================
	AcChain.GetEvent(0);

	double promptECutLow = Ac_promptECut[0], promptECutHigh = Ac_promptECut[1];	//MeVee
	double delayECutLow  = Ac_delayECut[0],  delayECutHigh  = Ac_delayECut[1];	
	
	double promptPSDCutLow = PSDCUTLOW, promptPSDCutHigh = PSDCUTHIGH;	//arb.
	double delayPSDCutLow  = PSDCUTLOW, delayPSDCutHigh = PSDCUTHIGH;

	double posCut = Ac_zCut;	//mm


	TFile *outFile = new TFile(Form("%s/Ac227Histograms.root",getenv("AD_AC227ANALYSIS_DATA_PLOTS")),"RECREATE");

//===============================================================
//Set up histograms
	printf("-------------------------------------------------------------------\n");
	printf("Initializing histograms \n");

	vector<TH1F*> vhSelectDt,        vhBGDt,        vhRnPoDt;
	vector<TH1F*> vhSelectPromptPSD, vhBGPromptPSD, vhPromptPSD;
	vector<TH1F*> vhSelectDelayPSD,  vhBGDelayPSD,  vhDelayPSD;
	vector<TH1F*> vhSelectPromptE,   vhBGPromptE,   vhPromptE;
	vector<TH1F*> vhSelectDelayE,    vhBGDelayE,    vhDelayE;
	vector<TH1F*> vhSelectDz,        vhBGDz,        vhRnPoDz;
	vector<TH1F*> vhSelectDelayPos,  vhBGDelayPos,  vhDelayPos;    	
	vector<TH2F*> vhSelectPSDvsZ,    vhBGPSDvsZ,    vhRnPoPSDvsZ;
	vector<TH2F*> vhSelectEvsZ,      vhBGEvsZ,      vhRnPoEvsZ;

	double POLIFETIME = 2.57;	//[ms]
	double TIMECUT = 0.5, TIMEWINDOW = 5.0*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;	//[ms]

	int numDtBins = 100;
	double selectDtMin = TIMECUT, selectDtMax = TIMEWINDOW;	//[ms]
	double BGDtMin = selectDtMin + TIMEOFFSET, BGDtMax = selectDtMax + TIMEOFFSET;

	int numPSDBins = 100;
	double PSDMin = 0.1, PSDMax = 0.4;	//[arb]

	int numEBins = 100;
	double EMin = 0.4, EMax = 1.2;			//[MeVee]

	int numDzBins = 100;
	double dzMin = -400.0, dzMax = 400.0;		//[mm]

	int numPosBins = 100;
	double posMin = -1200, posMax = 1200;	//[mm]

	for(int i=0;i<NUMCELLS;i++){
		TH1F* hNewSelectDt = new TH1F(Form("hSelectDt_%i",i),Form("Cell %i: Dt of Selection Events;dt [ms];Counts",i),numDtBins,selectDtMin,selectDtMax);
		TH1F* hNewBGDt     = new TH1F(Form("hBGDt_%i",i),Form("Cell %i: Dt of BG Events;dt [ms];Counts",i),numDtBins,BGDtMin,BGDtMax);
		vhSelectDt.push_back(hNewSelectDt);
		vhBGDt.push_back(hNewBGDt);

		TH1F* hNewSelectPromptPSD = new TH1F(Form("hSelectPromptPSD_%i",i),Form("Cell %i: PSD of Selection Prompt Events;PSD [arb];Counts",i),numPSDBins,PSDMin,PSDMax);
		TH1F* hNewBGPromptPSD     = new TH1F(Form("hBGPromptPSD_%i",i),Form("Cell %i: PSD of BG Prompt Events;PSD [arb];Counts",i),numPSDBins,PSDMin,PSDMax);
		TH1F* hNewSelectDelayPSD  = new TH1F(Form("hSelectDelayPSD_%i",i),Form("Cell %i: PSD of Selection Delay Events;PSD [arb];Counts",i),numPSDBins,PSDMin,PSDMax);
		TH1F* hNewBGDelayPSD      = new TH1F(Form("hBGDelayPSD_%i",i),Form("Cell %i: PSD of BG Delay Events;PSD [arb];Counts",i),numPSDBins,PSDMin,PSDMax);
		vhSelectPromptPSD.push_back(hNewSelectPromptPSD);
		vhBGPromptPSD.push_back(hNewBGPromptPSD);
		vhSelectDelayPSD.push_back(hNewSelectDelayPSD);
		vhBGDelayPSD.push_back(hNewBGDelayPSD);
	
		TH1F* hNewSelectPromptE = new TH1F(Form("hSelectPromptE_%i",i),Form("Cell %i: Energy of Selection Prompt Events;Energy [MeVee];Counts",i),numEBins,EMin,EMax);
		TH1F* hNewBGPromptE     = new TH1F(Form("hBGPromptE_%i",i),Form("Cell %i: Energy of BG Prompt Events;Energy [MeVee];Counts",i),numEBins,EMin,EMax);
		TH1F* hNewSelectDelayE  = new TH1F(Form("hSelectDelayE_%i",i),Form("Cell %i: Energy of Selection Delay Events;Energy [MeVee];Counts",i),numEBins,EMin,EMax);
		TH1F* hNewBGDelayE      = new TH1F(Form("hBGDelayE_%i",i),Form("Cell %i: Energy of BG Delay Events;Energy [MeVee];Counts",i),numEBins,EMin,EMax);
		vhSelectPromptE.push_back(hNewSelectPromptE);
		vhBGPromptE.push_back(hNewBGPromptE);
		vhSelectDelayE.push_back(hNewSelectDelayE);
		vhBGDelayE.push_back(hNewBGDelayE);

		TH1F* hNewSelectDz = new TH1F(Form("hSelectDz_%i",i),Form("Cell %i: Dz of Selection Events;dz [mm];Counts",i),numDzBins,dzMin,dzMax);
		TH1F* hNewBGDz     = new TH1F(Form("hBGDz_%i",i),Form("Cell %i: Dz of BG Events;dz [mm];Counts",i),numDzBins,dzMin,dzMax);
		vhSelectDz.push_back(hNewSelectDz);
		vhBGDz.push_back(hNewBGDz);

		TH1F* hNewSelectDelayPos = new TH1F(Form("hSelectDelayPos_%i",i),Form("Cell %i: Position of Selection Delay Events;z [mm];Counts",i),numPosBins,posMin,posMax);
		TH1F* hNewBGDelayPos     = new TH1F(Form("hBGDelayPos_%i",i),Form("Cell %i: Position of BG Delay Events;z [mm];Counts",i),numPosBins,posMin,posMax);
		vhSelectDelayPos.push_back(hNewSelectDelayPos);
		vhBGDelayPos.push_back(hNewBGDelayPos);

		TH2F* hNewSelectPSDvsZ = new TH2F(Form("hSelectPSDvsZ_%i",i),Form("Cell %i: PSD vs Position Selection Events;z [mm];PSD [arb]",i),200,posMin,posMax,200,PSDMin,PSDMax);
		TH2F* hNewBGPSDvsZ = new TH2F(Form("hBGPSDvsZ_%i",i),Form("Cell %i: PSD vs Position BG Events;z [mm];PSD [arb]",i),200,posMin,posMax,200,PSDMin,PSDMax);
		vhSelectPSDvsZ.push_back(hNewSelectPSDvsZ);
		vhBGPSDvsZ.push_back(hNewBGPSDvsZ);

		TH2F* hNewSelectEvsZ = new TH2F(Form("hSelectEvsZ_%i",i),Form("Cell %i: Energy vs Position Selection Events;z [mm];E [MeVee]",i),200,posMin,posMax,200,EMin,EMax);
		TH2F* hNewBGEvsZ = new TH2F(Form("hBGEvsZ_%i",i),Form("Cell %i: Energy vs Position BG Events;z [mm];E [MeVee]",i),200,posMin,posMax,200,EMin,EMax);
		vhSelectEvsZ.push_back(hNewSelectEvsZ);
		vhBGEvsZ.push_back(hNewBGEvsZ);
		

	}

	TH1F* hSelectDt_AllCells = new TH1F("hSelectDt_AllCells","All Cells: Dt of Selection Events;dt [ms];Counts",numDtBins,selectDtMin,selectDtMax);
	TH1F* hBGDt_AllCells     = new TH1F("hBGDt_AllCells","All Cells: Dt of BG Events;dt [ms];Counts",numDtBins,BGDtMin,BGDtMax);

	TH1F* hSelectPromptPSD_AllCells = new TH1F("hSelectPromptPSD_AllCells","All Cells: PSD of Selection Prompt Events;PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	TH1F* hBGPromptPSD_AllCells     = new TH1F("hBGPromptPSD_AllCells","All Cells: PSD of BG Prompt Events;PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	TH1F* hSelectDelayPSD_AllCells  = new TH1F("hSelectDelayPSD_AllCells","All Cells: PSD of Selection Delay Events;PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);
	TH1F* hBGDelayPSD_AllCells      = new TH1F("hBGDelayPSD_AllCells","All Cells: PSD of BG Delay Events;PSD [arb];Counts",numPSDBins,PSDMin,PSDMax);

	TH1F* hSelectPromptE_AllCells = new TH1F("hSelectPromptE_AllCells","All Cells: Energy of Selection Prompt Events;Energy [MeVee];Counts",numEBins,EMin,EMax);
	TH1F* hBGPromptE_AllCells     = new TH1F("hBGPromptE_AllCells","All Cells: Energy of BG Prompt Events;Energy [MeVee];Counts",numEBins,EMin,EMax);
	TH1F* hSelectDelayE_AllCells  = new TH1F("hSelectDelayE_AllCells","All Cells: Energy of Selection Delay Events;Energy [MeVee];Counts",numEBins,EMin,EMax);
	TH1F* hBGDelayE_AllCells      = new TH1F("hBGDelayE_AllCells","All Cells: Energy of BG Delay Events;Energy [MeVee];Counts",numEBins,EMin,EMax);

	TH1F* hSelectDz_AllCells = new TH1F("hSelectDz_AllCells","All Cells: Dz of Selection Events;dz [mm];Counts",numDzBins,dzMin,dzMax);
	TH1F* hBGDz_AllCells     = new TH1F("hBGDz_AllCells","All Cells: Dz of BG Events;dz [mm];Counts",numDzBins,dzMin,dzMax);

	TH1F* hSelectDelayPos_AllCells = new TH1F("hSelectDelayPos_AllCells","All Cells: Position of Selection Delay Events;z [mm];Counts",numPosBins,posMin,posMax);
	TH1F* hBGDelayPos_AllCells     = new TH1F("hBGDelayPos_AllCells","All Cells: Position of BG Delay Events;z [mm];Counts",numPosBins,posMin,posMax);

	printf("Histograms created \n");


//===============================================================
//Loop through trees and fill histograms
	printf("-------------------------------------------------------------------\n");
	printf("Looping over trees to fill histograms \n");
	printf("Starting at idx: %i \n",IDX);

	double CONVERTnsTOms = 1e-6;	//1*10^(-6) ms/ns
	double LargestDt = 0.0;
	int segNum; 
	double selectDt, BGDt;
	float  selectDz, BGDz;

	double firstTime = 0.0;

	double livetime = 0.0;	//ms
	double prevTime = 0.0;

	AcChain.GetEvent(IDX);
	double tstamp = Ac_tstamp;

	int countNumTrees = 0;

	TH2F* hSelectPSDvsEn = new TH2F("hSelectPSDvsEn","Selected Events: PSD vs. Energy;Energy [MeVee];PSD [arb]",150,promptECutLow,delayECutHigh,150,promptPSDCutLow,promptPSDCutHigh);
	TH2F* hBGPSDvsEn = new TH2F("hBGPSDvsEn","BG Events: PSD vs. Energy;Energy [MeVee];PSD [arb]",150,promptECutLow,delayECutHigh,150,promptPSDCutLow,promptPSDCutHigh);
	TH2F* hSelectDelayEnVsPromptEn = new TH2F("hSelectDelayEnVsPromptEn","Selected Events: Delay vs. Prompt Energy;Energy [MeVee];Energy [MeVee]",150,promptECutLow,delayECutHigh,150,promptECutLow,delayECutHigh);
	TH2F* hBGDelayEnVsPromptEn = new TH2F("hBGDelayEnVsPromptEn","BG Events: Delay vs. Prompt Energy;Energy [MeVee];Energy [MeVee]",150,promptECutLow,delayECutHigh,150,promptECutLow,delayECutHigh);
	
	TH1F* hRnPoSelectSeg = new TH1F("hRnPoSelectSeg","RnPo Event Segment;Segment;Counts",154,0,154);
	TH1F* hRnPoBGSeg = new TH1F("hRnPoBGSeg","RnPo Event Segment;Segment;Counts",154,0,154);

	double promptTime, promptPSD;

	for(int i=IDX;i<nAcEvents;i++){
		if(i%10000==0) printf("Event: %d \n",i);
		
		AcChain.GetEvent(i);

		promptPSD = Ac_PSD[0];
		if(promptPSD<PSDCUTLOW || promptPSD>PSDCUTHIGH || Ac_seg[0]>200) goto skipEntry;

		promptTime = Ac_t[0]*CONVERTnsTOms;

		//if starting to look at a new tree calculate livetime for previous tree
		if(promptTime<prevTime){
			livetime = livetime + prevTime;
			countNumTrees = countNumTrees + 1;

			if(countNumTrees==NUMTREES){
				IDX = i;
				printf("Finished bin. Number of trees: %d \n",countNumTrees);
				break;
			}
		}

		prevTime = promptTime;	

		segNum = Ac_seg[0];
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), segNum) != end(ExcludeCellArr);


		//check if there is a selected delay event
		if(Ac_evt[1]!=0 && !exclude && Ac_PSD[1]>PSDCUTLOW && Ac_PSD[1]<PSDCUTHIGH && abs(Ac_z[1])<ZFIDCUT){
			segNum = Ac_seg[1];
			hRnPoSelectSeg->Fill(segNum);

			selectDt = (Ac_t[1] - Ac_t[0])*CONVERTnsTOms;
			vhSelectDt[segNum]->Fill(selectDt);
			hSelectDt_AllCells->Fill(selectDt);		
			
			vhSelectPromptPSD[segNum]->Fill(Ac_PSD[0]);	
			vhSelectDelayPSD[segNum]->Fill(Ac_PSD[1]);
			hSelectPromptPSD_AllCells->Fill(Ac_PSD[0]);
			hSelectDelayPSD_AllCells->Fill(Ac_PSD[1]);

			vhSelectPromptE[segNum]->Fill(Ac_E[0]);
			vhSelectDelayE[segNum]->Fill(Ac_E[1]);
			hSelectPromptE_AllCells->Fill(Ac_E[0]);
			hSelectDelayE_AllCells->Fill(Ac_E[1]);

			selectDz = Ac_z[1] - Ac_z[0];
			vhSelectDz[segNum]->Fill(selectDz);
			hSelectDz_AllCells->Fill(selectDz);

			vhSelectDelayPos[segNum]->Fill(Ac_z[1]);
			hSelectDelayPos_AllCells->Fill(Ac_z[1]);

			hSelectPSDvsEn->Fill(Ac_E[0],Ac_PSD[0]);
			hSelectPSDvsEn->Fill(Ac_E[1],Ac_PSD[1]);

			hSelectDelayEnVsPromptEn->Fill(Ac_E[0],Ac_E[1]);

			vhSelectPSDvsZ[segNum]->Fill(Ac_z[0],Ac_PSD[0]);
			vhSelectPSDvsZ[segNum]->Fill(Ac_z[1],Ac_PSD[1]);

			vhSelectEvsZ[segNum]->Fill(Ac_z[0],Ac_E[0]);
			vhSelectEvsZ[segNum]->Fill(Ac_z[1],Ac_E[1]);
		}

		//check if there is a background delay event
		if(Ac_evt[2]!=0 && !exclude && Ac_PSD[2]>PSDCUTLOW && Ac_PSD[2]<PSDCUTHIGH && abs(Ac_z[2])<ZFIDCUT){
			segNum = Ac_seg[2];	
			hRnPoBGSeg->Fill(segNum);	

			BGDt = (Ac_t[2] - Ac_t[0])*CONVERTnsTOms;
			vhBGDt[segNum]->Fill(BGDt);
			hBGDt_AllCells->Fill(BGDt);

			vhBGPromptPSD[segNum]->Fill(Ac_PSD[0]);
			vhBGDelayPSD[segNum]->Fill(Ac_PSD[2]);
			hBGPromptPSD_AllCells->Fill(Ac_PSD[0]);
			hBGDelayPSD_AllCells->Fill(Ac_PSD[2]);

			vhBGPromptE[segNum]->Fill(Ac_E[0]);
			vhBGDelayE[segNum]->Fill(Ac_E[2]);
			hBGPromptE_AllCells->Fill(Ac_E[0]);
			hBGDelayE_AllCells->Fill(Ac_E[2]);
		
			BGDz = Ac_z[2] - Ac_z[0];
			vhBGDz[segNum]->Fill(BGDz);
			hBGDz_AllCells->Fill(BGDz);
	
			vhBGDelayPos[segNum]->Fill(Ac_z[2]);
			hBGDelayPos_AllCells->Fill(Ac_z[2]);

			hBGPSDvsEn->Fill(Ac_E[0],Ac_PSD[0]);
			hBGPSDvsEn->Fill(Ac_E[2],Ac_PSD[2]);

			hBGDelayEnVsPromptEn->Fill(Ac_E[0],Ac_E[2]);

			vhBGPSDvsZ[segNum]->Fill(Ac_z[0],Ac_PSD[0]);
			vhBGPSDvsZ[segNum]->Fill(Ac_z[2],Ac_PSD[2]);

			vhBGEvsZ[segNum]->Fill(Ac_z[0],Ac_E[0]);
			vhBGEvsZ[segNum]->Fill(Ac_z[2],Ac_E[2]);
		}

		skipEntry:

		IDX = i;
	}	//end for loop over AcChain

	printf("Histograms filled \n");
	printf("Finished at idx: %i \n",IDX);

//===============================================================
//Setup variables for results for each cell

	double dtBinWidth = (selectDtMax - selectDtMin)/(double)numDtBins; 
	string cfDtExp = Form("[0]*exp(-x/[1])*(%f/[1])",dtBinWidth);

	TF1 *fDtExp 		= nullptr; 
	TF1 *fPromptPSDGaus = nullptr;
	TF1 *fDelayPSDGaus  = nullptr; 
	TF1 *fGaus0         = nullptr;
	TF1 *fGaus1         = nullptr;
	TF1 *fPromptETot    = nullptr;
	TF1 *fDelayEGaus    = nullptr;
	TF1 *fDzGaus        = nullptr;

	TH1F *hRnPoDt 	 = nullptr;
	TH1F *hPromptPSD = nullptr;
	TH1F *hDelayPSD  = nullptr;		
	TH1F *hPromptE   = nullptr;
	TH1F *hDelayE    = nullptr;
	TH1F *hRnPoDz    = nullptr;

	TH1F *hRnPoDt_AllCells    = nullptr;
	TH1F *hPromptPSD_AllCells = nullptr;
	TH1F *hDelayPSD_AllCells  = nullptr;
	TH1F *hPromptE_AllCells   = nullptr;
	TH1F *hDelayE_AllCells    = nullptr;
	TH1F *hRnPoDz_AllCells    = nullptr;

	double par[6];

	double rate,    N,    lifetime,    promptPSDEff,    delayPSDEff,     promptEnEff,    delayEnEff,    posEff,    totEff;
	double rateErr, NErr, lifetimeErr, promptPSDEffErr, delayPSDEffErr,  promptEnEffErr, delayEnEffErr, posEffErr, totEffErr;

	double PoEnMean, 	PoEnSigma,    RnPSDMean,    PoPSDMean,    RnPoDzMean,    RnPoDzSigma,    PoPosMean,    PoPosSigma;
	double PoEnMeanErr, PoEnSigmaErr, RnPSDMeanErr, PoPSDMeanErr, RnPoDzMeanErr, RnPoDzSigmaErr, PoPosMeanErr, PoPosSigmaErr;

	vector<double> vRate,    vN,    vLifetime,    vPromptPSDEff,    vDelayPSDEff,    vPromptEnEff,    vDelayEnEff,    vPosEff,    vTotEff;
	vector<double> vRateErr, vNErr, vLifetimeErr, vPromptPSDEffErr, vDelayPSDEffErr, vPromptEnEffErr, vDelayEnEffErr, vPosEffErr, vTotEffErr; 

	vector<double> vPoEnMean,    vPoEnSigma,    vRnPSDMean,    vPoPSDMean,    vRnPoDzMean,    vRnPoDzSigma,    vPoPosMean,    vPoPosSigma;
	vector<double> vPoEnMeanErr, vPoEnSigmaErr, vRnPSDMeanErr, vPoPSDMeanErr, vRnPoDzMeanErr, vRnPoDzSigmaErr, vPoPosMeanErr, vPoPosSigmaErr;

	double rateAllCells,    NAllCells,    lifetimeAllCells,    promptPSDEffAllCells,    delayPSDEffAllCells,     promptEnEffAllCells,    delayEnEffAllCells,    posEffAllCells,    totEffAllCells;
	double rateErrAllCells, NErrAllCells, lifetimeErrAllCells, promptPSDEffErrAllCells, delayPSDEffErrAllCells,  promptEnEffErrAllCells, delayEnEffErrAllCells, posEffErrAllCells, totEffErrAllCells;

	double PoEnMeanAllCells,    PoEnSigmaAllCells,    RnPSDMeanAllCells,    PoPSDMeanAllCells,    RnPoDzMeanAllCells,    RnPoDzSigmaAllCells,    PoPosMeanAllCells,    PoPosSigmaAllCells;
	double PoEnMeanErrAllCells, PoEnSigmaErrAllCells, RnPSDMeanErrAllCells, PoPSDMeanErrAllCells, RnPoDzMeanErrAllCells, RnPoDzSigmaErrAllCells, PoPosMeanErrAllCells, PoPosSigmaErrAllCells;

//===============================================================

	printf("-------------------------------------------------------------------\n");
	printf("Subtracting histograms \n");

//===============================================================
//Set up tree for rate results

	TFile rateTreeFile("ADAcRate_ZFidCut.root","RECREATE");
	TTree rateTree("AcRate","AcRate");

	int brSegment;
	double brRate, brRateErr;

	TBranch *bSegment = rateTree.Branch("seg",&brSegment,"seg/I");
	TBranch *bRate = rateTree.Branch("rate",&brRate,"rate/D");
	TBranch *bRateErr = rateTree.Branch("rateErr",&brRateErr,"rateErr/D");

	for(int i=0;i<NUMCELLS;i++){		
		rate = 0.0,    N = 0.0,    lifetime = 0.0,    promptPSDEff = 0.0,    delayPSDEff = 0.0,     promptEnEff = 0.0,    delayEnEff = 0.0,    posEff = 0.0,    totEff = 0.0,    PoEnMean = 0.0, PoPSDMean = 0.0;
		rateErr = 0.0, NErr = 0.0, lifetimeErr = 0.0, promptPSDEffErr = 0.0, delayPSDEffErr = 0.0,  promptEnEffErr = 0.0, delayEnEffErr = 0.0, posEffErr = 0.0, totEffErr = 0.0, PoEnMeanErr = 0.0, PoPSDMeanErr = 0.0;

//		if(PLOTFLAG==1){

		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(!exclude){
		//--------------------------------------------------------------------------
		//Subtract dt histograms
		hRnPoDt = (TH1F*)vhSelectDt[i]->Clone();
		hRnPoDt->SetNameTitle(Form("RnPoDt_%i",i),Form("Cell %i: Dt of RnPo Events",i));

		vhBGDt[i]->GetXaxis()->SetLimits(selectDtMin,selectDtMax);	

		hRnPoDt->Add(vhBGDt[i],-1);

		//Set errros on RnPoDt histogram
		int hRnPoDt_size = hRnPoDt->GetSize();
		for(int j=0;j<hRnPoDt_size;j++){
			double selectDtErr = vhSelectDt[i]->GetBinError(j);
			double BGDtErr     = vhBGDt[i]->GetBinError(j);
			double RnPoDtErr   = sqrt(pow(selectDtErr,2) + pow(BGDtErr,2));
			hRnPoDt->SetBinError(j,RnPoDtErr);
		}

		vhRnPoDt.push_back(hRnPoDt);

		fDtExp = new TF1("fDtExp",cfDtExp.c_str(),selectDtMin,selectDtMax);

		fDtExp->SetParName(0,"N");
		fDtExp->SetParName(1,"PoLifetime");
		fDtExp->SetParameter(1,2.57);

		hRnPoDt->Fit(fDtExp,"RQ");

		//--------------------------------------------------------------------------
		//Subtract prompt and delay PSD histograms
		hPromptPSD = (TH1F*)vhSelectPromptPSD[i]->Clone();
		hPromptPSD->SetNameTitle(Form("hPromptPSD_%i",i),Form("Cell %i: PSD of Prompt Events",i));
	
		hPromptPSD->Sumw2();	
		hPromptPSD->Add(vhBGPromptPSD[i],-1);	

		vhPromptPSD.push_back(hPromptPSD);			
	

		hDelayPSD = (TH1F*)vhSelectDelayPSD[i]->Clone();
		hDelayPSD->SetNameTitle(Form("hDelayPSD_%i",i),Form("Cell %i: PSD of Delay Events",i));
		
		hDelayPSD->Sumw2();
		hDelayPSD->Add(vhBGDelayPSD[i],-1);

		vhDelayPSD.push_back(hDelayPSD);

		//Fit prompt and delay PSD histograms
		fPromptPSDGaus = new TF1("fPromptPSDGaus","gaus",PSDMin,PSDMax);
		hPromptPSD->Fit(fPromptPSDGaus,"RQ");

		fDelayPSDGaus = new TF1("fDelayPSDGaus","gaus",PSDMin,PSDMax);
		hDelayPSD->Fit(fDelayPSDGaus,"RQ");

		//--------------------------------------------------------------------------
		//Subtract prompt and delay energy histograms
		hPromptE = (TH1F*)vhSelectPromptE[i]->Clone();
		hPromptE->SetNameTitle(Form("hPromptE_%i",i),Form("Cell %i: Energy of Prompt Events",i));
		
		hPromptE->Sumw2();
		hPromptE->Add(vhBGPromptE[i],-1);

		vhPromptE.push_back(hPromptE);

		
		hDelayE = (TH1F*)vhSelectDelayE[i]->Clone();
		hDelayE->SetNameTitle(Form("hDelayE_%i",i),Form("Cell %i: Energy of Delay Events",i));
	
		hDelayE->Sumw2();
		hDelayE->Add(vhBGDelayE[i],-1);

		vhDelayE.push_back(hDelayE);	

		//Fit prompt and delay energy histograms
		fGaus0 = new TF1("fGaus0","gaus",0.6,0.8);
		fGaus1 = new TF1("fGaus1","gaus",0.8,1.0);
		fGaus1->SetParLimits(1,0.6,0.9);
		fPromptETot = new TF1("fPromptETot","gaus(0) + gaus(3)",0.5,1.1);

		hPromptE->Fit(fGaus0,"RQ0");
		hPromptE->Fit(fGaus1,"RQ0+");

		fGaus0->GetParameters(&par[0]);
		fGaus1->GetParameters(&par[3]);
		fPromptETot->SetParameters(par);
		fPromptETot->SetParLimits(4,0.6,0.9);

		hPromptE->Fit(fPromptETot,"RQ0+");
		fPromptETot->SetRange(EMin,EMax);

		fDelayEGaus = new TF1("fDelayEGaus","gaus",EMin,EMax);
		hDelayE->Fit(fDelayEGaus,"RQ");
	
		//--------------------------------------------------------------------------
		//Subtract prompt-delay and prompt-BG dz histograms
		hRnPoDz = (TH1F*)vhSelectDz[i]->Clone();
		hRnPoDz->SetNameTitle(Form("hRnPoDz_%i",i),Form("Cell %i: Dz of RnPo Events",i));
		
		hRnPoDz->Sumw2();
		hRnPoDz->Add(vhBGDz[i],-1);

		vhRnPoDz.push_back(hRnPoDz);	

		//Fit RnPo dz histogram
		fDzGaus = new TF1("fDzGaus","gaus",dzMin,dzMax);
		hRnPoDz->Fit(fDzGaus,"RQ");

		//--------------------------------------------------------------------------
		//Subtract prompt-delay and prompt-BG position histograms
		TH1F *hPoPos = (TH1F*)vhSelectDelayPos[i]->Clone();
		hPoPos->SetName("hPoPos");
		hPoPos->Sumw2();
		hPoPos->Add(vhBGDelayPos[i],-1);	

		//--------------------------------------------------------------------------
		//Subtract PSD vs Z and E vs Z plots
		TH2F* hRnPoPSDvsZ = (TH2F*)vhSelectPSDvsZ[i]->Clone();
		hRnPoPSDvsZ->SetNameTitle(Form("hRnPoPSDvsZ_%i",i),Form("Cell %i: RnPo PSD vs Z",i));
		hRnPoPSDvsZ->Add(vhBGPSDvsZ[i],-1);
		vhRnPoPSDvsZ.push_back(hRnPoPSDvsZ);

		TH2F* hRnPoEvsZ = (TH2F*)vhSelectEvsZ[i]->Clone();
		hRnPoEvsZ->SetNameTitle(Form("hRnPoEvsZ_%i",i),Form("Cell %i: RnPo Energy vs Z",i));
		hRnPoEvsZ->Add(vhBGEvsZ[i],-1);
		vhRnPoEvsZ.push_back(hRnPoEvsZ);
		//--------------------------------------------------------------------------
		//Calculate efficiencies and Ac227 rate

		//Prompt PSD Eff
		promptPSDEff = fPromptPSDGaus->Integral(promptPSDCutLow,promptPSDCutHigh)/fPromptPSDGaus->Integral(PSDMin,PSDMax);	
		promptPSDEffErr = sqrt((promptPSDEff*(1.0-promptPSDEff))/hPromptPSD->GetEntries());

		//DelayPSDEff
		delayPSDEff = fDelayPSDGaus->Integral(delayPSDCutLow,delayPSDCutHigh)/fDelayPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1.0-delayPSDEff))/hDelayPSD->GetEntries());

		//Prompt Energy Eff
		promptEnEff = fPromptETot->Integral(promptECutLow,promptECutHigh)/fPromptETot->Integral(EMin,EMax);
		promptEnEffErr = sqrt((promptEnEff*(1.0-promptEnEff))/hPromptE->GetEntries());

		//Delay Energy Eff
		delayEnEff = fDelayEGaus->Integral(delayECutLow,delayECutHigh)/fDelayEGaus->Integral(EMin,EMax);
		delayEnEffErr = sqrt((delayEnEff*(1.0-delayEnEff))/hDelayE->GetEntries());

		//Position Eff
		posEff = fDzGaus->Integral(-posCut,posCut)/fDzGaus->Integral(dzMin,dzMax);
		posEffErr = sqrt((posEff*(1.0-posEff))/hRnPoDz->GetEntries());


		totEff = promptPSDEff*delayPSDEff*promptEnEff*delayEnEff*posEff;
		totEffErr = totEff*sqrt(pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(posEffErr/posEff,2)); 

		//==================================================	
		N = fDtExp->GetParameter(0);
		NErr = fDtExp->GetParError(0);
		lifetime = fDtExp->GetParameter(1);
		lifetimeErr = fDtExp->GetParError(1);

		rate = N/(livetime*totEff);
		rateErr = rate*sqrt(pow(NErr/N,2) + pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(posEffErr/posEff,2)); 

		PoEnMean = fDelayEGaus->GetParameter(1);
		PoEnMeanErr = fDelayEGaus->GetParError(1);

		PoEnSigma = fDelayEGaus->GetParameter(2);
		PoEnSigmaErr = fDelayEGaus->GetParError(2);

		RnPSDMean = fPromptPSDGaus->GetParameter(1);
		RnPSDMeanErr = fPromptPSDGaus->GetParError(1);

		PoPSDMean = fDelayPSDGaus->GetParameter(1);
		PoPSDMeanErr = fDelayPSDGaus->GetParError(1);

		RnPoDzMean = fDzGaus->GetParameter(1);
		RnPoDzMeanErr = fDzGaus->GetParError(1);

		RnPoDzSigma = fDzGaus->GetParameter(2);
		RnPoDzSigmaErr = fDzGaus->GetParError(2);
		
		PoPosMean = hPoPos->GetMean();
		PoPosMeanErr = hPoPos->GetMeanError();

		PoPosSigma = hPoPos->GetRMS();
		PoPosSigmaErr = hPoPos->GetRMSError();

		//==================================================	

		//Fill rate results tree	
		brSegment = i;
		brRate = rate;
		brRateErr = rateErr;
		
		rateTree.Fill();
		
		//--------------------------------------------------------------------------
		//Draw histograms
		gStyle->SetOptStat(11);
		gStyle->SetOptFit(1111);

		TLegend *leg = nullptr;

		//if(i==110 || i==109 || i==111 || i==124 || i==96 || i==123 || i==125 || i==95 || i==97){
		if(i<300){
		TCanvas *cRnPoDt = new TCanvas(Form("cRnPoDt_%i",i),"cRnPoDt",1);
		vhSelectDt[i]->SetLineColor(kBlue);
		vhSelectDt[i]->SetTitle(Form("Cell %i: Dt of RnPo Events",i));
		vhSelectDt[i]->Draw();
		vhBGDt[i]->SetLineColor(6);
		vhBGDt[i]->Draw("sames");
		hRnPoDt->SetLineColor(kBlack);
		hRnPoDt->SetMarkerColor(kBlack);
		cRnPoDt->Update();
		hRnPoDt->Draw("sames");
		leg = new TLegend(0.78,0.45,0.98,0.60);
		leg->AddEntry(vhSelectDt[i],"Selection","l");
		leg->AddEntry(vhBGDt[i],"Background","l");
		leg->AddEntry(hRnPoDt,"RnPo","l");
		leg->AddEntry(fDtExp,"Fit","l");
		leg->Draw();

		TCanvas *cRnPoPSD = new TCanvas(Form("cRnPoPSD_%i",i),"cRnPoPSD",1);
		vhSelectDelayPSD[i]->SetTitle(Form("Cell %i: PSD of RnPo Events",i));
		vhSelectDelayPSD[i]->SetLineColor(6);	
		vhSelectDelayPSD[i]->Draw();
		vhBGDelayPSD[i]->SetLineColor(6);
		vhBGDelayPSD[i]->SetLineStyle(6);
		vhBGDelayPSD[i]->Draw("same");
		hDelayPSD->SetMarkerColor(6);
		hDelayPSD->SetLineColor(6);
		hDelayPSD->SetLineWidth(2);
		hDelayPSD->Draw("sames");
		gPad->Update();
		TPaveStats *std = (TPaveStats*)hDelayPSD->FindObject("stats");
		std->SetX1NDC(0.7); 
		std->SetX2NDC(0.98); 
		std->SetY1NDC(0.68); 
		std->SetY2NDC(0.44); 
		vhSelectPromptPSD[i]->SetLineColor(4);
		vhSelectPromptPSD[i]->Draw("same");
		vhBGPromptPSD[i]->SetLineColor(4);
		vhBGPromptPSD[i]->SetLineStyle(2);
		vhBGPromptPSD[i]->Draw("same");
		hPromptPSD->SetMarkerColor(kBlue);
		hPromptPSD->SetLineColor(kBlue);
		hPromptPSD->SetLineWidth(2);
		hPromptPSD->Draw("sames");
		gPad->Update();	
		TPaveStats *stp = (TPaveStats*)hPromptPSD->FindObject("stats");
		stp->SetX1NDC(0.7); 
		stp->SetX2NDC(0.98); 
		stp->SetY1NDC(0.7); 
		stp->SetY2NDC(0.94); 
		leg = new TLegend(0.14,0.67,0.36,0.89);
		leg->AddEntry(vhSelectPromptPSD[i],"Selection Prompt","l");
		leg->AddEntry(vhBGPromptPSD[i],"BG Prompt","l");
		leg->AddEntry(hPromptPSD,"Prompt","l");
		leg->AddEntry(vhSelectDelayPSD[i],"Selection Delay","l");
		leg->AddEntry(vhBGDelayPSD[i],"BG Delay","l");
		leg->AddEntry(hDelayPSD,"Delay","l");
		leg->Draw();

		TCanvas *cRnPoE = new TCanvas(Form("cRnPoE_%i",i),"cRnPoE",1);
		vhSelectDelayE[i]->SetTitle(Form("Cell %i: Energy of RnPo Events",i));
		vhSelectDelayE[i]->SetLineColor(6);
		vhSelectDelayE[i]->Draw();
		vhBGDelayE[i]->SetLineColor(6);
		vhBGDelayE[i]->SetLineStyle(2);
		vhBGDelayE[i]->SetLineWidth(2);
		vhBGDelayE[i]->Draw("same");
		hDelayE->SetLineColor(6);
		hDelayE->SetMarkerColor(6);
		hDelayE->SetLineWidth(2);
		hDelayE->Draw("sames");
		gPad->Update();
		std = (TPaveStats*)hDelayE->FindObject("stats");
		std->SetX1NDC(0.7); 
		std->SetX2NDC(0.98); 
		std->SetY1NDC(0.68); 
		std->SetY2NDC(0.44); 
		vhSelectPromptE[i]->SetLineColor(4);
		vhSelectPromptE[i]->Draw("same");
		vhBGPromptE[i]->SetLineColor(4);
		vhBGPromptE[i]->SetLineStyle(2);
		vhBGPromptE[i]->SetLineWidth(2);
		vhBGPromptE[i]->Draw("same");
		hPromptE->SetLineColor(4);
		hPromptE->SetMarkerColor(4);
		hPromptE->SetLineWidth(2);
		hPromptE->Draw("sames");
		fPromptETot->Draw("sames");
		fGaus0->SetLineColor(8);
		fGaus0->SetRange(EMin,EMax);
		fGaus0->SetLineStyle(2);
		fGaus0->Draw("same");
		fGaus1->SetLineColor(8);
		fGaus1->Draw("same");
		fGaus1->SetRange(EMin,EMax);
		fGaus1->SetLineStyle(2);
		
		gPad->Update();	
		stp = (TPaveStats*)hPromptE->FindObject("stats");
		stp->SetX1NDC(0.7); 
		stp->SetX2NDC(0.98); 
		stp->SetY1NDC(0.7); 
		stp->SetY2NDC(0.94); 
		leg = new TLegend(0.14,0.67,0.36,0.89);
		leg->AddEntry(vhSelectPromptE[i],"Selection Prompt","l");
		leg->AddEntry(vhBGPromptE[i],"BG Prompt","l");
		leg->AddEntry(hPromptE,"Prompt","l");
		leg->AddEntry(vhSelectDelayE[i],"Selection Delay","l");
		leg->AddEntry(vhBGDelayE[i],"BG Delay","l");
		leg->AddEntry(hDelayE,"Delay","l");
		leg->Draw();

		TCanvas *cRnPoDz = new TCanvas(Form("cRnPoDz_%i",i),"cRnPoDz",1);
		vhSelectDz[i]->SetTitle(Form("Cell %i: Dz of RnPo Events",i));
		vhSelectDz[i]->SetLineColor(4);
		vhSelectDz[i]->Draw();
		vhBGDz[i]->SetLineColor(6);
		vhBGDz[i]->Draw("same");
		hRnPoDz->SetLineColor(1);
		hRnPoDz->Draw("same");
		leg = new TLegend(0.78,0.45,0.98,0.60);
		leg->AddEntry(vhSelectDz[i],"Selection","l");
		leg->AddEntry(vhBGDz[i],"Background","l");
		leg->AddEntry(hRnPoDz,"RnPo","l");
		leg->AddEntry(fDzGaus,"Fit","l");
		leg->Draw();


		setup_PROSPECT_style();
	    gROOT->ForceStyle();

		TCanvas *cAllPSDvsZ = new TCanvas(Form("cPSDvsZ_%i",i),"cPSDvsZ",1);
		TH2F* hAllPSDvsZ = (TH2F*)vhSelectPSDvsZ[i]->Clone();
		hAllPSDvsZ->Add(vhBGPSDvsZ[i]);
		hAllPSDvsZ->Draw("colz");

		TCanvas *cAllEvsZ = new TCanvas(Form("cEvsZ_%i",i),"cEvsZ",1);
		TH2F* hAllEvsZ = (TH2F*)vhSelectEvsZ[i]->Clone();
		hAllEvsZ->Add(vhBGEvsZ[i]);
		hAllEvsZ->Draw("colz");

		TCanvas *cRnPoPSDvsZ = new TCanvas(Form("cRnPoPSDvsZ_%i",i),"cRnPoPSDvsZ",1);
		hRnPoPSDvsZ->Draw("colz");		
	
		TCanvas *cRnPoEvsZ = new TCanvas(Form("cRnPoEvsZ_%i",i),"cRnPoEvsZ",1);
		hRnPoEvsZ->Draw("colz");

		cRnPoDt->SaveAs(Form("%s/Cell%i_RnPoDt_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));	
		cRnPoPSD->SaveAs(Form("%s/Cell%i_RnPoPSD_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cRnPoE->SaveAs(Form("%s/Cell%i_RnPoE_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cRnPoDz->SaveAs(Form("%s/Cell%i_RnPoDz_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cAllPSDvsZ->SaveAs(Form("%s/Cell%i_AllPSDvsZ_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cAllEvsZ->SaveAs(Form("%s/Cell%i_AllEvsZ_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cRnPoPSDvsZ->SaveAs(Form("%s/Cell%i_RnPoPSDvsZ_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));
		cRnPoEvsZ->SaveAs(Form("%s/Cell%i_RnPoEvsZ_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),i,nLoop));

		delete cRnPoDt;
		delete cRnPoPSD;
		delete cRnPoE;
		delete cRnPoDz;
		delete cAllPSDvsZ;
		delete cAllEvsZ;
		delete cRnPoPSDvsZ;
		delete cRnPoEvsZ;
		}

		}

//		}	//if(PLOTFLAG==1)
		vRate.push_back(rate);
		vRateErr.push_back(rateErr);
		vN.push_back(N);
		vNErr.push_back(NErr);
		vLifetime.push_back(lifetime);
		vLifetimeErr.push_back(lifetimeErr);
		vPromptPSDEff.push_back(promptPSDEff);
		vPromptPSDEffErr.push_back(promptPSDEffErr);
		vDelayPSDEff.push_back(delayPSDEff);
		vDelayPSDEffErr.push_back(delayPSDEffErr);
		vPromptEnEff.push_back(promptEnEff);
		vPromptEnEffErr.push_back(promptEnEffErr);
		vDelayEnEff.push_back(delayEnEff);
		vDelayEnEffErr.push_back(delayEnEffErr);
		vPosEff.push_back(posEff);
		vPosEffErr.push_back(posEffErr);
		vTotEff.push_back(totEff);
		vTotEffErr.push_back(totEffErr);

		vPoEnMean.push_back(PoEnMean);
		vPoEnMeanErr.push_back(PoEnMeanErr);		
		vPoEnSigma.push_back(PoEnSigma);
		vPoEnSigmaErr.push_back(PoEnSigmaErr);
		vRnPSDMean.push_back(RnPSDMean);
		vRnPSDMeanErr.push_back(RnPSDMeanErr);
		vPoPSDMean.push_back(PoPSDMean);
		vPoPSDMeanErr.push_back(PoPSDMeanErr);		
		vRnPoDzMean.push_back(RnPoDzMean);
		vRnPoDzMeanErr.push_back(RnPoDzMeanErr);
		vRnPoDzSigma.push_back(RnPoDzSigma);
		vRnPoDzSigmaErr.push_back(RnPoDzSigmaErr);
		vPoPosMean.push_back(PoPosMean);
		vPoPosMeanErr.push_back(PoPosMeanErr);
		vPoPosSigma.push_back(PoPosSigma);
		vPoPosSigmaErr.push_back(PoPosSigmaErr);

	}

//===============================================================
//Find rate in total detector

	//--------------------------------------------------------------------------
	//Subtract dt histograms
	hRnPoDt_AllCells = (TH1F*)hSelectDt_AllCells->Clone();
	hRnPoDt_AllCells->SetNameTitle("RnPoDt_AllCells","All Cells: Dt of RnPo Events");

	hBGDt_AllCells->GetXaxis()->SetLimits(selectDtMin,selectDtMax);

	hRnPoDt_AllCells->Add(hBGDt_AllCells,-1);

	int hRnPoDt_AllCells_size = hRnPoDt_AllCells->GetSize();
	for(int j=0;j<hRnPoDt_AllCells_size;j++){
		double selectDtErr = hSelectDt_AllCells->GetBinError(j);
		double BGDtErr     = hBGDt_AllCells->GetBinError(j);
		double RnPoDtErr   = sqrt(pow(selectDtErr,2) + pow(BGDtErr,2));
		hRnPoDt_AllCells->SetBinError(j,RnPoDtErr);
	}

	fDtExp = new TF1("fDtExp",cfDtExp.c_str(),selectDtMin,selectDtMax);
	fDtExp->SetParName(0,"N");
	fDtExp->SetParName(1,"PoLifetime");
	fDtExp->SetParameter(1,2.57);

	hRnPoDt_AllCells->Fit(fDtExp,"RQ");

	//--------------------------------------------------------------------------
	//Subtract prompt and delay PSD histograms
	hPromptPSD_AllCells = (TH1F*)hSelectPromptPSD_AllCells->Clone();
	hPromptPSD_AllCells->SetNameTitle("hPromptPSD_AllCells","All Cells: PSD of Prompt Events");
	
	hPromptPSD_AllCells->Sumw2();	
	hPromptPSD_AllCells->Add(hBGPromptPSD_AllCells,-1);	


	hDelayPSD_AllCells = (TH1F*)hSelectDelayPSD_AllCells->Clone();
	hDelayPSD_AllCells->SetNameTitle("hDelayPSD_AllCells","All Cells: PSD of Delay Events");
		
	hDelayPSD_AllCells->Sumw2();
	hDelayPSD_AllCells->Add(hBGDelayPSD_AllCells,-1);

	//Fit prompt and delay PSD histograms
	fPromptPSDGaus = new TF1("fPromptPSDGaus","gaus",PSDMin,PSDMax);
	hPromptPSD_AllCells->Fit(fPromptPSDGaus,"RQ");

	fDelayPSDGaus = new TF1("fDelayPSDGaus","gaus",PSDMin,PSDMax);
	hDelayPSD_AllCells->Fit(fDelayPSDGaus,"RQ");
	
	//--------------------------------------------------------------------------
	//Subtract prompt and delay energy histograms
	hPromptE_AllCells = (TH1F*)hSelectPromptE_AllCells->Clone();
	hPromptE_AllCells->SetNameTitle("hPromptE_AllCells","All Cells: Energy of Prompt Events");
		
	hPromptE_AllCells->Sumw2();
	hPromptE_AllCells->Add(hBGPromptE_AllCells,-1);

		
	hDelayE_AllCells = (TH1F*)hSelectDelayE_AllCells->Clone();
	hDelayE_AllCells->SetNameTitle("hDelayE_AllCells","All Cells: Energy of Delay Events");
	
	hDelayE_AllCells->Sumw2();
	hDelayE_AllCells->Add(hBGDelayE_AllCells,-1);

	double promptEFitMin = promptECutLow+0.1, promptEFitMax = promptECutHigh-0.1;
	double delayEFitMin = delayECutLow+0.1, delayEFitMax = delayECutHigh-0.1;

	delayEFitMin = 0.6; delayEFitMax = 1.1;

	//Fit prompt and delay energy histograms
	fGaus0 = new TF1("fGaus0","gaus",0.55,0.8);
	fGaus1 = new TF1("fGaus1","gaus",0.82,1.1);
	fPromptETot = new TF1("fPromptETot","gaus(0) + gaus(3)",0.55,1.1);
	
	hPromptE_AllCells->Fit(fGaus0,"RQ0");
	hPromptE_AllCells->Fit(fGaus1,"RQ0+");
	
	fGaus0->GetParameters(&par[0]);
	fGaus1->GetParameters(&par[3]);
	fPromptETot->SetParameters(par);

	hPromptE_AllCells->Fit(fPromptETot,"RQ0+");
	fPromptETot->SetRange(EMin,EMax);

	fDelayEGaus = new TF1("fDelayEGaus","gaus",delayEFitMin,delayEFitMax);
	hDelayE_AllCells->Fit(fDelayEGaus,"RQ");
	fDelayEGaus->SetRange(EMin,EMax);

	//--------------------------------------------------------------------------
	//Subtract prompt-delay and prompt-BG dz histograms
	hRnPoDz_AllCells = (TH1F*)hSelectDz_AllCells->Clone();
	hRnPoDz_AllCells->SetNameTitle("hRnPoDz_AllCells","All Cells: Dz of RnPo Events");
		
	hRnPoDz_AllCells->Sumw2();
	hRnPoDz_AllCells->Add(hBGDz_AllCells,-1);

	//Fit RnPo dz histogram
	fDzGaus = new TF1("fDzGaus","gaus",dzMin,dzMax);
	hRnPoDz_AllCells->Fit(fDzGaus,"RQ");

	//--------------------------------------------------------------------------
	//Subtract prompt-delay and prompt-BG position histograms
	TH1F *hPoPos_AllCells = (TH1F*)hSelectDelayPos_AllCells->Clone();
	hPoPos_AllCells->SetName("hPoPos_AllCells");
	hPoPos_AllCells->Sumw2();
	hPoPos_AllCells->Add(hBGDelayPos_AllCells,-1);	

	//--------------------------------------------------------------------------
	//Calculate efficiencies and Ac227 rate

	//Prompt PSD Eff
	promptPSDEffAllCells = fPromptPSDGaus->Integral(promptPSDCutLow,promptPSDCutHigh)/fPromptPSDGaus->Integral(PSDMin,PSDMax);	
	promptPSDEffErrAllCells = sqrt((promptPSDEffAllCells*(1.0-promptPSDEffAllCells))/hPromptPSD_AllCells->GetEntries());

	//DelayPSDEff
	delayPSDEffAllCells = fDelayPSDGaus->Integral(delayPSDCutLow,delayPSDCutHigh)/fDelayPSDGaus->Integral(PSDMin,PSDMax);
	delayPSDEffErrAllCells = sqrt((delayPSDEffAllCells*(1.0-delayPSDEffAllCells))/hDelayPSD_AllCells->GetEntries());

	//Prompt Energy Eff
	promptEnEffAllCells = fPromptETot->Integral(promptECutLow,promptECutHigh)/fPromptETot->Integral(EMin,EMax);
	promptEnEffErrAllCells = sqrt((promptEnEffAllCells*(1.0-promptEnEffAllCells))/hPromptE_AllCells->GetEntries());

cout<<"Energy eff: "<<promptEnEffAllCells<<" +/- "<<promptEnEffErrAllCells<<endl;

	//Delay Energy Eff
	delayEnEffAllCells = fDelayEGaus->Integral(delayECutLow,delayECutHigh)/fDelayEGaus->Integral(EMin,EMax);
	delayEnEffErrAllCells = sqrt((delayEnEffAllCells*(1.0-delayEnEffAllCells))/hDelayE_AllCells->GetEntries());

	//Position Eff
	posEffAllCells = fDzGaus->Integral(-posCut,posCut)/fDzGaus->Integral(dzMin,dzMax);
	posEffErrAllCells = sqrt((posEffAllCells*(1.0-posEffAllCells))/hRnPoDz_AllCells->GetEntries());

	totEffAllCells = promptPSDEffAllCells*delayPSDEffAllCells*promptEnEffAllCells*delayEnEffAllCells*posEffAllCells;
	totEffErrAllCells = totEff*sqrt(pow(promptPSDEffErrAllCells/promptPSDEffAllCells,2) + pow(delayPSDEffErrAllCells/delayPSDEffAllCells,2) + pow(promptEnEffErrAllCells/promptEnEffAllCells,2) + pow(delayEnEffErrAllCells/delayEnEffAllCells,2) + pow(posEffErrAllCells/posEffAllCells,2)); 

	//==================================================	
	NAllCells = fDtExp->GetParameter(0);
	NErrAllCells = fDtExp->GetParError(0);
	lifetimeAllCells = fDtExp->GetParameter(1);
	lifetimeErrAllCells = fDtExp->GetParError(1);

	rateAllCells = NAllCells/(livetime*totEffAllCells);
	rateErrAllCells = rateAllCells*sqrt(pow(NErrAllCells/NAllCells,2) + pow(promptPSDEffErrAllCells/promptPSDEffAllCells,2) + pow(delayPSDEffErrAllCells/delayPSDEffAllCells,2) + pow(promptEnEffErrAllCells/promptEnEffAllCells,2) + pow(delayEnEffErrAllCells/delayEnEffAllCells,2) + pow(posEffErrAllCells/posEffAllCells,2)); 

	PoEnMeanAllCells = fDelayEGaus->GetParameter(1);
	PoEnMeanErrAllCells = fDelayEGaus->GetParError(1);

	PoEnSigmaAllCells = fDelayEGaus->GetParameter(2);
	PoEnSigmaErrAllCells = fDelayEGaus->GetParError(2);

	RnPSDMeanAllCells = fPromptPSDGaus->GetParameter(2);
	RnPSDMeanErrAllCells = fPromptPSDGaus->GetParError(2);	

	PoPSDMeanAllCells = fDelayPSDGaus->GetParameter(1);
	PoPSDMeanErrAllCells = fDelayPSDGaus->GetParError(1);

	RnPoDzMeanAllCells = fDzGaus->GetParameter(1);
	RnPoDzMeanErrAllCells = fDzGaus->GetParError(1);

	RnPoDzSigmaAllCells = fDzGaus->GetParameter(2);
	RnPoDzSigmaErrAllCells = fDzGaus->GetParError(2);

	PoPosMeanAllCells = hPoPos_AllCells->GetMean();
	PoPosMeanErrAllCells = hPoPos_AllCells->GetMeanError();

	PoPosSigmaAllCells = hPoPos_AllCells->GetRMS();
	PoPosSigmaErrAllCells = hPoPos_AllCells->GetRMSError();

	//==================================================	

	//Fill rate results tree	
	brSegment = 155;
	brRate = rateAllCells;
	brRateErr = rateErrAllCells;

	rateTree.Fill();

	printf("==================================================================================\n");
	printf("RATE TOTAL DETECTOR: %f +/- %f Hz\n",rateAllCells*1000,rateErrAllCells*1000);
	printf("TOTAL LIVETIME: %f hrs\n",livetime/(3.6e6));
	printf("==================================================================================\n");


	vRate.push_back(rateAllCells);
	vRateErr.push_back(rateErrAllCells);
	vN.push_back(NAllCells);
	vNErr.push_back(NErrAllCells);
	vLifetime.push_back(lifetimeAllCells);
	vLifetimeErr.push_back(lifetimeErrAllCells);
	vPromptPSDEff.push_back(promptPSDEffAllCells);
	vPromptPSDEffErr.push_back(promptPSDEffErrAllCells);
	vDelayPSDEff.push_back(delayPSDEffAllCells);
	vDelayPSDEffErr.push_back(delayPSDEffErrAllCells);
	vPromptEnEff.push_back(promptEnEffAllCells);
	vPromptEnEffErr.push_back(promptEnEffErrAllCells);
	vDelayEnEff.push_back(delayEnEffAllCells);
	vDelayEnEffErr.push_back(delayEnEffErrAllCells);
	vPosEff.push_back(posEffAllCells);
	vPosEffErr.push_back(posEffErrAllCells);
	vTotEff.push_back(totEffAllCells);
	vTotEffErr.push_back(totEffErrAllCells);

	vPoEnMean.push_back(PoEnMeanAllCells);
	vPoEnMeanErr.push_back(PoEnMeanErrAllCells);		
	vPoEnSigma.push_back(PoEnSigmaAllCells);
	vPoEnSigmaErr.push_back(PoEnSigmaErrAllCells);
	vRnPSDMean.push_back(RnPSDMeanAllCells);
	vRnPSDMeanErr.push_back(RnPSDMeanErrAllCells);
	vPoPSDMean.push_back(PoPSDMeanAllCells);
	vPoPSDMeanErr.push_back(PoPSDMeanErrAllCells);		
	vRnPoDzMean.push_back(RnPoDzMeanAllCells);
	vRnPoDzMeanErr.push_back(RnPoDzMeanErrAllCells);
	vRnPoDzSigma.push_back(RnPoDzSigmaAllCells);
	vRnPoDzSigmaErr.push_back(RnPoDzSigmaErrAllCells);
	vPoPosMean.push_back(PoPosMeanAllCells);
	vPoPosMeanErr.push_back(PoPosMeanErrAllCells);
	vPoPosSigma.push_back(PoPosSigmaAllCells);
	vPoPosSigmaErr.push_back(PoPosSigmaErrAllCells);

	if(PLOTFLAG==1){

	//--------------------------------------------------------------------------
	//Plot histograms

	gStyle->SetOptStat(11);
	gStyle->SetOptFit(1111);

	TH2F* hRnPoPSDvsEn = (TH2F*)hSelectPSDvsEn->Clone();
	hRnPoPSDvsEn->SetName("hRnPoPSDvsEn");
	hRnPoPSDvsEn->Add(hBGPSDvsEn,-1);

	TCanvas *cSelectPSDvsEn = new TCanvas("cSelectPSDvsEn","cSelectPSDvsEn",1);
	hSelectPSDvsEn->SetTitle("Selection Events PSD vs. Energy;Energy [MeVee];PSD [arb]");
	hSelectPSDvsEn->Draw("colz");

	TCanvas *cBGPSDvsEn = new TCanvas("cBGPSDvsEn","cBGPSDvsEn",1);
	hBGPSDvsEn->SetTitle("BG Events PSD vs. Energy;Energy [MeVee];PSD [arb]");
	hBGPSDvsEn->Draw("colz");

	TCanvas *cRnPoPSDvsEn = new TCanvas("cRnPoPSDvsEn","cRnPoPSDvsEn",1);
	hRnPoPSDvsEn->SetTitle("RnPo PSD vs. Energy;Energy [MeVee];PSD [arb]");
	hRnPoPSDvsEn->Draw("colz");	

	TH2F* hPoEnVsRnEn = (TH2F*)hSelectDelayEnVsPromptEn->Clone();
	hPoEnVsRnEn->SetName("hPoEnVsRnEn");
	hPoEnVsRnEn->Add(hBGDelayEnVsPromptEn,-1);

	TCanvas *cPoEnVsRnEn = new TCanvas("cPoEnVsRnEn","cPoEnVsRnEn",1200,1200);
	cPoEnVsRnEn->SetRightMargin(0.13);
	cPoEnVsRnEn->SetLeftMargin(0.13);
	cPoEnVsRnEn->SetTopMargin(0.13);
	cPoEnVsRnEn->SetBottomMargin(0.13);
	hPoEnVsRnEn->SetTitle("Delay vs. Prompt Energy;Prompt Energy [MeVee];Delay Energy [MeVee]");
	hPoEnVsRnEn->Draw("colz");

	TCanvas *cRnPoDtAll = new TCanvas("cRnPoDtAll","cRnPoDtAll",1);
	gPad->SetGrid();
	cRnPoDtAll->SetLogy();
	hSelectDt_AllCells->SetMinimum(1);
	hSelectDt_AllCells->SetMaximum(1e5);
	hSelectDt_AllCells->SetLineColor(kBlue);
	hSelectDt_AllCells->SetTitle("All Cells: Dt of RnPo Events");
	hSelectDt_AllCells->GetXaxis()->SetRangeUser(0,selectDtMax);
	hSelectDt_AllCells->Draw();
	hBGDt_AllCells->SetLineColor(6);
	hBGDt_AllCells->Draw("sames");
	hRnPoDt_AllCells->SetLineColor(kBlack);
	hRnPoDt_AllCells->SetMarkerColor(kBlack);
	hRnPoDt_AllCells->Draw("sames");
	TLegend	*leg = new TLegend(0.78,0.20,0.98,0.35);
	leg->AddEntry(hSelectDt_AllCells,"Selection","l");
	leg->AddEntry(hBGDt_AllCells,"Background","l");
	leg->AddEntry(hRnPoDt_AllCells,"RnPo","l");
	leg->AddEntry(fDtExp,"Fit","l");
	leg->Draw();

	TCanvas *cRnPoPSDAll = new TCanvas("cRnPoPSDAll","cRnPoPSDAll",1);
	hSelectDelayPSD_AllCells->SetTitle("All Cells: PSD of RnPo Events");
	hSelectDelayPSD_AllCells->SetLineColor(6);	
	hSelectDelayPSD_AllCells->Draw();
	hBGDelayPSD_AllCells->SetLineColor(6);
	hBGDelayPSD_AllCells->SetLineStyle(6);
	hBGDelayPSD_AllCells->Draw("same");
	hDelayPSD_AllCells->SetMarkerColor(6);
	hDelayPSD_AllCells->SetLineColor(6);
	hDelayPSD_AllCells->SetLineWidth(2);
	hDelayPSD_AllCells->Draw("sames");
	gPad->Update();
	TPaveStats *std = (TPaveStats*)hDelayPSD_AllCells->FindObject("stats");
	std->SetX1NDC(0.7); 
	std->SetX2NDC(0.98); 
	std->SetY1NDC(0.68); 
	std->SetY2NDC(0.44); 
	hSelectPromptPSD_AllCells->SetLineColor(4);
	hSelectPromptPSD_AllCells->Draw("same");
	hBGPromptPSD_AllCells->SetLineColor(4);
	hBGPromptPSD_AllCells->SetLineStyle(2);
	hBGPromptPSD_AllCells->Draw("same");
	hPromptPSD_AllCells->SetMarkerColor(kBlue);
	hPromptPSD_AllCells->SetLineColor(kBlue);
	hPromptPSD_AllCells->SetLineWidth(2);
	hPromptPSD_AllCells->Draw("sames");
	gPad->Update();	
	TPaveStats *stp = (TPaveStats*)hPromptPSD_AllCells->FindObject("stats");
	stp->SetX1NDC(0.7); 
	stp->SetX2NDC(0.98); 
	stp->SetY1NDC(0.7); 
	stp->SetY2NDC(0.94); 
	leg = new TLegend(0.7,0.2,0.95,0.4);
	leg->AddEntry(hSelectPromptPSD_AllCells,"Selection Prompt","l");
	leg->AddEntry(hBGPromptPSD_AllCells,"BG Prompt","l");
	leg->AddEntry(hPromptPSD_AllCells,"Prompt","l");
	leg->AddEntry(hSelectDelayPSD_AllCells,"Selection Delay","l");
	leg->AddEntry(hBGDelayPSD_AllCells,"BG Delay","l");
	leg->AddEntry(hDelayPSD_AllCells,"Delay","l");
	leg->Draw();

	TCanvas *cRnPoEAll = new TCanvas("cRnPoEAll","cRnPoEAll",1);
	hSelectDelayE_AllCells->SetTitle("All Cells: Energy of RnPo Events");
	hSelectDelayE_AllCells->SetLineColor(6);
	hSelectDelayE_AllCells->Draw();
	hBGDelayE_AllCells->SetLineColor(6);
	hBGDelayE_AllCells->SetLineStyle(2);
	hBGDelayE_AllCells->SetLineWidth(2);
	hBGDelayE_AllCells->Draw("same");
	hDelayE_AllCells->SetLineColor(6);
	hDelayE_AllCells->SetMarkerColor(6);
	hDelayE_AllCells->SetLineWidth(2);
	hDelayE_AllCells->Draw("sames");
	fDelayEGaus->Draw("sames");
	gPad->Update();
	std = (TPaveStats*)hDelayE_AllCells->FindObject("stats");
	std->SetX1NDC(0.7); 
	std->SetX2NDC(0.98); 
	std->SetY1NDC(0.68); 
	std->SetY2NDC(0.44); 
	hSelectPromptE_AllCells->SetLineColor(4);
	hSelectPromptE_AllCells->Draw("same");
	hBGPromptE_AllCells->SetLineColor(4);
	hBGPromptE_AllCells->SetLineStyle(2);
	hBGPromptE_AllCells->SetLineWidth(2);
	hBGPromptE_AllCells->Draw("same");
	hPromptE_AllCells->SetLineColor(4);
	hPromptE_AllCells->SetMarkerColor(4);
	hPromptE_AllCells->SetLineWidth(2);
	hPromptE_AllCells->Draw("sames");
	fPromptETot->Draw("sames");

	fGaus0->SetLineColor(8);	
	fGaus0->Draw("sames");
	fGaus1->SetLineColor(kBlack);
	fGaus1->Draw("sames");

	gPad->Update();	
	stp = (TPaveStats*)hPromptE_AllCells->FindObject("stats");
	stp->SetX1NDC(0.7); 
	stp->SetX2NDC(0.98); 
	stp->SetY1NDC(0.7); 
	stp->SetY2NDC(0.94); 
	leg = new TLegend(0.14,0.67,0.36,0.89);
	leg->AddEntry(hSelectPromptE_AllCells,"Selection Prompt","l");
	leg->AddEntry(hBGPromptE_AllCells,"BG Prompt","l");
	leg->AddEntry(hPromptE_AllCells,"Prompt","l");
	leg->AddEntry(hSelectDelayE_AllCells,"Selection Delay","l");
	leg->AddEntry(hBGDelayE_AllCells,"BG Delay","l");
	leg->AddEntry(hDelayE_AllCells,"Delay","l");
	leg->Draw();

	TCanvas *cRnPoDzAll = new TCanvas("cRnPoDzAll","cRnPoDzAll",1);
	hSelectDz_AllCells->SetTitle("All Cells: Dz of RnPo Events");
	hSelectDz_AllCells->SetLineColor(4);
	hSelectDz_AllCells->Draw();
	hBGDz_AllCells->SetLineColor(6);
	hBGDz_AllCells->Draw("same");
	hRnPoDz_AllCells->SetLineColor(1);
	hRnPoDz_AllCells->Draw("same");
	leg = new TLegend(0.78,0.45,0.98,0.60);
	leg->AddEntry(hSelectDz_AllCells,"Selection","l");
	leg->AddEntry(hBGDz_AllCells,"Background","l");
	leg->AddEntry(hRnPoDz_AllCells,"RnPo","l");
	leg->AddEntry(fDzGaus,"Fit","l");
	leg->Draw();

	TCanvas *cRnPoPosAll = new TCanvas("cRnPoPosAll","cRnPoPosAll",1);
	hSelectDelayPos_AllCells->SetTitle("All Cells: Position of Po Events");
	hSelectDelayPos_AllCells->SetLineColor(4);
	hSelectDelayPos_AllCells->Draw();
	hBGDelayPos_AllCells->SetLineColor(6);
	hBGDelayPos_AllCells->Draw("same");
	hPoPos_AllCells->SetLineColor(1);
	hPoPos_AllCells->Draw("same");
	leg = new TLegend(0.78,0.45,0.98,0.60);
	leg->AddEntry(hSelectDelayPos_AllCells,"Selection","l");
	leg->AddEntry(hBGDelayPos_AllCells,"Background","l");
	leg->AddEntry(hPoPos_AllCells,"RnPo","l");
	leg->Draw();

	TH1F* hRnPoSeg = (TH1F*)hRnPoSelectSeg->Clone();
	hRnPoSeg->SetName("hRnPoSeg");
	hRnPoSeg->Sumw2();
	hRnPoSeg->Add(hRnPoBGSeg,-1);

	TCanvas *cRnPoSeg = new TCanvas("cRnPoSeg","cRnPoSeg",1);
/*
	hRnPoSelectSeg->SetLineColor(4);
	hRnPoSelectSeg->SetMaximum(9000);
	hRnPoSelectSeg->Draw();
	hRnPoBGSeg->SetLineColor(6);
	hRnPoBGSeg->Draw("sames");
*/
	hRnPoSeg->SetLineColor(1);
	//hRnPoSeg->Draw("sames");
//	hRnPoSeg->SetMaximum(9000);
	hRnPoSeg->SetLineColor(kBlue);
	hRnPoSeg->Draw("HIST");
	hRnPoSeg->Fit("pol0");
	hRnPoSeg->GetListOfFunctions()->FindObject("pol0")->Draw("same");
/*
	leg = new TLegend(0.78,0.40,0.98,0.55);
	leg->AddEntry(hRnPoSelectSeg,"Selection","l");
	leg->AddEntry(hRnPoBGSeg,"Background","l");
	leg->AddEntry(hRnPoSeg,"RnPo","l");
	leg->Draw();
*/
	cSelectPSDvsEn->SaveAs(Form("%s/SelectPSDvsEn_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cBGPSDvsEn->SaveAs(Form("%s/BGPSDvsEn_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cRnPoPSDvsEn->SaveAs(Form("%s/RnPoPSDvsEn_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cPoEnVsRnEn->SaveAs(Form("%s/PoEnVsRnEn_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cRnPoDtAll->SaveAs(Form("%s/AllCells_RnPoDt_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cRnPoPSDAll->SaveAs(Form("%s/AllCells_RnPoPSD_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoEAll->SaveAs(Form("%s/AllCells_RnPoE_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoDzAll->SaveAs(Form("%s/AllCells_RnPoDz_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoPosAll->SaveAs(Form("%s/AllCells_RnPoPos_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoSeg->SaveAs(Form("%s/RnPoSeg_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));

	TFile *saveHistsFile = new TFile("AD_PhysNeutrino_Hists_noZFid_TotDetector.root","RECREATE");
	hRnPoPSDvsEn->Write();
	hPoEnVsRnEn->Write();
	hSelectDt_AllCells->Write();
	hBGDt_AllCells->Write();
	hRnPoDt_AllCells->Write();
	hPromptPSD_AllCells->Write();
	hDelayPSD_AllCells->Write();		
	hPromptE_AllCells->Write();
	hDelayE_AllCells->Write();
	hRnPoDz_AllCells->Write();
	hPoPos_AllCells->Write();
	hRnPoSeg->Write();
	saveHistsFile->Close();	

	delete cRnPoPSDvsEn;
	delete cPoEnVsRnEn;
	delete cRnPoDtAll;
	delete cRnPoPSDAll;
	delete cRnPoEAll;
	delete cRnPoDzAll;
	delete cRnPoPosAll;
	delete cRnPoSeg;
	}


	//--------------------------------------------------------------------------
	printf("Histograms subtracted \n");
	printf("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+ \n");
	printf("================================================================================================= \n");
	outFile->Write();
	outFile->Close();

	vector<vector<double>> vvRate, vvN, vvLifetime, vvPSDEff, vvEnEff, vvPosEff, vvTotEff;
	vector<vector<double>> vvPoEnMean, vvPoEnSigma, vvRnPSDMean, vvPoPSDMean, vvRnPoDzMean, vvRnPoDzSigma, vvPoPosMean, vvPoPosSigma;

	vvRate.push_back(vRate);
	vvRate.push_back(vRateErr);
	vvN.push_back(vN);
	vvN.push_back(vNErr);
	vvLifetime.push_back(vLifetime);
	vvLifetime.push_back(vLifetimeErr);
	vvPSDEff.push_back(vPromptPSDEff);
	vvPSDEff.push_back(vPromptPSDEffErr);
	vvPSDEff.push_back(vDelayPSDEff);
	vvPSDEff.push_back(vDelayPSDEffErr);
	vvEnEff.push_back(vPromptEnEff);
	vvEnEff.push_back(vPromptEnEffErr);
	vvEnEff.push_back(vDelayEnEff);
	vvEnEff.push_back(vDelayEnEffErr);
	vvPosEff.push_back(vPosEff);
	vvPosEff.push_back(vPosEffErr);
	vvTotEff.push_back(vTotEff);
	vvTotEff.push_back(vTotEffErr);

	vvPoEnMean.push_back(vPoEnMean);
	vvPoEnMean.push_back(vPoEnMeanErr);
	vvPoEnSigma.push_back(vPoEnSigma);
	vvPoEnSigma.push_back(vPoEnSigmaErr);
	vvRnPSDMean.push_back(vRnPSDMean);
	vvRnPSDMean.push_back(vRnPSDMeanErr);
	vvPoPSDMean.push_back(vPoPSDMean);
	vvPoPSDMean.push_back(vPoPSDMeanErr);
	vvRnPoDzMean.push_back(vRnPoDzMean);
	vvRnPoDzMean.push_back(vRnPoDzMeanErr);
	vvRnPoDzSigma.push_back(vRnPoDzSigma);
	vvRnPoDzSigma.push_back(vRnPoDzSigmaErr);
	vvPoPosMean.push_back(vPoPosMean);	
	vvPoPosMean.push_back(vPoPosMeanErr);	
	vvPoPosSigma.push_back(vPoPosSigma);	
	vvPoPosSigma.push_back(vPoPosSigmaErr);	

	if(PLOTFLAG==1){
	rateTreeFile.cd();
	//rateTree.Write("",TObject::kOverwrite);
	rateTree.Write();
	rateTreeFile.Close();
	}

	return make_tuple(IDX, tstamp, vvRate, vvN, vvLifetime, vvPSDEff, vvEnEff, vvPosEff, vvTotEff, vvPoEnMean, vvPoEnSigma, vvRnPSDMean, vvPoPSDMean, vvRnPoDzMean, vvRnPoDzSigma, vvPoPosMean, vvPoPosSigma);
}

void PlotResults(const int NUMCELLS, int NUMTREES, int PLOTFLAG, double PSDCUTLOW, double PSDCUTHIGH){

	const int NUMEXCLUDECELLS = 31;
	int ExcludeCellArr[NUMEXCLUDECELLS] = {0,1,2,3,4,5,6,9,10,11,12,13,18,21,23,24,27,32,34,40,44,52,68,79,86,102,115,122,127,130,139};
	bool exclude;

	const int NUMETCELLS = 46, NUMEXCLUDEETCELLS = 13 ;
	int ETCellArr[NUMETCELLS] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,27,28,41,42,55,56,69,70,83,84,97,98,111,112,125,126,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153};
	bool ET;

	int IDX = 0;

	//--------------------------------------------------------------------------
	//Figure number of bins for histograms based on time of all trees and chosen time segment
	TChain AcChain("TAc");

	AcChain.Add(Form("%s/WetCommissioning/Series015_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series016_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series017_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series018_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series019_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series020_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series021_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series022_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set0.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set1.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees_Set2.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180316_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series001_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180417_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180420_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));

	AcChain.Add(Form("%s/180501_ReactorOn/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series001_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series002_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180501_ReactorOn/Series004_AcTrees.root",getenv("P2X_ANALYZED")));

	Double_t Ac_tstamp;
	Double_t Ac_t[3];
	AcChain.SetBranchAddress("Ac_tstamp", &Ac_tstamp);
	AcChain.SetBranchAddress("Ac_t", Ac_t);	//[s] in epoch time

	Int_t nAcEvents = AcChain.GetEntries();

	double prevTime = 0.0;

	double firstTime = 0.0;
	double lastTime  = 0.0;

	int countNumTrees = 0;
	
	for(int i=0;i<nAcEvents;i++){
		AcChain.GetEvent(i);

		if(i%1000000==0) printf("Event: %d \n",i);
		double t = Ac_t[0]*(1e-6);	//ms
		lastTime = Ac_tstamp;

		if(t<prevTime){
			countNumTrees = countNumTrees + 1;
		}
	
		prevTime = t;
	}

	countNumTrees = countNumTrees + 1;
	int numBins = ceil((double)countNumTrees/(double)NUMTREES);	

	if(PLOTFLAG==2)	numBins = 57;
	printf("Number of trees: %d Number of bins: %d \n",countNumTrees,numBins);

	//--------------------------------------------------------------------------
	//Set up histograms	
	int hNumBins = 70;
	double min = 0.0, max = NUMCELLS;

	TH2F* hRatePerCell = new TH2F("hRatePerCell","Rate Per Cell;;;Rate [mHz]",14,0,14,11,0,11);
	TH2F* hRelRatePerCell = new TH2F("hRelRatePerCell","Relative Rate Per Cell;;;Relative Rate",14,0,14,11,0,11);

	//--------------------------------------------------------------------------
	//Set up graphs
	double x[numBins], y[numBins];
	double xErr[numBins], yErr[numBins];

	TGraphErrors *grRate_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grN_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grLifetime_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);

	TGraphErrors *grPromptEnEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grDelayEnEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPromptPSDEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grDelayPSDEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPosEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);

	TGraphErrors *grTotEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoEnMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoEnSigma_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRnPoDzSigma_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoPosMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoPosSigma_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);

	int NumActiveCells = NUMCELLS - NUMEXCLUDECELLS;
	double xPerCell[NumActiveCells], xErrPerCell[NumActiveCells];

	vector<TGraphErrors*> vgrRatePerCell, vgrRelRatePerCell, vgrPromptPSDEffPerCell, vgrDelayPSDEffPerCell, vgrPromptEnEffPerCell, vgrDelayEnEffPerCell, vgrPosEffPerCell, vgrEffPerCell;
	vector<TGraphErrors*> vgrPoEnMeanPerCell, vgrPoEnSigmaPerCell, vgrRelPoEnMeanPerCell, vgrRelPoEnSigmaPerCell, vgrRnPSDMeanPerCell, vgrPoPSDMeanPerCell, vgrRnPoDzMeanPerCell, vgrRnPoDzSigmaPerCell, vgrPoPosMeanPerCell, vgrPoPosSigmaPerCell, vgrRelRnPoDzMeanPerCell;

	int NumActiveCellsET = NUMETCELLS - NUMEXCLUDEETCELLS;
	double xPerCellET[NumActiveCellsET], xErrPerCellET[NumActiveCellsET];

	vector<TGraphErrors*> vgrRatePerCellET, vgrRelRatePerCellET, vgrPromptPSDEffPerCellET, vgrDelayPSDEffPerCellET, vgrPromptEnEffPerCellET, vgrDelayEnEffPerCellET, vgrPosEffPerCellET, vgrEffPerCellET;
	vector<TGraphErrors*> vgrPoEnMeanPerCellET, vgrPoEnSigmaPerCellET, vgrRelPoEnMeanPerCellET, vgrRelPoEnSigmaPerCellET, vgrRnPSDMeanPerCellET, vgrPoPSDMeanPerCellET, vgrRnPoDzMeanPerCellET, vgrRnPoDzSigmaPerCellET, vgrPoPosMeanPerCellET, vgrPoPosSigmaPerCellET;

	TH1F* hgrRate_AllCells = new TH1F("hgrRate_AllCells","Rate for all cells;Time;Rate [Hz]",numBins,0,numBins);

	TH1F* hPoEnMean_PerCell = new TH1F("hPoEnMean_PerCell","Po Energy Mean Per Cell;Cell;Energy [keV]",154,0,154);
	TH1F* hPoEnSigmaPerCell = new TH1F("hPoEnSigmaPerCell","Po Energy Sigma Per Cell;Cell;Energy [keV]",154,0,154);
	TH1F* hPoPosMeanPerCell = new TH1F("hPoPosMeanPerCell","Po Position Mean Per Cell;Cell;Position [mm]",154,0,154);
	TH1F* hPoPosSigmaPerCell = new TH1F("hPoPosSigmaPerCell","Po Position RMS Per Cell;Cell;Position [mm]",154,0,154);
	TH1F* hRnPoDzMeanPerCell = new TH1F("hRnPoDzMeanPerCell","RnPo dZ Mean Per Cell;Cell;dZ [mm]",154,0,154);
	TH1F* hRnPoDzSigmaPerCell = new TH1F("hRnPoDzSigmaPerCell","RnPo dZ Sigma Per Cell;Cell;dZ [mm]",154,0,154);

	//============================================================================================================
	for(int i=0;i<numBins;i++){
		TGraphErrors *grRatePerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrRatePerCell.push_back(grRatePerCell_new);

		TGraphErrors *grRelRatePerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrRelRatePerCell.push_back(grRelRatePerCell_new);

		TGraphErrors *grPromptPSDEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrPromptPSDEffPerCell.push_back(grPromptPSDEffPerCell_new);

		TGraphErrors *grDelayPSDEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrDelayPSDEffPerCell.push_back(grDelayPSDEffPerCell_new);

		TGraphErrors *grPromptEnEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrPromptEnEffPerCell.push_back(grPromptEnEffPerCell_new);

		TGraphErrors *grDelayEnEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrDelayEnEffPerCell.push_back(grDelayEnEffPerCell_new);

		TGraphErrors *grPosEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrPosEffPerCell.push_back(grPosEffPerCell_new);

		TGraphErrors *grEffPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrEffPerCell.push_back(grEffPerCell_new);

		TGraphErrors *grPoEnMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrPoEnMeanPerCell.push_back(grPoEnMeanPerCell_new);

		TGraphErrors *grPoEnSigmaPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrPoEnSigmaPerCell.push_back(grPoEnSigmaPerCell_new);

		// The mean of the Po Energy distribution in keV, relative to cell 76
		TGraphErrors *grRelPoEnMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrRelPoEnMeanPerCell.push_back(grRelPoEnMeanPerCell_new);

		TGraphErrors *grRelPoEnSigmaPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);
		vgrRelPoEnSigmaPerCell.push_back(grRelPoEnSigmaPerCell_new);

		TGraphErrors *grRnPSDMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrRnPSDMeanPerCell.push_back(grRnPSDMeanPerCell_new);

		TGraphErrors *grPoPSDMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrPoPSDMeanPerCell.push_back(grPoPSDMeanPerCell_new);

		TGraphErrors *grRnPoDzMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrRnPoDzMeanPerCell.push_back(grRnPoDzMeanPerCell_new);

		TGraphErrors *grRnPoDzSigmaPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrRnPoDzSigmaPerCell.push_back(grRnPoDzSigmaPerCell_new);

		TGraphErrors *grPoPosMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrPoPosMeanPerCell.push_back(grPoPosMeanPerCell_new);

		TGraphErrors *grPoPosSigmaPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrPoPosSigmaPerCell.push_back(grPoPosSigmaPerCell_new);

		TGraphErrors *grRelRnPoDzMeanPerCell_new = new TGraphErrors(NumActiveCells,xPerCell,y,xErrPerCell,yErr);	
		vgrRelRnPoDzMeanPerCell.push_back(grRelRnPoDzMeanPerCell_new);


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
		//Graphs for ET cells

		TGraphErrors *grRatePerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrRatePerCellET.push_back(grRatePerCellET_new);

		TGraphErrors *grRelRatePerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrRelRatePerCellET.push_back(grRelRatePerCellET_new);

		TGraphErrors *grPromptPSDEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrPromptPSDEffPerCellET.push_back(grPromptPSDEffPerCellET_new);

		TGraphErrors *grDelayPSDEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrDelayPSDEffPerCellET.push_back(grDelayPSDEffPerCellET_new);

		TGraphErrors *grPromptEnEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrPromptEnEffPerCellET.push_back(grPromptEnEffPerCellET_new);

		TGraphErrors *grDelayEnEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrDelayEnEffPerCellET.push_back(grDelayEnEffPerCellET_new);

		TGraphErrors *grPosEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrPosEffPerCellET.push_back(grPosEffPerCellET_new);

		TGraphErrors *grEffPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrEffPerCellET.push_back(grEffPerCellET_new);

		TGraphErrors *grPoEnMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrPoEnMeanPerCellET.push_back(grPoEnMeanPerCellET_new);

		TGraphErrors *grPoEnSigmaPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrPoEnSigmaPerCellET.push_back(grPoEnSigmaPerCellET_new);

		// The mean of the Po Energy distribution in keV, relative to cell 76
		TGraphErrors *grRelPoEnMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrRelPoEnMeanPerCellET.push_back(grRelPoEnMeanPerCellET_new);

		TGraphErrors *grRelPoEnSigmaPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);
		vgrRelPoEnSigmaPerCellET.push_back(grRelPoEnSigmaPerCellET_new);

		TGraphErrors *grRnPSDMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrRnPSDMeanPerCellET.push_back(grRnPSDMeanPerCellET_new);

		TGraphErrors *grPoPSDMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrPoPSDMeanPerCellET.push_back(grPoPSDMeanPerCellET_new);

		TGraphErrors *grRnPoDzMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrRnPoDzMeanPerCellET.push_back(grRnPoDzMeanPerCellET_new);

		TGraphErrors *grRnPoDzSigmaPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrRnPoDzSigmaPerCellET.push_back(grRnPoDzSigmaPerCellET_new);

		TGraphErrors *grPoPosMeanPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrPoPosMeanPerCellET.push_back(grPoPosMeanPerCellET_new);

		TGraphErrors *grPoPosSigmaPerCellET_new = new TGraphErrors(NumActiveCellsET,xPerCellET,y,xErrPerCellET,yErr);	
		vgrPoPosSigmaPerCellET.push_back(grPoPosSigmaPerCellET_new);
	}

	//--------------------------------------------------------------------------
	//Fill histograms
	vector<vector<double>> vRate, vN, vLifetime, vPSDEff, vEnEff, vPosEff, vTotEff;
	vector<vector<double>> vPoEnMean, vPoEnSigma, vRnPSDMean, vPoPSDMean, vRnPoDzMean, vRnPoDzSigma, vPoPosMean, vPoPosSigma;
	double tstamp;

	int n=0;
	while(IDX>=0){
		TH1F* hRate = new TH1F("hRate","Rate Per Cell;Rate [mHz];Counts",hNumBins,3,3.5);
		TH1F* hN = new TH1F("hN","N;N [arb];Counts",hNumBins,6300,7000);
		TH1F* hLifetime = new TH1F("hLifetime","Lifetime;#tau [ms];Counts",hNumBins,2.3,2.9);
		TH1F* hTotEff = new TH1F("hTotEff","Total Efficiency;Efficiency;Counts",hNumBins,0.94,1.0);
		TH1F* hPoEnMeanPerCell = new TH1F("hPoEnMeanPerCell","Po Energy Mean Per Cell;Energy [MeVee];Counts",hNumBins,0.79,0.83);
		TH1F* hDzMeanPerCell = new TH1F("hDzMeanPerCell","Dz Mean Per Cell;#Deltaz [mm];Counts",hNumBins,-5,5);
cout<<"GETTING RESULTS: bin "<<n<<endl;		

//		if(PLOTFLAG==2 && n==(numBins-1)) NUMTREES = NUMTREES+2;

		auto results = ReadAc227Trees(NUMCELLS, IDX, NUMTREES, n, PLOTFLAG, PSDCUTLOW, PSDCUTHIGH);
		IDX 	  	 = get<0>(results);
		tstamp    	 = get<1>(results);

		vRate 	  	 = get<2>(results);
		vN        	 = get<3>(results);
		vLifetime    = get<4>(results);
		vPSDEff      = get<5>(results);
		vEnEff       = get<6>(results);
		vPosEff      = get<7>(results);
		vTotEff      = get<8>(results);

		vPoEnMean    = get<9>(results);	
		vPoEnSigma   = get<10>(results);
		vRnPSDMean   = get<11>(results);
		vPoPSDMean   = get<12>(results);	
		vRnPoDzMean  = get<13>(results);
		vRnPoDzSigma = get<14>(results);
		vPoPosMean   = get<15>(results);
		vPoPosSigma  = get<16>(results);

		if(IDX<0) break;

		int relIdx = 140;
		double rateCellRel = vRate[0][relIdx];

		double totalRate = vRate[0][NUMCELLS];
		double totalRateErr = vRate[1][NUMCELLS];

		int grPt = 0, grPtET = 0;

		double centerCellPoEnMean = vPoEnMean[0][76], centerCellPoEnMeanErr = vPoEnMean[1][76];
		double centerCellPoEnSigma = vPoEnSigma[0][76], centerCellPoEnSigmaErr = vPoEnSigma[1][76];

		if(PLOTFLAG==1){
		for(int i=0;i<NUMCELLS;i++){

			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
			if(!exclude){
			
			int biny = i+1;

			hRate->Fill(vRate[0][i]*1e6);
			hN->Fill(vN[0][i]);
			hLifetime->Fill(vLifetime[0][i]);
			hTotEff->Fill(vTotEff[0][i]);
			hPoEnMeanPerCell->Fill(vPoEnMean[0][i]);
			hDzMeanPerCell->Fill(vRnPoDzMean[0][i]);

			//Heat maps of rate	
			int binx = i%14 + 1;
			biny = (i/14) + 1;
			hRatePerCell->SetBinContent(binx,biny,vRate[0][i]*1e6);

			double relRate = vRate[0][i]/rateCellRel;
			hRelRatePerCell->SetBinContent(binx,biny,relRate);

			//Graphs of variables per cell
			vgrRatePerCell[n]->SetPoint(grPt,i,vRate[0][i]*1e6);
			vgrRatePerCell[n]->SetPointError(grPt,0,vRate[1][i]*1e6);

			vgrPoEnMeanPerCell[n]->SetPoint(grPt,i,vPoEnMean[0][i]*1000);
			vgrPoEnMeanPerCell[n]->SetPointError(grPt,0,vPoEnMean[1][i]*1000);

			hPoEnMean_PerCell->SetBinContent(i,vPoEnMean[0][i]*1000);
			hPoEnMean_PerCell->SetBinError(i,vPoEnMean[1][i]*1000);
			
			vgrPoEnSigmaPerCell[n]->SetPoint(grPt,i,vPoEnSigma[0][i]*1000);
			vgrPoEnSigmaPerCell[n]->SetPointError(grPt,0,vPoEnSigma[1][i]*1000);

			hPoEnSigmaPerCell->SetBinContent(i,vPoEnSigma[0][i]*1000);
			hPoEnSigmaPerCell->SetBinError(i,vPoEnSigma[1][i]*1000);

			vgrPromptPSDEffPerCell[n]->SetPoint(grPt,i,vPSDEff[0][i]);
			vgrPromptPSDEffPerCell[n]->SetPointError(grPt,0,vPSDEff[1][i]);

			vgrDelayPSDEffPerCell[n]->SetPoint(grPt,i,vPSDEff[2][i]);
			vgrDelayPSDEffPerCell[n]->SetPointError(grPt,0,vPSDEff[3][i]);

			vgrPromptEnEffPerCell[n]->SetPoint(grPt,i,vEnEff[0][i]);
			vgrPromptEnEffPerCell[n]->SetPointError(grPt,0,vEnEff[1][i]);

			vgrDelayEnEffPerCell[n]->SetPoint(grPt,i,vEnEff[2][i]);
			vgrDelayEnEffPerCell[n]->SetPointError(grPt,0,vEnEff[3][i]);

			vgrPosEffPerCell[n]->SetPoint(grPt,i,vPosEff[0][i]);
			vgrPosEffPerCell[n]->SetPointError(grPt,0,vPosEff[1][i]);

			vgrEffPerCell[n]->SetPoint(grPt,i,vTotEff[0][i]);
			vgrEffPerCell[n]->SetPointError(grPt,0,vTotEff[1][i]);

			vgrRnPSDMeanPerCell[n]->SetPoint(grPt,i,vRnPSDMean[0][i]);
			vgrRnPSDMeanPerCell[n]->SetPointError(grPt,0,vRnPSDMean[1][i]);
	
			vgrPoPSDMeanPerCell[n]->SetPoint(grPt,i,vPoPSDMean[0][i]);
			vgrPoPSDMeanPerCell[n]->SetPointError(grPt,0,vPoPSDMean[1][i]);

			vgrRnPoDzMeanPerCell[n]->SetPoint(grPt,i,vRnPoDzMean[0][i]);
			vgrRnPoDzMeanPerCell[n]->SetPointError(grPt,0,vRnPoDzMean[1][i]);

			hRnPoDzMeanPerCell->SetBinContent(i,vRnPoDzMean[0][i]);
			hRnPoDzMeanPerCell->SetBinError(i,vRnPoDzMean[1][i]);

			vgrRnPoDzSigmaPerCell[n]->SetPoint(grPt,i,vRnPoDzSigma[0][i]);
			vgrRnPoDzSigmaPerCell[n]->SetPointError(grPt,0,vRnPoDzSigma[1][i]);

			hRnPoDzSigmaPerCell->SetBinContent(i,vRnPoDzSigma[0][i]);
			hRnPoDzSigmaPerCell->SetBinError(i,vRnPoDzSigma[1][i]);

			vgrPoPosMeanPerCell[n]->SetPoint(grPt,i,vPoPosMean[0][i]);
			vgrPoPosMeanPerCell[n]->SetPointError(grPt,0,vPoPosMean[1][i]);

			hPoPosMeanPerCell->SetBinContent(i,vPoPosMean[0][i]);
			hPoPosMeanPerCell->SetBinError(i,vPoPosMean[1][i]);

			vgrPoPosSigmaPerCell[n]->SetPoint(grPt,i,vPoPosSigma[0][i]);
			vgrPoPosSigmaPerCell[n]->SetPointError(grPt,0,vPoPosSigma[1][i]);

			hPoPosSigmaPerCell->SetBinContent(i,vPoPosSigma[0][i]);
			hPoPosSigmaPerCell->SetBinError(i,vPoPosSigma[1][i]);

			grPt++;

				ET = find(begin(ETCellArr), end(ETCellArr), i) != end(ETCellArr);
				if(ET){
				//Graphs of variables per cell
				vgrRatePerCellET[n]->SetPoint(grPtET,i,vRate[0][i]*1e6);
				vgrRatePerCellET[n]->SetPointError(grPtET,0,vRate[1][i]*1e6);

				vgrPoEnMeanPerCellET[n]->SetPoint(grPtET,i,vPoEnMean[0][i]*1000);
				vgrPoEnMeanPerCellET[n]->SetPointError(grPtET,0,vPoEnMean[1][i]*1000);
			
				vgrPoEnSigmaPerCellET[n]->SetPoint(grPtET,i,vPoEnSigma[0][i]*1000);
				vgrPoEnSigmaPerCellET[n]->SetPointError(grPtET,0,vPoEnSigma[1][i]*1000);

				vgrPromptPSDEffPerCellET[n]->SetPoint(grPtET,i,vPSDEff[0][i]);
				vgrPromptPSDEffPerCellET[n]->SetPointError(grPtET,0,vPSDEff[1][i]);

				vgrDelayPSDEffPerCellET[n]->SetPoint(grPtET,i,vPSDEff[2][i]);
				vgrDelayPSDEffPerCellET[n]->SetPointError(grPtET,0,vPSDEff[3][i]);

				vgrPromptEnEffPerCellET[n]->SetPoint(grPtET,i,vEnEff[0][i]);
				vgrPromptEnEffPerCellET[n]->SetPointError(grPtET,0,vEnEff[1][i]);

				vgrDelayEnEffPerCellET[n]->SetPoint(grPtET,i,vEnEff[2][i]);
				vgrDelayEnEffPerCellET[n]->SetPointError(grPtET,0,vEnEff[3][i]);

				vgrPosEffPerCellET[n]->SetPoint(grPtET,i,vPosEff[0][i]);
				vgrPosEffPerCellET[n]->SetPointError(grPtET,0,vPosEff[1][i]);

				vgrEffPerCellET[n]->SetPoint(grPtET,i,vTotEff[0][i]);
				vgrEffPerCellET[n]->SetPointError(grPtET,0,vTotEff[1][i]);

				vgrRnPSDMeanPerCellET[n]->SetPoint(grPtET,i,vRnPSDMean[0][i]);
				vgrRnPSDMeanPerCellET[n]->SetPointError(grPtET,0,vRnPSDMean[1][i]);
	
				vgrPoPSDMeanPerCellET[n]->SetPoint(grPtET,i,vPoPSDMean[0][i]);
				vgrPoPSDMeanPerCellET[n]->SetPointError(grPtET,0,vPoPSDMean[1][i]);

				vgrRnPoDzMeanPerCellET[n]->SetPoint(grPtET,i,vRnPoDzMean[0][i]);
				vgrRnPoDzMeanPerCellET[n]->SetPointError(grPtET,0,vRnPoDzMean[1][i]);

				vgrRnPoDzSigmaPerCellET[n]->SetPoint(grPtET,i,vRnPoDzSigma[0][i]);
				vgrRnPoDzSigmaPerCellET[n]->SetPointError(grPtET,0,vRnPoDzSigma[1][i]);

				vgrPoPosMeanPerCellET[n]->SetPoint(grPtET,i,vPoPosMean[0][i]);
				vgrPoPosMeanPerCellET[n]->SetPointError(grPtET,0,vPoPosMean[1][i]);

				vgrPoPosSigmaPerCellET[n]->SetPoint(grPtET,i,vPoPosSigma[0][i]);
				vgrPoPosSigmaPerCellET[n]->SetPointError(grPtET,0,vPoPosSigma[1][i]);
			
				grPtET++;
				}
			}	
		}

		} //if(PLOTFLAG==1)
		grRate_AllCells->SetPoint(n,tstamp,vRate[0][NUMCELLS]*1000);	//convert to Hz
		grRate_AllCells->SetPointError(n,0,vRate[1][NUMCELLS]*1000);

		hgrRate_AllCells->SetBinContent(n,vRate[0][NUMCELLS]*1000);
		hgrRate_AllCells->SetBinError(n,vRate[1][NUMCELLS]*1000);
	
		grN_AllCells->SetPoint(n,tstamp,vN[0][NUMCELLS]);
		grN_AllCells->SetPointError(n,0,vN[1][NUMCELLS]);
	
		grLifetime_AllCells->SetPoint(n,tstamp,vLifetime[0][NUMCELLS]);
		grLifetime_AllCells->SetPointError(n,0,vLifetime[1][NUMCELLS]);
cout<<"Filling effciency graphs"<<endl;

		grPromptEnEff_AllCells->SetPoint(n,tstamp,vEnEff[0][NUMCELLS]);
		grPromptEnEff_AllCells->SetPointError(n,0,vEnEff[1][NUMCELLS]);

		grDelayEnEff_AllCells->SetPoint(n,tstamp,vEnEff[2][NUMCELLS]);
		grDelayEnEff_AllCells->SetPointError(n,0,vEnEff[3][NUMCELLS]);

		grPromptPSDEff_AllCells->SetPoint(n,tstamp,vPSDEff[0][NUMCELLS]);
		grPromptPSDEff_AllCells->SetPointError(n,0,vPSDEff[1][NUMCELLS]);

		grDelayPSDEff_AllCells->SetPoint(n,tstamp,vPSDEff[2][NUMCELLS]);
		grDelayPSDEff_AllCells->SetPointError(n,0,vPSDEff[3][NUMCELLS]);

		grPosEff_AllCells->SetPoint(n,tstamp,vPosEff[0][NUMCELLS]);
		grPosEff_AllCells->SetPointError(n,0,vPosEff[1][NUMCELLS]);

		grTotEff_AllCells->SetPoint(n,tstamp,vTotEff[0][NUMCELLS]);
		grTotEff_AllCells->SetPointError(n,0,vTotEff[1][NUMCELLS]);	
	
		grPoEnMean_AllCells->SetPoint(n,tstamp,vPoEnMean[0][NUMCELLS]);
		grPoEnMean_AllCells->SetPointError(n,0,vPoEnMean[1][NUMCELLS]);	
	
		grPoEnSigma_AllCells->SetPoint(n,tstamp,vPoEnSigma[0][NUMCELLS]);
		grPoEnSigma_AllCells->SetPointError(n,0,vPoEnSigma[1][NUMCELLS]);	

		grPoPSDMean_AllCells->SetPoint(n,tstamp,vPoPSDMean[0][NUMCELLS]);
		grPoPSDMean_AllCells->SetPointError(n,0,vPoPSDMean[1][NUMCELLS]);	

		grRnPoDzMean_AllCells->SetPoint(n,tstamp,vRnPoDzMean[0][NUMCELLS]);
		grRnPoDzMean_AllCells->SetPointError(n,0,vRnPoDzMean[1][NUMCELLS]);

		grRnPoDzSigma_AllCells->SetPoint(n,tstamp,vRnPoDzSigma[0][NUMCELLS]);
		grRnPoDzSigma_AllCells->SetPointError(n,0,vRnPoDzSigma[1][NUMCELLS]);

		grPoPosMean_AllCells->SetPoint(n,tstamp,vPoPosMean[0][NUMCELLS]);
		grPoPosMean_AllCells->SetPointError(n,0,vPoPosMean[1][NUMCELLS]);

		grPoPosSigma_AllCells->SetPoint(n,tstamp,vPoPosSigma[0][NUMCELLS]);
		grPoPosSigma_AllCells->SetPointError(n,0,vPoPosSigma[1][NUMCELLS]);

		//=================================================================================================
		//Plot histograms 
		
		if(PLOTFLAG==1){
		TCanvas *cRate = new TCanvas("cRate","Rate Per Cell",1);
		hRate->Draw();
		hRate->Fit("gaus");
		cRate->SaveAs(Form("%s/HistRatePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		TCanvas *cN = new TCanvas("cN","N Per Cell",1);
		hN->Draw();
		cN->SaveAs(Form("%s/HistNPerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));
	
		TCanvas *cLifetime = new TCanvas("cLifetime","Lifetime Per Cell",1);
		hLifetime->Draw();
		cLifetime->SaveAs(Form("%s/HistLifetimePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		TCanvas *cTotEff = new TCanvas("cTotEff","TotEff Per Cell",1);
		hTotEff->Draw();
		cTotEff->SaveAs(Form("%s/HistTotEffPerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		TCanvas *cPoEnMeanPerCell = new TCanvas("cPoEnMeanPerCell","Po Energy Mean Per Cell",1);
		hPoEnMeanPerCell->Draw();
		cPoEnMeanPerCell->SaveAs(Form("%s/HistPoEnergyMeanPerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		TCanvas *cDzMeanPerCell = new TCanvas("cDzMeanPerCell","Po Energy Mean Per Cell",1);
		hDzMeanPerCell->Draw();
		hDzMeanPerCell->Fit("gaus");
		cDzMeanPerCell->SaveAs(Form("%s/HistDzMeanPerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		gStyle->SetPaintTextFormat("2.2f");
//		gStyle->SetOptStat(0);
		TCanvas *cRatePerCell = new TCanvas("cRatePerCell","Rate per Cell",1);
		gPad->SetRightMargin(0.16);
		hRatePerCell->SetMinimum(3.0);
		hRatePerCell->Draw("colz && text");
		cRatePerCell->SaveAs(Form("%s/RatePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));
	
		TCanvas *cRelRatePerCell = new TCanvas("cRelRatePerCell","Relative Rate per Cell",1);
		gPad->SetRightMargin(0.16);
		hRelRatePerCell->SetMinimum(0.95);
		hRelRatePerCell->Draw("colz && text");
		cRelRatePerCell->SaveAs(Form("%s/RelRatePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		TFile *f = new TFile("AD_PhysNeutrino_RateHist_noZFid_TotDetector.root","RECREATE");
		hRatePerCell->Write();
		f->Close();

		delete cRate;
		delete cN;
		delete cLifetime;
		delete cTotEff;
		delete cPoEnMeanPerCell;
		delete cRatePerCell;
		delete cRelRatePerCell;
		delete cDzMeanPerCell;

		}
		delete hRate;
		delete hN;
		delete hLifetime;
		delete hTotEff;
		delete hPoEnMeanPerCell;
		delete hDzMeanPerCell;

		n = n + 1;
	}


	//--------------------------------------------------------------------------
	//Plot Graphs
	if(PLOTFLAG==2){
	gStyle->SetOptStat(11);
	gStyle->SetOptFit(1111);
cout<<"Make canvas"<<endl;
	TCanvas *cRateAllCells = new TCanvas("cRateAllCells","Rate All Cells",1);
	gPad->SetGrid();
	grRate_AllCells->SetTitle("Total Rate of Detector");
	grRate_AllCells->GetYaxis()->SetTitle("Rate (Hz)");
	grRate_AllCells->GetYaxis()->SetTitleOffset(1.3);
	grRate_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grRate_AllCells->SetMarkerStyle(20);
	grRate_AllCells->SetMarkerColor(kBlue);
	grRate_AllCells->SetLineColor(kBlue);
	grRate_AllCells->Draw("AP");
	grRate_AllCells->Fit("pol0");
	grRate_AllCells->GetXaxis()->SetTimeDisplay(1);
	grRate_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cRateAllCells->SaveAs(Form("%s/Rate_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cNAllCells = new TCanvas("cNAllCells","N All Cells",1);
	gPad->SetGrid();
	grN_AllCells->SetTitle("Total N of Detector");
	grN_AllCells->GetYaxis()->SetTitle("N");
	grN_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grN_AllCells->SetMarkerStyle(20);
	grN_AllCells->SetMarkerColor(kBlue);
	grN_AllCells->SetLineColor(kBlue);
	grN_AllCells->Draw("AP");
	grN_AllCells->Fit("pol0");
	grN_AllCells->GetXaxis()->SetTimeDisplay(1);
	grN_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cNAllCells->SaveAs(Form("%s/N_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cLifetimeAllCells = new TCanvas("cLifetimeAllCells","Lifetime All Cells",1);
	gPad->SetGrid();
	grLifetime_AllCells->SetTitle("Total Po-215 Lifetime of Detector");
	grLifetime_AllCells->GetYaxis()->SetTitle("Lifetime (ms)");
	grLifetime_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grLifetime_AllCells->SetMarkerStyle(20);
	grLifetime_AllCells->SetMarkerColor(kBlue);
	grLifetime_AllCells->SetLineColor(kBlue);
	grLifetime_AllCells->Draw("AP");
	grLifetime_AllCells->Fit("pol0");
	grLifetime_AllCells->GetXaxis()->SetTimeDisplay(1);
	grLifetime_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cLifetimeAllCells->SaveAs(Form("%s/Lifetime_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cAllEffAllCells = new TCanvas("cAllEffAllCells","All Efficiency All Cells",1);
	gPad->SetGrid();
	grTotEff_AllCells->SetTitle("Cut Efficiencies;;Eff");
	grTotEff_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grTotEff_AllCells->GetXaxis()->SetTimeDisplay(1);
	grTotEff_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	grTotEff_AllCells->SetMarkerStyle(34);
	grTotEff_AllCells->SetMarkerColor(kBlack);
	grTotEff_AllCells->SetLineColor(kBlack);
	grTotEff_AllCells->Draw("AP");

	grPromptEnEff_AllCells->SetMarkerStyle(20);
	grPromptEnEff_AllCells->SetMarkerColor(kBlue);
	grPromptEnEff_AllCells->SetLineColor(kBlue);
	grPromptEnEff_AllCells->Draw("P");
	grDelayEnEff_AllCells->SetMarkerStyle(20);
	grDelayEnEff_AllCells->SetMarkerColor(kMagenta);
	grDelayEnEff_AllCells->SetLineColor(kMagenta);
	grDelayEnEff_AllCells->Draw("P");

	grPromptPSDEff_AllCells->SetMarkerStyle(21);
	grPromptPSDEff_AllCells->SetMarkerColor(kBlue);
	grPromptPSDEff_AllCells->SetLineColor(kBlue);
	grPromptPSDEff_AllCells->Draw("P");
	grDelayPSDEff_AllCells->SetMarkerStyle(21);
	grDelayPSDEff_AllCells->SetMarkerColor(kMagenta);
	grDelayPSDEff_AllCells->SetLineColor(kMagenta);
	grDelayPSDEff_AllCells->Draw("P");

	grPosEff_AllCells->SetMarkerStyle(22);
	grPosEff_AllCells->SetMarkerColor(8);
	grPosEff_AllCells->SetLineColor(8);
	grPosEff_AllCells->Draw("P");

	TLegend *leg = new TLegend(.78,.78,.98,.98);
	leg->AddEntry(grPromptPSDEff_AllCells,"Prompt PSD","p");
	leg->AddEntry(grDelayPSDEff_AllCells,"Delay PSD","p");
	leg->AddEntry(grPromptEnEff_AllCells,"Prompt Energy","p");
	leg->AddEntry(grDelayEnEff_AllCells,"Delay Energy","p");
	leg->AddEntry(grPosEff_AllCells,"Position","p");
	leg->AddEntry(grTotEff_AllCells,"Total","p");
	leg->Draw();
	cAllEffAllCells->SaveAs(Form("%s/grAllEffAllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cTotEffAllCells = new TCanvas("cTotEffAllCells","Tot Eff All Cells",1);
	gPad->SetGrid();
	grTotEff_AllCells->SetTitle("Total Efficiency of Detector");
	grTotEff_AllCells->GetYaxis()->SetTitle("Eff");
	grTotEff_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grTotEff_AllCells->SetMarkerStyle(20);
	grTotEff_AllCells->SetMarkerColor(kBlue);
	grTotEff_AllCells->SetLineColor(kBlue);
	grTotEff_AllCells->GetYaxis()->SetRangeUser(0.95,1.02);
	grTotEff_AllCells->Draw("AP");
	grTotEff_AllCells->Fit("pol0");
	grTotEff_AllCells->GetXaxis()->SetTimeDisplay(1);
	grTotEff_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cTotEffAllCells->SaveAs(Form("%s/TotEff_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoEnMeanAllCells = new TCanvas("cPoEnMeanAllCells","Po Mean All Cells",1);
	gPad->SetGrid();
	grPoEnMean_AllCells->SetTitle("Total Po-215 Energy Mean of Detector");
	grPoEnMean_AllCells->GetYaxis()->SetTitle("Energy [MeVee]");
	grPoEnMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoEnMean_AllCells->GetYaxis()->SetTitleOffset(1.43);
	grPoEnMean_AllCells->SetMarkerStyle(20);
	grPoEnMean_AllCells->SetMarkerColor(kBlue);
	grPoEnMean_AllCells->SetLineColor(kBlue);
	grPoEnMean_AllCells->Draw("AP");
	grPoEnMean_AllCells->Fit("pol0");
	grPoEnMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoEnMean_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cPoEnMeanAllCells->SaveAs(Form("%s/PoEnMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoEnSigmaAllCells = new TCanvas("cPoEnSigmaAllCells","Po Sigma All Cells",1);
	gPad->SetGrid();
	grPoEnSigma_AllCells->SetTitle("Total Po-215 Energy Sigma of Detector");
	grPoEnSigma_AllCells->GetYaxis()->SetTitle("Energy [MeVee]");
	grPoEnSigma_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoEnSigma_AllCells->GetYaxis()->SetTitleOffset(1.43);
	grPoEnSigma_AllCells->SetMarkerStyle(20);
	grPoEnSigma_AllCells->SetMarkerColor(kBlue);
	grPoEnSigma_AllCells->SetLineColor(kBlue);
	grPoEnSigma_AllCells->Draw("AP");
	grPoEnSigma_AllCells->Fit("pol1");
	grPoEnSigma_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoEnSigma_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cPoEnSigmaAllCells->SaveAs(Form("%s/PoEnSigma_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPSDMeanAllCells = new TCanvas("cPoPSDMeanAllCells","Po PSD Mean All Cells",1);
	gPad->SetGrid();
	grPoPSDMean_AllCells->SetTitle("Total Po-215 PSD Mean of Detector");
	grPoPSDMean_AllCells->GetYaxis()->SetTitle("PSD");
	grPoPSDMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoPSDMean_AllCells->SetMarkerStyle(20);
	grPoPSDMean_AllCells->SetMarkerColor(kBlue);
	grPoPSDMean_AllCells->SetLineColor(kBlue);
	grPoPSDMean_AllCells->Draw("AP");
	grPoPSDMean_AllCells->Fit("pol1");
	grPoPSDMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoPSDMean_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cPoPSDMeanAllCells->SaveAs(Form("%s/PoPSDMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRnPoDzMeanAllCells = new TCanvas("cRnPoDzMeanAllCells","RnPo Dz Mean All Cells",1);
	gPad->SetGrid();
	grRnPoDzMean_AllCells->SetTitle("Total RnPo Dz Mean of Detector");
	grRnPoDzMean_AllCells->GetYaxis()->SetTitle("Dz [mm]");
	grRnPoDzMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grRnPoDzMean_AllCells->SetMarkerStyle(20);
	grRnPoDzMean_AllCells->SetMarkerColor(kBlue);
	grRnPoDzMean_AllCells->SetLineColor(kBlue);
	grRnPoDzMean_AllCells->Draw("AP");
	grRnPoDzMean_AllCells->Fit("pol0");
	grRnPoDzMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grRnPoDzMean_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cRnPoDzMeanAllCells->SaveAs(Form("%s/RnPoDzMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRnPoDzSigmaAllCells = new TCanvas("cRnPoDzSigmaAllCells","RnPo Dz Sigma All Cells",1);
	gPad->SetGrid();
	grRnPoDzSigma_AllCells->SetTitle("Total RnPo Dz Sigma of Detector");
	grRnPoDzSigma_AllCells->GetYaxis()->SetTitle("Dz [mm]");
	grRnPoDzSigma_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grRnPoDzSigma_AllCells->SetMarkerStyle(20);
	grRnPoDzSigma_AllCells->SetMarkerColor(kBlue);
	grRnPoDzSigma_AllCells->SetLineColor(kBlue);
	grRnPoDzSigma_AllCells->Draw("AP");
	grRnPoDzSigma_AllCells->Fit("pol0");
	grRnPoDzSigma_AllCells->GetXaxis()->SetTimeDisplay(1);
	grRnPoDzSigma_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cRnPoDzSigmaAllCells->SaveAs(Form("%s/RnPoDzSigma_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPosMeanAllCells = new TCanvas("cPoPosMeanAllCells","RnPo Dz Mean All Cells",1);
	gPad->SetGrid();
	grPoPosMean_AllCells->SetTitle("Total Po Position Mean of Detector");
	grPoPosMean_AllCells->GetYaxis()->SetTitle("Mean");
	grPoPosMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoPosMean_AllCells->SetMarkerStyle(20);
	grPoPosMean_AllCells->SetMarkerColor(kBlue);
	grPoPosMean_AllCells->SetLineColor(kBlue);
	grPoPosMean_AllCells->Draw("AP");
	grPoPosMean_AllCells->Fit("pol0");
	grPoPosMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoPosMean_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cPoPosMeanAllCells->SaveAs(Form("%s/PoPosMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPosSigmaAllCells = new TCanvas("cPoPosSigmaAllCells","RnPo Dz Mean All Cells",1);
	gPad->SetGrid();
	grPoPosSigma_AllCells->SetTitle("Total Po Position Sigma of Detector");
	grPoPosSigma_AllCells->GetYaxis()->SetTitle("Sigma");
	grPoPosSigma_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoPosSigma_AllCells->SetMarkerStyle(20);
	grPoPosSigma_AllCells->SetMarkerColor(kBlue);
	grPoPosSigma_AllCells->SetLineColor(kBlue);
	grPoPosSigma_AllCells->Draw("AP");
	grPoPosSigma_AllCells->Fit("pol0");
	grPoPosSigma_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoPosSigma_AllCells->GetXaxis()->SetTimeFormat("%m/%d");
	cPoPosSigmaAllCells->SaveAs(Form("%s/PoPosSigma_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	//==============================================================
	//create relative graphs
	double x[numBins], y[numBins];
	double xErr[numBins], yErr[numBins];

	TGraphErrors *grRelRate_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRelPoEnMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRelPoEnSigma_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRelPoPosRMS_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grRelRnPoDzMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);

	TF1 *fRatePol = new TF1("fRelRatePol","pol0");
	TF1 *fEMeanPol = new TF1("fEMeanPol","pol0");
	TF1 *fESigmaPol = new TF1("fESigmaPol","pol0");
	TF1 *fPosRMSPol = new TF1("fPosRMSPol","pol0");
	TF1 *fDzMeanPol = new TF1("fDzMeanPol","pol0");

	grRate_AllCells->Fit(fRatePol);
	grPoEnMean_AllCells->Fit(fEMeanPol);
	grPoEnSigma_AllCells->Fit(fESigmaPol);
	grPoPosSigma_AllCells->Fit(fPosRMSPol);
	grRnPoDzMean_AllCells->Fit(fDzMeanPol);

	double avgRate = fRatePol->GetParameter(0), avgRateErr = fRatePol->GetParError(0);
	double avgEMean = fEMeanPol->GetParameter(0), avgEMeanErr = fEMeanPol->GetParError(0);
	double avgESigma = fESigmaPol->GetParameter(0), avgESigmaErr = fESigmaPol->GetParError(0);
	double avgPosRMS = fPosRMSPol->GetParameter(0), avgPosRMSErr = fPosRMSPol->GetParError(0);
	double avgDzMean = fDzMeanPol->GetParameter(0), avgDzMeanErr = fDzMeanPol->GetParError(0);

	double relRate, relRateErr;
	double relEMean, relEMeanErr;
	double relESigma, relESigmaErr;
	double relPosRMS, relPosRMSErr;
	double relDzMean, relDzMeanErr;

	double grx, gry, gryErr;
	for(int i=0;i<numBins;i++){
		grRate_AllCells->GetPoint(i,grx,gry);
		gryErr = grRate_AllCells->GetErrorY(i);

		relRate = gry/avgRate;
		relRateErr = relRate * sqrt(pow(gryErr/gry,2) + pow(avgRateErr/avgRate,2));

		grRelRate_AllCells->SetPoint(i,grx,relRate);
		grRelRate_AllCells->SetPointError(i,0,relRateErr);

		//--------------------------------------------------
		grPoEnMean_AllCells->GetPoint(i,grx,gry);
		gryErr = grPoEnMean_AllCells->GetErrorY(i);
	
		relEMean = gry/avgEMean;
		relEMeanErr = relEMean * sqrt(pow(gryErr/gry,2) + pow(avgEMeanErr/avgEMean,2));	

		grRelPoEnMean_AllCells->SetPoint(i,grx,relEMean);
		grRelPoEnMean_AllCells->SetPointError(i,0,relEMeanErr);
	
		//--------------------------------------------------
		grPoEnSigma_AllCells->GetPoint(i,grx,gry);
		gryErr = grPoEnSigma_AllCells->GetErrorY(i);
	
		relESigma = gry/avgESigma;
		relESigmaErr = relESigma * sqrt(pow(gryErr/gry,2) + pow(avgESigmaErr/avgESigma,2));	

		grRelPoEnSigma_AllCells->SetPoint(i,grx,relESigma);
		grRelPoEnSigma_AllCells->SetPointError(i,0,relESigmaErr);

		//--------------------------------------------------
		grPoPosSigma_AllCells->GetPoint(i,grx,gry);
		gryErr = grPoPosSigma_AllCells->GetErrorY(i);
	
		relPosRMS = gry/avgPosRMS;
		relPosRMSErr = relPosRMS * sqrt(pow(gryErr/gry,2) + pow(avgPosRMSErr/avgPosRMS,2));	

		grRelPoPosRMS_AllCells->SetPoint(i,grx,relPosRMS);
		grRelPoPosRMS_AllCells->SetPointError(i,0,relPosRMSErr);

		//--------------------------------------------------
		grRnPoDzMean_AllCells->GetPoint(i,grx,gry);
		gryErr = grRnPoDzMean_AllCells->GetErrorY(i);
	
		relDzMean = gry/avgDzMean;
		relDzMeanErr = relDzMean * sqrt(pow(gryErr/gry,2) + pow(avgDzMeanErr/avgDzMean,2));	

		grRelRnPoDzMean_AllCells->SetPoint(i,grx,relDzMean);
		grRelRnPoDzMean_AllCells->SetPointError(i,0,relDzMeanErr);
	}


	
	//==============================================================
	TFile *saveGraphsFile = new TFile("AD_NeutrinoTGraphs_noZFidCut.root","RECREATE");
	hgrRate_AllCells->Write();
	grRate_AllCells->Write("grRate");
	grPromptEnEff_AllCells->Write("grPromptEnEff");
	grDelayEnEff_AllCells->Write("grDelayEnEff");
	grPromptPSDEff_AllCells->Write("grPromptPSDEff");
	grDelayPSDEff_AllCells->Write("grDelayPSDEff");
	grPosEff_AllCells->Write("grPosEff");
	grTotEff_AllCells->Write("grTotEff");
	grPoEnMean_AllCells->Write("grPoEnMean");
	grPoEnSigma_AllCells->Write("grPoEnSigma");
	grRnPoDzMean_AllCells->Write("grRnPoDzMean");
	grRnPoDzSigma_AllCells->Write("grRnPoDzSigma");
	grPoPosMean_AllCells->Write("grPoPosMean");
	grPoPosSigma_AllCells->Write("grPoPosSigma");

	grRelRate_AllCells->Write("grRelRate");
	grRelPoEnMean_AllCells->Write("grRelPoEnMean");
	grRelPoEnSigma_AllCells->Write("grRelPoEnSigma");
	grRelPoPosRMS_AllCells->Write("grRelPoPosRMS");
	grRelRnPoDzMean_AllCells->Write("grRelRnPoDzMean");

	saveGraphsFile->Close(); 

	}

	//==========================================================================================
	//Plots of variables per cell

	if(PLOTFLAG==1){	
	TF1 *fRatePol = new TF1("fRatePol","pol0");

	TCanvas *cRatePerCell = new TCanvas("cRatePerCell","Rate Per Cell",1);
	gPad->SetGrid();
	vgrRatePerCell[0]->SetTitle("Rate Per Cell;Cell;Rate [mHz]");
	vgrRatePerCell[0]->SetMarkerStyle(20);
	vgrRatePerCell[0]->SetMarkerColor(kBlue); 
	vgrRatePerCell[0]->SetLineColor(kBlue); 
	vgrRatePerCell[0]->Draw("AP");
	vgrRatePerCellET[0]->SetMarkerStyle(20);
	vgrRatePerCellET[0]->SetMarkerColor(kRed); 
	vgrRatePerCellET[0]->SetLineColor(kRed); 
	vgrRatePerCellET[0]->Draw("P");
	vgrRatePerCell[0]->Fit(fRatePol);
	cRatePerCell->SaveAs(Form("%s/grRatePerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRatePerCell->SaveAs(Form("%s/RatePerCell.C",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	double avgRate = fRatePol->GetParameter(0), avgRateErr = fRatePol->GetParError(0);
	int grPt = 0, grPtET = 0;
	for(int i=0;i<NUMCELLS;i++){
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(!exclude){
		double relRate = (vRate[0][i]*1e6)/avgRate;
		double relRateErr = relRate * sqrt(pow(vRate[1][i]/vRate[0][i],2)+pow(avgRateErr/avgRate,2)); 
		vgrRelRatePerCell[0]->SetPoint(grPt,i,relRate);
		vgrRelRatePerCell[0]->SetPointError(grPt,0,relRateErr);
	
		grPt++;
			ET = find(begin(ETCellArr), end(ETCellArr), i) != end(ETCellArr);
			if(ET){
			vgrRelRatePerCellET[0]->SetPoint(grPtET,i,relRate);
			vgrRelRatePerCellET[0]->SetPointError(grPtET,0,relRateErr);

			grPtET++;
			}
		}		
	}

	TCanvas *cRelRatePerCell = new TCanvas("cRelRatePerCell","Relative Rate Per Cell",1);
	gPad->SetGrid();
	vgrRelRatePerCell[0]->SetTitle("Rate Per Cell/Average Rate Per Cell;Cell;Relative Rate");
	vgrRelRatePerCell[0]->SetMarkerStyle(20);
	vgrRelRatePerCell[0]->SetMarkerColor(kBlue);
	vgrRelRatePerCell[0]->SetLineColor(kBlue);
	vgrRelRatePerCell[0]->Draw("AP");
	vgrRelRatePerCellET[0]->SetMarkerStyle(20);
	vgrRelRatePerCellET[0]->SetMarkerColor(kRed);
	vgrRelRatePerCellET[0]->SetLineColor(kRed);
	vgrRelRatePerCellET[0]->Draw("P");
	vgrRelRatePerCell[0]->Fit("pol0");
	cRelRatePerCell->SaveAs(Form("%s/grRelRatePerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRelRatePerCell->SaveAs(Form("%s/RelRatePerCell.C",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cAllEffPerCell = new TCanvas("cAllEffPerCell","AllEfficiency Per Cell",1);
	gPad->SetGrid();
	vgrPromptPSDEffPerCell[0]->SetTitle("Cut Efficiency Per Cell;Cell;Eff");
	vgrPromptPSDEffPerCell[0]->GetYaxis()->SetRangeUser(0.94,1.02);

	vgrPromptPSDEffPerCell[0]->SetMarkerStyle(20);
	vgrPromptPSDEffPerCell[0]->SetMarkerColor(kBlue);
	vgrPromptPSDEffPerCell[0]->SetLineColor(kBlue);
	vgrPromptPSDEffPerCell[0]->Draw("AP");
	vgrDelayPSDEffPerCell[0]->SetMarkerStyle(20);
	vgrDelayPSDEffPerCell[0]->SetMarkerColor(kMagenta);
	vgrDelayPSDEffPerCell[0]->SetLineColor(kMagenta);
	vgrDelayPSDEffPerCell[0]->Draw("P");

	vgrPromptPSDEffPerCellET[0]->SetMarkerStyle(24);
	vgrPromptPSDEffPerCellET[0]->SetLineColor(kBlue);
	vgrPromptPSDEffPerCellET[0]->Draw("P");
	vgrDelayPSDEffPerCellET[0]->SetMarkerStyle(24);
	vgrDelayPSDEffPerCellET[0]->SetLineColor(kMagenta);
	vgrDelayPSDEffPerCellET[0]->Draw("P");

	vgrPromptEnEffPerCell[0]->SetMarkerStyle(21);
	vgrPromptEnEffPerCell[0]->SetMarkerColor(kBlue);
	vgrPromptEnEffPerCell[0]->SetLineColor(kBlue);
	vgrPromptEnEffPerCell[0]->Draw("P");
	vgrDelayEnEffPerCell[0]->SetMarkerStyle(21);
	vgrDelayEnEffPerCell[0]->SetMarkerColor(kMagenta);
	vgrDelayEnEffPerCell[0]->SetLineColor(kMagenta);
	vgrDelayEnEffPerCell[0]->Draw("P");

	vgrPromptEnEffPerCellET[0]->SetMarkerStyle(25);
	vgrPromptEnEffPerCellET[0]->SetLineColor(kBlue);
	vgrPromptEnEffPerCellET[0]->Draw("P");
	vgrDelayEnEffPerCellET[0]->SetMarkerStyle(25);
	vgrDelayEnEffPerCellET[0]->SetLineColor(kMagenta);
	vgrDelayEnEffPerCellET[0]->Draw("P");

	vgrPosEffPerCell[0]->SetMarkerStyle(20);
	vgrPosEffPerCell[0]->SetMarkerColor(8);
	vgrPosEffPerCell[0]->SetLineColor(8);
	vgrPosEffPerCell[0]->Draw("P");

	vgrPosEffPerCellET[0]->SetMarkerStyle(24);
	vgrPosEffPerCellET[0]->SetLineColor(8);
	vgrPosEffPerCellET[0]->Draw("P");

	vgrEffPerCell[0]->SetMarkerStyle(20);
	vgrEffPerCell[0]->SetMarkerColor(kBlack);
	vgrEffPerCell[0]->SetLineColor(kBlack);
	vgrEffPerCell[0]->Draw("P");

	vgrEffPerCellET[0]->SetMarkerStyle(24);
	vgrEffPerCellET[0]->SetLineColor(kBlack);
	vgrEffPerCellET[0]->Draw("P");

	TLegend *leg = new TLegend(0.78,0.78,0.98,0.98);
	leg->AddEntry(vgrPromptPSDEffPerCell[0],"Prompt PSD","p");
	leg->AddEntry(vgrDelayPSDEffPerCell[0],"Delay PSD","p");
	leg->AddEntry(vgrPromptEnEffPerCell[0],"Prompt Energy","p");
	leg->AddEntry(vgrDelayEnEffPerCell[0],"Delay Energy","p");
	leg->AddEntry(vgrPosEffPerCell[0],"Position","p");
	leg->AddEntry(vgrEffPerCell[0],"Total","p");
	leg->Draw();
	cAllEffPerCell->SaveAs(Form("%s/grAllEffPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cEffPerCell = new TCanvas("cEffPerCell","Efficiency Per Cell",1);
	gPad->SetGrid();
	vgrEffPerCell[0]->SetTitle("Efficiency Per Cell;Cell;Eff");
	vgrEffPerCell[0]->GetYaxis()->SetRangeUser(0.94,1.02);
	vgrEffPerCell[0]->SetMarkerStyle(20);
	vgrEffPerCell[0]->SetMarkerColor(kBlue);
	vgrEffPerCell[0]->SetLineColor(kBlue);
	vgrEffPerCell[0]->Draw("AP");
	vgrEffPerCellET[0]->SetMarkerStyle(20);
	vgrEffPerCellET[0]->SetMarkerColor(kRed);
	vgrEffPerCellET[0]->SetLineColor(kRed);
	vgrEffPerCellET[0]->Draw("P");
	vgrEffPerCell[0]->Fit("pol0");
	cEffPerCell->SaveAs(Form("%s/grEffPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoEnMeanPerCell = new TCanvas("cPoEnMeanPerCell","Po Energy Mean Per Cell",1);
	gPad->SetGrid();
	vgrPoEnMeanPerCell[0]->SetTitle("Po Energy Mean Per Cell;Cell;Energy [keV]");
	vgrPoEnMeanPerCell[0]->SetMarkerStyle(20);
	vgrPoEnMeanPerCell[0]->SetMarkerColor(kBlue);
	vgrPoEnMeanPerCell[0]->SetLineColor(kBlue);
	vgrPoEnMeanPerCell[0]->Draw("AP");
	vgrPoEnMeanPerCellET[0]->SetMarkerStyle(20);
	vgrPoEnMeanPerCellET[0]->SetMarkerColor(kRed);
	vgrPoEnMeanPerCellET[0]->SetLineColor(kRed);
	vgrPoEnMeanPerCellET[0]->Draw("P");
	vgrPoEnMeanPerCell[0]->Fit("pol0");
	cPoEnMeanPerCell->SaveAs(Form("%s/PoEnMeanPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoEnSigmaPerCell = new TCanvas("cPoEnSigmaPerCell","Po Energy Sigma Per Cell",1);
	gPad->SetGrid();
	vgrPoEnSigmaPerCell[0]->SetTitle("Po Energy Sigma Per Cell;Cell;Sigma [keV]");
	vgrPoEnSigmaPerCell[0]->SetMarkerStyle(20);
	vgrPoEnSigmaPerCell[0]->SetMarkerColor(kBlue);
	vgrPoEnSigmaPerCell[0]->SetLineColor(kBlue);
	vgrPoEnSigmaPerCell[0]->Draw("AP");
	vgrPoEnSigmaPerCellET[0]->SetMarkerStyle(20);
	vgrPoEnSigmaPerCellET[0]->SetMarkerColor(kRed);
	vgrPoEnSigmaPerCellET[0]->SetLineColor(kRed);
	vgrPoEnSigmaPerCellET[0]->Draw("P");
	vgrPoEnSigmaPerCell[0]->Fit("pol0");
	cPoEnSigmaPerCell->SaveAs(Form("%s/PoEnSigmaPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRnPoPSDMeanPerCell = new TCanvas("cRnPoPSDMeanPerCell","RnPo PSD Mean Per Cell",1);
	gPad->SetGrid();
	vgrRnPSDMeanPerCell[0]->SetTitle("RnPo PSD Mean Per Cell;Cell;PSD [arb]");
	vgrRnPSDMeanPerCell[0]->SetMarkerStyle(20);
	vgrRnPSDMeanPerCell[0]->SetMarkerColor(kBlue);
	vgrRnPSDMeanPerCell[0]->SetLineColor(kBlue);
	vgrRnPSDMeanPerCell[0]->Draw("AP");
	vgrPoPSDMeanPerCell[0]->SetMarkerStyle(20);
	vgrPoPSDMeanPerCell[0]->SetMarkerColor(kMagenta);
	vgrPoPSDMeanPerCell[0]->SetLineColor(kMagenta);
	vgrPoPSDMeanPerCell[0]->Draw("P");
	cRnPoPSDMeanPerCell->SaveAs(Form("%s/RnPoPSDMeanPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRnPoDzMeanPerCell = new TCanvas("cRnPoDzMeanPerCell","RnPo Dz Mean Per Cell",1);
	gPad->SetGrid();
	vgrRnPoDzMeanPerCell[0]->SetTitle("RnPo dz Mean Per Cell;Cell;dz [mm]");
	vgrRnPoDzMeanPerCell[0]->SetMarkerStyle(20);
	vgrRnPoDzMeanPerCell[0]->SetMarkerColor(kBlue);
	vgrRnPoDzMeanPerCell[0]->SetLineColor(kBlue);
	vgrRnPoDzMeanPerCell[0]->Draw("AP");
	vgrRnPoDzMeanPerCellET[0]->SetMarkerStyle(20);
	vgrRnPoDzMeanPerCellET[0]->SetMarkerColor(kRed);
	vgrRnPoDzMeanPerCellET[0]->SetLineColor(kRed);
	vgrRnPoDzMeanPerCellET[0]->Draw("P");
	vgrRnPoDzMeanPerCell[0]->Fit("pol0");
	cRnPoDzMeanPerCell->SaveAs(Form("%s/RnPoDzMeanPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRnPoDzSigmaPerCell = new TCanvas("cRnPoDzSigmaPerCell","RnPo Dz Sigma Per Cell",1);
	gPad->SetGrid();
	vgrRnPoDzSigmaPerCell[0]->SetTitle("RnPo dz Sigma Per Cell;Cell;Sigma [mm]");
	vgrRnPoDzSigmaPerCell[0]->SetMarkerStyle(20);
	vgrRnPoDzSigmaPerCell[0]->SetMarkerColor(kBlue);
	vgrRnPoDzSigmaPerCell[0]->SetLineColor(kBlue);
	vgrRnPoDzSigmaPerCell[0]->Draw("AP");
	vgrRnPoDzSigmaPerCellET[0]->SetMarkerStyle(20);
	vgrRnPoDzSigmaPerCellET[0]->SetMarkerColor(kRed);
	vgrRnPoDzSigmaPerCellET[0]->SetLineColor(kRed);
	vgrRnPoDzSigmaPerCellET[0]->Draw("P");
	vgrRnPoDzSigmaPerCell[0]->Fit("pol0");
	cRnPoDzSigmaPerCell->SaveAs(Form("%s/RnPoDzSigmaPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPosMeanPerCell = new TCanvas("cPoPosMeanPerCell","Po Position Mean Per Cell",1);
	gPad->SetGrid();
	vgrPoPosMeanPerCell[0]->SetTitle("Po Position Mean Per Cell;Cell;Mean [mm]");
	vgrPoPosMeanPerCell[0]->SetMarkerStyle(20);
	vgrPoPosMeanPerCell[0]->SetMarkerColor(kBlue);
	vgrPoPosMeanPerCell[0]->SetLineColor(kBlue);
	vgrPoPosMeanPerCell[0]->Draw("AP");
	vgrPoPosMeanPerCellET[0]->SetMarkerStyle(20);
	vgrPoPosMeanPerCellET[0]->SetMarkerColor(kRed);
	vgrPoPosMeanPerCellET[0]->SetLineColor(kRed);
	vgrPoPosMeanPerCellET[0]->Draw("P");
	vgrPoPosMeanPerCell[0]->Fit("pol0");
	cPoPosMeanPerCell->SaveAs(Form("%s/PoPosMeanPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPosSigmaPerCell = new TCanvas("cPoPosSigmaPerCell","Po Position RMS Per Cell",1);
	gPad->SetGrid();
	vgrPoPosSigmaPerCell[0]->SetTitle("Po Position RMS Per Cell;Cell;RMS [mm]");
	vgrPoPosSigmaPerCell[0]->SetMarkerStyle(20);
	vgrPoPosSigmaPerCell[0]->SetMarkerColor(kBlue);
	vgrPoPosSigmaPerCell[0]->SetLineColor(kBlue);
	vgrPoPosSigmaPerCell[0]->Draw("AP");
	vgrPoPosSigmaPerCellET[0]->SetMarkerStyle(20);
	vgrPoPosSigmaPerCellET[0]->SetMarkerColor(kRed);
	vgrPoPosSigmaPerCellET[0]->SetLineColor(kRed);
	vgrPoPosSigmaPerCellET[0]->Draw("P");
	vgrPoPosSigmaPerCell[0]->Fit("pol0");
	cPoPosSigmaPerCell->SaveAs(Form("%s/PoPosSigmaPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	//========================================================================================================
	//Create graphs relative to the mean

	TF1 *fEnMeanPol = new TF1("fEnMeanPol","pol0");
	TF1 *fEnSigmaPol = new TF1("fEnSigmaPol","pol0");
	TF1 *fDzMeanPol = new TF1("fDzMeanPol","pol0");
	
	vgrPoEnMeanPerCell[0]->Fit(fEnMeanPol);
	vgrPoEnSigmaPerCell[0]->Fit(fEnSigmaPol);
	vgrRnPoDzMeanPerCell[0]->Fit(fDzMeanPol);

	double avgPoEnMean = fEnMeanPol->GetParameter(0)/1000.0, avgPoEnMeanErr = fEnMeanPol->GetParError(0)/1000.0;
	double avgPoEnSigma = fEnSigmaPol->GetParameter(0)/1000.0, avgPoEnSigmaErr = fEnSigmaPol->GetParError(0)/1000.0;		
	double avgRnPoDzMean = fDzMeanPol->GetParameter(0), avgRnPoDzMeanErr = fDzMeanPol->GetParError(0);

	double relPoEnMean, relPoEnMeanErr, relPoEnSigma, relPoEnSigmaErr, relRnPoDzMean, relRnPoDzMeanErr;

	grPt = 0;
	for(int i=0;i<NUMCELLS;i++){
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(!exclude){

			relPoEnMean = vPoEnMean[0][i]/avgPoEnMean;
			relPoEnMeanErr = relPoEnMean * sqrt(pow(vPoEnMean[1][i]/vPoEnMean[0][i],2) + pow(avgPoEnMeanErr/avgPoEnMean,2));

			vgrRelPoEnMeanPerCell[0]->SetPoint(grPt,i,relPoEnMean);
			vgrRelPoEnMeanPerCell[0]->SetPointError(grPt,0,relPoEnMeanErr);
		
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			relPoEnSigma = vPoEnSigma[0][i]/avgPoEnSigma;
			relPoEnSigmaErr = relPoEnSigma * sqrt(pow(vPoEnSigma[1][i]/vPoEnSigma[0][i],2) + pow(avgPoEnSigmaErr/avgPoEnSigma,2));

			vgrRelPoEnSigmaPerCell[0]->SetPoint(grPt,i,relPoEnSigma);
			vgrRelPoEnSigmaPerCell[0]->SetPointError(grPt,0,relPoEnSigmaErr);

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			relRnPoDzMean = vRnPoDzMean[0][i]/avgRnPoDzMean;
			relRnPoDzMeanErr = relRnPoDzMean * sqrt(pow(vRnPoDzMean[1][i]/vRnPoDzMean[0][i],2) + pow(avgRnPoDzMeanErr/avgRnPoDzMean,2));				

			vgrRelRnPoDzMeanPerCell[0]->SetPoint(grPt,i,relRnPoDzMean);
			vgrRelRnPoDzMeanPerCell[0]->SetPointError(grPt,0,relRnPoDzMeanErr);

			grPt++;
		}
	}


	vgrRelPoEnMeanPerCell[0]->SetTitle("Relative Po Energy Mean Per Cell;Cell;Energy [MeV]");
	vgrRelPoEnSigmaPerCell[0]->SetTitle("Relative Po Energy Sigma Per Cell;Cell;Sigma [MeV]");

/*
	TCanvas *cRelPoEnMeanPerCell = new TCanvas("cRelPoEnMeanPerCell","Po Energy Mean Per Cell",1);
	gPad->SetGrid();
	vgrRelPoEnMeanPerCell[0]->SetTitle("Relative Po Energy Mean Per Cell;Cell;Energy [keV]");
	vgrRelPoEnMeanPerCell[0]->SetMarkerStyle(20);
	vgrRelPoEnMeanPerCell[0]->SetMarkerColor(kBlue);
	vgrRelPoEnMeanPerCell[0]->SetLineColor(kBlue);
	vgrRelPoEnMeanPerCell[0]->Draw("AP");
	vgrRelPoEnMeanPerCellET[0]->SetMarkerStyle(20);
	vgrRelPoEnMeanPerCellET[0]->SetMarkerColor(kRed);
	vgrRelPoEnMeanPerCellET[0]->SetLineColor(kRed);
	vgrRelPoEnMeanPerCellET[0]->Draw("P");
	vgrRelPoEnMeanPerCell[0]->Fit("pol0");
	cRelPoEnMeanPerCell->SaveAs(Form("%s/RelPoEnMeanPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cRelPoEnSigmaPerCell = new TCanvas("cRelPoEnSigmaPerCell","Po Energy Sigma Per Cell",1);
	gPad->SetGrid();
	vgrRelPoEnSigmaPerCell[0]->SetTitle("Relative Po Energy Sigma Per Cell;Cell;Sigma [keV]");
	vgrRelPoEnSigmaPerCell[0]->SetMarkerStyle(20);
	vgrRelPoEnSigmaPerCell[0]->SetMarkerColor(kBlue);
	vgrRelPoEnSigmaPerCell[0]->SetLineColor(kBlue);
	vgrRelPoEnSigmaPerCell[0]->Draw("AP");
	vgrRelPoEnSigmaPerCellET[0]->SetMarkerStyle(20);
	vgrRelPoEnSigmaPerCellET[0]->SetMarkerColor(kRed);
	vgrRelPoEnSigmaPerCellET[0]->SetLineColor(kRed);
	vgrRelPoEnSigmaPerCellET[0]->Draw("P");
	vgrRelPoEnSigmaPerCell[0]->Fit("pol0");
	cRelPoEnSigmaPerCell->SaveAs(Form("%s/RelPoEnSigmaPerCell.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
*/



	TFile *fOutPerCell = new TFile("AD_PhysNeutrino_TGraphs_noZFid_TotDetector.root","RECREATE");
	vgrRatePerCell[0]->Write("grRatePerCell");
	vgrRatePerCellET[0]->Write("grRatePerCellET");
	vgrRelRatePerCell[0]->Write("grRelRatePerCell");
	vgrRelRatePerCellET[0]->Write("grRelRatePerCellET");
	vgrPromptPSDEffPerCell[0]->Write("grPromptPSDEff");
	vgrDelayPSDEffPerCell[0]->Write("grDelayPSDEff");
	vgrPromptEnEffPerCell[0]->Write("grPromptEnEff");
	vgrDelayEnEffPerCell[0]->Write("grDelayEnEff");
	vgrPosEffPerCell[0]->Write("grPosEff");			
	vgrEffPerCell[0]->Write("grTotEff");
	vgrEffPerCellET[0]->Write("grTotEffET");
	vgrPoEnMeanPerCell[0]->Write("grPoEnMean");
	vgrPoEnMeanPerCellET[0]->Write("grPoEnMeanET");
	vgrPoEnSigmaPerCell[0]->Write("grPoEnSigma");		
	vgrPoEnSigmaPerCellET[0]->Write("grPoEnSigmaET");		
	vgrRelPoEnMeanPerCellET[0]->Write("grRelPoEnMeanET");
	vgrRelPoEnSigmaPerCellET[0]->Write("grRelPoEnSigmaET");
	vgrRnPSDMeanPerCell[0]->Write("grRnPSDMean");
	vgrRnPSDMeanPerCellET[0]->Write("grRnPSDMeanET");
	vgrRnPoDzMeanPerCell[0]->Write("grRnPoDzMean");
	vgrRnPoDzMeanPerCellET[0]->Write("grRnPoDzMeanET");
	vgrRnPoDzSigmaPerCell[0]->Write("grRnPoDzSigma");
	vgrRnPoDzSigmaPerCellET[0]->Write("grRnPoDzSigmaET");
	vgrPoPosMeanPerCell[0]->Write("grPoPosMean");	
	vgrPoPosMeanPerCellET[0]->Write("grPoPosMeanET");	
	vgrPoPosSigmaPerCell[0]->Write("grPoPosSigma");
	vgrPoPosSigmaPerCellET[0]->Write("grPoPosSigmaET");

	vgrRelPoEnMeanPerCell[0]->Write("grRelPoEnMean");
	vgrRelPoEnSigmaPerCell[0]->Write("grRelPoEnSigma");
	vgrRelRnPoDzMeanPerCell[0]->Write("grRelRnPoDzMean");

	hPoEnMean_PerCell->Write();
	hPoEnSigmaPerCell->Write();
	hPoPosMeanPerCell->Write();
	hPoPosSigmaPerCell->Write();
	hRnPoDzMeanPerCell->Write();
	hRnPoDzSigmaPerCell->Write();
	fOutPerCell->Close();

	}

}
