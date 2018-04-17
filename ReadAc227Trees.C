#include <fstream>

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

tuple<int, double, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>>ReadAc227Trees(string TreeDir, const int NUMCELLS, int IDX, int NUMTREES, int nLoop, int PLOTFLAG){

	int ExcludeCellArr[28] = {0,1,2,3,4,5,6,9,11,12,13,18,21,23,24,27,32,40,44,68,73,79,102,107,122,127,130,139};
	bool exclude;

//===============================================================
//Set variables for Ac-227 tree
	
	printf("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+ \n");
	printf("Number of cells: %i \n",NUMCELLS);

	TChain AcChain("TAc");
	AcChain.Add(Form("%s/WetCommissioning/Series009_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series012_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series013_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series014_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series015_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series016_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series017_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series018_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series019_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series020_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series021_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series022_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Rampdown/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series001_AcTrees.root",getenv("P2X_ANALYZED")));

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
		return make_tuple(IDX, ts, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull);
	}
	else if((IDX+1)>=nAcEvents){ 
		printf("FINISHED ANALYSIS! :D \n");
		IDX = -1;
		double ts = 0.0;
		vector<double> vNull;
		vNull.push_back(0.0);
		vector<vector<double>> vvNull;
		vvNull.push_back(vNull); 
		return make_tuple(IDX, ts, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull, vvNull);
	}

//===============================================================
	AcChain.GetEvent(0);

	double promptECutLow = Ac_promptECut[0], promptECutHigh = Ac_promptECut[1];	//MeVee
	double delayECutLow  = Ac_delayECut[0],  delayECutHigh  = Ac_delayECut[1];	
	
	double promptPSDCutLow = Ac_promptPSDCut[0], promptPSDCutHigh = Ac_promptPSDCut[1];	//arb.
	double delayPSDCutLow  = Ac_delayPSDCut[0],  delayPSDCutHigh  = Ac_delayPSDCut[1];

	double posCut = Ac_zCut;	//mm


	TFile *outFile = new TFile("Ac227Histograms.root","RECREATE");

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

	double POLIFETIME = 2.57;	//[ms]
	double TIMECUT = 0.5, TIMEWINDOW = 4.5*POLIFETIME, TIMEOFFSET = 10.0*POLIFETIME;	//[ms]

	int numDtBins = 100;
	double selectDtMin = TIMECUT, selectDtMax = TIMEWINDOW;	//[ms]
	double BGDtMin = selectDtMin + TIMEOFFSET, BGDtMax = selectDtMax + TIMEOFFSET;

	int numPSDBins = 100;
	double PSDMin = 0.15, PSDMax = 0.45;	//[arb]

	int numEBins = 100;
	double EMin = 0.2, EMax = 1.4;			//[MeVee]

	int numDzBins = 100;
	double dzMin = -400.0, dzMax = 400.0;		//[mm]

	int numPosBins = 100;
	double posMin = -1200, posMax = 1200;	//[mm]

	for(int i=0;i<NUMCELLS;i++){
		exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
		if(!exclude){
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
		}
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

	for(int i=IDX;i<nAcEvents;i++){
		AcChain.GetEvent(i);

		double promptTime = Ac_t[0]*CONVERTnsTOms;

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
		if(Ac_evt[1]!=0 && !exclude){
			segNum = Ac_seg[1];

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
		}

		//check if there is a background delay event
		if(Ac_evt[2]!=0 && !exclude){
			segNum = Ac_seg[2];	
		
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
		}


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
	TF1 *fPromptEGaus   = nullptr;
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

	double rate,    N,    lifetime,    promptPSDEff,    delayPSDEff,     promptEnEff,    delayEnEff,    posEff,    totEff,    PoMean, 	 PoPSDMean;
	double rateErr, NErr, lifetimeErr, promptPSDEffErr, delayPSDEffErr,  promptEnEffErr, delayEnEffErr, posEffErr, totEffErr, PoMeanErr, PoPSDMeanErr;

	vector<double> vRate,    vN,    vLifetime,    vPromptPSDEff,    vDelayPSDEff,    vPromptEnEff,    vDelayEnEff,    vPosEff,    vTotEff,   vPoMean, vPoPSDMean;
	vector<double> vRateErr, vNErr, vLifetimeErr, vPromptPSDEffErr, vDelayPSDEffErr, vPromptEnEffErr, vDelayEnEffErr, vPosEffErr, vTotEffErr, vPoMeanErr, vPoPSDMeanErr; 

	double rateAllCells,    NAllCells,    lifetimeAllCells,    promptPSDEffAllCells,    delayPSDEffAllCells,     promptEnEffAllCells,    delayEnEffAllCells,    posEffAllCells,    totEffAllCells,    PoMeanAllCells, PoPSDMeanAllCells;
	double rateErrAllCells, NErrAllCells, lifetimeErrAllCells, promptPSDEffErrAllCells, delayPSDEffErrAllCells,  promptEnEffErrAllCells, delayEnEffErrAllCells, posEffErrAllCells, totEffErrAllCells, PoMeanErrAllCells, PoPSDMeanErrAllCells;

//===============================================================

	printf("-------------------------------------------------------------------\n");
	printf("Subtracting histograms \n");

	for(int i=0;i<NUMCELLS;i++){		
		rate = 0.0,    N = 0.0,    lifetime = 0.0,    promptPSDEff = 0.0,    delayPSDEff = 0.0,     promptEnEff = 0.0,    delayEnEff = 0.0,    posEff = 0.0,    totEff = 0.0,    PoMean = 0.0, PoPSDMean = 0.0;
		rateErr = 0.0, NErr = 0.0, lifetimeErr = 0.0, promptPSDEffErr = 0.0, delayPSDEffErr = 0.0,  promptEnEffErr = 0.0, delayEnEffErr = 0.0, posEffErr = 0.0, totEffErr = 0.0, PoMeanErr = 0.0, PoPSDMeanErr = 0.0;

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
		fPromptEGaus = new TF1("fPromptEGaus","gaus",EMin,EMax);		
		hPromptE->Fit(fPromptEGaus,"RQ");

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
		//Calculate efficiencies and Ac227 rate

		//Prompt PSD Eff
		promptPSDEff = fPromptPSDGaus->Integral(promptPSDCutLow,promptPSDCutHigh)/fPromptPSDGaus->Integral(PSDMin,PSDMax);	
		promptPSDEffErr = sqrt((promptPSDEff*(1.0-promptPSDEff))/hPromptPSD->GetEntries());

		//DelayPSDEff
		delayPSDEff = fDelayPSDGaus->Integral(delayPSDCutLow,delayPSDCutHigh)/fDelayPSDGaus->Integral(PSDMin,PSDMax);
		delayPSDEffErr = sqrt((delayPSDEff*(1.0-delayPSDEff))/hDelayPSD->GetEntries());

		//Prompt Energy Eff
		promptEnEff = fPromptEGaus->Integral(promptECutLow,promptECutHigh)/fPromptEGaus->Integral(EMin,EMax);
		promptEnEffErr = sqrt((promptEnEff*(1.0-promptEnEff))/hPromptE->GetEntries());

		//Delay Energy Eff
		delayEnEff = fDelayEGaus->Integral(delayECutLow,delayECutHigh)/fDelayEGaus->Integral(EMin,EMax);
		delayEnEffErr = sqrt((delayEnEff*(1.0-delayEnEff))/hDelayE->GetEntries());

		//Position Eff
		posEff = fDzGaus->Integral(-posCut,posCut)/fDzGaus->Integral(dzMin,dzMax);
		posEffErr = sqrt((posEff*(1.0-posEff))/hRnPoDz->GetEntries());


		totEff = promptPSDEff*delayPSDEff*promptEnEff*delayEnEff*posEff;
		totEffErr = totEff*sqrt(pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(posEffErr/posEff,2)); 

		PoMean = fDelayEGaus->GetParameter(1);
		PoMeanErr = fDelayEGaus->GetParError(1);

		PoPSDMean = fDelayPSDGaus->GetParameter(1);
		PoPSDMeanErr = fDelayPSDGaus->GetParError(1);
	
		N = fDtExp->GetParameter(0);
		NErr = fDtExp->GetParError(0);
		lifetime = fDtExp->GetParameter(1);
		lifetimeErr = fDtExp->GetParError(1);

		rate = N/(livetime*totEff);
		rateErr = rate*sqrt(pow(NErr/N,2) + pow(promptPSDEffErr/promptPSDEff,2) + pow(delayPSDEffErr/delayPSDEff,2) + pow(promptEnEffErr/promptEnEff,2) + pow(delayEnEffErr/delayEnEff,2) + pow(posEffErr/posEff,2)); 

		//--------------------------------------------------------------------------
		//Draw histograms
		gStyle->SetOptStat(11);
		gStyle->SetOptFit(1111);

		TLegend *leg = nullptr;

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

/*
		cRnPoDt->SaveAs(Form("Cell%i_RnPoDt_%i.pdf",i,nLoop));	
		cRnPoPSD->SaveAs(Form("Cell%i_RnPoPSD_%i.pdf",i,nLoop));
		cRnPoE->SaveAs(Form("Cell%i_RnPoE_%i.pdf",i,nLoop));
		cRnPoDz->SaveAs(Form("Cell%i_RnPoDz_%i.pdf",i,nLoop));
*/

		delete cRnPoDt;
		delete cRnPoPSD;
		delete cRnPoE;
		delete cRnPoDz;
		}

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
		vPoMean.push_back(PoMean);
		vPoMeanErr.push_back(PoMeanErr);		
		vPoPSDMean.push_back(PoPSDMean);
		vPoPSDMeanErr.push_back(PoPSDMeanErr);		

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
	//fPromptEGaus = new TF1("fPromptEGaus","gaus",EMin,EMax);		
	fPromptEGaus = new TF1("fPromptEGaus","gaus",promptEFitMin,promptEFitMax);		
	hPromptE_AllCells->Fit(fPromptEGaus,"RQ");
	fPromptEGaus->SetRange(EMin,EMax);

	//fDelayEGaus = new TF1("fDelayEGaus","gaus",EMin,EMax);
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
	//Calculate efficiencies and Ac227 rate

	//Prompt PSD Eff
	promptPSDEffAllCells = fPromptPSDGaus->Integral(promptPSDCutLow,promptPSDCutHigh)/fPromptPSDGaus->Integral(PSDMin,PSDMax);	
	promptPSDEffErrAllCells = sqrt((promptPSDEffAllCells*(1.0-promptPSDEffAllCells))/hPromptPSD_AllCells->GetEntries());

	//DelayPSDEff
	delayPSDEffAllCells = fDelayPSDGaus->Integral(delayPSDCutLow,delayPSDCutHigh)/fDelayPSDGaus->Integral(PSDMin,PSDMax);
	delayPSDEffErrAllCells = sqrt((delayPSDEffAllCells*(1.0-delayPSDEffAllCells))/hDelayPSD_AllCells->GetEntries());

	//Prompt Energy Eff
	promptEnEffAllCells = fPromptEGaus->Integral(promptECutLow,promptECutHigh)/fPromptEGaus->Integral(EMin,EMax);
	promptEnEffErrAllCells = sqrt((promptEnEffAllCells*(1.0-promptEnEffAllCells))/hPromptE_AllCells->GetEntries());

	//Delay Energy Eff
	delayEnEffAllCells = fDelayEGaus->Integral(delayECutLow,delayECutHigh)/fDelayEGaus->Integral(EMin,EMax);
	delayEnEffErrAllCells = sqrt((delayEnEffAllCells*(1.0-delayEnEffAllCells))/hDelayE_AllCells->GetEntries());

	//Position Eff
	posEffAllCells = fDzGaus->Integral(-posCut,posCut)/fDzGaus->Integral(dzMin,dzMax);
	posEffErrAllCells = sqrt((posEffAllCells*(1.0-posEffAllCells))/hRnPoDz_AllCells->GetEntries());

	totEffAllCells = promptPSDEffAllCells*delayPSDEffAllCells*promptEnEffAllCells*delayEnEffAllCells*posEffAllCells;
	totEffErrAllCells = totEff*sqrt(pow(promptPSDEffErrAllCells/promptPSDEffAllCells,2) + pow(delayPSDEffErrAllCells/delayPSDEffAllCells,2) + pow(promptEnEffErrAllCells/promptEnEffAllCells,2) + pow(delayEnEffErrAllCells/delayEnEffAllCells,2) + pow(posEffErrAllCells/posEffAllCells,2)); 

	PoMeanAllCells = fDelayEGaus->GetParameter(1);
	PoMeanErrAllCells = fDelayEGaus->GetParError(1);

	PoPSDMeanAllCells = fDelayPSDGaus->GetParameter(1);
	PoPSDMeanErrAllCells = fDelayPSDGaus->GetParError(1);
	
	NAllCells = fDtExp->GetParameter(0);
	NErrAllCells = fDtExp->GetParError(0);
	lifetimeAllCells = fDtExp->GetParameter(1);
	lifetimeErrAllCells = fDtExp->GetParError(1);

	rateAllCells = NAllCells/(livetime*totEffAllCells);
	rateErrAllCells = rateAllCells*sqrt(pow(NErrAllCells/NAllCells,2) + pow(promptPSDEffErrAllCells/promptPSDEffAllCells,2) + pow(delayPSDEffErrAllCells/delayPSDEffAllCells,2) + pow(promptEnEffErrAllCells/promptEnEffAllCells,2) + pow(delayEnEffErrAllCells/delayEnEffAllCells,2) + pow(posEffErrAllCells/posEffAllCells,2)); 

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
	vPoMean.push_back(PoMeanAllCells);
	vPoMeanErr.push_back(PoMeanErrAllCells);		
	vPoPSDMean.push_back(PoPSDMeanAllCells);
	vPoPSDMeanErr.push_back(PoPSDMeanErrAllCells);		


	if(PLOTFLAG==1){

	//--------------------------------------------------------------------------
	//Plot histograms

	gStyle->SetOptStat(11);
	gStyle->SetOptFit(1111);

	TCanvas *cRnPoDtAll = new TCanvas("cRnPoDtAll","cRnPoDtAll",1);
	gPad->SetGrid();
	cRnPoDtAll->SetLogy();
	hSelectDt_AllCells->SetMinimum(1);
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
	//cRnPoPSDAll->SetLogy();
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
	fPromptEGaus->Draw("sames");
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


	cRnPoDtAll->SaveAs(Form("%s/AllCells_RnPoDt_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));	
	cRnPoPSDAll->SaveAs(Form("%s/AllCells_RnPoPSD_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoEAll->SaveAs(Form("%s/AllCells_RnPoE_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));
	cRnPoDzAll->SaveAs(Form("%s/AllCells_RnPoDz_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),nLoop));


	delete cRnPoDtAll;
	delete cRnPoPSDAll;
	delete cRnPoEAll;
	delete cRnPoDzAll;

	}

	//--------------------------------------------------------------------------
	printf("Histograms subtracted \n");
	printf("+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+ \n");
	printf("================================================================================================= \n");
	outFile->Write();

	vector<vector<double>> vvRate, vvN, vvLifetime, vvPSDEff, vvEnEff, vvPosEff, vvTotEff, vvPoMean, vvPoPSDMean;
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
	vvPoMean.push_back(vPoMean);
	vvPoMean.push_back(vPoMeanErr);
	vvPoPSDMean.push_back(vPoPSDMean);
	vvPoPSDMean.push_back(vPoPSDMeanErr);


	return make_tuple(IDX, tstamp, vvRate, vvN, vvLifetime, vvPSDEff, vvEnEff, vvPosEff, vvTotEff, vvPoMean, vvPoPSDMean);
}

void PlotResults(string TreeDir, const int NUMCELLS, int NUMTREES, int PLOTFLAG){

	int ExcludeCellArr[28] = {0,1,2,3,4,5,6,9,11,12,13,18,21,23,24,27,32,40,44,68,73,79,102,107,122,127,130,139};
	bool exclude;

	int IDX = 0;

	//--------------------------------------------------------------------------
	//Figure number of bins for histograms based on time of all trees and chosen time segment
	TChain AcChain("TAc");
	AcChain.Add(Form("%s/WetCommissioning/Series009_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series012_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series013_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series014_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series015_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series016_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series017_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series018_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series019_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series020_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series021_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series022_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/WetCommissioning/Series023_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Rampdown/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series000_AcTrees.root",getenv("P2X_ANALYZED")));
	AcChain.Add(Form("%s/180316_Background/Series001_AcTrees.root",getenv("P2X_ANALYZED")));

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

		if(i%10000==0) printf("Event: %d \n",i);
		double t = Ac_t[0]*(1e-6);	//ms
		lastTime = Ac_tstamp;

		if(t<prevTime){
			countNumTrees = countNumTrees + 1;
		}
	
		prevTime = t;
	}

	countNumTrees = countNumTrees + 1;
	
	int numBins = ceil((double)countNumTrees/(double)NUMTREES);

	printf("Number of trees: %d Number of bins: %d \n",countNumTrees,numBins);

	//--------------------------------------------------------------------------
	//Set up histograms	
	double min = 0.0, max = NUMCELLS;
	
	TH2F* hRate = new TH2F("hRate","Rate Per Cell;;Cell;Rate [mHz]",numBins,firstTime,lastTime,NUMCELLS,min,max);
	TH2F* hRateErr = new TH2F("hRateErr","Rate Error Per Cell;;Cell;Rate [mHz]",numBins,firstTime,lastTime,NUMCELLS,min,max);
	
	TH2F* hN = new TH2F("hN","N;;Cell;N",numBins,firstTime,lastTime,NUMCELLS,min,max);
	TH2F* hNErr = new TH2F("hNErr","N Error;;Cell;#deltaN",numBins,firstTime,lastTime,NUMCELLS,min,max);

	TH2F* hLifetime = new TH2F("hLifetime","Lifetime;;Cell;#tau [ms]",numBins,firstTime,lastTime,NUMCELLS,min,max);
	TH2F* hLifetimeErr = new TH2F("hLifetimeErr","Lifetime Error;;Cell;#delta#tau [ms]",numBins,firstTime,lastTime,NUMCELLS,min,max);

	TH2F* hTotEff = new TH2F("hTotEff","Total Efficiency;;Cell;Eff",numBins,firstTime,lastTime,NUMCELLS,min,max);
	TH2F* hTotEffErr = new TH2F("hTotEffErr","Error on Total Efficiency;;Cell;#deltaEff",numBins,firstTime,lastTime,NUMCELLS,min,max);

	TH2F* hPoMean = new TH2F("hPoMean","Po-215 Energy Mean;;Cell;Energy [MeVee]",numBins,firstTime,lastTime,NUMCELLS,min,max);
	TH2F* hPoPSDMean = new TH2F("hPoPSDMean","Po-215 PSD Mean;;Cell;Energy [MeVee]",numBins,firstTime,lastTime,NUMCELLS,min,max);

	TH2F* hRatePerCell = new TH2F("hRatePerCell","Rate Per Cell;;;Rate [mHz]",14,0,14,11,0,11);
	TH2F* hRelRatePerCell = new TH2F("hRelRatePerCell","Relative Rate Per Cell;;;Relative Rate",14,0,14,11,0,11);

	
	//--------------------------------------------------------------------------
	//Set up graphs
	double x[numBins], y[numBins];
	double xErr[numBins], yErr[numBins];

	TGraphErrors *grRate_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grN_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grLifetime_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grTotEff_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);
	TGraphErrors *grPoPSDMean_AllCells = new TGraphErrors(numBins,x,y,xErr,yErr);

	vector<TGraphErrors*> vRateGraphs;
	for(int i=0;i<NUMCELLS;i++){
		TGraphErrors *grRate_new = new TGraphErrors(numBins,x,y,xErr,yErr);
		vRateGraphs.push_back(grRate_new);
	}

	//--------------------------------------------------------------------------
	//Fill histograms
	vector<vector<double>> vRate, vN, vLifetime, vTotEff, vPoMean, vPoPSDMean;
	double tstamp;

	int n=1;
	while(IDX>=0){
		auto results = ReadAc227Trees(TreeDir, NUMCELLS, IDX, NUMTREES, n-1, PLOTFLAG);
		IDX 	  = get<0>(results);
		tstamp    = get<1>(results);
		vRate 	  = get<2>(results);
		vN        = get<3>(results);
		vLifetime = get<4>(results);
		vTotEff   = get<8>(results);
		vPoMean   = get<9>(results);	
		vPoPSDMean   = get<10>(results);	

		if(IDX<0) break;

		int relIdx = 140;
		double rateCellRel = vRate[0][relIdx];

		for(int i=0;i<NUMCELLS;i++){

			exclude = find(begin(ExcludeCellArr), end(ExcludeCellArr), i) != end(ExcludeCellArr);
			if(!exclude){
			
			int biny = i+1;

			hRate->SetBinContent(n,biny,vRate[0][i]*1e6);	//convert to mHz
			hRateErr->SetBinContent(n,biny,vRate[1][i]*1e6);

			hN->SetBinContent(n,biny,vN[0][i]);
			hNErr->SetBinContent(n,biny,vN[1][i]);

			hLifetime->SetBinContent(n,biny,vLifetime[0][i]);
			hLifetimeErr->SetBinContent(n,biny,vLifetime[1][i]);
		
			hTotEff->SetBinContent(n,biny,vTotEff[0][i]);
			hTotEffErr->SetBinContent(n,biny,vTotEff[1][i]);
		
			hPoMean->SetBinContent(n,biny,vPoMean[0][i]);			
			hPoPSDMean->SetBinContent(n,biny,vPoPSDMean[0][i]);			
		
			int binx = i%14 + 1;
			biny = (i/14) + 1;
			hRatePerCell->SetBinContent(binx,biny,vRate[0][i]*1e6);

			double relRate = vRate[0][i]/rateCellRel;
			hRelRatePerCell->SetBinContent(binx,biny,relRate);

			int pt = n-1;
			vRateGraphs[i]->SetPoint(pt,tstamp,vRate[0][i]*1e6);
			vRateGraphs[i]->SetPointError(pt,0,vRate[1][i]*1e6);		

			}	
		}

		gStyle->SetPaintTextFormat("2.2f");
		gStyle->SetOptStat(0);
		TCanvas *cRatePerCell = new TCanvas("cRatePerCell","Rate per Cell",1);
		gPad->SetRightMargin(0.16);
		hRatePerCell->Draw("colz && text");
		cRatePerCell->SaveAs(Form("%s/RatePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));
	
		TCanvas *cRelRatePerCell = new TCanvas("cRelRatePerCell","Relative Rate per Cell",1);
		gPad->SetRightMargin(0.16);
		hRelRatePerCell->Draw("colz && text");
		cRelRatePerCell->SaveAs(Form("%s/RelRatePerCell_%i.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS"),n));

		delete cRatePerCell;
		delete cRelRatePerCell;

		int pt = n-1;
		grRate_AllCells->SetPoint(pt,tstamp,vRate[0][NUMCELLS]*1000);	//convert to Hz
		grRate_AllCells->SetPointError(pt,0,vRate[1][NUMCELLS]*1000);

		grN_AllCells->SetPoint(pt,tstamp,vN[0][NUMCELLS]);
		grN_AllCells->SetPointError(pt,0,vN[1][NUMCELLS]);
	
		grLifetime_AllCells->SetPoint(pt,tstamp,vLifetime[0][NUMCELLS]);
		grLifetime_AllCells->SetPointError(pt,0,vLifetime[1][NUMCELLS]);

		grTotEff_AllCells->SetPoint(pt,tstamp,vTotEff[0][NUMCELLS]);
		grTotEff_AllCells->SetPointError(pt,0,vTotEff[1][NUMCELLS]);	
	
		grPoMean_AllCells->SetPoint(pt,tstamp,vPoMean[0][NUMCELLS]);
		grPoMean_AllCells->SetPointError(pt,0,vPoMean[1][NUMCELLS]);	

		grPoPSDMean_AllCells->SetPoint(pt,tstamp,vPoPSDMean[0][NUMCELLS]);
		grPoPSDMean_AllCells->SetPointError(pt,0,vPoPSDMean[1][NUMCELLS]);	

		n = n + 1;
	}

	//--------------------------------------------------------------------------
	//Plot Histograms
	gStyle->SetOptStat(0);
	gStyle->SetPalette(kTemperatureMap);

	int numCellDivisions = (NUMCELLS/11) + 100*2 + 10000*0;
	//--------------------------------------------------------------------------
	//Plot Graphs

	TCanvas *cRateAllCells = new TCanvas("cRateAllCells","Rate All Cells",1);
	grRate_AllCells->SetTitle("Total Rate of Detector");
	grRate_AllCells->GetYaxis()->SetTitle("Rate (Hz)");
	grRate_AllCells->GetYaxis()->SetTitleOffset(1.3);
	grRate_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grRate_AllCells->SetMarkerStyle(20);
	grRate_AllCells->SetMarkerColor(kBlue);
	grRate_AllCells->SetLineColor(kBlue);
	grRate_AllCells->Draw("AP");
	grRate_AllCells->GetXaxis()->SetTimeDisplay(1);
	grRate_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cRateAllCells->SaveAs(Form("%s/Rate_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cNAllCells = new TCanvas("cNAllCells","N All Cells",1);
	grN_AllCells->SetTitle("Total N of Detector");
	grN_AllCells->GetYaxis()->SetTitle("N");
	grN_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grN_AllCells->SetMarkerStyle(20);
	grN_AllCells->SetMarkerColor(kBlue);
	grN_AllCells->SetLineColor(kBlue);
	grN_AllCells->Draw("AP");
	grN_AllCells->GetXaxis()->SetTimeDisplay(1);
	grN_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cNAllCells->SaveAs(Form("%s/N_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cLifetimeAllCells = new TCanvas("cLifetimeAllCells","Lifetime All Cells",1);
	grLifetime_AllCells->SetTitle("Total Po-215 Lifetime of Detector");
	grLifetime_AllCells->GetYaxis()->SetTitle("Lifetime (ms)");
	grLifetime_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grLifetime_AllCells->SetMarkerStyle(20);
	grLifetime_AllCells->SetMarkerColor(kBlue);
	grLifetime_AllCells->SetLineColor(kBlue);
	grLifetime_AllCells->Draw("AP");
	grLifetime_AllCells->GetXaxis()->SetTimeDisplay(1);
	grLifetime_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cLifetimeAllCells->SaveAs(Form("%s/Lifetime_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cTotEffAllCells = new TCanvas("cTotEffAllCells","Tot Eff All Cells",1);
	grTotEff_AllCells->SetTitle("Total Efficiency of Detector");
	grTotEff_AllCells->GetYaxis()->SetTitle("Eff");
	grTotEff_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grTotEff_AllCells->SetMarkerStyle(20);
	grTotEff_AllCells->SetMarkerColor(kBlue);
	grTotEff_AllCells->SetLineColor(kBlue);
	grTotEff_AllCells->Draw("AP");
	grTotEff_AllCells->GetXaxis()->SetTimeDisplay(1);
	grTotEff_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cTotEffAllCells->SaveAs(Form("%s/TotEff_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoMeanAllCells = new TCanvas("cPoMeanAllCells","Po Mean All Cells",1);
	grPoMean_AllCells->SetTitle("Total Po-215 Mean of Detector");
	grPoMean_AllCells->GetYaxis()->SetTitle("N");
	grPoMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoMean_AllCells->SetMarkerStyle(20);
	grPoMean_AllCells->SetMarkerColor(kBlue);
	grPoMean_AllCells->SetLineColor(kBlue);
	grPoMean_AllCells->Draw("AP");
	grPoMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoMean_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cPoMeanAllCells->SaveAs(Form("%s/PoMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	TCanvas *cPoPSDMeanAllCells = new TCanvas("cPoPSDMeanAllCells","Po PSD Mean All Cells",1);
	grPoPSDMean_AllCells->SetTitle("Total Po-215 PSD Mean of Detector");
	grPoPSDMean_AllCells->GetYaxis()->SetTitle("N");
	grPoPSDMean_AllCells->GetXaxis()->SetLabelOffset(0.028);
	grPoPSDMean_AllCells->SetMarkerStyle(20);
	grPoPSDMean_AllCells->SetMarkerColor(kBlue);
	grPoPSDMean_AllCells->SetLineColor(kBlue);
	grPoPSDMean_AllCells->Draw("AP");
	grPoPSDMean_AllCells->GetXaxis()->SetTimeDisplay(1);
	grPoPSDMean_AllCells->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	cPoPSDMeanAllCells->SaveAs(Form("%s/PoPSDMean_AllCells.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));

	//--------------------------------------------------------------------------
	//Plot rate per cell by row and column
		
	vector<int> vColor;
	vColor.push_back(1);
	vColor.push_back(2);
	vColor.push_back(3);
	vColor.push_back(4);
	vColor.push_back(41);
	vColor.push_back(6);
	vColor.push_back(7);
	vColor.push_back(1);
	vColor.push_back(2);
	vColor.push_back(3);
	vColor.push_back(4);
	vColor.push_back(41);
	vColor.push_back(6);
	vColor.push_back(7);

	vector<int> vMarker;
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(20);
	vMarker.push_back(22);
	vMarker.push_back(22);
	vMarker.push_back(22);
	vMarker.push_back(22);
	vMarker.push_back(22);
	vMarker.push_back(22);
	vMarker.push_back(22);

	vector<int> vLine;
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(1);
	vLine.push_back(9);
	vLine.push_back(9);
	vLine.push_back(9);
	vLine.push_back(9);
	vLine.push_back(9);
	vLine.push_back(9);
	vLine.push_back(9);

	double minRate = 2, maxRate = 3.5;
	int cellIdx = 0;
	double markerSize = 2.0;

	TCanvas *cRateRow0_1 = new TCanvas("cRateRow0_1","Rate Row 0-1",1800,1800);
	cRateRow0_1->Divide(1,2);

	cRateRow0_1->cd(2);
	gPad->SetRightMargin(0.2);
	int styleIdx = 0;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 0;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	TLegend *leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	cRateRow0_1->cd(1);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 1;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();


	TCanvas *cRateRow2_3 = new TCanvas("cRateRow2_3","Rate Row 2-3",1800,1800);
	cRateRow2_3->Divide(1,2);

	cRateRow2_3->cd(2);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 2;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	cRateRow2_3->cd(1);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 3;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	TCanvas *cRateRow4_5 = new TCanvas("cRateRow4_5","Rate Row 4-5",1800,1800);
	cRateRow4_5->Divide(1,2);

	cRateRow4_5->cd(2);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 4;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	cRateRow4_5->cd(1);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 5;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	TCanvas *cRateRow6_7 = new TCanvas("cRateRow6_7","Rate Row 6-7",1800,1800);
	cRateRow6_7->Divide(1,2);

	cRateRow6_7->cd(2);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 6;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	cRateRow6_7->cd(1);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 7;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	TCanvas *cRateRow8_9 = new TCanvas("cRateRow8_9","Rate Row 8-9",1800,1800);
	cRateRow8_9->Divide(1,2);

	cRateRow8_9->cd(2);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 8;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	cRateRow8_9->cd(1);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 9;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	TCanvas *cRateRow10 = new TCanvas("cRateRow10","Rate Row 10",1800,900);
	cellIdx = cellIdx + 14;
	styleIdx = 0;
	gPad->SetRightMargin(0.2);
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Row 10;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(2);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;

	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	for(int i=cellIdx+1;i<cellIdx+14;i++){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(2);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateRow0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
	}
	leg->Draw();

	//--------------------------------------------------------------------------
	cellIdx = 0;

	TCanvas *cRateCol0_1 = new TCanvas("cRateCol0_1","Rate Column 0-1",3600,1200);
	cRateCol0_1->Divide(2,1);

	cRateCol0_1->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 0;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	int i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol0_1->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 1;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol0_1->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol2_3 = new TCanvas("cRateCol2_3","Rate Column 2-3",3600,1200);
	cRateCol2_3->Divide(2,1);

	cRateCol2_3->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 2;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol2_3->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol2_3->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 3;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol2_3->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol4_5 = new TCanvas("cRateCol4_5","Rate Column 4-5",3600,1200);
	cRateCol4_5->Divide(2,1);

	cRateCol4_5->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 4;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol4_5->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol4_5->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 5;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol4_5->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol6_7 = new TCanvas("cRateCol6_7","Rate Column 6-7",3600,1200);
	cRateCol6_7->Divide(2,1);

	cRateCol6_7->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 6;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol6_7->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol6_7->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 7;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol6_7->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol8_9 = new TCanvas("cRateCol8_9","Rate Column 8-9",3600,1200);
	cRateCol8_9->Divide(2,1);

	cRateCol8_9->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 8;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol8_9->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol8_9->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 9;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol8_9->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol10_11 = new TCanvas("cRateCol10_11","Rate Column 10-11",3600,1200);
	cRateCol10_11->Divide(2,1);

	cRateCol10_11->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 10;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol10_11->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol10_11->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 11;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol10_11->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	TCanvas *cRateCol12_13 = new TCanvas("cRateCol12_13","Rate Column 12-13",3600,1200);
	cRateCol12_13->Divide(2,1);

	cRateCol12_13->cd(1);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 12;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol12_13->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();
	
	cRateCol12_13->cd(2);
	gPad->SetRightMargin(0.2);
	styleIdx = 1;
	cellIdx++;
	vRateGraphs[cellIdx]->SetTitle("Rate Per Cell: Column 13;;Rate (mHz)");
	vRateGraphs[cellIdx]->SetMarkerStyle(vMarker[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetMarkerSize(markerSize);
	vRateGraphs[cellIdx]->SetLineColor(vColor[styleIdx]);
	vRateGraphs[cellIdx]->SetLineWidth(1);
	vRateGraphs[cellIdx]->SetLineStyle(vLine[styleIdx]);
	vRateGraphs[cellIdx]->Draw("ALP");
	vRateGraphs[cellIdx]->GetYaxis()->SetRangeUser(minRate,maxRate);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeDisplay(1);
	vRateGraphs[cellIdx]->GetXaxis()->SetTimeFormat("#splitline{%m/%d}{%H:%M}");
	vRateGraphs[cellIdx]->GetXaxis()->SetLabelOffset(0.028);
	styleIdx++;
		
	leg = new TLegend(0.85,0.1,1.0,0.9);
	leg->AddEntry(vRateGraphs[cellIdx],Form("Cell %d",cellIdx),"pl");
	i = cellIdx + 14;
	while(i<154){
		vRateGraphs[i]->SetMarkerStyle(vMarker[styleIdx]);
		vRateGraphs[i]->SetMarkerColor(vColor[styleIdx]);	
		vRateGraphs[i]->SetMarkerSize(markerSize);	
		vRateGraphs[i]->SetLineColor(vColor[styleIdx]);
		vRateGraphs[i]->SetLineWidth(1);
		vRateGraphs[i]->SetLineStyle(vLine[styleIdx]);
		vRateGraphs[i]->Draw("LP");
		cRateCol12_13->Update();
		styleIdx++;
		leg->AddEntry(vRateGraphs[i],Form("Cell %d",i),"pl");
		i = i + 14;
	}
	leg->Draw();

	//=========================================================================================================
	if(PLOTFLAG==1){
	cRateRow0_1->SaveAs(Form("%s/RatePerCell_Row0and1.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateRow2_3->SaveAs(Form("%s/RatePerCell_Row2and3.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateRow4_5->SaveAs(Form("%s/RatePerCell_Row4and5.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateRow6_7->SaveAs(Form("%s/RatePerCell_Row6and7.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateRow8_9->SaveAs(Form("%s/RatePerCell_Row8and9.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateRow10->SaveAs(Form("%s/RatePerCell_Row10.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol0_1->SaveAs(Form("%s/RatePerCell_Col0and1.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol2_3->SaveAs(Form("%s/RatePerCell_Col2and3.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol4_5->SaveAs(Form("%s/RatePerCell_Col4and5.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol6_7->SaveAs(Form("%s/RatePerCell_Col6and7.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol8_9->SaveAs(Form("%s/RatePerCell_Col8and9.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol10_11->SaveAs(Form("%s/RatePerCell_Col10and11.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	cRateCol12_13->SaveAs(Form("%s/RatePerCell_Col12and13.pdf",getenv("AD_AC227ANALYSIS_DATA_PLOTS")));
	}
}
