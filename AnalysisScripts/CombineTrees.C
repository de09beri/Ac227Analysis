//Combine all trees of a series into one tree file

#include <fstream>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TExec.h"
#include <sstream>
#include "TChain.h"
#include "TString.h"

using namespace std;

void CombineTrees(int DataType){

	TChain AcChain("Ac227TreePlugin/TAc");

	string line;
	ifstream file(Form("%s/filePath.txt",getenv("AD_AC227ANALYSIS_SCRIPTS")));

	int numline = 0;
	int series = 0;

	if(file.is_open()){
		while(getline(file,line)){
			cout<<line<<endl;

			if(numline==0){
				series = atoi(line.c_str());
			}
			else{
				AcChain.Add(line.c_str());
			}

			numline++;
		}

		file.close();
	}

	TTree *T = AcChain.CloneTree();

	if(DataType==0){
		TFile *fout = TFile::Open(Form("%s/WetCommissioning/Series%03d_AcTrees.root",getenv("P2X_ANALYZED"),series),"RECREATE");
		T->Write();
		fout->Close();
	}
	if(DataType==1){ 
		TFile *fout = TFile::Open(Form("%s/180316_Rampdown/Series%03d_AcTrees.root",getenv("P2X_ANALYZED"),series),"RECREATE");
		T->Write();
		fout->Close();
	}
	if(DataType==2){ 
		TFile *fout = TFile::Open(Form("%s/180316_Background/Series%03d_AcTrees.root",getenv("P2X_ANALYZED"),series),"RECREATE");
		T->Write();
		fout->Close();

	}
}
