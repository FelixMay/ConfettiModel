//---------------------------------------------------------------------------
#include "Forest.h"
#include <time.h>
#include <string>
using std::string;

#include <fstream>
using std::ifstream;

#include <iostream>
using std::cout;
using std::cin;
using std::endl;

//---------------------------------------------------------------------------
int StringToInt(std::string S)
{
	int result;
	std::istringstream is(S);
	is>>result;
	return result;
}

int main(int argc, char* argv[])
{
	string SettingsFileName = "Settings1";

	//uniform initialization in C++11
	string ParaFileName{"Para"};
	//string ParaFileName{"ParaGS1"};

	string InFileID{"1"};
	string OutFileID{"1"};

	if (argc == 2){
		ParaFileName = argv[1];
	}

	if (argc == 3){
      SettingsFileName = argv[1];
		ParaFileName = argv[2];
	}

	if (argc == 4){
      SettingsFileName = argv[1];
		ParaFileName = argv[2];
		OutFileID = InFileID = argv[3];
	}

	if (argc == 5){
      SettingsFileName = argv[1];
		ParaFileName = argv[2];
		InFileID = argv[3];
		OutFileID = argv[4];
	}

   SettingsFileName = "Input/" + SettingsFileName + ".txt";
   CModelSettings* pSettings = new CModelSettings();
	pSettings->ReadSettings(SettingsFileName);

	ParaFileName = "Input/"+ParaFileName + InFileID +".csv";

	ifstream InFile;
	InFile.open(ParaFileName.c_str());

	cout<<"Settings-File:\t"<<SettingsFileName<<endl;
	cout<<"Parameter-File:\t"<<ParaFileName<<endl;

	cout<<"Generations:\t"<<pSettings->nGen<<endl;
	cout<<"Replicates:\t"<<pSettings->nRep<<endl;

	time_t start, end;

	start = time(0);

	string line1;
	getline(InFile,line1);

	CPara* pPara = new CPara();

	int seed = start + StringToInt(OutFileID);
   //int seed = 99;

	CForest* pForest = new CForest(seed, pSettings);
	pForest->FileOpen(OutFileID);

	int isim;
	double jm;
	char delim;

	if (InFile.good()) {

		//read first parameter set
		InFile>>isim>>delim;
		InFile>>pPara->theta>>delim;
		InFile>>jm>>delim;
      pPara->Jm = static_cast<int>(jm);
		InFile>>pPara->metaSR>>delim;
		InFile>>pPara->metaCV>>delim;
		InFile>>pPara->m>>delim;
		InFile>>pPara->r_max>>delim;
		InFile>>pPara->aRec>>delim;
		//InFile>>pPara->aSurv>>delim;
		//InFile>>pPara->bSurv>>delim;
		InFile>>pPara->m_dm_spec>>delim;
		InFile>>pPara->sd_dm_spec>>delim;
		InFile>>pPara->m_JCspec>>delim;
		InFile>>pPara->cv_JCspec>>delim;
		InFile>>pPara->niche_breadth>>delim;
		InFile>>pPara->n_hills;

		while (InFile.good()) {

			cout<<"Sim "<<isim<<endl;
			pForest->pPars = pPara;

			//pForest->initSpecies();
			//pForest->writeInteractMat(isim);
			//pForest->writeSpecies(isim);

			//run simulations
			for (int irep=1; irep <= pSettings->nRep; ++irep) {

				cout<<"  Rep "<< irep <<endl;

            pForest->initSpecies();
            //pForest->writeInteractMat(isim, irep);

            pForest->initTrees();

				pForest->OneRun(isim, irep);

				//pForest->GetPPA();
				//pForest->writeSpecies(isim, irep);

				pForest->clearTrees();
				pForest->clearSpecies();
			}  // end irep

			//pForest->clearSpecies();

			//try to read new parameter set
			InFile>>isim>>delim;
			InFile>>pPara->theta>>delim;
			InFile>>jm>>delim;
         pPara->Jm = static_cast<int>(jm);
			InFile>>pPara->metaSR>>delim;
         InFile>>pPara->metaCV>>delim;
			InFile>>pPara->m>>delim;
			InFile>>pPara->r_max>>delim;
			InFile>>pPara->aRec>>delim;
			//InFile>>pPara->aSurv>>delim;
			//InFile>>pPara->bSurv>>delim;
			InFile>>pPara->m_dm_spec>>delim;
			InFile>>pPara->sd_dm_spec>>delim;
			InFile>>pPara->m_JCspec>>delim;
			InFile>>pPara->cv_JCspec>>delim;
			InFile>>pPara->niche_breadth>>delim;
			InFile>>pPara->n_hills;

		}
	}
	else cout<<"Error SimFile"<<endl;

	delete pForest;
	delete pPara;
	delete pSettings;

	end = time(0);

	cout<<"\nRuntime: "<<end - start<<" seconds"<<endl;

	cin.ignore();

	return 0;
}
//---------------------------------------------------------------------------
