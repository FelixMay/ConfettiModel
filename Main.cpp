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
	string SettingsFileName = "Input\\Settings1.txt";

	//string ParaFileName = "Para";
	//string FileLabel = "1";

	//uniform initialization in C++11
	string ParaFileName{"Para"};
	string FileLabel{"1"};

	if (argc == 2){
		ParaFileName = argv[1];
	}

	if (argc == 3){
		ParaFileName = argv[1];
		FileLabel = argv[2];
	}

	CModelSettings* pSettings = new CModelSettings();

	pSettings->ReadSettings(SettingsFileName);

	ParaFileName = "Input\\"+ParaFileName + FileLabel +".txt";

	ifstream InFile;
	InFile.open(ParaFileName.c_str());

	cout<<"Settings-File:\t"<<SettingsFileName<<endl;
	cout<<"Parameter-File:\t"<<ParaFileName<<endl;
	cout<<"Sim-Label:\t"<<FileLabel<<endl;
	cout<<"Generations:\t"<<pSettings->nGen<<endl;
	cout<<"Replicates:\t"<<pSettings->nRep<<endl;


	time_t start, end;

	start = time(0);

	string line1;
	getline(InFile,line1);

	CPara* pPara = new CPara();

	int seed = (int) start;
   //int seed = 99;

//	else if (forest == "Sin") {
//		xmax = 500;
//		ymax = 500;
//		map_cell_size = 500.0/26.0;
//		map_file_name = "InOut\\HabitatMapSinharaja_Ruwan1.txt";
//		rel_dens_file_name = "InOut\\RelativeDensitySinharaja_Ruwan_n50.txt";
//		n_hab_types = 5;
//	}
//	else cout<<"Error Forest Name"<<endl;

	CForest* pForest = new CForest(seed, pSettings);
	pForest->FileOpen(FileLabel);

	int isim;

	if (InFile.good()) {

		//read first parameter set
		InFile>>isim;
		InFile>>pPara->theta;
		InFile>>pPara->metaSR;
		InFile>>pPara->metaCV;
		InFile>>pPara->m;
		InFile>>pPara->r_max;
		InFile>>pPara->aRec;
		InFile>>pPara->aHab;
		InFile>>pPara->aSurv;
		InFile>>pPara->bSurv;
		InFile>>pPara->m_dm_spec;
		InFile>>pPara->sd_dm_spec;
		InFile>>pPara->m_JCspec;
		InFile>>pPara->sd_JCspec;

		while (InFile.good()) {

			cout<<"Sim "<<isim<<endl;
			pForest->pPars = pPara;

			//run simulations
			for (int irep=1; irep <= pSettings->nRep; ++irep) {

				cout<<"  Rep "<<irep<<endl;
				pForest->OneRun(isim, irep);
				pForest->ClearForest();
			}  // end irep

			//try to read new parameter set
			InFile>>isim;
			InFile>>pPara->theta;
			InFile>>pPara->metaSR;
         InFile>>pPara->metaCV;
			InFile>>pPara->m;
			InFile>>pPara->r_max;
			InFile>>pPara->aRec;
			InFile>>pPara->aHab;
			InFile>>pPara->aSurv;
			InFile>>pPara->bSurv;
			InFile>>pPara->m_dm_spec;
			InFile>>pPara->sd_dm_spec;
			InFile>>pPara->m_JCspec;
			InFile>>pPara->sd_JCspec;
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
