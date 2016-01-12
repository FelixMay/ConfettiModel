//---------------------------------------------------------------------------
#include "Forest.h"
#include <time.h>
#include <string>
using namespace std;

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
	int NRep = 1;  //number replicates
	int NGen = 100;  //number generations (# complete turnover of community)
	bool StepsOut = false;
	bool R_Mode = false;

	string forest = "BCI";
	//string forest = "Sin";
	string SimFileName = "Para";
	string FileLabel = "1";

	if (argc == 2){
		SimFileName = argv[1];
	}

	if (argc == 3){
		SimFileName = argv[1];
		FileLabel = argv[2];
	}

	if (argc == 4){
		SimFileName = argv[1];
		FileLabel = argv[2];
		NGen = StringToInt(argv[3]);
	}

	if (argc == 5){
		SimFileName = argv[1];
		FileLabel = argv[2];
		NGen = StringToInt(argv[3]);
		NRep = StringToInt(argv[4]);
	}

	if (argc == 6){
		SimFileName = argv[1];
		FileLabel = argv[2];
		NGen = StringToInt(argv[3]);
		NRep = StringToInt(argv[4]);
		forest = argv[5];
	}

	SimFileName = "InOut\\"+SimFileName + FileLabel +".txt";
	//SimFileName = "InOut\\"+SimFileName +".txt";

	ifstream InFile;
	InFile.open(SimFileName.c_str());

	cout<<"Input-File:\t"<<SimFileName<<endl;
	cout<<"Sim-Label:\t"<<FileLabel<<endl;
	cout<<"Generations:\t"<<NGen<<endl;
	cout<<"Replicates:\t"<<NRep<<endl;
	cout<<"Forest plot:\t"<<forest<<endl<<endl;

	time_t start, end;

	start = time(0);

	string line1;
	getline(InFile,line1);

	CPara* pPara = new CPara();

	//int seed = (int) start;
	int seed = 99;

	double xmax, ymax, map_cell_size;
	char* map_file_name;
	char* rel_dens_file_name;
	int n_hab_types;

	if (forest == "BCI"){
		xmax = 1000;
		ymax = 500;
		map_cell_size = 20.0;
      //map_file_name = "InOut\\HabitatMapBCI_Harms2001.txt";
		map_file_name = "InOut\\HabitatMapBCI_Harms_nomixed.txt";
		//map_file_name = "InOut\\HabitatMapBCI_fake.txt";

		rel_dens_file_name = "InOut\\RelativeDensityBCI_harms_nomixed_n50.txt";
		n_hab_types = 6;
	}
	else if (forest == "Sin") {
		xmax = 500;
		ymax = 500;
		map_cell_size = 500.0/26.0;
		map_file_name = "InOut\\HabitatMapSinharaja_Ruwan1.txt";
		rel_dens_file_name = "InOut\\RelativeDensitySinharaja_Ruwan_n50.txt";
		n_hab_types = 5;
	}
	else cout<<"Error Forest Name"<<endl;

	CForest* pForest = new CForest(seed, xmax, ymax,
								   map_cell_size,
								   map_file_name,
								   rel_dens_file_name,
								   n_hab_types);
	pForest->FileOpen(FileLabel);

	int isim;

	if (InFile.good()) {

		//read first parameter set
		InFile>>isim;
		InFile>>pPara->NTrees;
		InFile>>pPara->Jmeta;
		InFile>>pPara->theta;
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

		//pPara->GetMuSigma();

		while (InFile.good()) {

			cout<<"Sim "<<isim<<endl;
			pForest->Pars = pPara;

			//run simulations
			for (int irep=1; irep <= NRep; ++irep) {

				cout<<"  Rep "<<irep<<endl;
				pForest->OneRun(isim, irep, NGen, StepsOut, R_Mode);
				pForest->ClearForest();
			}  // end irep

			//try to read new parameter set
			InFile>>isim;
			InFile>>pPara->NTrees;
			InFile>>pPara->Jmeta;
			InFile>>pPara->theta;
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

			//pPara->GetMuSigma();
		}
	}
	else cout<<"Error SimFile"<<endl;

	delete pForest;
	delete pPara;

	end = time(0);

	cout<<"\nRuntime: "<<end - start<<" seconds"<<endl;

	cin.ignore();

	return 0;
}
//---------------------------------------------------------------------------
