//---------------------------------------------------------------------------


#include <R.h>
#include "Forest.h"
#include <fstream>
using std::string;

extern "C" {

void PredFDP(int* seed,
             char** settings_file,
             //parameters
             double* theta,
			    double* m,
			    double* rmax,
			    double* a_rec,
			    double* a_hab,
			    double* a_surv,
			    double* b_surv,
			    double* m_dm_spec,
			    double* sd_dm_spec,
			    double* m_JC,
			    double* sd_JC,
			    //summary statistics,
			    int* NSpec,
			    double* Shannon,
			    int* Abund,
			    int* nClassSAD,
			    int* SAD,
			    double* AnnualMort,
			    int* BD5,
			    int* BDtotal,
			    double* PCF,
			    double* PropCon,
			    double* SARq,
			    //species specific patterns
			    //double* MetaRelAbund,
			    //double* muDisp,
			    //double* sigmaDisp,
			    //double* JCfac,
			    double* Kcon10,
			    double* Kcon50,
             //double* mNNdist,
			    double* Khet10,
			    double* Khet50
			    // tree specific measures
			    //double* gx,
			    //double* gy,
			    //int* sp
			    )
{

	CModelSettings* pSettings = new CModelSettings();

	//string SettingsFileName = "Input\\Settings_R1.txt";
	string SettingsFileName(settings_file[0]);

//	std::ofstream Test;
//	Test.open("Test.txt");
//	Test<<SettingsFileName;
//	Test<<settings_file[0]<<"\n";
// Test.close();

   pSettings->ReadSettings(SettingsFileName);

	int metaSR = 0;
	double metaCV = 0.0;

	//create parameter set
	CPara* pPara = new CPara(theta[0],
                            metaSR,
                            metaCV,
                            m[0],
                            rmax[0],
                            a_rec[0],
                            a_hab[0],
                            a_surv[0],
                            b_surv[0],
                            m_dm_spec[0],
                            sd_dm_spec[0],
                            m_JC[0],
                            sd_JC[0]
                           );

   CForest* pForest = new CForest(seed[0], pSettings);

   pForest->pPars = pPara;

	//may override values in Settings file
	pSettings->steps_out = false;
	pSettings->R_mode = true;

	pForest->FileOpen("1");

	pForest->OneRun(1,1);

	NSpec[0] = pForest->SpecAbund.size();
	Shannon[0] = pForest->GetShannon();
	BD5[0] = pForest->BD_5years;
	BDtotal[0] = pForest->BD_total;

	//Annual mortality assuming that one model time step equals five years
	AnnualMort[0] = (log(pForest->NTrees) - log(pForest->NTrees - pForest->BD_5years))/5.0 * 100.0;

	pForest->GetSAD();
	for (int iclass = 0; iclass < nClassSAD[0]; ++iclass)
		SAD[iclass] = pForest->SAD[iclass];

	pForest->GetPPA();
	pForest->GetSARq();

	//pForest->WriteOutput(1,1,1);

	std::map<int,int>::iterator spec_it1;
	std::map<int,int>::iterator spec_it2;

	int i1 = 0, i2 = 0, i3 = 0;

	int i10m = floor(10.0/pSettings->bw2) - 1;
	int i50m = floor(50.0/pSettings->bw2) - 1;

	for (spec_it1 = pForest->SpecAbund.begin(); spec_it1!=pForest->SpecAbund.end(); ++spec_it1){

		Abund[i1] = spec_it1->second;
		++i1;

		if (spec_it1->second >= pSettings->minAbund){

			Kcon10[i2] = pForest->KSpecIJ[spec_it1->first][spec_it1->first][i10m];
			Kcon50[i2] = pForest->KSpecIJ[spec_it1->first][spec_it1->first][i50m];
			++i2;

			for (spec_it2 = pForest->SpecAbund.begin(); spec_it2!=pForest->SpecAbund.end(); ++spec_it2){
				if (spec_it2->second >= pSettings->minAbund){
					if (spec_it1->first != spec_it2->first){

						Khet10[i3] = pForest->KSpecIJ[spec_it1->first][spec_it2->first][i10m];
						Khet50[i3] = pForest->KSpecIJ[spec_it1->first][spec_it2->first][i50m];

						//Dhet10[i3] = pForest->DSpecIJ[spec_it1->first][spec_it2->first][1];
						//Dhet50[i3] = pForest->DSpecIJ[spec_it1->first][spec_it2->first][9];

						//AbundSpec2[i3] = spec_it2->second;

						//xPOD50[i3] = pForest->xPOD[spec_it1->first][spec_it2->first];
						//NNDist[i3] = pForest->NNDistSpecIJ[spec_it1->first][spec_it2->first];

						++i3;
					}
				}
			}
		}
	}

	for (int ir=0; ir<pForest->nBins1; ++ir){
		PCF[ir] = pForest->PCF_all[ir];
		PropCon[ir] = pForest->PropCon[ir];
		//SAR1[ir] = pForest->SAR[ir];
	}

	for (int ir=0; ir<pForest->SARq_n; ++ir)
		SARq[ir] = pForest->SARq_m[ir];


	pForest->ClearForest();

   delete pForest;
   delete pPara;
   delete pSettings;
}

}  //end extern "C"

//------------------------------------------------------------------------------
