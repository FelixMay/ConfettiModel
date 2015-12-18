//---------------------------------------------------------------------------


#include <R.h>
#include "Forest.h"
using namespace std;

extern "C" {

void PredFDP(int* nGen,
			    int* seed,
             //model settings
			    double* xmax,
			    double* ymax,
			    int* ntrees,
			    int* jmeta,
			    double* map_cell_size,
			    char** map_file,
			    char** rel_dens_file,
			    int* n_hab_types,
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
			    double* SAR1,
			    double* SAR2_m,
			    double* SAR2_sd,
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
	//create parameter set
	CPara* pPara = new CPara(ntrees[0],
							 jmeta[0],
							 theta[0],
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

	CForest* pForest = new CForest(seed[0],
								   xmax[0],ymax[0],
								   map_cell_size[0],
								   map_file[0],
								   rel_dens_file[0],
								   n_hab_types[0]);

	pForest->Pars = pPara;

	bool StepsOut = false;
	bool R_Mode = true;

	//pForest->FileOpen("1");

	pForest->OneRun(1,1,nGen[0],StepsOut,R_Mode);

	NSpec[0] = pForest->SpecAbund.size();
	Shannon[0] = pForest->GetShannon();
	BD5[0] = pForest->BD_5years;
	BDtotal[0] = pForest->BD_total;
	AnnualMort[0] = (log(pForest->NTrees) - log(pForest->NTrees - pForest->BD_5years))/5.0 * 100.0;

	pForest->GetSAD();
	for (int iclass = 0; iclass < nClassSAD[0]; ++iclass)
		SAD[iclass] = pForest->SAD[iclass];


	pForest->GetPPA();
	pForest->GetSAR2();

	//pForest->WriteOutput(1,1,1);

	map<int,int>::iterator spec_it1;
	map<int,int>::iterator spec_it2;

	int i1 = 0, i2 = 0, i3 = 0;
	int minAbund = pForest->minAbund;

	for (spec_it1 = pForest->SpecAbund.begin(); spec_it1!=pForest->SpecAbund.end(); ++spec_it1){

		Abund[i1] = spec_it1->second;
		++i1;

		if (spec_it1->second >= minAbund){

			Kcon10[i2] = pForest->KSpecIJ[spec_it1->first][spec_it1->first][1];
			Kcon50[i2] = pForest->KSpecIJ[spec_it1->first][spec_it1->first][9];
			++i2;

			for (spec_it2 = pForest->SpecAbund.begin(); spec_it2!=pForest->SpecAbund.end(); ++spec_it2){
				if (spec_it2->second >= minAbund){
					if (spec_it1->first != spec_it2->first){

						Khet10[i3] = pForest->KSpecIJ[spec_it1->first][spec_it2->first][1];
						Khet50[i3] = pForest->KSpecIJ[spec_it1->first][spec_it2->first][9];

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

	for (int ir=0; ir<pForest->nBins1; ir++){
		PCF[ir] = pForest->PCF_all[ir];
		PropCon[ir] = pForest->PropCon[ir];
		SAR1[ir] = pForest->SAR[ir];
	}

	for (int ir=0; ir<pForest->SAR2_n; ir++){
		SAR2_m[ir] = pForest->SAR2_m[ir];
		SAR2_sd[ir] = pForest->SAR2_sd[ir];
	}



	pForest->ClearForest();

	delete pForest;
	delete pPara;
}

}  //end extern "C"

//------------------------------------------------------------------------------
