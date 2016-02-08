// ---------------------------------------------------------------------------
#include "Forest.h"
#include <time.h>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

// ---------------------------------------------------------------------------
CForest::CForest(int seed, double xmax, double ymax, double map_cell,
				     char* map_file, char* rel_dens_file, int n_hab_types){

	// Landscape - forest plot
	Xmax = xmax;
	Ymax = ymax;

	//Read habitat map file
	MapCellSize = map_cell;
	MapXcells = Xmax/MapCellSize;
	MapYcells = Ymax/MapCellSize;
	nHabTypes = n_hab_types;

	Map = new int*[MapXcells];
	for (int x = 0; x < MapXcells; x++)
		Map[x] = new int [MapYcells];

	ifstream MapFile;
	MapFile.open(map_file);

	int hab_type = -99;
	for (int x = 0; x < MapXcells; x++){
		for (int y = 0; y < MapYcells; y++){
			MapFile>>hab_type;
			Map[x][y] = hab_type;
		}
	}
	MapFile.close();

	//Calculate habitat proportions in neighborhood
   int smooth_habitat = 1; //9 neighbors

   HabitatProp = new double**[MapXcells];
   for (int x = 0; x < MapXcells; x++){
		HabitatProp[x] = new double*[MapYcells];
		for (int y = 0; y < MapYcells; y++)
		   HabitatProp[x][y] = new double [nHabTypes];
   }

   //init
   for (int x = 0; x < MapXcells; x++)
		for (int y = 0; y < MapYcells; y++)
         for (int ihab = 0; ihab < nHabTypes; ihab++)
            HabitatProp[x][y][ihab] = 0.0;

   //calculate
   int* hab_prop = new int[nHabTypes];
   int ncells;

   int x2, y2;

   for (int x = 0; x < MapXcells; x++){
      for (int y = 0; y < MapYcells; y++){

         for (int ihab=0; ihab < nHabTypes; ihab++)
            hab_prop[ihab] = 0;

         // loop over neighboring cells
         ncells = 0;
         for (int dx = -smooth_habitat; dx <= smooth_habitat; dx++) {
            for (int dy = -smooth_habitat; dy <= smooth_habitat; dy++) {

               //boundary condition
               x2 = x + dx;
               y2 = y + dy;
               if (((x2 >= 0) && (x2 < MapXcells)) && ((y2 >= 0) && (y2 < MapYcells))){
                  hab_type = Map[x2][y2];
                  if (hab_type < nHabTypes)
                     hab_prop[hab_type]++;
                  ncells++;
               }

            } // for dy
         } // for dx

         for (int ihab=0; ihab < nHabTypes; ihab++)
            HabitatProp[x][y][ihab] = (double) hab_prop[ihab]/ncells;

      } //for y
   } // for x

   delete[] hab_prop;

   /*
   ofstream SmoothMap;
   SmoothMap.open("InOut\\SmoothMap.txt");

   SmoothMap<<"X\tY\t";
   for (int ihab=0; ihab < nHabTypes; ihab++)
      SmoothMap<<"Hab"<<ihab<<"\t";
   SmoothMap<<endl;

   // write out realized habitat map
   for (int x = 0; x < MapXcells; x++){
      for (int y = 0; y < MapYcells; y++){
         SmoothMap<<x*MapCellSize<<"\t"
                  <<y*MapCellSize<<"\t";
         for (int ihab=0; ihab < nHabTypes; ihab++)
            SmoothMap<<HabitatProp[x][y][ihab]<<"\t";
         SmoothMap<<"\n";
      }
   }

   SmoothMap.close();
   */

	// read species relative densities calculated from the data
	ifstream RelHabDensFile;
	RelHabDensFile.open(rel_dens_file);

	vector<double> rel_dens_spec(nHabTypes,0.0);
	string dummy;

	if (RelHabDensFile.good()) {

		getline(RelHabDensFile,dummy);
		RelHabDensFile>>dummy;
		for (int ihab=0; ihab < nHabTypes; ihab++)
			RelHabDensFile>>rel_dens_spec[ihab];
		RelHabDensData.push_back(rel_dens_spec);

		while (RelHabDensFile.good()){
			RelHabDensFile>>dummy;
			for (int ihab=0; ihab < nHabTypes; ihab++)
				RelHabDensFile>>rel_dens_spec[ihab];
			RelHabDensData.push_back(rel_dens_spec);
		}
	}
	RelHabDensFile.close();

	// Grid
	CellSize = 5; // grid cell size for recruitment grid

	XCells = (int)Xmax / CellSize;
	YCells = (int)Ymax / CellSize;

	//Grid = new (CCell*[XCells]);
	Grid = new CCell*[XCells];
	for (int x = 0; x < XCells; x++)
		//Grid[x] = new (CCell[YCells]);
		Grid[x] = new CCell[YCells];


	// Random number generators
	// int seed = (int) time(0);            // random seed
	// RandGen1 = new CRandomMersenneA(seed);

	RandGen1 = new CRandomMersenne(seed);
	RandGen2 = new StochasticLib1(seed);

	// Point pattern variables all trees
	Rmax1 = 50.0;
	BW1 = 1.0;
	nBins1 = (int)ceil(Rmax1 / BW1);

	rvec1 = new double[nBins1];
	for (int i = 0; i < nBins1; i++)
		rvec1[i] = BW1 * i + 0.5 * BW1;

	CountAll = new int[nBins1];
	CountCon = new int[nBins1];
	PropCon = new double[nBins1];

	PCF_all = new double[nBins1];
	NennerPCF = new double[nBins1];
	SAR = new double[nBins1];

	int dmin;
	if (Xmax <= Ymax) dmin = floor(Xmax);
	else              dmin = floor(Ymax);

	SAR2_scales.clear();
	for (int i = 1; i <= dmin; i++)
		if (dmin % i == 0)
			SAR2_scales.push_back(i*i);

	//double scales[] = {1,4,25,100,400,625,2500,10000,15625,62500,250000,500000};

	if (Xmax > Ymax) SAR2_scales.push_back(Xmax*Ymax);

	//if (Xmax < 1000.0) SAR2_n = (int)sizeof(scales) / sizeof(double) - 1;
	//else               SAR2_n =  (int)sizeof(scales) / sizeof(double);

	//for (int i = 0; i < SAR2_n; ++i)
	//	SAR2_scales.push_back(scales[i]);
	SAR2_n = SAR2_scales.size();

	SAR2_m = new double[SAR2_n];
	SAR2_sd = new double[SAR2_n];

	// bivariate point pattern variables
	minAbund = 50;

	Rmax2 = 50; // needs to be smaller than Rmax1!
	BW2 = 5.0;
	nBins2 = (int) ceil(Rmax2 / BW2);

	rvec2 = new double[nBins2];
	for (int i = 0; i < nBins2; i++)
		rvec2[i] = BW2 * i + 0.5 * BW2;

	// Input-Output
	isim = 1;
	irep = 1;

	//TestFile.close();
}

// ---------------------------------------------------------------------------
CForest::~CForest() {
	for (TreeIterV itree = TreeList.begin(); itree != TreeList.end(); ++itree)
		delete(*itree);
	TreeList.clear();

	CumRelAbundMeta.clear();
	SpecAbund.clear();
	SAR2_scales.clear();
	SpecPars.clear();

	for (int x = 0; x < MapXcells; x++)
		delete[] Map[x];
	delete[] Map;

   //  Deallocate 3D array
   for(int x = 0; x < MapXcells; x++){
      for(int y = 0; y < MapYcells; y++)
         delete[] HabitatProp[x][y];
      delete[] HabitatProp[x];
   }
   delete[] HabitatProp;

	for (int x = 0; x < XCells; x++)
		delete[] Grid[x];
	delete[] Grid;

	delete RandGen1;
	delete RandGen2;

	delete[] rvec1;
	delete[] rvec2;

	delete[] CountAll;
	delete[] CountCon;
	delete[] PropCon;

	delete[] PCF_all;
	delete[] NennerPCF;
	delete[] SAR;

	delete[] SAR2_m;
	delete[] SAR2_sd;

	// File
	DivFile.close();
	SAD_File.close();
	AbundFile.close();
	PCF_File.close();
	PropConFile.close();
	SAR1_File.close();
	SAR2_File.close();
	Lf20_File.close();
	Kcon20_File.close();
	Khet20_File.close();

	// CrossPCF_File.close();
	// xPOD_File.close();

	// TestFile.close();

}

//void CForest::SetPars(CPara* pars) {
//}

// ---------------------------------------------------------------------------
void CForest::FileOpen(string label) {
	string FileName;
	string FileNameEnd = label + ".csv";

	FileName = "InOut\\Diversity" + FileNameEnd;
	DivFile.open(FileName.c_str(), ios::in); {
		if (DivFile.good()) {
			DivFile.close();
			DivFile.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			DivFile.clear();
			DivFile.open(FileName.c_str(), ios::out);
			DivFile << "SimNr; RepNr; Step; NSpec; BDevents; Shannon" << endl;
		}
	}


	FileName = "InOut\\Abund" + FileNameEnd;
	AbundFile.open(FileName.c_str(), ios::in); {
		if (AbundFile.good()) {
			AbundFile.close();
			AbundFile.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			AbundFile.clear();
			AbundFile.open(FileName.c_str(), ios::out);
		}
	}

	FileName = "InOut\\SAD" + FileNameEnd;
	SAD_File.open(FileName.c_str(), ios::in); {
		if (SAD_File.good()) {
			SAD_File.close();
			SAD_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			SAD_File.clear();
			SAD_File.open(FileName.c_str(), ios::out);

			SAD_File << "A1; A2_3; A4_7; A8_15; A16_31; A32_63; A64_127; A128_255; "
						<< "A256_511; A512_1023; A1024_2047; A2048-Inf"
						<< endl;
			// SAD_File<<"A1; A2; A3_4; A5_8; A9_16; A17_32; A33_64; A65_128;"
			// <<"A129_256; A257_512; A513_1024; A1025_2048; A2049-Inf;"<<endl;
		}
	}

	FileName = "InOut\\PCFall" + FileNameEnd;
	PCF_File.open(FileName.c_str(), ios::in); {
		if (PCF_File.good()) {
			PCF_File.close();
			PCF_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			PCF_File.clear();
			PCF_File.open(FileName.c_str(), ios::out);
			for (int ibin1 = 0; ibin1 < (nBins1-1); ibin1++)
				PCF_File << rvec1[ibin1] << "; ";
			PCF_File << rvec1[nBins1-1] << endl;
		}
	}

	FileName = "InOut\\PropCon" + FileNameEnd;
	PropConFile.open(FileName.c_str(), ios::in); {
		if (PropConFile.good()) {
			PropConFile.close();
			PropConFile.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			PropConFile.clear();
			PropConFile.open(FileName.c_str(), ios::out);
			for (int ibin1 = 0; ibin1 < (nBins1-1); ibin1++)
				PropConFile << rvec1[ibin1] << "; ";
			PropConFile << rvec1[nBins1-1] <<endl;
		}
	}

	FileName = "InOut\\SAR1_" + FileNameEnd;
	SAR1_File.open(FileName.c_str(), ios::in); {
		if (SAR1_File.good()) {
			SAR1_File.close();
			SAR1_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			SAR1_File.clear();
			SAR1_File.open(FileName.c_str(), ios::out);
			for (int ibin1 = 0; ibin1 < (nBins1-1); ibin1++)
				SAR1_File << rvec1[ibin1] << "; ";
			SAR1_File << rvec1[nBins1-1] <<endl;
		}
	}

	FileName = "InOut\\SAR2_" + FileNameEnd;
	SAR2_File.open(FileName.c_str(), ios::in); {
		if (SAR2_File.good()) {
			SAR2_File.close();
			SAR2_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			SAR2_File.clear();
			SAR2_File.open(FileName.c_str(), ios::out);
			for (int i = 0; i < SAR2_n; i++)
				SAR2_File << "m" << SAR2_scales[i] << "; ";
			for (int i = 0; i < (SAR2_n-1); i++)
				SAR2_File << "sd" << SAR2_scales[i] << "; ";
			SAR2_File << SAR2_scales[SAR2_n-1] <<endl;
		}
	}

	/*
	FileName = "InOut\\Lf20_" + FileNameEnd;
	Lf20_File.open(FileName.c_str(), ios::in); {
		if (Lf20_File.good()) {
			Lf20_File.close();
			Lf20_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			Lf20_File.clear();
			Lf20_File.open(FileName.c_str(), ios::out);
		}
	}
	*/

	FileName = "InOut\\Kcon20_" + FileNameEnd;
	Kcon20_File.open(FileName.c_str(), ios::in); {
		if (Kcon20_File.good()) {
			Kcon20_File.close();
			Kcon20_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			Kcon20_File.clear();
			Kcon20_File.open(FileName.c_str(), ios::out);
		}
	}

	FileName = "InOut\\Khet20_" + FileNameEnd;
	Khet20_File.open(FileName.c_str(), ios::in); {
		if (Khet20_File.good()) {
			Khet20_File.close();
			Khet20_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			Khet20_File.clear();
			Khet20_File.open(FileName.c_str(), ios::out);
		}
	}

	/*
	FileName = "InOut\\CrossPCF_"+FileNameEnd;
	CrossPCF_File.open(FileName.c_str(),ios::in);
	{
	if(CrossPCF_File.good()){
	CrossPCF_File.close();
	CrossPCF_File.open(FileName.c_str(),ios::out|ios::app);
	}
	else {
	CrossPCF_File.clear();
	CrossPCF_File.open(FileName.c_str(),ios::out);
	CrossPCF_File<<"SpecID1; SpecID2; ";
	for (int ibin2=0; ibin2 < nBins2; ibin2++) CrossPCF_File<<rvec2[ibin2]<<";";
	CrossPCF_File<<"\n";
	}
	}
	 */

	/*
	FileName = "InOut\\xPOD_" + FileNameEnd;
	xPOD_File.open(FileName.c_str(), ios::in); {
		if (xPOD_File.good()) {
			xPOD_File.close();
			xPOD_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			xPOD_File.clear();
			xPOD_File.open(FileName.c_str(), ios::out);
		}
	}

	FileName = "InOut\\NNDist_" + FileNameEnd;
	NNDist_File.open(FileName.c_str(), ios::in); {
		if (NNDist_File.good()) {
			NNDist_File.close();
			NNDist_File.open(FileName.c_str(), ios::out | ios::app);
		}
		else {
			NNDist_File.clear();
			NNDist_File.open(FileName.c_str(), ios::out);
		}
	}
	*/
}

// ---------------------------------------------------------------------------
void CForest::ClearForest() {
	for (TreeIterV itree = TreeList.begin(); itree != TreeList.end(); ++itree)
		delete (*itree);
	TreeList.clear();

	for (unsigned int ispec = 0; ispec < SpecMax; ispec++)
		delete[] InteractMat[ispec];
   delete[] InteractMat;

	/*
	for (TreeIterV itree = RecruitList.begin(); itree != RecruitList.end(); ++itree)
	delete (*itree);
	RecruitList.clear();
	 */

	SpecAbund.clear();
	CumRelAbundMeta.clear();
	// SpecPars.clear();
}

// ---------------------------------------------------------------------------
inline double CForest::Distance(double x1, double y1, double x2, double y2) {
	return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

// ---------------------------------------------------------------------------
/*
inline double CForest::DistSquare(double x1, double y1, double x2, double y2)
{
return ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}
 */

// ---------------------------------------------------------------------------
inline void CForest::PeriodBound(double& xx, double& yy) {
	xx = Xmax * (xx / Xmax - floor(xx / Xmax));
	yy = Ymax * (yy / Ymax - floor(yy / Ymax));
}

// ---------------------------------------------------------------------------
inline void CForest::BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax) {
	xx = xx % Xmax;
	if (xx < 0)
		xx = Xmax + xx;

	yy = yy % Ymax;
	if (yy < 0)
		yy = Ymax + yy;
}

// ---------------------------------------------------------------------------
inline void CForest::BoundGrid(int& xx, int& yy, double& xb, double& yb) {
	if (xx >= XCells) {
		xx = xx - XCells;
		xb = +Xmax;
	}
	else if (xx < 0) {
		xx = xx + XCells;
		xb = -Xmax;
	}

	if (yy >= YCells) {
		yy = yy - YCells;
		yb = +Ymax;
	}
	else if (yy < 0) {
		yy = yy + YCells;
		yb = -Ymax;
	}
}

// ---------------------------------------------------------------------------
inline int CForest::GetRandSpec() {
	unsigned int ispec = 0;
	double r1 = RandGen1->Random();
	bool choose = false;

	while (choose == false) {
		if (r1 <= CumRelAbundMeta[ispec])
			choose = true;
		// if (r1 <= CumProbImmi[ispec]) choose = true;
		else
			++ispec;
	}

	return(ispec);
}

// ---------------------------------------------------------------------------
vector<double>CForest::SeqConstruct(unsigned int J, double theta, double m) {
	unsigned int a = 0, s = 0;
	unsigned int index;

	double I, R1, R2, x, y;

	vector<int>SpecLabelInd(J, 0);
	vector<int>AncLabelInd(J, 0);
	vector<int>SpecLabelAnc(J, 0);

	if (m < 1.0)
		I = m * (J - 1) / (1 - m);

	for (unsigned int j = 1; j <= J; ++j) {
		if (m < 1.0)
			R1 = (double)I / (I + j - 1);
		else
			R1 = 1.0;

		x = RandGen1->Random();
		if (x <= R1) { // new immigrating ancestor
			++a;
			AncLabelInd[j - 1] = a;
			R2 = (double)theta / (theta + a - 1);
			y = RandGen1->Random();
			if (y <= R2) { // immigrant of new species
				++s;
				SpecLabelInd[j - 1] = s;
				SpecLabelAnc[a - 1] = s;
			}
			else { // immigrant of previous species
				index = RandGen1->IRandom(0, a - 2);
				SpecLabelInd[j - 1] = SpecLabelAnc[index];
				SpecLabelAnc[a - 1] = SpecLabelAnc[index];
			}
		}
		else { // offspring of local community
			index = RandGen1->IRandom(0, j - 2);
			AncLabelInd[j - 1] = AncLabelInd[index];
			SpecLabelInd[j - 1] = SpecLabelInd[index];
		}
		// check = SpecLabelInd[j-1];
	}

	unsigned int Smax = s;

	// convert to abundance
	vector<int>Abund(Smax, 0);
	for (unsigned int j = 0; j < J; ++j) {
		++Abund[SpecLabelInd[j] - 1];
	}

	// convert to relative and cumulative rel. abund.
	vector<double>RelAbund(Smax, 0);
	vector<double>CumRelAbund(Smax, 0);

	for (unsigned int is = 0; is < Smax; ++is) {
		RelAbund[is] = (double)Abund[is] / J;
		if (is == 0)
			CumRelAbund[is] = RelAbund[is];
		else
			CumRelAbund[is] = CumRelAbund[is - 1] + RelAbund[is];
	}

	return(CumRelAbund);
}

// ---------------------------------------------------------------------------
void CForest::Init() {

	CTree* pTree1;

	int iX1, iY1;

	BD_5years = 0;
	BD_total = 0;

	// Init Grid Cells
	for (iX1 = 0; iX1 < XCells; iX1++)
		for (iY1 = 0; iY1 < YCells; iY1++)
			Grid[iX1][iY1].InitCell(iX1, iY1,1);

	//grid_steps = (int) ceil(2.0 * Pars->r_max / CellSize);
	grid_steps = (int) ceil(Pars->r_max / CellSize);

	// Init metacommunity
	CumRelAbundMeta.clear();
	CumRelAbundMeta = SeqConstruct(Pars->Jmeta, Pars->theta);
	SpecMax = CumRelAbundMeta.size();

	double meanD_spec, sdD_spec, JC_spec;

	//Init species parameters
	double sigma_dm = sqrt(log(1.0 + (Pars->sd_dm_spec*Pars->sd_dm_spec)
												 /(Pars->m_dm_spec*Pars->m_dm_spec)));
	double mu_dm = log(Pars->m_dm_spec) - 0.5 * sigma_dm*sigma_dm;

	double sigma_JC = sqrt(log(1.0 + (Pars->sd_JCspec *Pars->sd_JCspec)
												 /(Pars->m_JCspec*Pars->m_JCspec)));
	double mu_JC = log(Pars->m_JCspec) - 0.5 * sigma_JC*sigma_JC;

	//Construct species interaction matrix
	InteractMat = new double*[SpecMax];
   for (unsigned int ispec = 0; ispec < SpecMax; ispec++)
      InteractMat[ispec] = new double[SpecMax];

   //Initialize interaction matrix
   for (unsigned int ispec=0; ispec<SpecMax; ispec++)
      for (unsigned int jspec = 0; jspec<SpecMax; jspec++)
         InteractMat[ispec][jspec] = 1.0;

	//double relabund;
	for (unsigned int ispec=0; ispec<SpecMax; ispec++){

		//DISPERSAL PARAMETERS

		//negative relationship between dispersal distance and
		//metacommunity abundance
      //parameters from Muller Landau 2008

		//if (ispec == 0) relabund = CumRelAbundMeta[0];
		//else  relabund = CumRelAbundMeta[ispec] - CumRelAbundMeta[ispec-1];

		//meanD_spec = exp(1.85 - 0.23 * log(relabund) + RandGen2->Normal(0,0.7));
		//sdD_spec = meanD_spec;

		// difference among species
		meanD_spec = exp(RandGen2->Normal(mu_dm,sigma_dm));
		sdD_spec = meanD_spec;

		//no difference among species
		//meanD_spec = Pars->m_dm_spec;
		//sdD_spec = Pars->sd_dm_spec;

		// NCDD - NEGATIVE CONSPECIFIC DENSITY DEPENDENCE

		//if (Pars->sd_JCspec == 0.0) JC_spec =  Pars->m_JCspec;
		//else
		//JC_spec = 1.0 + exp(RandGen2->Normal(mu_JC,sigma_JC));  //log-normal distribution

		//JC_spec = exp(RandGen2->Normal(mu_JC,sigma_JC));
		JC_spec = RandGen2->Normal(Pars->m_JCspec,Pars->sd_JCspec); //normal distribution
		if (JC_spec < 1.0) //exclude positive density dependence
         JC_spec = 1.0;

		InteractMat[ispec][ispec] = JC_spec;

		//relationship between NCDD and metacommunity abundance
		//JC_spec = 1.0 + exp(Pars->m_JCspec + Pars->sd_JCspec * (log10(relabund) + 2.0));
		//JC_spec = 1.0 + Pars->m_JCspec * (log10(relabund) + Pars->sd_dm_spec);

		//Random species interactions
		//for (unsigned int jspec = 0; jspec < SpecMax; jspec++){
         //InteractMat[ispec][jspec] = RandGen2->Normal(5.0,2.0);
         //InteractMat[ispec][jspec] = RandGen1->Random()*10;

         //if (InteractMat[ispec][jspec]<0)
         //   InteractMat[ispec][jspec] = 0;
		//}


		// HABITAT ASSOCIATIONS
		int irand = RandGen1->IRandom(0, RelHabDensData.size() - 1);
		SpecPars[ispec] = CSpecPara(meanD_spec,sdD_spec,JC_spec,RelHabDensData[irand]);
	}

	// Init immigration rate
	if (Pars->m < 0.0)
		m = (2 * Xmax + 2 * Ymax) * Pars->m_dm_spec / (Pi * Xmax * Ymax);
	else
		m = Pars->m;

	// Init Trees
	double x, y;
	int SpecID;

	//for (unsigned int i=0; i < SpecMax; i++) SpecAbund[i] = 0;
	SpecAbund.clear();

	// random init
	TreeID = 0;
	NTrees = Pars->NTrees;
	for (int i = 0; i < NTrees; i++) {
		x = RandGen1->Random() * Xmax;
		y = RandGen1->Random() * Ymax;
		SpecID = GetRandSpec();

		pTree1 = new CTree(TreeID, x, y, SpecID, Pars->r_max, Pars->bSurv);
		TreeList.push_back(pTree1);
		++TreeID;

		iX1 = (int) floor(x / CellSize);
		iY1 = (int) floor(y / CellSize);
		Grid[iX1][iY1].TreeList.push_back(pTree1);

		++SpecAbund[pTree1->SpecID];
	}

	/*
	ifstream BCI_File("InOut\\bci2005_10cm.txt");

	string line;
	getline(BCI_File,line);

	int tag;
	string SpecStr, GrForm;
	double dbh;

	do {
	BCI_File>>tag;
	BCI_File>>x;
	BCI_File>>y;
	BCI_File>>SpecStr;
	BCI_File>>SpecID;
	BCI_File>>dbh;
	BCI_File>>GrForm;

	if (!BCI_File.eof()){
	pTree1 = new CTree(x,y,SpecID,Pars->rTreeMax,0);
	TreeList.push_back(pTree1);

	iX1 = (int) floor(pTree1->X/CellSize);
	iY1 = (int) floor(pTree1->Y/CellSize);
	Grid[iX1][iY1].TreeList.push_back(pTree1);

	++SpecAbund[pTree1->SpecID];

	}
	} while(!BCI_File.eof());

	BCI_File.close();
	NTrees = TreeList.size();
	 */

   /*
   int iX2, iY2;
   CTree* pTree2;

   double d12, dmin;
	double xbound, ybound;

	// Init competition and recruitment
	for (iX1 = 0; iX1 < XCells; iX1++) {
		for (iY1 = 0; iY1 < YCells; iY1++) {

			// loop over trees in cell 1
			for (TreeIterL itree1 = Grid[iX1][iY1].TreeList.begin();
				itree1 != Grid[iX1][iY1].TreeList.end(); itree1++) {

				pTree1 = *itree1;

				// loop over neighboring cells
				for (int dx = -grid_steps; dx <= grid_steps; dx++) {
					for (int dy = -grid_steps; dy <= grid_steps; dy++) {

						iX2 = iX1 + dx;
						iY2 = iY1 + dy;

						xbound = 0.0;
						ybound = 0.0;

						// periodic boundary
						BoundGrid(iX2, iY2, xbound, ybound);

						// loop over trees in cell 2
						for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
							itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

							// update competition ------------------------------
							pTree2 = *itree2;
							if (pTree1 != pTree2) {
								d12 = Distance(pTree1->X, pTree1->Y,
									pTree2->X + xbound, pTree2->Y + ybound);

								// Distance based
								dmin = pTree1->R;

								if (d12 < dmin) { // do trees overlap?

									if (d12<0.0001) d12 = 0.0001;


//									if (pTree1->SpecID == pTree2->SpecID)
//										pTree1->NCI = pTree1->NCI + (SpecPars[pTree1->SpecID].JCfac
//																		  * (1.0 - d12/Pars->r_max));
//
//									else
//										pTree1->NCI = pTree1->NCI + (1.0 - d12/Pars->r_max);


                           //pTree1->NCI = pTree1->NCI + InteractMat[pTree1->SpecID][pTree2->SpecID]*(1.0 - d12/Pars->r_max);
                           //pTree2->NCI = pTree2->NCI + InteractMat[pTree2->SpecID][pTree1->SpecID]*(1.0 - d12/Pars->r_max);

								} // if overlap

							} // if tree1 != tree2
							// ------------------------------------------------

						} // end tree 2

					} // end dy
				} // end dx

			} // end tree 1

		} // grid y
	} // grid x

	// calculate mortality for each tree;
	for (TreeIterV itree = TreeList.begin(); itree != TreeList.end(); ++itree)
		(*itree)->GetPSurv2(Pars->aSurv, Pars->bSurv);
   */
}

// ---------------------------------------------------------------------------
double CForest::GetProbRecruit(double x1, double y1, unsigned int spec_id)
{
	int iX1 = (int)floor(x1/CellSize);
	int iY1 = (int)floor(y1/CellSize);

	int iX2, iY2;
	double xbound, ybound, d12;
	CTree* pTree2;

	// calculate NCI
	double NCI = 0.0;

	for (int dxi = -grid_steps; dxi <= grid_steps; dxi++) {
		for (int dyi = -grid_steps; dyi <= grid_steps; dyi++) {

			iX2 = iX1 + dxi;
			iY2 = iY1 + dyi;

			xbound = 0.0;
			ybound = 0.0;

			// torus boundary
			BoundGrid(iX2, iY2, xbound, ybound);

			// loop over trees in cell 2
			if (Grid[iX2][iY2].TreeList.size() > 0) {
				for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
					itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

					pTree2 = (*itree2);
					d12 = Distance(x1, y1, pTree2->X + xbound,pTree2->Y + ybound);

					//distance based
					if (d12 < Pars->r_max) { // do trees overlap?
						if (d12<0.0001) d12 = 0.0001;

							/*
							if (spec_id == pTree2->SpecID){
								NCI = NCI + (SpecPars[spec_id].JCfac
										  * (1.0 - d12/Pars->r_max));
							}
							else {
								NCI = NCI + (1.0 - d12/Pars->r_max);
							}
							*/

							NCI = NCI + InteractMat[spec_id][pTree2->SpecID]*(1.0 - d12/Pars->r_max);

						} // if overlap
				} // end tree 2
			} // if TreeList > 0

		} // end dx
	} // end dy

	double prob_rec1;
	if (NCI == 0) prob_rec1 = 1.0;
	else prob_rec1 = 1.0 - NCI/(Pars->aRec + NCI);

	// Habitat effect
	iX1 = (int) floor(x1 / MapCellSize);
	iY1 = (int) floor(y1 / MapCellSize);


   //calculate average recruitment probability in habitat here
	//int hab_type =  Map[iX1][iY1];
	//double rel_dens =  SpecPars[spec_id].RelHabDens[hab_type];

	double rel_dens = 0;
	for (int ihab = 0; ihab < nHabTypes; ihab++)
      rel_dens += HabitatProp[iX1][iY1][ihab] * SpecPars[spec_id].RelHabDens[ihab];

	double prob_rec2;
	if (Pars->aHab > 0.0001)
		prob_rec2 = prob_rec1*pow(rel_dens,Pars->aHab);
   else prob_rec2 = prob_rec1;

	return(prob_rec2);
}

// ---------------------------------------------------------------------------
// add tree to competition neighborhood and to recruitment grid
void CForest::AddTree(CTree* pTree1)
{
	// add tree grid and overlap cells
	int iX1 = (int)floor(pTree1->X / CellSize);
	int iY1 = (int)floor(pTree1->Y / CellSize);

	Grid[iX1][iY1].TreeList.push_back(pTree1);

   /*
   double d12, dmin, xbound, ybound;
	int iX2, iY2;

	CTree* pTree2;

	pTree1->NCI = 0.0;

	for (int dxi = -grid_steps; dxi <= grid_steps; dxi++) {
		for (int dyi = -grid_steps; dyi <= grid_steps; dyi++) {

			iX2 = iX1 + dxi;
			iY2 = iY1 + dyi;

			xbound = 0.0;
			ybound = 0.0;

			// torus boundary
			BoundGrid(iX2, iY2, xbound, ybound);

			// loop over trees in cell 2
			if (Grid[iX2][iY2].TreeList.size() > 0) {
				for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
					itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

					pTree2 = (*itree2);
					if (pTree1 != pTree2) {
						d12 = Distance(pTree1->X, pTree1->Y, pTree2->X + xbound, pTree2->Y + ybound);

						//distance based
						dmin = pTree1->R;

						if (d12 < dmin) { // do trees overlap?

							if (d12<0.0001) d12 = 0.0001;


//							if (pTree1->SpecID == pTree2->SpecID){
//								pTree1->NCI = pTree1->NCI + (SpecPars[pTree1->SpecID].JCfac
//																  * (1.0 - d12/Pars->r_max));
//								pTree2->NCI = pTree2->NCI + (SpecPars[pTree2->SpecID].JCfac
//																  * (1.0 - d12/Pars->r_max));
//							}
//							else {
//								pTree1->NCI = pTree1->NCI + (1.0 - d12/Pars->r_max);
//								pTree2->NCI = pTree2->NCI + (1.0 - d12/Pars->r_max);
//							}

							pTree1->NCI = pTree1->NCI + InteractMat[pTree1->SpecID][pTree2->SpecID]*(1.0 - d12/Pars->r_max);
							pTree2->NCI = pTree2->NCI + InteractMat[pTree2->SpecID][pTree1->SpecID]*(1.0 - d12/Pars->r_max);

							pTree2->GetPSurv2(Pars->aSurv, Pars->bSurv);

						} // if overlap
					} // if tree1 != tree2
				} // end tree 2
			}

		} // end dx
	} // end dy

	pTree1->GetPSurv2(Pars->aSurv, Pars->bSurv);

*/

}

// ---------------------------------------------------------------------------
// remove tree from competition neighborhood and from recruitment grid
void CForest::RemoveTree(CTree* pTree1) {

	// remove tree from grid
	int iX1 = (int)floor(pTree1->X / CellSize);
	int iY1 = (int)floor(pTree1->Y / CellSize);

	// remove tree from competition grid ----------------------------------------
	// and update competition indices of neighbor
	/*
	int iX2, iY2;
	double xbound, ybound;
	double dmin, d12;

	CTree* pTree2;

	for (int diX = -grid_steps; diX <= grid_steps; diX++) {
		for (int diY = -grid_steps; diY <= grid_steps; diY++) {

			iX2 = iX1 + diX;
			iY2 = iY1 + diY;

			xbound = 0.0;
			ybound = 0.0;

			BoundGrid(iX2, iY2, xbound, ybound);

			// loop over trees in cell 2
			if (Grid[iX2][iY2].TreeList.size() > 0) {
				for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
					itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

					pTree2 = *itree2;
					if (pTree1 != pTree2) {
						d12 = Distance(pTree1->X, pTree1->Y, pTree2->X + xbound,
							pTree2->Y + ybound);

						// distance based
						dmin = pTree1->R;

						if (d12 < dmin) { // do trees overlap?

							if (d12<0.0001) d12 = 0.0001;

//							if (pTree1->SpecID == pTree2->SpecID)
//								pTree2->NCI = pTree2->NCI - (SpecPars[pTree2->SpecID].JCfac
//																  * (1.0 - d12/Pars->r_max));
//
//							else
//								pTree2->NCI = pTree2->NCI - (1.0 - d12/Pars->r_max);
//								//pTree2->NCI = pTree2->NCI - 1.0/pow(d12,Pars->beta_r);
//

							pTree2->NCI = pTree2->NCI - InteractMat[pTree2->SpecID][pTree1->SpecID]*(1.0 - d12/Pars->r_max);

							if (pTree2->NCI < 0.0)
								pTree2->NCI = 0.0;

							pTree2->GetPSurv2(Pars->aSurv, Pars->bSurv);

						} // if overlap
					} // if tree1 != tree2
				} // end tree 2
			}
		} // end dx
	} // end dy

*/
	pTree1->NCI = 0;

	Grid[iX1][iY1].TreeList.remove(pTree1);
}

// ---------------------------------------------------------------------------
void CForest::GetNewXY(double &x1, double &y1, int idspec) {

	double x0 = x1;
	double y0 = y1;

	if (Pars->m_dm_spec < 999.0) {

		double r_dist = exp(RandGen2->Normal(SpecPars[idspec].muDisp, SpecPars[idspec].sigmaDisp));
		//double r_dist = exp(RandGen2->Normal(Pars->muDisp, Pars->sigmaDisp));
		double r_angle = RandGen1->Random() * 2.0 * Pi;

		// TestFile<<Pars->muDisp<<"\t"<<Pars->sigmaDisp<<endl;
		// TestFile<<r_dist<<"\t"<<r_angle<<endl;

		x1 = x0 + cos(r_angle) * r_dist;
		y1 = y0 + sin(r_angle) * r_dist;

		PeriodBound(x1, y1);
	}

	else {
		// global dispersal
		x1 = RandGen1->Random() * Xmax; // random coordinates
		y1 = RandGen1->Random() * Ymax;
	}

	/*
	dist = Distance(x0,y0,x1,y1);

	//reflecting boundary condition - draw again
	if ((x1 < 0) || (y1 < 0) || (x1 >= Xmax) || (y1 >= Ymax) || (dist >= MaxDist)){
	do {

	//r_dist = -Pars->mDisp*log(RandGen1->Random());
	r_dist = exp(RandGen2->Normal(d_mu,d_sigma));
	r_angle = RandGen1->Random()*2.0*Pi;

	x1 = x0 + cos(r_angle)*r_dist;
	y1 = y0 + sin(r_angle)*r_dist;

	dist = Distance(x0,y0,x1,y1);
	}
	while ((x1 < 0) || (y1 < 0) || (x1 >= Xmax) || (y1 >= Ymax) || (dist >= MaxDist));
	}
	 */
}

// ---------------------------------------------------------------------------
bool CForest::BirthDeathAsync() {

	CTree* pTree;

	int ID_Spec;
	double xnew, ynew;
	CTree* pMotherTree;

	bool recruit = false;
	bool stop = false;

	//double rand1;
	double ProbRecruit;
	// int iX, iY;

	// BD_5years = 0;

	// choose random tree´and kill this tree
	pTree = TreeList[RandGen1->IRandom(0, TreeList.size() - 1)];

	double pkill = 1.0 - pTree->pSurv;
	//double pkill = 1.0 - Pars->bSurv;

	if (RandGen1->Random() < pkill) {
		//kill = true;
		++BD_5years;
		++BD_total;

		RemoveTree(pTree);

		if (SpecAbund[pTree->SpecID] > 0) --SpecAbund[pTree->SpecID];
		if (SpecAbund[pTree->SpecID] == 0)  SpecAbund.erase(pTree->SpecID);

		// recruitment and speciation
      recruit = false;
		int ntrials = 0;
		const int max_trials = 999;

		// immigration from metacommunity pool
      if (RandGen1->Random() < m) {
         ID_Spec = GetRandSpec();

         do {
            xnew = RandGen1->Random() * Xmax; // random coordinates
            ynew = RandGen1->Random() * Ymax;

            ProbRecruit = GetProbRecruit(xnew,ynew,ID_Spec);

            ++ntrials;
            if (RandGen1->Random() < ProbRecruit)
               recruit = true;
         } while ((recruit == false) && (ntrials < max_trials));
      }

      // mother tree within the plot - local recruitment
      else {
         do {
            pMotherTree = TreeList[RandGen1->IRandom(0, TreeList.size() - 1)];
            ID_Spec = pMotherTree->SpecID;

            xnew = pMotherTree->X;
            ynew = pMotherTree->Y;
            GetNewXY(xnew, ynew, ID_Spec);

            ProbRecruit = GetProbRecruit(xnew,ynew,ID_Spec);

            ++ntrials;
            if (RandGen1->Random() < ProbRecruit)
               recruit = true;

         } while ((recruit == false) && (ntrials < max_trials));
      }

		if (recruit == true) {
			++TreeID;
			pTree->TreeID = TreeID;
			pTree->X = xnew;
			pTree->Y = ynew;
			pTree->SpecID = ID_Spec;
			++SpecAbund[pTree->SpecID];
			AddTree(pTree);
		}

		else {
			stop = true;
		}
		//cout<<ntrials1<<"\t"<<ntrials2<<endl;
	}

	return(stop);
}

// ---------------------------------------------------------------------------
void CForest::OneRun(int isim, int irep, int ngen, bool steps_out, bool r_mode)
{
	Init();

	if (steps_out){
		GetPPA();
		GetSAR2();
		WriteOutput(BD_total,isim,irep);
		//WriteTreesTime(0);
	}

	BD_max = ngen * TreeList.size();

	int count = 0;
	bool stoprun = false;

	int nspec = GetSAD();

	//while (nspec > 1){
	while (BD_total < BD_max) {

		// one loop representing five years = NTrees test for survival/mortality
		BD_5years = 0;
		for (int i = 0; i < NTrees; ++i) {
			stoprun = BirthDeathAsync();
			if (stoprun == true) break;
		}

		if (stoprun == true) break;

		count++;

		// if (kill && steps_out){
		if (count % 200  == 0) {
			if (steps_out) {
				cout << "     Step " << BD_total << endl;
				cout << "     Species "<<nspec<<"\n";
				GetPPA();
				GetSAR2();
				WriteOutput(BD_total, isim, irep);
				// WriteTrees();
				//WriteTreesTime(count);
			}
		}

		nspec = GetSAD();

	}  // while BD_total < BD_max

	unsigned int nCensusOut = 0;

	/*
	for (int icensus=0; icensus < nCensusOut; ++icensus) {
		// one loop representing five years = NTrees test for survival/mortality
		BD_5years = 0;
		for (int i = 0; i < NTrees; ++i) {
			stoprun = BirthDeathAsync();
			if (stoprun == true) break;
		}

		if (stoprun == true) break;

		count++;

		cout << "     Step " << BD_total << endl;
		GetPPA();
		GetSAR2();
		WriteOutput(BD_total, isim, irep);
		WriteTreesTime(count);
	}
	*/

	if ((!steps_out && !r_mode) && (nCensusOut==0)){
		GetPPA();
		GetSAR2();
		WriteOutput(BD_total,isim,irep);
		//WriteTreesTime(1);
		//WriteSpecPar();
	}
}

// ---------------------------------------------------------------------------
void CForest::WriteOutput(int istep, int isim = 1, int irep = 1) {

	// Community data  ----------------------------------------------
	double Shannon = GetShannon();
	int NSpec = GetSAD();

	DivFile << isim << "; " << irep << "; " << istep << "; " << NSpec
		 << "; " << BD_5years << "; " << Shannon << endl;

	for (int i = 0; i < (MaxSAD-1); ++i)
		SAD_File << SAD[i] << "; ";
	SAD_File << SAD[MaxSAD-1] << endl;

	for (int ibin1 = 0; ibin1 < (nBins1-1); ++ibin1) {
		PCF_File << PCF_all[ibin1] << "; ";
		PropConFile << PropCon[ibin1] << "; ";
		SAR1_File << SAR[ibin1] << "; ";
	}
	PCF_File << PCF_all[nBins1-1] << endl;
	PropConFile << PropCon[nBins1-1] << endl;
	SAR1_File << SAR[nBins1-1] << endl;

	// Species data ------------------------------------------------------------

	map<int, int>::iterator spec_it1, spec_it2;

	for (spec_it1 = SpecAbund.begin(); spec_it1!=SpecAbund.end(); ++spec_it1){

		AbundFile << spec_it1->second <<"; ";

		if (spec_it1->second >= minAbund){

			//Kf20_File<<Kf[spec_it1->first]<<";";
			//Lf20_File<<Lf[spec_it1->first]<<";";

			for (spec_it2 = SpecAbund.begin(); spec_it2!=SpecAbund.end(); ++spec_it2){
				if (spec_it2->second >= minAbund){
					if (spec_it1->first == spec_it2->first)
						Kcon20_File<<KSpecIJ[spec_it1->first][spec_it2->first][3]<<"; ";
					else {
						Khet20_File<<KSpecIJ[spec_it1->first][spec_it2->first][3]<<"; ";
						//xPOD_File<<xPOD[spec_it1->first][spec_it2->first]<<"; ";
						//NNDist_File<<NNDistSpecIJ[spec_it1->first][spec_it2->first]<<"; ";
					}

					//CrossPCF_File<<spec_it1->first<<"; "<<spec_it2->first<<"; ";

					//for (int ibin2=0; ibin2  <nBins2; ++ibin2)
					//	CrossPCF_File<<gSpecIJ[spec_it1->first][spec_it2->first][ibin2]<<"; ";
					//CrossPCF_File<<endl;

				}  // if spec2 >= minAbund
			} // for spec 2
		} // if spec1 >= minAbund


	}  // for spec1

	AbundFile<<endl;

	Kcon20_File<<endl;
	Khet20_File<<endl;

//	xPOD_File<<endl;
//	NNDist_File<<endl;

	for (int iscale = 0; iscale < SAR2_n; ++iscale)
		SAR2_File << SAR2_m[iscale] << "; ";
	for (int iscale = 0; iscale < (SAR2_n - 1); ++iscale)
		SAR2_File << SAR2_sd[iscale] << "; ";
	SAR2_File <<  SAR2_sd[SAR2_n-1] << endl;
}

// ---------------------------------------------------------------------------
void CForest::WriteSpecPar() {
	ofstream OutFile1("InOut\\SpecOut1.txt");

	OutFile1<< "SpecID\t"
			<< "MetaAbund\t"
			<< "muDisp\t"
			<< "sigmaDisp\t"
			<< "JC\t";
	for (int ihab=0; ihab < nHabTypes; ihab++)
		OutFile1 << "RelDens"<<ihab<<"\t";

	OutFile1<< "LocalAbund\t"
		   << "Kcon20"
		   << endl;

   ofstream OutFile2("InOut\\InteractMat1.txt");

	//map<int, CSpecPara>::iterator spec_it1;

	double relabund, Kcon20;

	//for (spec_it1 = SpecPars.begin(); spec_it1!=SpecPars.end(); ++spec_it1){
	for (unsigned int i=0; i < SpecMax; ++i) {

		if (i==0) relabund = CumRelAbundMeta[i];
		else      relabund = CumRelAbundMeta[i] -  CumRelAbundMeta[i-1];

		if (SpecAbund[i] >= minAbund) Kcon20 = KSpecIJ[i][i][3];
		else                          Kcon20 = 0.0;


		OutFile1<< i <<"\t"
				<<relabund <<"\t"
				<<SpecPars[i].muDisp<<"\t"
				<<SpecPars[i].sigmaDisp<<"\t"
				<<SpecPars[i].JCfac<<"\t";

		for (int ihab=0; ihab < nHabTypes; ihab++)
			OutFile1<<SpecPars[i].RelHabDens[ihab]<<"\t";

		OutFile1<<SpecAbund[i]<<"\t"
				<<Kcon20
				<<endl;

		for (unsigned int j=0; j < (SpecMax-1); ++j)
         OutFile2<<InteractMat[i][j]<<"\t";
      OutFile2<<InteractMat[i][SpecMax-1]<<endl;

	}
	OutFile1.close();

}

// ---------------------------------------------------------------------------
/*
void CForest::WriteTrees() {
	ofstream OutFile1("InOut\\TreesOut1.txt");
	OutFile1 << "Nr\t"
	// <<"TreeID\t"
	<< "X\t" << "Y\t" << "Radius\t" << "SpecID\t" << "Azoi\t" << "Acomp\t" <<
		 "Comp\t" << "State" << endl;

	CTree* pTree;

	for (unsigned int i = 0; i < TreeList.size(); i++) {

		pTree = TreeList[i];

		OutFile1 << i + 1 << "\t"
			  // <<pTree->TreeID<<"\t"
				 << pTree->X << "\t"
				 << pTree->Y << "\t"
				 << pTree->R << "\t"
				 << pTree->SpecID << "\t"
				 << pTree->pSurv << "\n";
	}

	OutFile1.close();

//	ofstream OutFile2("InOut\\Grid1.txt");
//	for (int x = 0; x < XCells; x++) {
//	for (int y = 0; y < YCells; y++) {
//	OutFile2<<Grid[x][y].nTreesOverlap<<"\t";
//	//OutFile2<<Grid[x][y].TreeList.size()<<"\t";
//	}
//	OutFile2<<"\n";
//	}
//
//	OutFile2.close();


//	ofstream OutFile2("InOut\\ISAR1.csv");
//	for (unsigned int i = 0; i < TreeList.size(); i++) {
//	for (int rbin = 0; rbin < nBins1; rbin++){
//	pTree = TreeList[i];
//	OutFile2 << pTree->ISAR[rbin]<<";";
//	}
//	OutFile2 << "\n";
//	}
//	OutFile2.close();
//}
*/

// ---------------------------------------------------------------------------
void CForest::WriteTreesTime(int tstep) {
	string FileName = "InOut\\TreesOut" + IntToString(tstep) + ".csv";
	ofstream OutFile1(FileName.c_str());
	OutFile1 << "Nr;"
			 <<"TreeID;"
			 << "X;"
			 << "Y;"
			 //<< "Radius;"
			 << "SpecID"
			 //<< "Azoi;"
			 //<< "Acomp;"
			 //<< "Psurv;"
			 //<< "State"
			 << endl;

	CTree* pTree;

	// int ID_SpecMax = GetSpecMaxAbund();

	for (unsigned int i = 0; i < TreeList.size(); i++) {

		pTree = TreeList[i];

		OutFile1 << i << ";"
				 <<pTree->TreeID<<";"
				 << pTree->X<< ";"
				 << pTree->Y << ";"
				 //<< pTree->R << ";"
				 <<pTree->SpecID
				 //<< pTree->pSurv << ";"
				 //<< pTree->State2
				 << endl;
	}

	OutFile1.close();
}

// ---------------------------------------------------------------------------
double CForest::GetShannon() {
	// map<int,int>::iterator spec_it;

	double shannon = 0, relabund;
	double sum = 0.0;

	int ntrees = 0;
	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		ntrees = ntrees + spec_it->second;
	}

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		if (spec_it->second > 0) {
			relabund = (double) spec_it->second / ntrees;
			shannon += log(relabund) * relabund;
			sum += relabund;
		}
	}

	return(-shannon);
}

// ---------------------------------------------------------------------------
int CForest::GetSAD() {
	// map<int,int>::iterator spec_it;

	double Abund;
	int log2Abund;
	int nspec = 0;


	for (int i = 0; i < MaxSAD; ++i)
		SAD[i] = 0;

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		Abund = spec_it->second;
		if (Abund > 0) {
			++nspec;
			log2Abund = (int)floor(log(Abund) / log(2.0));
			// log2Abund = (int) ceil(log(Abund)/log(2.0));
			if (log2Abund >= MaxSAD)
				++SAD[MaxSAD - 1];
			else
				++SAD[log2Abund];
		}
	}

	return(nspec);
}

// Point pattern functions ---------------------------------------------------

// ---------------------------------------------------------------------------
// wird benutzt um Anteil Kreisbogen in W und Anteil Kreis in W auszurechnen
inline void CForest::fKK(double D1, double D2, double r, double* ra) {
	double D3, W;
	double S1, S2;
	if (D1 < r) {
		D3 = sqrt(r * r - D1 * D1);
		if (D3 < D2) {
			if (D1 == 0)
				S1 = Pi;
			else
				S1 = atan(D2 / D1);
			if (D1 == 0)
				S2 = Pi;
			else
				S2 = atan(D3 / D1);
			W = S1 - S2;
			ra[0] = ra[0] + 0.5 * (D1 * D3 + r * r * W);
			ra[1] = ra[1] + W * r;
		}
		else {
			ra[0] = ra[0] + 0.5 * D1 * D2;
		}
	}
	else {
		if (D1 == 0)
			W = Pi;
		else
			W = atan(D2 / D1);
		ra[0] = ra[0] + 0.5 * r * r * W;
		ra[1] = ra[1] + W * r;
	}
}

// einmal global definieren: ra:array[0..1] of real

// ------------------------------------------------------------------------------
double CForest::FracRingInWin2(double x, double y, double R) {
	double result;
	double a = Xmax, b = Ymax;
	double ra[2];

	if (R == 0)
		result = 1.0;
	else {
		/*
		if ( Xmax < Ymax )
		{
		a = Xmax;
		b = Ymax;
		}
		else
		{
		a = Ymax;
		b = Xmax;
		}
		 */
		ra[0] = 0;
		ra[1] = 0;
		fKK(a - x, b - y, R, ra);
		fKK(a - x, y, R, ra);
		fKK(y, a - x, R, ra);
		fKK(y, x, R, ra);
		fKK(x, y, R, ra);
		fKK(x, b - y, R, ra);
		fKK(b - y, x, R, ra);
		fKK(b - y, a - x, R, ra);
		result = ra[1] / (2 * Pi * R);
	}
	return result;
}

// ------------------------------------------------------------------------------
double CForest::FracCircleInWin2(double x, double y, double R) {
	double result;
	double a = Xmax, b = Ymax;
	double ra[2];

	if (R == 0)
		result = 1.0;
	else {
		/*
		if ( Xmax < Ymax )
		{
		a = Xmax;
		b = Ymax;
		}
		else
		{
		a = Ymax;
		b = Xmax;
		}
		 */
		ra[0] = 0;
		ra[1] = 0;
		fKK(a - x, b - y, R, ra);
		fKK(a - x, y, R, ra);
		fKK(y, a - x, R, ra);
		fKK(y, x, R, ra);
		fKK(x, y, R, ra);
		fKK(x, b - y, R, ra);
		fKK(b - y, x, R, ra);
		fKK(b - y, a - x, R, ra);
		result = ra[0] / (Pi * R * R);
	}
	return result;
}

// ---------------------------------------------------------------------------
// calculates the fraction of the circle with radius r and center x,y in Window
inline double CForest::FracRingInWin(double x, double y, double R) {
	// see Diggle 2003, pg. 51 and Goreaud & Pelissier 1999. J. Veg. Sci.
	double d1 = min(x, Xmax - x);
	double d2 = min(y, Ymax - y);
	double d3 = sqrt(d1 * d1 + d2 * d2);

	double alpha_out;

	if (R < d3) // corner outside the circle
		alpha_out = 2 * acos(min(d1, R) / R) + 2 * acos(min(d2, R) / R);
	else // corner inside
		alpha_out = Pi * 0.5 + acos(d1 / R) + acos(d2 / R);

	return 1.0 - (alpha_out / (2 * Pi));
}

// ---------------------------------------------------------------------------
void CForest::GetPPA() {
	CTree* pTree1;
	CTree* pTree2;
	CSAR_Point* pPoint1;

	// clear all containers
	Mff.clear();
	Mfo.clear();

	nSpecIJ.clear();
	gSpecIJ.clear();
	KSpecIJ.clear();
	ARingI.clear();
	ACircleI.clear();

	for (TreeIterV itree1 = TreeList.begin(); itree1 != TreeList.end(); itree1++) {
		pTree1 = *itree1;
		pTree1->NNDistSpec.clear();
	}

	// generate central points for SAR
	list<CSAR_Point*>PointList;
	double x, y;

	/* random
	int nPoints = 900;
	for (int ip=0; ip<nPoints; ++ip){
	x = RandGen1->Random()*(Xmax-2*Rmax1) + Rmax1;
	y = RandGen1->Random()*(Ymax-2*Rmax1) + Rmax1;
	pPoint1 = new CSAR_Point(x,y);
	PointList.push_back(pPoint1);
	}
	 */

	///*
	// regular grid
	int nPoints = 0;
	double dgrid = 10.0;

	for (x = Rmax1; x < (Xmax - Rmax1); x = x + dgrid) {
		for (y = Rmax1; y < (Ymax - Rmax1); y = y + dgrid) {
			pPoint1 = new CSAR_Point(x, y);
			PointList.push_back(pPoint1);
			nPoints++;
		}
	}
	// */
	// ----------------------
	int Rgrid = (int)ceil(Rmax1 / CellSize);

	map<int, int>::iterator spec_it1;
	map<int, int>::iterator spec_it2;

	for (spec_it1 = SpecAbund.begin(); spec_it1 != SpecAbund.end(); ++spec_it1) {

		if (spec_it1->second >= minAbund) {

			// NennerSpecK[spec_it1->first] = 0.0;

			Mff[spec_it1->first] = 0;
			Mfo[spec_it1->first] = 0;

			//Lf[spec_it1->first] = 0.0;

			for (int ibin2 = 0; ibin2 < nBins2; ibin2++) {

				// Annuli and circle areas
				ARingI[spec_it1->first].push_back(0.0);
				ACircleI[spec_it1->first].push_back(0.0);

				// species-pair-counts
				for (spec_it2 = SpecAbund.begin();
					  spec_it2 != SpecAbund.end(); ++spec_it2)
				{
					if (spec_it2->second >= minAbund) {
						nSpecIJ[spec_it1->first][spec_it2->first].push_back(0);
						gSpecIJ[spec_it1->first][spec_it2->first].push_back(0.0);
						KSpecIJ[spec_it1->first][spec_it2->first].push_back(0.0);
						//DSpecIJ[spec_it1->first][spec_it2->first].push_back(0.0);
					}
				}
			}

			/*
			for (spec_it2 = SpecAbund.begin();
				 spec_it2 != SpecAbund.end(); ++spec_it2)
			{
				if (spec_it2->second >= minAbund) {
					xPOD[spec_it1->first][spec_it2->first] = 0.0;
					NNDistSpecIJ[spec_it1->first][spec_it2->first] = Xmax + Ymax;
				}
			}
			*/
		}

		for (TreeIterV itree1 = TreeList.begin(); itree1 != TreeList.end(); itree1++) {
			pTree1 = *itree1;
			pTree1->NNDistSpec[spec_it1->first] = Xmax + Ymax;
		}

		for (PointIterL ipoint1 = PointList.begin(); ipoint1 != PointList.end();
			ipoint1++) {
			pPoint1 = *ipoint1;
			pPoint1->NNDistSpec[spec_it1->first] = Xmax + Ymax;
		}
	}

	for (int ibin1 = 0; ibin1 < nBins1; ibin1++) {
		CountAll[ibin1] = 0;
		CountCon[ibin1] = 0;
		PropCon[ibin1] = 0.0;
		NennerPCF[ibin1] = 0.0;
		PCF_all[ibin1] = 0.0;
		SAR[ibin1] = 0.0;
	}

	int ibin1, ibin2;
	int iX2, iY2;
	int Xcells2 = Xmax / CellSize;
	int Ycells2 = Ymax / CellSize;

	double xydist, xbound, ybound;
	double nndist;


	//main nested loop from tree i to tree 2
	for (int iX1 = 0; iX1 < Xcells2; ++iX1) {
		for (int iY1 = 0; iY1 < Ycells2; ++iY1) {

			// loop over trees in cell 1
			for (TreeIterL itree1 = Grid[iX1][iY1].TreeList.begin();
				itree1 != Grid[iX1][iY1].TreeList.end(); itree1++) {

				pTree1 = *itree1;

				// loop over neighboring cells up to maximum distance
				for (int diX = -Rgrid; diX <= Rgrid; diX++) {
					for (int diY = -Rgrid; diY <= Rgrid; diY++) {

						iX2 = iX1 + diX;
						iY2 = iY1 + diY;

						xbound = 0;
						ybound = 0;

						// torus edge correction
						// BoundRecGrid(x2,y2,xbound,ybound);

						// no edge correcion
						if ((iX2 < 0 || iY2 < 0) || (iX2 >= XCells) || (iY2 >= YCells))
							continue;

						// loop over trees in cell 2
						for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
							itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

							pTree2 = *itree2;

							if (pTree1 != pTree2) {
								xydist = Distance(pTree1->X, pTree1->Y,
									pTree2->X + xbound, pTree2->Y + ybound);
								if (xydist < Rmax1) {
									ibin1 = (int)floor(xydist / BW1);
									++CountAll[ibin1];

									nndist = pTree1->NNDistSpec[pTree2->SpecID];
									if (xydist < nndist)
										pTree1->NNDistSpec[pTree2->SpecID] = xydist;

									if (pTree1->SpecID == pTree2->SpecID)
										++CountCon[ibin1];
								} // if (d < Rmax1)

								// uni and bivariate patterns
								if ((xydist < Rmax2) &&
									(SpecAbund[pTree1->SpecID] >= minAbund)) {
									if (pTree1->SpecID == pTree2->SpecID)
										++Mff[pTree1->SpecID];
									else
										++Mfo[pTree1->SpecID];

									if (SpecAbund[pTree2->SpecID] >= minAbund) {
										ibin2 = (int)floor(xydist / BW2);
										++nSpecIJ[pTree1->SpecID][pTree2->SpecID][ibin2];
									}

								} // if (d < Rmax2)
							} // if tree1 != tree2

						} // end tree 2

					} // end dy
				} // end dx

			} // end for tree 1

		} // end for y1
	} // end for x1
	// end main loop

	// calculate Wiegand-Moloney edge correction
	for (int ibin1 = 0; ibin1 < nBins1; ibin1++) {
		for (unsigned int itree = 0; itree < TreeList.size(); itree++) {
			pTree1 = TreeList[itree];
				NennerPCF[ibin1] = NennerPCF[ibin1] +
					 (2.0 * Pi * rvec1[ibin1] * BW1) // area of annuli
				* FracRingInWin(pTree1->X, pTree1->Y, rvec1[ibin1]);

				// double ec1 = FracRingInWin(pTree1->X,pTree1->Y,rvec[rbin]);
				// double ec2 = FracRingInWin2(pTree1->X,pTree1->Y,rvec[rbin]);
		}

		PCF_all[ibin1] = (Xmax * Ymax / (NTrees - 1))
			 * CountAll[ibin1] / NennerPCF[ibin1];

		if (CountAll[ibin1] > 0)
			PropCon[ibin1] = (double)CountCon[ibin1] / CountAll[ibin1];
	}

	/*
	for (unsigned int itree=0; itree < TreeList.size(); itree++) {
	pTree1 = TreeList[itree];
	if ((pTree1->X >= x0 && pTree1->X < x1) && (pTree1->Y >= y0 && pTree1->Y < y1)){

	NennerSpecK[pTree1->SpecID] = NennerSpecK[pTree1->SpecID]
	+ FracCircleInWin2(pTree1->X,pTree1->Y,Rmax2);

	for (spec_it = SpecAbund.begin(); spec_it!=SpecAbund.end(); ++spec_it){

	nndist = pTree1->NNDistSpec[spec_it->first];
	if (nndist < Rmax1){
	rbin0 = (int) floor(nndist/BW1);
	for (rbin = rbin0; rbin < nBins1; ++rbin){
	//++pTree1->ISAR[rbin];
	++ISAR[rbin];
	}

	//if (nndist < Rmax2) ++ISAR_f[pTree1->SpecID];
	}
	}

	//	}
	//}

	//for (int rbin=0; rbin< nBins1; rbin++) ISAR[rbin] = (double) ISAR[rbin]/nTrees2;
	 */

	// calculate Lf 50
	/*
	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it) {
		if (spec_it->second >= minAbund) {

			// int abund1 = spec_it->second;
			// int count1 =  Mff[spec_it->first];
			// double nenner1 =  NennerSpecK[spec_it->first];

			// Kf[spec_it->first] = (double) (Xmax*Ymax/(spec_it->second - 1))
			// * Mff[spec_it->first]/NennerSpecK[spec_it->first];

			Lf[spec_it->first] = (double) Mff[spec_it->first] /
				 (Mff[spec_it->first] + Mfo[spec_it->first]);

			// Af[spec_it->first] = Mff[spec_it->first]/(spec_it->second)
			// - (spec_it->second)/(Xmax*Ymax)*Rmax2*Rmax2*Pi;
			// Mf[spec_it->first] = (double) Mff[spec_it->first]/(spec_it->second);
			// Mf[spec_it->first] = (double) Kf[spec_it->first] * (spec_it->second - 1)/(Xmax*Ymax);
			// ISAR_f[spec_it->first] = ISAR_f[spec_it->first]/spec_it->second;

		}
	}
	*/

	//calculate nearest neighbour distances between species
	// and nearest neighbour distribution for species pairs
   int ibin0;

	/*
	for (TreeIterV itree1 = TreeList.begin(); itree1 != TreeList.end(); itree1++) {
		pTree1 = *itree1;
		if (SpecAbund[pTree1->SpecID] >= minAbund) {
			for (spec_it2 = SpecAbund.begin(); spec_it2 != SpecAbund.end(); ++spec_it2) {
				if (spec_it2->second >= minAbund) {

					nndist = pTree1->NNDistSpec[spec_it2->first];
					if (nndist < Rmax2) {
						ibin0 = (int)floor(nndist / BW2);
						for (ibin2 = ibin0; ibin2 < nBins2; ++ibin2)
							++DSpecIJ[pTree1->SpecID][spec_it2->first][ibin2];
					}

					if (pTree1->NNDistSpec[spec_it2->first]
						 <	NNDistSpecIJ[pTree1->SpecID][spec_it2->first])
					{
						NNDistSpecIJ[pTree1->SpecID][spec_it2->first] = pTree1->NNDistSpec[spec_it2->first];
					}
				}
			}
		}
	}
	*/

	// edge correction for bivariate measures
	for (int ibin2 = 0; ibin2 < nBins2; ibin2++) {
		for (unsigned int itree = 0; itree < TreeList.size(); itree++) {
			pTree1 = TreeList[itree];
			if (SpecAbund[pTree1->SpecID] >= minAbund) {

				// double check2 = FracRingInWin2(pTree1->X,pTree1->Y,rvec2[ibin2]);
				// double check3 = FracCircleInWin2(pTree1->X,pTree1->Y,rvec2[ibin2]+0.5*BW2);

				// Annulus area
				ARingI[pTree1->SpecID][ibin2] = ARingI[pTree1->SpecID][ibin2] +
					 (2.0 * Pi * rvec2[ibin2] * BW2) * FracRingInWin(pTree1->X,
					pTree1->Y, rvec2[ibin2]);
				// Circle area
				ACircleI[pTree1->SpecID][ibin2] = ACircleI[pTree1->SpecID][ibin2]
					 + FracCircleInWin2(pTree1->X, pTree1->Y,
					rvec2[ibin2] + 0.5 * BW2);

			} // if Abund >= minAbund
		} // for trees
	} // for ibin

	// calculate bivariate K(r) and g(r)
	int nInd, nPoints2;

	for (spec_it1 = SpecAbund.begin(); spec_it1 != SpecAbund.end(); ++spec_it1) {
		if (spec_it1->second >= minAbund) {
			for (spec_it2 = SpecAbund.begin(); spec_it2 != SpecAbund.end();
				++spec_it2) {
				if (spec_it2->second >= minAbund) {

					for (int ibin2 = 0; ibin2 < nBins2; ++ibin2) {

						nInd = 0;
						for (int ibin3 = 0; ibin3 <= ibin2; ibin3++)
							nInd = nInd + nSpecIJ[spec_it1->first][spec_it2->first][ibin3];

						if (spec_it1->first != spec_it2->first)
							nPoints2 = spec_it2->second;
						else
							nPoints2 = spec_it2->second - 1;

						KSpecIJ[spec_it1->first][spec_it2->first][ibin2] =
							 (Xmax * Ymax / nPoints2) * nInd /
							 ACircleI[spec_it1->first][ibin2];

						gSpecIJ[spec_it1->first][spec_it2->first][ibin2] =
							 (Xmax * Ymax / nPoints2) *
							 nSpecIJ[spec_it1->first][spec_it2->first][ibin2] /
							 ARingI[spec_it1->first][ibin2];

						/*
						if (spec_it1->first != spec_it2->first) {
							if (gSpecIJ[spec_it1->first][spec_it2->first][ibin2] > 1.0e-6)
								xPOD[spec_it1->first][spec_it2->first] +=
								log(gSpecIJ[spec_it1->first][spec_it2->first][ibin2])* BW2;
						}

						DSpecIJ[spec_it1->first][spec_it2->first][ibin2] =
							DSpecIJ[spec_it1->first][spec_it2->first][ibin2] /
							spec_it1->second;
						*/
					}
				} // if Abund2 >= minAbund
			} // for Spec1
		} // if Abund1 >= minAbund
	} // for Spec1
	// */

	// calculate SAR with random points as centers
	// int iX1, iY1;

	// loop over sampling points 1
	for (PointIterL ipoint1 = PointList.begin(); ipoint1 != PointList.end();
		ipoint1++) {

		pPoint1 = *ipoint1;
		int iX1 = (int)pPoint1->X / CellSize;
		int iY1 = (int)pPoint1->Y / CellSize;

		// loop over neighboring cells up to maximum distance
		for (int diX = -Rgrid; diX <= Rgrid; diX++) {
			for (int diY = -Rgrid; diY <= Rgrid; diY++) {

				iX2 = iX1 + diX;
				iY2 = iY1 + diY;

				// no edge correcion
				if ((iX2 < 0 || iY2 < 0) || (iX2 >= XCells) || (iY2 >= YCells))
					continue;

				// loop over trees in cell 2
				for (TreeIterL itree2 = Grid[iX2][iY2].TreeList.begin();
					itree2 != Grid[iX2][iY2].TreeList.end(); itree2++) {

					pTree2 = *itree2;

					xydist = Distance(pPoint1->X, pPoint1->Y, pTree2->X, pTree2->Y);
					if (xydist < Rmax1) {
						nndist = pPoint1->NNDistSpec[pTree2->SpecID];
						if (xydist < nndist)
							pPoint1->NNDistSpec[pTree2->SpecID] = xydist;
					} // if (d < Rmax1)

				} // end tree 2

			} // end dy
		} // end dx

	} // end for point 1

	for (PointIterL ipoint1 = PointList.begin(); ipoint1 != PointList.end();
		ipoint1++) {

		pPoint1 = *ipoint1;

		for (spec_it1 = SpecAbund.begin(); spec_it1 != SpecAbund.end(); ++spec_it1) {
			nndist = pPoint1->NNDistSpec[spec_it1->first];
			if (nndist < Rmax1) {
				ibin0 = (int) floor(nndist / BW1);
				for (ibin1 = ibin0; ibin1 < nBins1; ++ibin1)
					++SAR[ibin1];
			}
		}
	}

	for (PointIterL ipoint1 = PointList.begin(); ipoint1 != PointList.end();
		ipoint1++) {
		pPoint1 = *ipoint1;
		pPoint1->NNDistSpec.clear();
		delete pPoint1;
	}

	for (int ibin1 = 0; ibin1 < nBins1; ibin1++)
		SAR[ibin1] = (double) SAR[ibin1] / nPoints;
}

// ------------------------------------------------------------------------------
void CForest::GetSRLocal(double sq_size, double& m_SR, double& sd_SR) {

	int Xsq = floor(Xmax / sq_size);
	int Ysq = floor(Ymax / sq_size);

	// 3D array

	CSpecSquare** SRgrid;

	//SRgrid = new (CSpecSquare*[Xsq]);
	SRgrid = new CSpecSquare*[Xsq];
	for (int iX = 0; iX < Xsq; iX++)
		//SRgrid[iX] = new (CSpecSquare[Ysq]);
		SRgrid[iX] = new CSpecSquare[Ysq];

	for (int iX = 0; iX < Xsq; ++iX)
		for (int iY = 0; iY < Ysq; ++iY)
			SRgrid[iX][iY].InitspecSquare(iX, iY);

	CTree* pTree1;

	int iX, iY;

	for (unsigned int itree = 0; itree < TreeList.size(); itree++) {
		pTree1 = TreeList[itree];
		iX = (int)floor(pTree1->X / sq_size);
		iY = (int)floor(pTree1->Y / sq_size);
		++SRgrid[iX][iY].Spec[pTree1->SpecID];
	}

	// calculate mean and sd
	int sum_SR = 0;
	for (int iX = 0; iX < Xsq; ++iX)
		for (int iY = 0; iY < Ysq; ++iY)
			sum_SR = sum_SR + SRgrid[iX][iY].Spec.size();

	m_SR = (double)sum_SR / (Xsq * Ysq);

	double SumSq = 0.0;
	for (int iX = 0; iX < Xsq; ++iX)
		for (int iY = 0; iY < Ysq; ++iY)
			SumSq = SumSq + (SRgrid[iX][iY].Spec.size() - m_SR) *
				 (SRgrid[iX][iY].Spec.size() - m_SR);

	sd_SR = sqrt(SumSq / (Xsq * Ysq - 1));

	// free storage
	for (int iX = 0; iX < Xsq; iX++)
		delete[] SRgrid[iX];
	delete[] SRgrid;

}

// ------------------------------------------------------------------------------
void CForest::GetSAR2() {
	for (int i = 0; i < (SAR2_n - 1); ++i) {
		GetSRLocal(sqrt(SAR2_scales[i]), SAR2_m[i], SAR2_sd[i]);
	}

	int nspec = 0;

	for (spec_it = SpecAbund.begin(); spec_it != SpecAbund.end(); ++spec_it)
		if (spec_it->second > 0)
			++nspec;

	SAR2_m[SAR2_n - 1] = nspec;
	SAR2_sd[SAR2_n - 1] = 0.0;
}

// ------------------------------------------------------------------------------
