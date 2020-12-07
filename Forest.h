//---------------------------------------------------------------------------

#ifndef ForestH
#define ForestH

#include "Tree.h"
#include "ParaSettings.h"
#include "randomc.h"
//#include "randoma\randoma.h"
#include "stocc.h"
#include <fstream>
#include <sstream>
#include <deque>

//---------------------------------------------------------------------------
class CForest
{
public:
	CPara *pPars;
	CModelSettings *pSettings;

	//Immigration rate
	double m;

	//Landscape
	double Xmax;
	double Ymax;

	//Trees
	int NTrees;
	std::vector<CTree*> TreeList;
	int64_t TreeID;

	//Species
	int SpecMax;
	std::vector<double> CumRelAbundMeta;
	//std::vector<double> CumProbImmi;

	std::map<int,int> SpecAbund;     // map with first --> key, second --> abund
	//std::map<int,int> SpecAbundWin;  // abundances only in sampling window
	std::map<int,CSpecPara> SpecPars;

	double **InteractMat;      //matrix with species-specific interaction coefficients

	//Runtime
	int64_t BD_max;

	//Grid
	int XCells;
	int YCells;
	int grid_steps;

	CCell **Grid;

	//HabitatMap
	double MapCellSize;
	int MapXcells;
	int MapYcells;
	int **Map;
	int nHabTypes;

	int smooth_habitat; //number of neighborhood grid cells for habitat proportion calculation
	                    // 1 ... 9 neighbors
	                    // 2 ... 25 neighbors

   double ***HabitatProp; //Habitat grid with proportions of nHabTypes habitats
                          //Dimensions MapXCells x MapYCells x nHabTypes

	//relative habitat densities from the field data
	std::vector<std::vector<double> > RelHabDensData;

	//random number generators
	CRandomMersenne *RandGen1;
	StochasticLib1 *RandGen2;

	//Point pattern variables
	int nBins1;
	int nBins2;
	double *rvec1;
	double *rvec2;

	int *CountAll;
	int *CountCon;
	double *PropCon;
	double *PCF_all;
	double *DenomPCF;
	double *SAR;

	// Point pattern variables single species
	std::map<int,int>::iterator spec_it;

	std::map <int,int> Mff;
	std::map <int,int> Mfo;
	//std::map <int,double> Lf;

	std::map<int, std::map<int, std::vector<int> > > nSpecIJ;  //count of individuals
	std::map<int, std::map<int, std::vector<double> > > gSpecIJ; // pair-correlation function
	std::map<int, std::map<int, std::vector<double> > > KSpecIJ;
	//std::map<int, std::map<int, std::vector<double> > > DSpecIJ; //nearest neighbour distribution function
																 //proportion of points of i with neighbour of j within r

	//std::map<int, std::map<int, double > > xPOD;    //integral of log(g_ij)
	//std::map<int, std::map<int, double > > NNDistSpecIJ;  //nearest neighbour distance

	std::map<int,std::vector<double> > ARingI;    // area of annuli
	std::map<int,std::vector<double> > ACircleI;  // area of circles

	// Variables for quadrat-based SAR
	std::vector<double> SARq_scales; // square side length
	int SARq_n;      //number of scales;
	double *SARq_m;  // mean species richness
	double *SARq_sd; // sd species richness

	//index variables
	int isim;
	int irep;

	//Files
	std::fstream AbundFile;
	std::fstream SAD_File;
	std::fstream DivFile;
	std::fstream PCF_File;
	std::fstream PropConFile;
	std::fstream SARq_File;   //quadrat-based SAR

	//std::fstream SARp_File; //point-based SAR
   //	std::fstream Lf20_File;
   //	std::fstream Kcon20_File;
   //	std::fstream Khet20_File;

	//std::fstream xPOD_File;
	//std::fstream CrossPCF_File;
	//std::fstream NNDist_File;

	//std::ofstream TestFile;

	//Output Variables
	int BD_5years;   //counter for birth-death events per loop
	int64_t BD_total;

	static const int MaxSAD = 12;
	int SAD[MaxSAD];      //Species abundance distributions as octave curve 2^0 - 2^11
	int SAD_Win[MaxSAD];  //SAD in sampling window

	//Private functions
	inline int GetRandSpec(); //draw a species from the species pool

	//inline double Overlap2(CTree* pTree1, CTree* pTree2, double d);
	inline double Distance(double x1, double y1, double x2, double y2);

	inline void PeriodBound(double& xx, double& yy);
	inline void BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax);
	inline void BoundGrid(int& xx, int& yy, double& xb, double& yb);

	inline double FracRingInWin(double x, double y, double R);
	inline void fKK( double D1, double D2, double r, double* ra);

	double FracRingInWin2( double x, double y, double R);
	double FracCircleInWin2( double x, double y, double R);

	void GetSRLocal(double sq_size, double& m_SR, double& sd_SR);

	double GetProbRecruit(double x0, double y0, int spec_id);

	void AddTree(CTree* pTree1);
	void RemoveTree(CTree* pTree1);

	//int GetSpecMaxAbund();
	void GetNewXY(double &x1, double &y1, int idspec);
	std::vector<double> SeqConstruct(int J1, double theta1, double m1 = 1.0); //log-series SAD
	std::vector<double> UniformSAD(int nSpecies);
	std::vector<double> LognormSAD(int nSpecies, int nIndividuals, double cv_abund);

//public:
	CForest(int seed, CModelSettings* pset);
	~CForest();

	void FileOpen(std::string label);
	void CreateHabitatMap();
	void ReadSADFile();
	void initTrees();
	void initSpecies();
	bool BirthDeathAsync();
	void UpdateTrees();
	void clearTrees();
	void clearSpecies();
	//void Loop1(int NSteps, bool output);
	void WriteTrees(int isim, int irep, int istep);
	void writeSpecies(int isim, int irep = 0);
	void writeInteractMat(int isim, int irep = 0);
	void WriteOutput(int isim, int irep);
	void GetPPA();
	double GetShannon();
	void GetDiversity(int &nspec, double &shannon, double &simpson);
	double getQueueCV(std::deque<double> &queue);
	int GetSAD();
	void GetSARq();
	void OneRun(int isim, int irep);

	std::string IntToString(int i)
	{
	  std::ostringstream os;
	  os<<i;
	  return os.str();
	}

};
//---------------------------------------------------------------------------
#endif
