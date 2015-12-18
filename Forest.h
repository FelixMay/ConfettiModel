//---------------------------------------------------------------------------

#ifndef ForestH
#define ForestH

#include "Tree.h"
#include "randomc.h"
//#include "randoma\randoma.h"
#include "stocc.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>

//---------------------------------------------------------------------------
class CSpecSquare
{
public:
	int iX;
	int iY;
	std::map<int,int> Spec;

	CSpecSquare(){
		iX = 0;
		iY = 0;
	};

	~CSpecSquare() {Spec.clear();};

	void InitspecSquare(int x, int y)
	{
		iX=x;
		iY=y;
		Spec.clear();
	};
};

//---------------------------------------------------------------------------
class CPara
{
public:
	int    NTrees;
	int    Jmeta;
	double theta;
	double m;
	double r_max;
	double aRec;
	double aHab;   // Habitat sensitivity: 0 ... no habitat effects,
                  //1 ... as in data
						// > 1 ... strong habitat effects
	double aSurv;
	double bSurv;
	double m_dm_spec;
	double sd_dm_spec;
	double m_JCspec;
	double sd_JCspec;  // factor to calculate heterospecific competition relative to conspecific competition

	CPara(){NTrees = 21000;
			Jmeta = 2e6;
			theta = 50;
			m     = 0.1;
			r_max = 20;
			aRec = 1000;
			aHab = 0.0;
			aSurv = 1000;
			bSurv = 0.9;
			m_dm_spec = 30;
			sd_dm_spec = 10;
			m_JCspec = 1.0;
			sd_JCspec = 0.0;
		   }

	CPara(int    ntrees,
			int    jm,
			double theta1,
			double m1,
			double rmax1,
			double a_rec,
			double a_hab,
			double a_surv,
			double b_surv,
			double m_disp,
			double sd_disp,
			double m_jc,
			double sd_jc
		  )
		  {NTrees = ntrees;
			Jmeta = jm;
			theta = theta1;
			m = m1;
			r_max = rmax1;
			aRec = a_rec;
			aHab = a_hab;
			aSurv = a_surv;
			bSurv = b_surv;
			m_dm_spec = m_disp;
			sd_dm_spec = sd_disp;
			m_JCspec = m_jc;
			sd_JCspec = sd_jc;
		  }

	~CPara(){};
};

//---------------------------------------------------------------------------
class CSpecPara
{
public:
	double muDisp;
	double sigmaDisp;
	double JCfac;

	std::vector<double> RelHabDens; //relative density in five habitat types

	CSpecPara(){};

	CSpecPara(double meanDisp,
			  double sdDisp,
			  double jc,
			  std::vector<double> rel_dens
			  )
	{
		sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
		muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;
		JCfac = jc;

		RelHabDens = rel_dens;
	};

	~CSpecPara(){};
};

//---------------------------------------------------------------------------
class CForest
{
public:
	CPara* Pars;

	//Immigration rate
	double m;

	//Landscape
	double Xmax;
	double Ymax;

	//Trees
	int NTrees;
	std::vector<CTree*> TreeList;
	unsigned int TreeID;

	//Species
	unsigned int SpecMax;
	std::vector<double> CumRelAbundMeta;
	//std::vector<double> CumProbImmi;

	std::map<int,int> SpecAbund;  // map with first --> key, second --> abund
	std::map<int,CSpecPara> SpecPars;

	double** InteractMat;      //matrix with species-specific interaction coefficients

	//Runtime
	int BD_max;

	//Grid
	int CellSize;
	int XCells;
	int YCells;
	int grid_steps;

	CCell** Grid;

	//HabitatMap
	double MapCellSize;
	int MapXcells;
	int MapYcells;
	int** Map;
	int nHabTypes;

	//relative habitat densities from the field data
	std::vector<std::vector<double> > RelHabDensData;

	//random number generators
	CRandomMersenne* RandGen1;
	StochasticLib1* RandGen2;

	//Point pattern variables all trees
	double Rmax1;
	double BW1;
	double* rvec1;
	int nBins1;

	int* CountAll;
	int* CountCon;
	double* PropCon;
	double* PCF_all;
	double* NennerPCF;
	double* SAR;

	// Point pattern variables single species
	int minAbund;  //abundance threshold for calculation
	               //of bivariate point patterns

	double Rmax2;  // radius for species-specific point patterns
	double BW2;
	double* rvec2;
	int nBins2;

	std::map<int,int>::iterator spec_it;

	//std::map <int,double> NennerSpecK;
	std::map <int,int> Mff;
	std::map <int,int> Mfo;
	std::map <int,double> Lf;

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
	std::vector<double> SAR2_scales; // square side length
	int SAR2_n;      //number of scales;
	double* SAR2_m;  // mean species richness
	double* SAR2_sd; // sd species richness

	//index variables
	int isim;
	int irep;

	//Files
	std::fstream AbundFile;
	std::fstream SAD_File;
	std::fstream DivFile;
	std::fstream PCF_File;
	std::fstream PropConFile;
	std::fstream Lf20_File;
	std::fstream Kcon20_File;
	std::fstream Khet20_File;
	std::fstream SAR1_File;
	std::fstream SAR2_File;
	//std::fstream xPOD_File;
	//std::fstream CrossPCF_File;
	//std::fstream NNDist_File;

	//std::ofstream TestFile;

	//Output Variables
	int BD_5years;   //counter for birth-death events per loop
	int BD_total;
	//double Shannon;
	static const int MaxSAD = 12;
	int SAD[MaxSAD];    //Species abundance distributions as octave curve 2^0 - 2^11

	//Private functions
	inline int GetRandSpec(); //draw a species from the species pool

	inline double Overlap2(CTree* pTree1, CTree* pTree2, double d);
	inline double Distance(double x1, double y1, double x2, double y2);

	inline void PeriodBound(double& xx, double& yy);
	inline void BoundIntGrid(int& xx, int& yy, int Xmax, int Ymax);
	inline void BoundGrid(int& xx, int& yy, double& xb, double& yb);

	inline double FracRingInWin(double x, double y, double R);
	inline void fKK( double D1, double D2, double r, double* ra);

	double FracRingInWin2( double x, double y, double R);
	double FracCircleInWin2( double x, double y, double R);

	void GetSRLocal(double sq_size, double& m_SR, double& sd_SR);

	double GetProbRecruit(double x0, double y0, unsigned int spec_id);

	void AddTree(CTree* pTree1);
	void RemoveTree(CTree* pTree1);

	//int GetSpecMaxAbund();
	void GetNewXY(double &x1, double &y1, int idspec);
	std::vector<double> SeqConstruct(unsigned int J1, double theta1, double m1 = 1.0);

//public:
	CForest(int seed,
			double xmax, double ymax,
			double map_cell,
			char* map_file,
			char* rel_dens_file,
			int n_hab_types
			);
	~CForest();

	//void SetPars(CPara* pars);
	void FileOpen(std::string label);
	void Init();
	void BirthDeathLoop();
	bool BirthDeathAsync();
	void UpdateTrees();
	void ClearForest();
	void Loop1(int NSteps, bool output);
	void WriteTrees();
	void WriteTreesTime(int tstep);
	void WriteSpecPar();
	void WriteOutput(int istep, int isim, int irep);
	void GetPPA();
	double GetShannon();
	int GetSAD();
	void GetSAR2();
	void OneRun(int isim, int irep, int ngen, bool steps_out, bool r_mode);

	std::string IntToString(int i)
	{
	  std::ostringstream os;
	  os<<i;
	  return os.str();
	}

};
//---------------------------------------------------------------------------
#endif
