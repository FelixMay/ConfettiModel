
// File with small helper classes
#ifndef ParaSettingsH
#define ParaSettingsH

#include <map>
#include <vector>
#include <string>
#include <cmath>

//---------------------------------------------------------------------------
class CSpecSquare
{
public:
	int iX = 0;
	int iY = 0;
	std::map<int,int> Spec;

	CSpecSquare(){};

	~CSpecSquare() {Spec.clear();};

	void InitspecSquare(int x, int y)
	{
		iX=x;
		iY=y;
		Spec.clear();
	};
};

//--------------------------------------------------------------------------------
//model settings that apply to several simulation runs and are usually not varied
//in parameter optimization
class CModelSettings
{
public:
   //technical settings
   int nRep = 1;  //number replicates
	int nGen = 100;  //number generations (# complete turnover of community)
	bool steps_out = false;
	bool R_mode = false;
	double cellSize = 5.0;  // size of neighborhood grid

   //metacommunity
   int metaSAD = 0;              //mode for metacommunity
                                 //0 ... logseries (theta)
                                 //1 ... uniform (metaSR)
                                 //2 ... lognormal (metaSR, metaCV)
                                 //3 ... read from file
                                 //other value... logseries (theta)
   std::string sad_file_name;    //file name
   int Jm = 2000000;             //metacommunity size in number of individuals

   //local community size - total extent of the simulated forest
   int nTrees = 21000;           //number of trees in local community
   double Xext = 1000.0;
   double Yext = 500.0;

   //point pattern output
   double rmax = 100.0;     //maximum neighborhood radius for point pattern
   double bw1 = 50.0;       //bandwidth for community level patterns
   double bw2 = 50.0;       //bandwidth for species level patterns
   int    minAbund = 50;    //abundance threshold for species level point patterns

//   //windows for sampling output
//   double xmin;
//   double xmax;
//   double ymin;
//   double ymax;

   //habitat map file
   bool habitat = false;
   double map_cell_size = 20;
   std::string map_file_name;
	std::string rel_dens_file_name;
	int n_hab_types;

   CModelSettings();
   ~CModelSettings();
   void ReadSettings(std::string file_name);
};

//---------------------------------------------------------------------------
//parameters that are varied in optimization approaches
class CPara
{
public:
	double theta = 50.0;
	int    metaSR = 400;  //species richness in case of uniform or log-normal metacommunity
	double metaCV = 1.0;  //cv of abundances for log-normal metacommunity
	double m = 0.1;
	double r_max = 10.0;
	double aRec = 0.005;
	double aHab = 1.0;   // Habitat sensitivity: 0 ... no habitat effects,
                        // 1 ... as in data
						      // > 1 ... strong habitat effects
	double aSurv = 999.0;
	double bSurv = 0.89;
	double m_dm_spec = 30.0;
	double sd_dm_spec = 0.0;
	double m_JCspec = 1.0; // factor to calculate heterospecific competition relative to conspecific competition
	double cv_JCspec = 0.0;
	double sd_JCspec = 0.0;

	CPara(){};
	CPara(double theta1,
         int    metaSR1,
         double metaCV1,
         double m1,
         double r_max1,
         double aRec1,
         double aHab1,
         double aSurv1,
         double bSurv1,
         double m_dm_spec1,
         double sd_dm_spec1,
         double m_JCspec1,
         double cv_JCspec1
        );
	~CPara(){};
};

//---------------------------------------------------------------------------
class CSpecPara
{
public:
	double meanDisp; //mean dispersal distance

	double muDisp;   //parameters of log-normal dispersal kernel
	double sigmaDisp;

	double comp_trait;

	//double pRec;     //recruitment probability without competition

	std::vector<double> RelHabDens; //relative density in habitat types

	CSpecPara(){};

	CSpecPara(double mDisp, double sdDisp, double trait_val) : meanDisp{mDisp}, comp_trait{trait_val}
	{
		sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
		muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;
		RelHabDens.clear();
	};

	~CSpecPara(){RelHabDens.clear();};
};

//---------------------------------------------------------------------------
#endif

