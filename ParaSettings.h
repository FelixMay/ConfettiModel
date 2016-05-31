
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

//--------------------------------------------------------------------------------
//model settings that apply to several simulation runs and are usually not varied
//in parameter optimization
class CModelSettings
{
public:
   //technical settings
   int nRep;  //number replicates
	int nGen;  //number generations (# complete turnover of community)
	bool steps_out;
	bool R_mode;
	double cellSize;  // size of neighborhood grid

   //metacommunity
   int metaSAD;                  //mode for metacommunity
                                 //0 ... logseries (theta)
                                 //1 ... uniform (metaSR)
                                 //2 ... lognormal (metaSR, metaCV)
                                 //3 ... read from file
                                 //other value... logseries (theta)
   std::string sad_file_name;    //file name
   int Jm;                       //metacommunity size in number of individuals

   //local community size - total extent of the simulated forest
   int nTrees;                   //number of trees in local community
   double Xext;
   double Yext;

   //point pattern output
   double rmax;      //maximum neighborhood radius for point pattern
   double bw1;       //bandwidth for community level patterns
   double bw2;       //bandwidth for species level patterns
   int    minAbund;  //abundance threshold for species level point patterns

   //windows for sampling output
   double xmin;
   double xmax;
   double ymin;
   double ymax;

   //habitat map file
   bool habitat;
   double map_cell_size;
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
	double theta;
	int    metaSR;  //species richness in case of uniform or log-normal metacommunity
	double metaCV;  //cv of abundances for log-normal metacommunity
	double m;
	double r_max;
	double aRec;
	double aHab;   // Habitat sensitivity: 0 ... no habitat effects,
                  // 1 ... as in data
						// > 1 ... strong habitat effects
	double aSurv;
	double bSurv;
	double m_dm_spec;
	double sd_dm_spec;
	double m_JCspec;
	double sd_JCspec;  // factor to calculate heterospecific competition relative to conspecific competition

	CPara();
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
         double sd_JCspec1
        );
	~CPara();
};

//---------------------------------------------------------------------------
class CSpecPara
{
public:
	double muDisp;
	double sigmaDisp;
	//double JCfac;

	std::vector<double> RelHabDens; //relative density in habitat types

	CSpecPara(){};

	CSpecPara(double meanDisp, double sdDisp)
	{
		sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(meanDisp*meanDisp)));
		muDisp = log(meanDisp) - 0.5 * sigmaDisp*sigmaDisp;
		//JCfac = jc;
		RelHabDens.clear();
	};

	~CSpecPara(){RelHabDens.clear();};
};

//---------------------------------------------------------------------------
#endif

