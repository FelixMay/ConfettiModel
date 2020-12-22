
// File with small helper classes
#ifndef ParaSettingsH
#define ParaSettingsH

#include <map>
#include <vector>
#include <string>
#include <cmath>

//---------------------------------------------------------------------------
//! \brief Helper class for SAR
//! This class is needed for the calculation of
//! quadrat-based species-area relationships (SAR)
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
//! Class for constant model settings
//!
//! Model settings that apply to several simulation runs and are usually not varied
//! in parameter optimization or brute-force simulations
class CModelSettings
{
public:
   //technical settings
   int nRep = 1;  //number replicates
	int nGen = 100;  //number generations (# complete turnover of community)
	bool steps_out = false;
	bool R_mode = false;
	double cellSize = 5.0;  // size of neighborhood grid in meters

   //metacommunity
   int metaSAD = 0;              //mode for metacommunity
                                 //0 ... logseries (theta, Jm)
                                 //1 ... uniform (metaSR)
                                 //2 ... lognormal (metaSR, metaCV)
                                 //3 ... read from file
                                 //other value... logseries (theta, Jm)
   //std::string sad_file_name;    //file name
   //int Jm = 2000000;             //metacommunity size in number of individuals

   //Default parameters
   //local community size - total extent of the simulated forest
   int nTrees = 21000;           //number of trees in local community
   double Xext = 1000.0;
   double Yext = 500.0;

   //point pattern output
   double rmax = 100.0;     //maximum neighborhood radius for point pattern
   double bw1 = 50.0;       //bandwidth for community level patterns
   double bw2 = 50.0;       //bandwidth for species level patterns
   int    minAbund = 50;    //abundance threshold for species level point patterns

   CModelSettings();
   ~CModelSettings();

   //! \brief Read simulation settings from file
   //!
   //! \param  file_name Name of the input file
   //!
   void ReadSettings(std::string file_name);
};

//---------------------------------------------------------------------------
//! \brief Class with model parameters
//!
//! Model parameters with ecological interpretation that
//! are varied in simulation sets or parameter optimization
//!
class CPara
{
public:
	double theta = 50.0;        //!> Fundamental biodiversity number following Hubbell 2001
	int    Jm    = 2000000;     //!> Metacommunity size in number of individuals
	int    metaSR = 400;        //!> Metacommunity species richness in case of uniform or log-normal metacommunity
	double metaCV = 1.0;        //!> Coefficient of Variation (CV) of abundances for log-normal metacommunity
	double m = 0.1;             //!> Immigration rate
	double r_max = 10.0;        //!> Neighbourhood radius for tree-tree interactions
	double aRec = 0.005;        //!> Model parameter for relationship between competition and recruitment probability
	//double aSurv = 999.0;       //!> Model parameter for relationship between competition and survival (not used in current version!)
	//double bSurv = 0.89;        //!> Survival rate without competition (not used in current version!)
	double m_dm_spec = 30.0;    //!> Mean dispersal distance of all species /(no interspecific variation in dispersal in the current model version)
	double sd_dm_spec = 0.0;    //!> Standard deviation of dispersal distance
	double m_JCspec = 1.0;      //!> Ratio of conspecific relative to heterospecific competition (= 1 when CNDD = HNDD)
	double cv_JCspec = 0.0;     //!> Coefficient of Variation among CNDD of species
	double sd_JCspec = 0.0;     //!> Standard deviation of interspecific variation in CNDD = mJC_spec * cv_JCspec
	double niche_breadth = 1.0; //!> Niche breadth (see Gravel et al. 2005 EcolLett, Eq. 3)
                               //Niche breadth i equal for all species in the current version, but could be species-specific
   int    n_hills = 1;         //!> Number of hills in the landscape in the cosinus-function landscapes

	//double sigma_comp = 0.0; // niche width for species competition, see Scheffer and van Nes (2006) PNAS Eq. 4

	CPara(){};
	CPara(double theta1,
         int Jm1,
         int    metaSR1,
         double metaCV1,
         double m1,
         double r_max1,
         double aRec1,
         //double aSurv1,
         //double bSurv1,
         double m_dm_spec1,
         double sd_dm_spec1,
         double m_JCspec1,
         double cv_JCspec1,
         double niche_breadth1,
         int    n_hills1
         //double sigma_comp1
        );
	~CPara(){};
};

//---------------------------------------------------------------------------
//! \brief Class with species-specific parameters
//!

class CSpecPara
{
public:
	//double meanDisp; //mean dispersal distance

	double muDisp;    //!< mu parameter of log-normal dispersal kernel
	double sigmaDisp; //!< sigma parameter of log-normal dispersal kernel

	double muEnvir;  //!< Species-specific environmental optimum, see Gravel et al. 2006 Ecology Letters, Eq. 3


	CSpecPara(){};

	CSpecPara(double mDisp, double sdDisp, double mu_opt) :  muEnvir(mu_opt)
	{
		sigmaDisp = sqrt(log(1.0 + (sdDisp*sdDisp)/(mDisp*mDisp)));
		muDisp = log(mDisp) - 0.5 * sigmaDisp*sigmaDisp;
	};

	~CSpecPara(){};
};

//---------------------------------------------------------------------------
#endif

