// ---------------------------------------------------------------------------
#include "ParaSettings.h"
#include <fstream>
#include <iostream>
#include <string>
using std::string;

// ---------------------------------------------------------------------------
CModelSettings::CModelSettings()
{
   //sad_file_name = "Input/MetacommunitySAD1.txt";

   map_file_name = "Input/HabitatMapBCI_Harms_nomixed.txt";
   rel_dens_file_name = "Input/RelativeDensityBCI_harms_nomixed_n50.txt";
   n_hab_types = 6;

//   //Sinharaja settings
//   map_cell_size = 500.0/26.0;//   map_file_name = "Input/HabitatMapBCI_Harms_nomixed.txt";
//   rel_dens_file_name = "Input/RelativeDensityBCI_harms_nomixed_n50.txt";
//   n_hab_types = 5;
}

// ---------------------------------------------------------------------------
CModelSettings::~CModelSettings(){};

// ---------------------------------------------------------------------------
void CModelSettings::ReadSettings(std::string file_name)
{
   std::ifstream InFile;
   std::string dummy;

   InFile.open(file_name.c_str());
   if (InFile.good()){
      InFile>>dummy>>nRep;
      InFile>>dummy>>nGen;
      InFile>>dummy>>steps_out;
      InFile>>dummy>>R_mode;
      InFile>>dummy>>cellSize;
      InFile>>dummy>>metaSAD;
      //InFile>>dummy>>sad_file_name;
      //InFile>>dummy>>Jm;
      InFile>>dummy>>nTrees;
      InFile>>dummy>>Xext;
      InFile>>dummy>>Yext;
      InFile>>dummy>>rmax;
      InFile>>dummy>>bw1;
      InFile>>dummy>>bw2;
      InFile>>dummy>>minAbund;
//      InFile>>dummy>>xmin;
//      InFile>>dummy>>xmax;
//      InFile>>dummy>>ymin;
//      InFile>>dummy>>ymax;
      InFile>>dummy>>habitat;
      if (habitat == true){
         InFile>>dummy>>map_cell_size;         InFile>>dummy>>map_file_name;
         InFile>>dummy>>rel_dens_file_name;
         InFile>>dummy>>n_hab_types;
      }
   }
   else std::cout<<"Error settings file"<<std::endl;
   InFile.close();
}

// ---------------------------------------------------------------------------
CPara::CPara(double theta1,
             int    Jm1,
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
             double cv_JCspec1,
             double sigma_comp1
             ) :
   theta{theta1},
   Jm{Jm1},
   metaSR{metaSR1},
	metaCV{metaCV1},
   m{m1},
   r_max{r_max1},
   aRec{aRec1},
   aHab{aHab1},
   aSurv{aSurv1},
   bSurv{bSurv1},
   m_dm_spec{m_dm_spec1},
   sd_dm_spec{sd_dm_spec},
   m_JCspec{m_JCspec1},
   cv_JCspec{cv_JCspec1},
   sd_JCspec{m_JCspec * cv_JCspec1},
   sigma_comp{sigma_comp1}
{
}

// ---------------------------------------------------------------------------

