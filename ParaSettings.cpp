// ---------------------------------------------------------------------------
#include "ParaSettings.h"
#include <fstream>
#include <iostream>
#include <string>
using std::string;

// ---------------------------------------------------------------------------
CModelSettings::CModelSettings()
{
   //technical settings
   nRep = 1;
   nGen = 100;
   steps_out = false;
   R_mode = false;
   cellSize = 5.0;
   //metacommunity
   metaSAD = 0;
   sad_file_name = "Input\\MetacommunitySAD1.txt";
   Jm = 2e6;
   //local community size
   nTrees = 21000;
   Xext = 1000.0;
   Yext = 500.0;
    //point pattern output
   rmax = 100.0;
   bw1 = 1.0;
   bw2 = 5.0;
   minAbund = 50;

   //windows for sampling output
   xmin = 0.0;
   xmax = 1000.0;
   ymin = 0.0;
   ymax = 500.0;
   //habitat map file
   habitat = true;
   //BCI settings
   map_cell_size = 20.0;   map_file_name = "Input\\HabitatMapBCI_Harms_nomixed.txt";
   rel_dens_file_name = "Input\\RelativeDensityBCI_harms_nomixed_n50.txt";
   n_hab_types = 6;
//   //Sinharaja settings
//   map_cell_size = 500.0/26.0;//   map_file_name = "Input\\HabitatMapBCI_Harms_nomixed.txt";
//   rel_dens_file_name = "Input\\RelativeDensityBCI_harms_nomixed_n50.txt";
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
      InFile>>dummy>>sad_file_name;
      InFile>>dummy>>Jm;
      InFile>>dummy>>nTrees;
      InFile>>dummy>>Xext;
      InFile>>dummy>>Yext;
      InFile>>dummy>>rmax;
      InFile>>dummy>>bw1;
      InFile>>dummy>>bw2;
      InFile>>dummy>>minAbund;
      InFile>>dummy>>xmin;
      InFile>>dummy>>xmax;
      InFile>>dummy>>ymin;
      InFile>>dummy>>ymax;
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
CPara::CPara()
{
   theta      = 50;
   metaSR     = 100;
	metaCV     = 1.0;
   m          = 0.1;
   r_max      = 10;
   aRec       = 1;
   aHab       = 0.0;
   aSurv      = 999;
   bSurv      = 0.9;
   m_dm_spec  = 30;
   sd_dm_spec = 0.0;
   m_JCspec    = 1.0;
   sd_JCspec   = 0.0;
}

CPara::CPara(double theta1,
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
             )
{
   theta      = theta1;
   metaSR     = metaSR1;
	metaCV     = metaCV1;
   m          = m1;
   r_max      = r_max1;
   aRec       = aRec1;
   aHab       = aHab1;
   aSurv      = aSurv1;
   bSurv      = bSurv1;
   m_dm_spec  = m_dm_spec1;
   sd_dm_spec = sd_dm_spec1;
   m_JCspec    = m_JCspec1;
   sd_JCspec   = sd_JCspec1;
}

// ---------------------------------------------------------------------------
CPara::~CPara(){};
