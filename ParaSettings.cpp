// ---------------------------------------------------------------------------
#include "ParaSettings.h"
#include <fstream>
#include <iostream>
#include <string>
using std::string;

// ---------------------------------------------------------------------------
CModelSettings::CModelSettings()
{}

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
             //double aSurv1,
             //double bSurv1,
             double m_dm_spec1,
             double sd_dm_spec1,
             double m_JCspec1,
             double cv_JCspec1,
             double niche_breadth1
             ) :
   theta{theta1},
   Jm{Jm1},
   metaSR{metaSR1},
   metaCV{metaCV1},
   m{m1},
   r_max{r_max1},
   aRec{aRec1},
   //aSurv{aSurv1},
   //bSurv{bSurv1},
   m_dm_spec{m_dm_spec1},
   sd_dm_spec{sd_dm_spec},
   m_JCspec{m_JCspec1},
   cv_JCspec{cv_JCspec1},
   sd_JCspec{m_JCspec * cv_JCspec1},
   niche_breadth{niche_breadth1}
{
}

// ---------------------------------------------------------------------------

