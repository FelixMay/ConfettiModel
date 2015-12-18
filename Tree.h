//---------------------------------------------------------------------------

#ifndef TreeH
#define TreeH

#include <math.h>
#include <vector>
#include <map>
#include "Cell.h"

//---------------------------------------------------------------------------
class CSAR_Point
{
public:
	double X;
	double Y;

	std::map<int,double> NNDistSpec; //nearest neighbour distance to each species

	CSAR_Point(double x, double y)
	{
		X = x; Y = y;
		NNDistSpec.clear();
	};

	~CSAR_Point()
	{
		NNDistSpec.clear();
	};
};

//---------------------------------------------------------------------------
class CTree
{
public:
	unsigned int TreeID;
	double X;
	double Y;
	unsigned int SpecID;

	double R; //ZOI Radius
	double NCI; // neighbourhood crowding index
	double pSurv;  //Survival probability

	std::list <CCell*> CellList; //list of grid cells overlapped by the tree
	std::map<int,double> NNDistSpec; //nearest neighbour distance to each species

	CTree(unsigned int id, double x, double y, int spec, double r, double p_surv)
	{
		TreeID = id;
		X = x; Y = y;
		SpecID = spec;
		R = r;
		pSurv = p_surv;
		NNDistSpec.clear();
	};

	~CTree()
	{
		CellList.clear();
		NNDistSpec.clear();
	};

	/*
	void GetPSurv(double a, double b)
	{
		if (A_zoi == 0.0) pSurv = b;
		else pSurv = b - b*(A_comp/A_zoi)/(a + A_comp/A_zoi);
	}
	*/

	void GetPSurv2(double a, double b)
	{
		if (NCI == 0.0) pSurv = b;
		else pSurv = b - b*(NCI)/(a + NCI);
	}
};


typedef std::vector<CTree*>::iterator TreeIterV;
typedef std::list<CTree*>::iterator TreeIterL;
typedef std::list<CSAR_Point*>::iterator PointIterL;

//---------------------------------------------------------------------------
#endif
