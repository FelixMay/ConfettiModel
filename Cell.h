//---------------------------------------------------------------------------

#ifndef CellH
#define CellH

#include <list>

//---------------------------------------------------------------------------
const double Pi = 3.14159265358979323846;

class CTree;

class CCell
{
public:
	int X;
	int Y;
	int HabitatType;

	//int nTreesOverlap; //number of trees overlapping the cell
	//double ProbRecruit;

	std::list<CTree*> TreeList; //list of trees with STEM in the cell

	CCell(){
		X = 0;
		Y = 0;
		HabitatType = - 99;

		//nTreesOverlap = 0;
		//ProbRecruit = 0;
	};

	~CCell() {TreeList.clear();};

	void InitCell(int x, int y, int hab_type)
	{
		X=x;
		Y=y;
		HabitatType = hab_type;
		//nTreesOverlap = 0;
		//ProbRecruit = 0;
		TreeList.clear();
	};

	/*
	double GetProbRecruit(double a)
	{
		if (nTreesOverlap == 0)
			return(1.0);
		else
			return (1.0 - nTreesOverlap/(a + nTreesOverlap));
	};
	*/
};

typedef std::list<CCell*>::iterator CellIterL;


//---------------------------------------------------------------------------
#endif
