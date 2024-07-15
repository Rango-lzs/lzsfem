#ifndef _INTER_POLATION_HH 
#define _INTER_POLATION_HH

#include "../algebra/Matrix.h"

class Element;

class Interpolation 
{
public:
	Interpolation(Element* elem);
	virtual ~Interpolation();

    /**
	 * get the shape function Ni(i =1 : nnode)
	 */
	virtual Vector<double> getN(double ksi, double eta) = 0;

    /**
	 * get the shape function derivation to element coordinate
	 * return matrix<2,nnode>
	 */
	virtual Matrix<double>  getdNds(double ksi, double eta) = 0;

    /**
	 * get the shape function derivation to natural coordinate
	 * return matrix<2,nnode>
	 */
	virtual Matrix<double> getdNdx(double r, double s) = 0;

	/**
	 *
	 *@return Matrix<2,2>
	 */
	virtual Matrix<double> jacob(double r, double s) = 0;

private:
	Matrix<double> getElemCoords();
	Element* elem_;
};

#endif //_INTER_POLATION_HH