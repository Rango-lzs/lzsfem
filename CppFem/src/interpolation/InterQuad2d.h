#ifndef _INTER_POLATION_QUAD_2D_HH 
#define _INTER_POLATION_QUAD_2D_HH

#include "Interpolation.h"
#include "../algebra/Matrix.h"

class Element;

class InterpQuad2d: public Interpolation
{
public:
	InterpQuad2d(Element* elem);

    /**
	 * get the shape function Ni(i =1 : nnode)
	 */
	Vector<double> getN(double ksi, double eta);

    /**
	 * get the shape function derivation to element coordinate
	 * return matrix<2,nnode>
	 */
	Matrix<double>  getdNds(double ksi, double eta);

    /**
	 * get the shape function derivation to natural coordinate
	 * return matrix<2,nnode>
	 */
	Matrix<double> getdNdx(double r, double s);

	/**
	 *
	 *@return Matrix<2,2>
	 */
	Matrix<double> jacob(double r, double s);

private:
	Matrix<double> getElemCoords();

	Element* elem_;
};

#endif