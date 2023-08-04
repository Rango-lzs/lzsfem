#include "InterQuad2d.h"
#include "../elements//Element.h"

InterpQuad2d::InterpQuad2d(Element* elem)
	:elem_(elem)
{

}

Vector<double> InterpQuad2d::getN(double ksi, double eta)
{
	return Vector<double> {(1 - ksi)* (1 - eta), (1 + ksi)* (1 - eta), (1 + ksi)* (1 + eta), (1 - ksi)* (1 + eta)};
}

Matrix<double> InterpQuad2d::getdNds(double ksi, double eta)
{
	Matrix<double> ret;
	ret.matrixResize(2, 4);
	ret(0, 0) = -(1 - eta);
	ret(0, 1) = (1 - eta);
	ret(0, 2) = (1 + eta);
	ret(0, 3) = -(1 + eta);

	ret(1, 0) = -(1 - ksi);
	ret(1, 1) = -(1 + ksi);
	ret(1, 2) = (1 + ksi);
	ret(1, 3) = (1 - eta);

	return ret;
}

//N(ksi,eta)
//x(ksi,eta) = N*xi ,dNdx related to (ksi,eta) and xi
Matrix<double> InterpQuad2d::getdNdx(double r, double s)
{
	// jacob(-1) * dNds
	Matrix<double> jac = jacob(r, s);
	Matrix<double> dnds = getdNds(r, s);
	return jac.inverse() * dnds;
}

Matrix<double> InterpQuad2d::jacob(double r, double s)
{
	//dNds*Coors  2:4 * 4:2
	Matrix<double> coords = getElemCoords();
	Matrix<double> dnds = getdNds(r, s);
	return dnds * coords;
}

Matrix<double> InterpQuad2d::getElemCoords()
{	
	Vector<Node*> nodes = elem_->getElementNodes();
	Matrix<double> ret(nodes.size(), 2);

	for (int i = 0; i < nodes.size(); ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			ret(i, j) = nodes[i]->getPosition()[j];
		}	
	}

	return ret;
}