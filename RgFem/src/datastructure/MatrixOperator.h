#pragma once
#include "femcore/fem_export.h"

// abstract base class for matrix operators, i.e. a class that can calculate a matrix-vector product
class FEM_EXPORT MatrixOperator
{
public:
	MatrixOperator() {}
	virtual ~MatrixOperator() {}

	// calculate the product Ax = y
	virtual bool mult_vector(double* x, double* y) = 0;
};
