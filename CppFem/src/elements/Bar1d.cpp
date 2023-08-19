/*****************************************************************//**
 * \file   Bar1d.cpp
 * \brief  
 * 
 * \author Leizs
 * \date   August 2023
 *********************************************************************/

#include "Bar1d.h"

Matrix<double>& Bar1d::stiffnessMatrix()
{
	Matrix<double> km(2, 2);
	double ea = 1e5;
	double length = 1.0;
	double one = 1.0;
	km(1, 1) = one;
	km(2, 2) = one;
	km(1, 2) = -one;
	km(2, 1) = -one;
	km *= ea / length;
}
