/*****************************************************************//**
 * \file   Bar1d.h
 * \brief  
 * 
 * \author Leizs
 * \date   August 2023
 *********************************************************************/
#ifndef BAR1D_H
#define BAR1D_H	

#include "Matrix.h"
#include "Element.h"

class Bar1d :public Element
{
public:
	Bar1d(Node* n1 = nullptr, Node* n2 = nullptr, Material* mat = nullptr);

	virtual void print()override;

	//This subroutine forms the stiffness matrix of a 1 - d "rod" element.
	virtual Matrix<double>& stiffnessMatrix() override;
	//Setters and getters
	//virtual void setNode(Node* n, int i)override;


private:

	virtual void calculateInternalForce()override;
	virtual void calculateBReducedIntegration()override;
	virtual void calculateJacobian(Matrix<double>& Jacobian)override;
	virtual void calculateG(Matrix<double>& G)override;
	virtual void assembleA(Matrix<double>& A, Matrix<double>& Jacobian)override;
	//DECLARATION OF THE JACOBIAN ELEMENT AND THE DERIVATIVES OF THE SHAPE FUNCTIONS

};
#endif //BAR1D_H