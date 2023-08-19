#ifndef QUAD2_H
#define QUAD2_H	

#include"Element.h"

class InterpQuad2d;

class Quad2 :public Element
{
public:
	Quad2(Node* n1 = nullptr, Node* n2 = nullptr, Node* n3 = nullptr, Node* n4 = nullptr, Material* mat = nullptr);

	virtual void print()override;
	/*!\brief Stiffness matrix (K) generation
	 *We assume that since we apply Gauss Integration on the matrices which are a fuction of ksi and epsilon, then we just have
	 *to implement the multiplication of all the matrices to get K, as seen on steps. Everything using a 1-point rule
	 *In Steps: 1.-calculate Jacobian and detJ
	 *			2.-calculate A matrix
	 *			3.-Calculate G
	 *			4.-calculate B=GA
	 *			5.-implement K=detJ*B_T*D*B
	 */
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

private:	
	InterpQuad2d* inter;
};
#endif