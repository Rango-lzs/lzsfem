#ifndef ELEMENT_H
#define ELEMENT_H

class Node;
class Material;
class Tensor;
class GaussIntegration;

class Element
{
	public:
		Element() = default;
		~Element() = default;

		//Setters and getters
		virtual std::string elementType() = 0;


	    virtual void setNode(Node* n, int i);
		virtual Vector<Node*> getElementNodes() = 0;

		void setMaterial(Material* mat);
		Material* getMaterial();

		virtual Matrix<double>& getMatrix() = 0;
		virtual Vector<double>& getElementSolutionVector() = 0;
		virtual Vector<double>& getInternalForce() = 0;

		Tensor& getStress();
		Tensor& getStrain();

		//FUNCTIONS THAT CALCULATE
		virtual void calculateMatrix()=0;
		void resizeElementSolutionVector(int n);
		void computeTensorialResults();
		void setNodalValues();
		void setNodalInternalForces();

		//OPERATOR OVERLOAD
		friend std::ostream& operator<<(std::ostream &out, Element& el);
		virtual void print() = 0;
		
	protected:
		//Stresses and strains
		Tensor* stress;/**<Default sized to 6 elements*/
		Tensor* strain;/**<Default sized to 6 elements*/
		Element(Material* mat=nullptr);
		Material* material;/**<Material mode, pointer of non-instantiable class*/

		GaussIntegration* g;
		Matrix<double> K;/**<Stiffness Matrix of the Element, calculation is derived-class-dependent*/
		Matrix<double> B;/**<B-operator matrix*/

		Vector<double> internalForce;/**<Internal force vector of the element*/
		Vector<Node*> nodes;/**<Stores nodes that compose the elment, size is type-dependant*/
		Vector<double> solution;/**<Stores values of solutions for each nodes and number of DOFs per node*/		
};
#endif