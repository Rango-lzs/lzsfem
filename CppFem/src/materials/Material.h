
#ifndef MATERIAL_H
#define MATERIAL_H

#include<iostream>
#include<string>
#include"../algebra/Matrix.h"
#include"../algebra/Vector.h"
#include"../tensors/Tensor.h"

/*
* Task: give the matrial constitutive matrix which related the the stress and strain status
* 
*/

class Material
{
	public:
		Material()  = default;
		//SETTERS AND GETTERS
		/*!\brief Returns Consitutive matrix
		 */
	    virtual Matrix<double>& getConstitutiveMatrix() = 0;
		/*!\brief Returns type of material for INSTANTIABLE class
		 */
		virtual std::string getType() = 0;
		
		/*!\brief Receives strain as vector and modifies tensors to 6x1 results
		 *DESCRIPTION: modifies strains and stresses tensors inside the element through reference. The stress is 
					   calculated inside this function as a temporary Vector object then it is resized to a correct 6x1.
					   The strains exist already and are just redistributed but not calculated. See ELEMENT.H method
					   COMPUTETENSORIALRESULTS for more information.
		 *@param[in] v element local strains vector
		 *@param[out] strains 6x1 strain tensor in Voigt notation calculated through \f$&epsilon;=B*v\f$
		 *@param[out] stresses 6x1 stress tensor in Voigt notation calculated through \f$&sigma;=C*&epsilon;\f$
		 */
		virtual void assembleTensors(Vector<double>& v, Tensor& strains, Tensor& stresses)=0;
		friend std::ostream& operator<<(std::ostream &out, Material& mat);
		
		//PLASTIC MODEL FUNCTIONS
	private:
		//These two are pure virtual in elastoplastic model
		//virtual void radialReturn()=0;
		//virtual double yieldFunction()=0;
		
	protected:
		virtual void setConstitutiveMatrix()=0;
	
	protected:
	
		bool plastic;//false upon creation
		std::string type;
		
		//ELASTIC PARAMETERS
		Matrix<double> C;/**<Constitutive matrix for \f$ &sigma;=C&epsilon;\f$*/
		Matrix<double> Cel=Matrix<double>(6,6);/**<Constitutive Matrix for Elastic behavior*/
		double mu;/**<Poisson's ratio */
		double E;/**<Young's modulus*/
		
		
		//PLASTIC PARAMETERS
		Matrix<double> Cep=Matrix<double>(6,6);/**<Constitutive Matrix for plastic behavior*/
		Tensor plasticStrain;
		double yieldStress;/**<ONLY PLASTIC: stress before plastic behavior &sigma;_y*/
		double plasticModulus;/**<ONLY LINEAR ISOTROPIC HARDENING: slope of yield stress increase*/
};
#endif