#ifndef ELASTOPLASTICMATERIAL
#define ELASTOPLASTICMATERIAL

#include"Material.h"
#include"../solver/descent/ConjugateGradientDescent.h"
#include"../solver/LinearIterativeSolver.h"
#include"../solver/Jacobi.h"

class ElastoPlasticMaterial:public Material
{
	protected:
		ElastoPlasticMaterial(double E, double mu, double yS, double pM);

		double getPoisson();
		double getModulus();
		double getYieldStress();
		double getPlasticModulus();
		Tensor& getPlasticStrain();
		bool isPlastic();
		
		//PLASTIC MODEL FUNCTIONS: these are to be used within assembleTensors function
	private:
		/*\brief Implements elastic predictor-plastic corrector algorithm
		 * The outcome of this algorithm must be the plastic strains and the updating of the stress state
		 param@[in]
		 */
		virtual Tensor radialReturn(Tensor& strains)=0;
		/*!\brief Returns the result of evaluating the yield function on a stress state
		 */
		virtual double yieldFunction(Tensor& stress)=0;
		
		
};
#endif