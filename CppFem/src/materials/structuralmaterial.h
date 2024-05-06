/*****************************************************************//**
 * \file   structuralmaterial.h
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#ifndef structuralmaterial_h
#define structuralmaterial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "mathfem.h"
#include "matconst.h"
#include "matstatus.h"
#include "valuemodetype.h"
#include <vector>

 ///@name Input fields for StructuralMaterial
 //@{
#define _IFT_StructuralMaterial_referencetemperature "referencetemperature"
#define _IFT_StructuralMaterial_talpha "talpha"
#define _IFT_StructuralMaterial_StressControl_stiffmode "scstiffmode"
#define _IFT_StructuralMaterial_StressControl_reltol "screltol"
#define _IFT_StructuralMaterial_StressControl_abstol "scabstol"
#define _IFT_StructuralMaterial_StressControl_maxiter "maxiter"
//@}

namespace fem
{
#define STRAIN_STEPS 10.0

	class GaussPoint;
	///@todo Update documentation
	/**
	 * Abstract base class for all "structural" constitutive models. It declares common  services provided
	 * by all structural material models. The implementation of these services is partly left on derived classes,
	 * which will implement constitutive model dependent part.
	 * Some general purpose services are implemented on this level. For details, how to store
	 * material model related history variables in integration points, see base class @ref Material documentation.
	 *
	 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
	 * Its capabilities can be examined using hasMaterialModeCapability  service.
	 * It is generally assumed, that results obtained from constitutive model services are according to
	 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
	 * services.
	 * Structural material introduces several stress/strain modes.
	 * Full and reduced formats of stress/strain vectors are introduced for convenience.
	 * The full format includes all components, even if they are zero due to stress/strain mode nature,
	 * but in the reduced format, only generally nonzero components are stored.
	 * (full format must used only if absolutely necessary, to avoid wasting of space. It is used
	 * by output routines to print results in general form). Methods for converting vectors between
	 * full and reduced format are provided.
	 *
	 * If in particular mode particular stress component is zero, the corresponding strain is not computed
	 * and not stored in reduced vector, and in full vector there is zero value on corresponding position.
	 * On the other hand, if some zero strain is imposed,
	 * On the other hand, if zero strain component is imposed this condition must be taken into account in geometrical
	 * relations (at element level), and corresponding component are included stress/strain reduced vectors.
	 *
	 * Structural material introduces following basic stress/strain modes
	 * - 3d state - all components of general stress/strain vector are generally nonzero.
	 *   General 3d strain vector has following components {sig_xx, sig_yy, sig_zz, tau_yz, tau_xz, tau_xy}
	 * - plane stress - sig_zz = tau_yz =  tau_xz = 0.
	 * - plane strain - eps_z = gamma_xz = gamma_yz = 0.
	 *   Note: as already described, if zero strain component is imposed
	 *   (Plane strain, ..) this condition must be taken into account in geometrical
	 *   relations, and corresponding component has to be included in reduced vector.
	 * - 1d uniaxial state - sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.
	 * - 2d beam layer - sigma_y=sigma_z=tau_zy=tau_xy = 0.
	 * - 3d shell layer, 2d plate layer - sigma_z = 0.
	 *
	 * Derived classes can of course extend those modes.
	 * Generally speaking, there are following major tasks, covered by declared services.
	 * - Computing real/second PK stress vector (tensor) at integration point for given strain increment and updating its
	 *   state (still temporary state, after overall equilibrium is reached).
	 * - Updating its state (final state), when equilibrium has been reached.
	 * - Returning its material stiffness (and/or flexibility) matrices for given material mode.
	 * - Storing/restoring its context to stream.
	 * - Returning its material properties.
	 *
	 * Structural material services should not be called directly by elements. Instead, they always should
	 * pass their requests to corresponding cross section model. Cross section performs all necessary integration over
	 * its volume and invokes material model services.
	 *
	 * @author almost everyone
	 * @author Jim Brouzoulis
	 * @author Mikael Öhman
	 */
	class FEM_EXPORT StructuralMaterial : public Material
	{
	private:
		/// stifness mode used in stress control
		MatResponseMode SCStiffMode = TangentStiffness; // tangent
		/// relative tolerance for stress control
		double SCRelTol = 1.e-3;  // relative stress control tolerance
		/// absolute stress tolerance for stress control
		double SCAbsTol = 1.e-12;
		/// maximum iterations for stress-control
		int SCMaxiter = 100000;
		// 材料的参考温度
		int mRefTemperature;

	public:
		/// Voigt index map
		static std::array< std::array< int, 3 >, 3 >vIindex;

		/// Symmetric Voigt index map
		static std::array< std::array< int, 3 >, 3 >svIndex;

		static int giveSymVI(int ind1, int ind2) { return svIndex[ind1 - 1][ind2 - 1]; }
		static int giveVI(int ind1, int ind2) { return vIindex[ind1 - 1][ind2 - 1]; }

		/**
		 * Constructor. Creates material with given number, belonging to given domain.
		 * @param n Material number.
		 * @param d Domain to which new material will belong.
		 */
		StructuralMaterial(int n, Domain* d);

		const char* giveClassName() const override { return "StructuralMaterial"; }
		void initializeFrom(InputRecord& ir) override;
		void giveInputRecord(DynamicInputRecord& input) override;

		/**
		 * Computes the stiffness matrix for giveRealStressVector of receiver in given integration point, respecting its history.
		 * The algorithm should use temporary or equilibrium  history variables stored in integration point status
		 * to compute and return required result.
		 * @param answer Contains result.
		 * @param mode Material response mode.
		 * @param gp Integration point.
		 * @param tStep Time step (most models are able to respond only when tStep is current time step).
		 */
		virtual void giveStiffnessMatrix(FloatMatrix& answer,
			MatResponseMode mode,
			GaussPoint* gp,
			TimeStep* tStep);

		/**
		 * Computes the real stress vector for given total strain and integration point.
		 * The total strain is defined as strain computed directly from displacement field at given time.
		 * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive driver.
		 * The service should use previously reached equilibrium history variables. Also
		 * it should update temporary history variables in status according to newly reached state.
		 * The temporary history variables are moved into equilibrium ones after global structure
		 * equilibrium has been reached by iteration process.
		 * @param answer Stress vector in reduced form. For large deformations it is treated as the second Piola-Kirchoff stress.
		 * @param gp Integration point.
		 * @param reducedStrain Strain vector in reduced form. For large deformations it is treated as the Green-Lagrange strain.
		 * @param tStep Current time step (most models are able to respond only when tStep is current time step).
		 */
		// 计算材料积分点应力，不同单元类型，计算方式不一样
		virtual void giveRealStressVector(FloatArray& answer, GaussPoint* gp,
			const FloatArray& reducedStrain, TimeStep* tStep);
		/// Default implementation relies on giveRealStressVector for second Piola-Kirchoff stress
		virtual FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 >& strain, GaussPoint* gp, TimeStep* tStep) const;
		/// Default implementation relies on giveRealStressVector_3d
		virtual FloatArrayF< 4 >giveRealStressVector_PlaneStrain(const FloatArrayF< 4 >& strain, GaussPoint* gp, TimeStep* tStep) const;
		/// Default implementation relies on giveRealStressVector_StressControl
		virtual FloatArrayF< 3 >giveRealStressVector_PlaneStress(const FloatArrayF< 3 >& reducedE, GaussPoint* gp, TimeStep* tStep) const;
		/// Default implementation relies on giveRealStressVector_StressControl
		virtual FloatArrayF< 1 >giveRealStressVector_1d(const FloatArrayF< 1 >& reducedE, GaussPoint* gp, TimeStep* tStep) const;

		/**
		 * @name Methods associated with the First PK stress tensor.
		 * Computes the first Piola-Kirchhoff stress vector for given total deformation gradient and integration point.
		 * The total deformation gradient is computed directly from displacement field at the given time step.
		 * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive
		 * driver.
		 * The service should use previously reached equilibrium history variables. Also
		 * it should update temporary history variables in status according to newly reached state.
		 * The temporary history variables are moved into equilibrium ones after global structure
		 * equilibrium has been reached by iteration process.
		 * The First Piola-Kirchhoff stress vector is computed in Total Lagrangian mode
		 *
		 * @param answer Contains result.
		 * @param gp Integration point.
		 * @param reducedF Deformation gradient in in reduced form.
		 * @param tStep Current time step (most models are able to respond only when tStep is current time step).
		 */
		 //@{
		 /// Default implementation relies on giveRealStressVector for second Piola-Kirchoff stress
		virtual FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 >& vF, GaussPoint* gp, TimeStep* tStep) const;
		/// Default implementation relies on giveFirstPKStressVector_3d
		virtual FloatArrayF< 5 >giveFirstPKStressVector_PlaneStrain(const FloatArrayF< 5 >& vF, GaussPoint* gp, TimeStep* tStep) const;
	
		/// Default implementation relies on giveFirstPKStressVector_StressControl
		virtual FloatArrayF< 4 >giveFirstPKStressVector_PlaneStress(const FloatArrayF< 4 >& vF, GaussPoint* gp, TimeStep* tStep) const;
		/// Default implementation relies on giveFirstPKStressVector_StressControl
		virtual FloatArrayF< 1 >giveFirstPKStressVector_1d(const FloatArrayF< 1 >& vF, GaussPoint* gp, TimeStep* tStep) const;
		//@}

		/**
		 * @name Methods associated with the Cauchy stress tensor.
		 * Computes the Cauchy stress vector for given increment of deformation gradient and given integration point.
		 * The increment of deformation gradient is computed directly from displacement field at the given time step
		 * and it is computed wrt configuration which was reached in the last step.
		 * The stress independent parts (temperature, eigenstrains) are subtracted in constitutive
		 * driver.
		 * The service should use previously reached equilibrium history variables. Also
		 * it should update temporary history variables in status according to newly reached state.
		 * The temporary history variables are moved into equilibrium ones after global structure
		 * equilibrium has been reached by iteration process.
		 * The Cauchy stress vector is computed in Updated Lagrangian mode
		 *
		 * @param answer Contains result.
		 * @param gp Integration point.
		 * @param reducedF Deformation gradient in in reduced form.
		 * @param tStep Current time step (most models are able to respond only when tStep is current time step).
		 */
		 //@{
		virtual void giveCauchyStressVector_3d(FloatArray& answer, GaussPoint* gp, const FloatArray& reducedF, TimeStep* tStep)
		{
			FEM_ERROR("not implemented ");
		}
		virtual void giveCauchyStressVector_PlaneStrain(FloatArray& answer, GaussPoint* gp, const FloatArray& reducedF, TimeStep* tStep)
		{
			FEM_ERROR("not implemented ");
		}
		virtual void giveCauchyStressVector_PlaneStress(FloatArray& answer, GaussPoint* gp, const FloatArray& reducedF, TimeStep* tStep)
		{
			FEM_ERROR("not implemented ");
		}
		virtual void giveCauchyStressVector_1d(FloatArray& answer, GaussPoint* gp, const FloatArray& reducedF, TimeStep* tStep)
		{
			FEM_ERROR("not implemented ");
		}
		//@}

		/**
		 * Computes full 3d material stiffness matrix at given integration point, time, respecting load history
		 * in integration point.
		 * @param answer Computed results.
		 * @param mode Material response mode.
		 * @param gp Integration point.
		 * @param tStep Time step (most models are able to respond only when tStep is current time step).
		 */
		virtual FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint* gp, TimeStep* tStep) const
		{
			FEM_ERROR("not implemented ");
		}

		/**
		 * Method for computing plane stress stiffness matrix of receiver.
		 * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
		 * reduces it to plane stress stiffness using reduce method described above.
		 * However, this reduction is quite time consuming and if it is possible,
		 * it is recommended to overload this method and provide direct method for computing
		 * particular stiffness matrix.
		 * @param answer Stiffness matrix.
		 * @param mmode Material response mode.
		 * @param gp Integration point, which load history is used.
		 * @param tStep Time step (most models are able to respond only when tStep is current time step).
		 */
		 //@{
		virtual FloatMatrixF< 3, 3 >givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;
		virtual FloatMatrixF< 4, 4 >givePlaneStressStiffnessMatrix_dPdF(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;

		virtual void givePlaneStressStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mmode, GaussPoint* gp,
			TimeStep* tStep);
		//@}

		/**
		 * Method for computing plane strain stiffness matrix of receiver.
		 * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
		 * reduces it to plane strain stiffness using reduce method described above.
		 * However, this reduction is quite time consuming and if it is possible,
		 * it is recommended to overload this method and provide direct method for computing
		 * particular stiffness matrix.
		 * Note: as already described, if zero strain component is imposed
		 * (Plane strain, ..) this condition must be taken into account in geometrical
		 * relations, and corresponding component has to be included in reduced vector.
		 * (So plane strain conditions are @f$ \epsilon_z = \gamma_{xz} = \gamma_{yz} = 0 @f$, but relations
		 * for @f$ \epsilon_z@f$ and @f$\sigma_z@f$ are included).
		 * @param answer Stiffness matrix.
		 * @param mmode Material response mode.
		 * @param gp Integration point, which load history is used.
		 * @param tStep Time step (most models are able to respond only when tStep is current time step).
		 */
		 //@{
		virtual FloatMatrixF< 4, 4 >givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;
		virtual FloatMatrixF< 5, 5 >givePlaneStrainStiffnessMatrix_dPdF(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;

		virtual void givePlaneStrainStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mmode, GaussPoint* gp,
			TimeStep* tStep);
		//@}

		/**
		 * Method for computing 1d stiffness matrix of receiver.
		 * Default implementation computes 3d stiffness matrix using give3dMaterialStiffnessMatrix and
		 * reduces it to 1d stiffness using reduce method described above.
		 * However, this reduction is quite time consuming and if it is possible,
		 * it is recommended to overload this method and provide direct method for computing
		 * particular stiffness matrix.
		 * @param mmode Material response mode.
		 * @param gp Integration point, which load history is used.
		 * @param tStep Time step (most models are able to respond only when tStep is current time step).
		 * @return Stiffness matrix.
		 */
		 //@{
		virtual FloatMatrixF< 1, 1 >give1dStressStiffMtrx(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;
		virtual FloatMatrixF< 1, 1 >give1dStressStiffnessMatrix_dPdF(MatResponseMode mmode, GaussPoint* gp, TimeStep* tStep) const;

		virtual void give1dStressStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mmode, GaussPoint* gp,
			TimeStep* tStep);
		//@}

		/**
		 * Computes equivalent of von Mises stress. Returns 0 if six stress components do not exist on the material.
		 * @param currentStress Stress vector given by 6 components.
		 */
		static double computeVonMisesStress(const FloatArray& currentStress);
		static double computeVonMisesStress_3D(const FloatArrayF< 6 >& stress);
		static double computeVonMisesStress_PlaneStress(const FloatArrayF< 3 >& stress);

		friend class CrossSection;
		friend class StructuralCrossSection;
		friend class SimpleCrossSection;
		friend class LayeredCrossSection;
	};
} // end namespace fem
#endif // structuralmaterial_h
