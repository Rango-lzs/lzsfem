/*****************************************************************//**
 * \file   structuralmaterial.cpp
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#include "structuralmaterial.h"
#include "domain.h"
#include "verbose.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "mathfem.h"
#include "engngm.h"
#include "field_manager.h"
#include "dynamic_input_record.h"

namespace fem
{
	StructuralMaterial::StructuralMaterial(int n, Domain* d) : Material(n, d) { }

	void StructuralMaterial::giveRealStressVector(FloatArray& answer, GaussPoint* gp, const FloatArray& reducedStrain, TimeStep* tStep)
	{
		///@todo Move this to StructuralCrossSection ?
		MaterialMode mode = gp->giveMaterialMode();
		if (mode == _3dMat) {
			answer = this->giveRealStressVector_3d(reducedStrain, gp, tStep);
		}
		else if (mode == _3dDegeneratedShell) {
			answer = this->giveRealStressVector_3d(reducedStrain, gp, tStep);
		}
		else if (mode == _PlaneStrain) {
			answer = this->giveRealStressVector_PlaneStrain(reducedStrain, gp, tStep);
		}
		else if (mode == _PlaneStress) {
			answer = this->giveRealStressVector_PlaneStress(reducedStrain, gp, tStep);
		}
		else if (mode == _1dMat) {
			answer = this->giveRealStressVector_1d(reducedStrain, gp, tStep);
		}
		else if (mode == _2dBeamLayer) {
			answer = this->giveRealStressVector_2dBeamLayer(reducedStrain, gp, tStep);
		}
		else if (mode == _PlateLayer) {
			answer = this->giveRealStressVector_PlateLayer(reducedStrain, gp, tStep);
		}
		else if (mode == _Fiber) {
			answer = this->giveRealStressVector_Fiber(reducedStrain, gp, tStep);
		}
		else if (mode == _2dPlateSubSoil) {
			answer = this->giveRealStressVector_2dPlateSubSoil(reducedStrain, gp, tStep);
		}
		else if (mode == _3dBeamSubSoil) {
			answer = this->giveRealStressVector_3dBeamSubSoil(reducedStrain, gp, tStep);
		}
	}


	FloatArrayF< 6 > StructuralMaterial::giveRealStressVector_3d(const FloatArrayF< 6 >& strain, GaussPoint* gp, TimeStep* tStep) const
	{
		FEM_ERROR("3d mode not supported");
	}


	FloatArrayF< 2 > StructuralMaterial::giveRealStressVector_Warping(const FloatArrayF< 2 >& reducedStrain, GaussPoint* gp, TimeStep* tStep) const
	{
		FEM_ERROR("Warping mode not supported");
	}


	FloatArrayF< 4 > StructuralMaterial::giveRealStressVector_PlaneStrain(const FloatArrayF< 4 >& strain, GaussPoint* gp, TimeStep* tStep) const
	{
		auto vS = this->giveRealStressVector_3d(assemble< 6 >(strain, { 0, 1, 2, 5 }), gp, tStep);
		return vS[{ 0, 1, 2, 5 }];
	}

	FloatArrayF< 3 > StructuralMaterial::giveRealStressVector_PlaneStress(const FloatArrayF< 3 >& reducedStrain, GaussPoint* gp, TimeStep* tStep) const
	{
		IntArray strainControl;
		StructuralMaterial::giveVoigtSymVectorMask(strainControl, _PlaneStress);
		return this->giveRealStressVector_StressControl(reducedStrain, strainControl, gp, tStep);
	}


	FloatArrayF< 1 >
		StructuralMaterial::giveRealStressVector_1d(const FloatArrayF< 1 >& reducedStrain, GaussPoint* gp, TimeStep* tStep) const
	{
		IntArray strainControl;
		StructuralMaterial::giveVoigtSymVectorMask(strainControl, _1dMat);
		return this->giveRealStressVector_StressControl(reducedStrain, strainControl, gp, tStep);
	}

	FloatArrayF< 9 >
		StructuralMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 >& vF, GaussPoint* gp, TimeStep* tStep) const
	{
		// Default implementation used if this method is not overloaded by the particular material model.
		// 1) Compute Green-Lagrange strain and call standard method for small strains.
		// 2) Treat stress as second Piola-Kirchhoff stress and convert to first Piola-Kirchhoff stress.
		// 3) Set state variables F, P

		auto F = from_voigt_form(vF);
		auto E = 0.5 * (Tdot(F, F) - eye< 3 >());
		auto vE = to_voigt_strain(E);
		auto vS = this->giveRealStressVector_3d(vE, gp, tStep);

		// Compute first PK stress from second PK stress
		auto status = static_cast<StructuralMaterialStatus*>(this->giveStatus(gp));
		auto S = from_voigt_stress(vS);
		auto P = dot(F, S);
		auto vP = to_voigt_form(P);
		status->letTempPVectorBe(vP);
		status->letTempFVectorBe(vF);

		return vP;
	}


	FloatArrayF< 5 >
		StructuralMaterial::giveFirstPKStressVector_PlaneStrain(const FloatArrayF< 5 >& vF, GaussPoint* gp, TimeStep* tStep) const
	{
		auto vP = this->giveFirstPKStressVector_3d(assemble< 9 >(vF, { 0, 1, 2, 5, 8 }), gp, tStep);
		return vP[{ 0, 1, 2, 5, 8 }];
	}

	FloatArrayF< 4 >
		StructuralMaterial::giveFirstPKStressVector_PlaneStress(const FloatArrayF< 4 >& reducedvF, GaussPoint* gp, TimeStep* tStep) const
	{
		IntArray FControl;
		StructuralMaterial::giveVoigtVectorMask(FControl, _PlaneStress);
		return this->giveRealStressVector_StressControl(reducedvF, FControl, gp, tStep);
	}


	FloatArrayF< 1 >
		StructuralMaterial::giveFirstPKStressVector_1d(const FloatArrayF< 1 >& reducedvF, GaussPoint* gp, TimeStep* tStep) const
	{
		IntArray FControl;
		StructuralMaterial::giveVoigtVectorMask(FControl, _PlaneStress);
		return this->giveRealStressVector_StressControl(reducedvF, FControl, gp, tStep);
	}

	void
		StructuralMaterial::giveStiffnessMatrix(FloatMatrix& answer,
			MatResponseMode rMode,
			GaussPoint* gp, TimeStep* tStep)
		//
		// Returns characteristic material stiffness matrix of the receiver
		//
	{
		MaterialMode mMode = gp->giveMaterialMode();
		switch (mMode) {
		case _3dMat:
			answer = this->give3dMaterialStiffnessMatrix(rMode, gp, tStep);
			break;
		case _PlaneStress:
			answer = this->givePlaneStressStiffMtrx(rMode, gp, tStep);
			break;
		case _PlaneStrain:
			answer = this->givePlaneStrainStiffMtrx(rMode, gp, tStep);
			break;
		case _1dMat:
			answer = this->give1dStressStiffMtrx(rMode, gp, tStep);
			break;
		case _PlateLayer:
			answer = this->givePlateLayerStiffMtrx(rMode, gp, tStep);
			break;
		case _2dBeamLayer:
			answer = this->give2dBeamLayerStiffMtrx(rMode, gp, tStep);
			break;
		case _Fiber:
			answer = this->giveFiberStiffMtrx(rMode, gp, tStep);
			break;
		case _Warping:
			answer = eye< 2 >();
			break;
		default:
			FEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode));
		}
	}



	FloatMatrixF< 5, 5 >
		StructuralMaterial::givePlaneStrainStiffnessMatrix_dPdF(MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix_dPdF(mode, gp, tStep);
		return m3d({ 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 });
	}


	FloatMatrixF< 4, 4 >
		StructuralMaterial::givePlaneStressStiffnessMatrix_dPdF(MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix_dPdF(mode, gp, tStep);
		auto c3d = inv(m3d);
		return inv(c3d({ 0, 1, 5, 8 }, { 0, 1, 5, 8 }));
	}


	FloatMatrixF< 1, 1 >
		StructuralMaterial::give1dStressStiffnessMatrix_dPdF(MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix(mode, gp, tStep);
		auto c3d = inv(m3d);
		return { 1. / c3d.at(1, 1) };
	}


	void
		StructuralMaterial::give3dMaterialStiffnessMatrix_dCde(FloatMatrix& answer,
			MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep)
	{
		///@todo what should be default implementaiton?
		FEM_ERROR("There is no default implementation");
	}


	void
		StructuralMaterial::givePlaneStressStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep)
	{
		FEM_ERROR("There is no default implementation");
	}


	void
		StructuralMaterial::givePlaneStrainStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep)
	{
		FEM_ERROR("There is no default implementation");
	}

	void
		StructuralMaterial::give1dStressStiffMtrx_dCde(FloatMatrix& answer,
			MatResponseMode mode,
			GaussPoint* gp, TimeStep* tStep)
	{
		FEM_ERROR("There is no default implementation");
	}

	FloatMatrixF< 3, 3 >
		StructuralMaterial::givePlaneStressStiffMtrx(MatResponseMode mode,
			GaussPoint* gp,
			TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix(mode, gp, tStep);
		auto c3d = inv(m3d);
		return inv(c3d({ 0, 1, 5 }, { 0, 1, 5 }));
	}

	FloatMatrixF< 4, 4 >
		StructuralMaterial::givePlaneStrainStiffMtrx(MatResponseMode mode,
			GaussPoint* gp,
			TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix(mode, gp, tStep);
		return m3d({ 0, 1, 2, 5 }, { 0, 1, 2, 5 });
	}

	FloatMatrixF< 1, 1 >
		StructuralMaterial::give1dStressStiffMtrx(MatResponseMode mode,
			GaussPoint* gp,
			TimeStep* tStep) const
	{
		auto m3d = this->give3dMaterialStiffnessMatrix(mode, gp, tStep);
		auto c3d = inv(m3d);
		return {
			1. / c3d.at(1, 1)
		};
	}

	double
		StructuralMaterial::computeVonMisesStress(const FloatArray& stress)
	{
		if (stress.giveSize() == 3) {
			return computeVonMisesStress_PlaneStress(stress);
		}
		else if (stress.giveSize() == 4) {
			// Plane strain
			double v1 = ((stress.at(1) - stress.at(2)) * (stress.at(1) - stress.at(2)));
			double v2 = ((stress.at(2) - stress.at(3)) * (stress.at(2) - stress.at(3)));
			double v3 = ((stress.at(3) - stress.at(1)) * (stress.at(3) - stress.at(1)));

			double J2 = (1. / 6.) * (v1 + v2 + v3) + stress.at(4) * stress.at(4);

			return sqrt(3 * J2);
		}
		else if (stress.giveSize() == 6) {
			return computeVonMisesStress_3D(stress);
		}
		else {
			return 0.0;
		}
	}

	double
		StructuralMaterial::computeVonMisesStress_PlaneStress(const FloatArrayF< 3 >& stress)
	{
		return sqrt(stress.at(1) * stress.at(1) + stress.at(2) * stress.at(2)
			- stress.at(1) * stress.at(2) + 3 * stress.at(3) * stress.at(3));
	}


	double
		StructuralMaterial::computeVonMisesStress_3D(const FloatArrayF< 6 >& stress)
	{
		double v1 = ((stress.at(1) - stress.at(2)) * (stress.at(1) - stress.at(2)));
		double v2 = ((stress.at(2) - stress.at(3)) * (stress.at(2) - stress.at(3)));
		double v3 = ((stress.at(3) - stress.at(1)) * (stress.at(3) - stress.at(1)));

		double J2 = (1. / 6.) * (v1 + v2 + v3) + stress.at(4) * stress.at(4) +
			stress.at(5) * stress.at(5) + stress.at(6) * stress.at(6);

		return sqrt(3 * J2);
	}

	void
		StructuralMaterial::initializeFrom(InputRecord& ir)
	{
		Material::initializeFrom(ir);

		mRefTemperature = 0.0;
		IR_GIVE_OPTIONAL_FIELD(ir, mRefTemperature, _IFT_StructuralMaterial_referencetemperature);

		double alpha = 0.0;
		IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_StructuralMaterial_talpha);
		if (!propertyDictionary.includes(tAlpha)) {
			//    if (alpha > 0.0 && !propertyDictionary.includes(tAlpha)) {
			// put isotropic thermal expansion coeff into dictionary, if provided
			// and not previosly defined
			propertyDictionary.add(tAlpha, alpha);
		}
		int stiffmode = 0;
		IR_GIVE_OPTIONAL_FIELD(ir, stiffmode, _IFT_StructuralMaterial_StressControl_stiffmode);
		this->SCStiffMode = (MatResponseMode)stiffmode;

		IR_GIVE_OPTIONAL_FIELD(ir, this->SCRelTol, _IFT_StructuralMaterial_StressControl_reltol);
		IR_GIVE_OPTIONAL_FIELD(ir, this->SCAbsTol, _IFT_StructuralMaterial_StressControl_abstol);
		IR_GIVE_OPTIONAL_FIELD(ir, this->SCMaxiter, _IFT_StructuralMaterial_StressControl_maxiter);
	}


	void
		StructuralMaterial::giveInputRecord(DynamicInputRecord& input)
	{
		Material::giveInputRecord(input);
		input.setField(this->referenceTemperature, _IFT_StructuralMaterial_referencetemperature);
		input.setField((int)this->SCStiffMode, _IFT_StructuralMaterial_StressControl_stiffmode);
		input.setField(this->SCRelTol, _IFT_StructuralMaterial_StressControl_reltol);
		input.setField(this->SCAbsTol, _IFT_StructuralMaterial_StressControl_abstol);
	}

} // end namespace oofem
