#include "FEBioModelBuilder.h"
#include "FEBioModel.h"
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEUDGHexDomain.h>
#include <FEBioMech/FEUT4Domain.h>
#include <FEBioMech/FESSIShellDomain.h>
#include <FEBioMech/RigidBC.h>
#include <FEBioMech/FERigidForce.h>
#include <FEBioMech/FEMechModel.h>

// In FEBio 3, the bulk modulus k must be defined at the top - level.
// However, this could break backward compatibility, so for older file version
// we apply this hack that collects the child moduli and assigns it to the top-level
void FixUncoupledMaterial(FEUncoupledMaterial* mat)
{
	double K = mat->m_K;
	for (int i = 0; i < mat->Properties(); ++i)
	{
		FEUncoupledMaterial* mati = dynamic_cast<FEUncoupledMaterial*>(mat->GetProperty(i));
		if (mati)
		{
			FixUncoupledMaterial(mati);
			K += mati->m_K;
			mati->m_K = 0.0;
		}
	}
	mat->m_K = K;
}

FEBioModelBuilder::FEBioModelBuilder(FEBioModel& fem) : FEModelBuilder(fem)
{

}

void FEBioModelBuilder::AddMaterial(FEMaterial* mat)
{
	FEModel& fem = GetFEModel();
	fem.AddMaterial(mat);

	// For uncoupled materials, we collect the bulk moduli of child materials
	// and assign it to the top-level material (this one)
	FEUncoupledMaterial* pucm = dynamic_cast<FEUncoupledMaterial*>(mat);
	if (pucm) FixUncoupledMaterial(pucm);
}

FEDomain* FEBioModelBuilder::CreateDomain(FE_Element_Spec espec, FEMaterial* mat)
{
	FEModel& fem = GetFEModel();

	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEDomain* pdom = febio.CreateDomain(espec, &fem.GetMesh(), mat);

	// Handle dome special cases
	// TODO: Find a better way of dealing with these special cases
	FEUDGHexDomain* udg = dynamic_cast<FEUDGHexDomain*>(pdom);
	if (udg)
	{
		udg->SetHourGlassParameter(m_udghex_hg);
	}

	FEUT4Domain* ut4 = dynamic_cast<FEUT4Domain*>(pdom);
	if (ut4)
	{
		ut4->SetUT4Parameters(m_ut4_alpha, m_ut4_bdev);
	}

	FESSIShellDomain* ssi = dynamic_cast<FESSIShellDomain*>(pdom);
	if (ssi) {
		ssi->m_bnodalnormals = espec.m_shell_norm_nodal;
	}

	return pdom;
}

//-----------------------------------------------------------------------------
void FEBioModelBuilder::AddRigidComponent(FEStepComponent* pmc)
{
	FEMechModel& fem = static_cast<FEMechModel&>(GetFEModel());

	AddComponent(pmc);

	FERigidFixedBC* prc = dynamic_cast<FERigidFixedBC*>(pmc);
	if (prc) { fem.AddRigidFixedBC(prc); return; }

	FERigidPrescribedBC* prf = dynamic_cast<FERigidPrescribedBC*>(pmc);
	if (prf) { fem.AddRigidPrescribedBC(prf); return; }

	FERigidIC* ric = dynamic_cast<FERigidIC*>(pmc);
	if (ric) { fem.AddRigidInitialCondition(ric); return; }

	FEModelLoad* pml = dynamic_cast<FEModelLoad*>(pmc);
	if (pml) { AddModelLoad(pml); return; }

	assert(false);
}
