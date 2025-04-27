#include "input/febio/FEBioModelBuilder.h"

FEBioModelBuilder::FEBioModelBuilder(FEModel& fem) : FEModelBuilder(fem)
{

}

void FEBioModelBuilder::AddMaterial(FEMaterial* mat)
{
	FEModel& fem = GetFEModel();
	fem.AddMaterial(mat);

	// For uncoupled materials, we collect the bulk moduli of child materials
	// and assign it to the top-level material (this one)
	/*FEUncoupledMaterial* pucm = dynamic_cast<FEUncoupledMaterial*>(mat);
	if (pucm) FixUncoupledMaterial(pucm);*/
}

FEDomain* FEBioModelBuilder::CreateDomain(FE_Element_Spec espec, FEMaterial* mat)
{
	FEModel& fem = GetFEModel();

	//FECoreKernel& febio = FECoreKernel::GetInstance();
	//FEDomain* pdom = febio.CreateDomain(espec, &fem.GetMesh(), mat);

	//// Handle dome special cases
	//// TODO: Find a better way of dealing with these special cases
	//FEUDGHexDomain* udg = dynamic_cast<FEUDGHexDomain*>(pdom);
	//if (udg)
	//{
	//	udg->SetHourGlassParameter(m_udghex_hg);
	//}

	//FEUT4Domain* ut4 = dynamic_cast<FEUT4Domain*>(pdom);
	//if (ut4)
	//{
	//	ut4->SetUT4Parameters(m_ut4_alpha, m_ut4_bdev);
	//}

	//FESSIShellDomain* ssi = dynamic_cast<FESSIShellDomain*>(pdom);
	//if (ssi) {
	//	ssi->m_bnodalnormals = espec.m_shell_norm_nodal;
	//}

	return nullptr;
}

//-----------------------------------------------------------------------------
void FEBioModelBuilder::AddRigidComponent(FEStepComponent* pmc)
{
	/*FEMechModel& fem = static_cast<FEMechModel&>(GetFEModel());

	AddComponent(pmc);

	FERigidFixedBC* prc = dynamic_cast<FERigidFixedBC*>(pmc);
	if (prc) { fem.AddRigidFixedBC(prc); return; }

	FERigidPrescribedBC* prf = dynamic_cast<FERigidPrescribedBC*>(pmc);
	if (prf) { fem.AddRigidPrescribedBC(prf); return; }

	FERigidIC* ric = dynamic_cast<FERigidIC*>(pmc);
	if (ric) { fem.AddRigidInitialCondition(ric); return; }

	FEModelLoad* pml = dynamic_cast<FEModelLoad*>(pmc);
	if (pml) { AddModelLoad(pml); return; }

	assert(false);*/
}
