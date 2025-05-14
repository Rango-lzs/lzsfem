#include "FESolidDomainFactory.h"
//#include "materials/FERigidMaterial.h"
//#include "FEUncoupledMaterial.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FELinearTrussDomain.h"
#include "femcore/Domain/FERigidSolidDomain.h"
#include "femcore/Domain/FERigidShellDomain.h"
#include "femcore/Domain/FERemodelingElasticDomain.h"
#include "femcore/Domain/FEUDGHexDomain.h"
#include "femcore/Domain/FEUT4Domain.h"
#include "femcore/Domain/FE3FieldElasticSolidDomain.h"
#include "femcore/Domain/FEDiscreteElasticDomain.h"
#include "materials/FERemodelingElasticMaterial.h"
#include "materials/FECore/FEDiscreteMaterial.h"
#include "materials/FEDiscreteElementMaterial.h"
#include "femcore/Domain/FESRIElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FEDomain* FESolidDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FEModel* pfem = pmat->GetFEModel();
	ElementCategory eclass = spec.eclass;
	ElementShape eshape = spec.eshape;
	ElementType etype = spec.etype;

	// this will store the domain we are going to allocate
	FEDomain* pd = nullptr;

	// try to allocate the domain
	if (eclass == FE_ELEM_SOLID)
	{
		if (dynamic_cast<FESolidMaterial*>(pmat) == nullptr) return nullptr;

		const char* sztype = nullptr;
		if      (dynamic_cast<FERigidMaterial*            >(pmat)) sztype = "rigid-solid";
		else if (dynamic_cast<FERemodelingElasticMaterial*>(pmat)) sztype = "remodeling-solid";
		else if (eshape == ET_HEX8)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else
			{
				if (etype == FE_HEX8G1) sztype = "udg-hex";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_TET10) || (eshape == ET_TET15) || (eshape == ET_TET20))
		{
			if ((etype == FE_TET10G8RI4) || (etype == FE_TET10G4RI1))
			{
				sztype = "sri-solid";
			}
			else
			{
				if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
				else sztype = "elastic-solid";
			}
		}
		else if ((eshape == ET_HEX20) || (eshape == ET_HEX27) || (eshape == ET_PYRA13))
		{
//			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			sztype = "elastic-solid";
		}
		else if (eshape == ET_TET4)
		{
			if (spec.m_but4) sztype = "ut4-solid";
			else sztype = "elastic-solid";
		}
		else if (eshape == ET_TET5)
		{
			sztype = "elastic-solid";
		}
		else if (eshape == ET_PENTA6)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else sztype = "elastic-solid";
		}
		else if (eshape == ET_PENTA15)
		{
			// three-field implementation for uncoupled materials
//          if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field_hex)) sztype = "three-field-solid";
			sztype = "elastic-solid";
		}
		else if (eshape == ET_PYRA5)
		{
			// three-field implementation for uncoupled materials
			if (dynamic_cast<FEUncoupledMaterial*>(pmat) && (spec.m_bthree_field)) sztype = "three-field-solid";
			else sztype = "elastic-solid";
		}

		if (sztype) pd = fecore_new<FESolidDomain>(sztype, pfem);
	}
	else if (eclass == FE_ELEM_SHELL)
	{
		if (dynamic_cast<FESolidMaterial*>(pmat) == nullptr) return nullptr;

		const char* sztype = nullptr;
		if (dynamic_cast<FERigidMaterial*>(pmat))
		{
			if (spec.m_shell_formulation == OLD_SHELL) sztype = "rigid-shell-old";
			else sztype = "rigid-shell";
		}
		else if (dynamic_cast<FESolidMaterial*>(pmat))
		{
			if ((eshape == ET_QUAD4) || (eshape == ET_TRI3) || (eshape == ET_QUAD8) || (eshape == ET_TRI6))
			{
				switch (spec.m_shell_formulation)
				{
				case NEW_SHELL:
					// three-field implementation for uncoupled materials
					if (dynamic_cast<FEUncoupledMaterial*>(pmat)) {
						if (spec.m_bthree_field) sztype = "three-field-shell";
						else if ((eshape == ET_QUAD4) && spec.m_bthree_field) sztype = "three-field-shell";
						else if ((eshape == ET_TRI3) && spec.m_bthree_field) sztype = "three-field-shell";
						else sztype = "elastic-shell";
					}
					else sztype = "elastic-shell";
					break;
				case OLD_SHELL: sztype = "elastic-shell-old"; break;
				case EAS_SHELL: sztype = "elastic-shell-eas"; break;
				case ANS_SHELL: sztype = "elastic-shell-ans"; break;
				default:
					return 0;
				}
			}
		}

		if (sztype) pd = fecore_new<FEShellDomain>(sztype, pfem);
	}
	else if (eclass == FE_ELEM_TRUSS)
	{
		const char* sztype = nullptr;
		if (dynamic_cast<FETrussMaterial*>(pmat))
			if (eshape == ET_TRUSS2) sztype = "linear-truss";

		if (sztype) pd = RANGO_NEW<FETrussDomain>(pfem, sztype);
	}
	else if (eclass == FE_ELEM_DISCRETE)
	{
		const char* sztype = nullptr;
		if (dynamic_cast<FEDiscreteMaterial*>(pmat))
		{
	//		if      (eclass == FE_ELEM_WIRE    ) sztype = "deformable-spring";
			if      (eclass == FE_ELEM_WIRE    ) sztype = "deformable-spring2";
			else if (eclass == FE_ELEM_DISCRETE)
			{
				if (dynamic_cast<FEDiscreteElasticMaterial*>(pmat)) sztype = "discrete";
			}
			else return 0;
		}

		if (sztype) pd = RANGO_NEW<FEDiscreteDomain>(sztype, pfem);
	}

	if (pd) pd->SetMaterial(pmat);
	return pd;
}
