#include "RgShellDomain.h"
#include "femcore/FEMesh.h"
#include "RgAssembler.h"

DEFINE_META_CLASS(RgShellDomain, RgDomain, "");

//-----------------------------------------------------------------------------
RgShellDomain::RgShellDomain(FEModel* fem) : RgDomain(FE_DOMAIN_SHELL, fem)
{
	m_assembler = nullptr;
}

//-----------------------------------------------------------------------------
void RgShellDomain::ForEachElement(std::function<void(FEElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}

//-----------------------------------------------------------------------------
void RgShellDomain::SetMaterial(FEMaterial* pmat)
{
	RgDomain::SetMaterial(pmat);
}

//-----------------------------------------------------------------------------
bool RgShellDomain::Init()
{
	return RgDomain::Init();
}

//-----------------------------------------------------------------------------
void RgShellDomain::Reset()
{
}

//-----------------------------------------------------------------------------
void RgShellDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
}

//-----------------------------------------------------------------------------
void RgShellDomain::Activate()
{
}

//-----------------------------------------------------------------------------
FEMaterial* RgShellDomain::GetMaterial()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
void RgShellDomain::UnpackLM(FEElement& el, std::vector<int>& lm)
{
}

//-----------------------------------------------------------------------------
void RgShellDomain::Update(const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
void RgShellDomain::InternalForces(FEGlobalVector& R)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InternalForces(R);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgShellDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForce(R, bf);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgShellDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InertialForces(R, F);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgShellDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->StiffnessMatrix(LS);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgShellDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->MassMatrix(LS, scale);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgShellDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForceStiffness(LS, bf);
		return;
	}
}