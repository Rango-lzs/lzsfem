#include "RgSolidDomain.h"
#include "femcore/FEMesh.h"
#include "RgAssembler.h"

DEFINE_META_CLASS(RgSolidDomain, RgDomain, "");

//-----------------------------------------------------------------------------
RgSolidDomain::RgSolidDomain(FEModel* fem) : RgDomain(fem)
{
	m_assembler = nullptr;
}

//-----------------------------------------------------------------------------
void RgSolidDomain::ForEachElement(std::function<void(RgElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}

//-----------------------------------------------------------------------------
void RgSolidDomain::SetMaterial(RgMaterial* pmat)
{
	RgDomain::SetMaterial(pmat);
}

//-----------------------------------------------------------------------------
bool RgSolidDomain::Init()
{
	return RgDomain::Init();
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Reset()
{
}

//-----------------------------------------------------------------------------
void RgSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Activate()
{
}

//-----------------------------------------------------------------------------
RgMaterial* RgSolidDomain::GetMaterial()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
void RgSolidDomain::UnpackLM(RgElement& el, std::vector<int>& lm)
{
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Update(const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
void RgSolidDomain::InternalForces(FEGlobalVector& R)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InternalForces(R);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForce(R, bf);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InertialForces(R, F);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->StiffnessMatrix(LS);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->MassMatrix(LS, scale);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForceStiffness(LS, bf);
		return;
	}
}