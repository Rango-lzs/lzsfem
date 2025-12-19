#include "RgBeamDomain.h"
#include "femcore/FEMesh.h"
#include "RgAssembler.h"

DEFINE_META_CLASS(RgBeamDomain, RgDomain, "");

//-----------------------------------------------------------------------------
RgBeamDomain::RgBeamDomain(FEModel* fem) : RgDomain(FE_DOMAIN_BEAM, fem)
{
	m_assembler = nullptr;
}

//-----------------------------------------------------------------------------
void RgBeamDomain::ForEachElement(std::function<void(FEElement& el)> f)
{
	int NE = Elements();
	for (int i = 0; i < NE; ++i) f(ElementRef(i));
}

//-----------------------------------------------------------------------------
void RgBeamDomain::SetMaterial(FEMaterial* pmat)
{
	RgDomain::SetMaterial(pmat);
}

//-----------------------------------------------------------------------------
bool RgBeamDomain::Init()
{
	return RgDomain::Init();
}

//-----------------------------------------------------------------------------
void RgBeamDomain::Reset()
{
}

//-----------------------------------------------------------------------------
void RgBeamDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
}

//-----------------------------------------------------------------------------
void RgBeamDomain::Activate()
{
}

//-----------------------------------------------------------------------------
FEMaterial* RgBeamDomain::GetMaterial()
{
	return nullptr;
}

//-----------------------------------------------------------------------------
void RgBeamDomain::UnpackLM(FEElement& el, std::vector<int>& lm)
{
}

//-----------------------------------------------------------------------------
void RgBeamDomain::Update(const FETimeInfo& tp)
{
}

//-----------------------------------------------------------------------------
void RgBeamDomain::InternalForces(FEGlobalVector& R)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InternalForces(R);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgBeamDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForce(R, bf);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgBeamDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InertialForces(R, F);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgBeamDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->StiffnessMatrix(LS);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgBeamDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->MassMatrix(LS, scale);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgBeamDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForceStiffness(LS, bf);
		return;
	}
}