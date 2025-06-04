#include "femcore/FEFixedDisplacement.h"
#include "femcore/FEParam.h"

BEGIN_PARAM_DEFINE(FEFixedDisplacement, FEFixedBC)
	ADD_PARAMETER(m_dofx, "x_dof")->setLongName("x-displacement");
	ADD_PARAMETER(m_dofy, "y_dof")->setLongName("y-displacement");
	ADD_PARAMETER(m_dofz, "z_dof")->setLongName("z-displacement");
END_PARAM_DEFINE();

FEFixedDisplacement::FEFixedDisplacement() : FEFixedBC()
{
	m_dofx = false;
	m_dofy = false;
	m_dofz = false;
}

bool FEFixedDisplacement::Init()
{
	std::vector<int> dofs;
	if (m_dofx) dofs.push_back(GetDOFIndex("x"));
	if (m_dofy) dofs.push_back(GetDOFIndex("y"));
	if (m_dofz) dofs.push_back(GetDOFIndex("z"));

	SetDOFList(dofs);
	return FEFixedBC::Init();
}
