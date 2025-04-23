#include "femcore/FEFixedRotation.h"
#include "femcore/FEParam.h"

BEGIN_PARAM_DEFINE(FEFixedRotation, FEFixedBC)
	ADD_PARAMETER(m_dof[0], "u_dof")->setLongName("x-rotation");
	ADD_PARAMETER(m_dof[1], "v_dof")->setLongName("y-rotation");
	ADD_PARAMETER(m_dof[2], "w_dof")->setLongName("z-rotation");
END_PARAM_DEFINE();

FEFixedRotation::FEFixedRotation(FEModel* fem) : FEFixedBC(fem)
{
	m_dof[0] = false;
	m_dof[1] = false;
	m_dof[2] = false;
}

bool FEFixedRotation::Init()
{
	std::vector<int> dofs;
	if (m_dof[0]) dofs.push_back(GetDOFIndex("u"));
	if (m_dof[1]) dofs.push_back(GetDOFIndex("v"));
	if (m_dof[2]) dofs.push_back(GetDOFIndex("w"));
	SetDOFList(dofs);

	return FEFixedBC::Init();
}

