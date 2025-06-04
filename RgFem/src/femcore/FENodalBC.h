#pragma once

#include "femcore/FEBoundaryCondition.h"

class FENodeSet;

class FEM_EXPORT FENodalBC : public FEBoundaryCondition
{
	DECLARE_META_CLASS(FENodalBC, FEBoundaryCondition);

public:
	FENodalBC();

	// Set the node set
	void SetNodeSet(FENodeSet* nodeSet);

	// Get the node set
	FENodeSet* GetNodeSet();

	void Serialize(DumpStream& ar) override;

private:
	FENodeSet* m_nodeSet;
};
