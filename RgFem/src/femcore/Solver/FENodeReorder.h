#pragma once

#include "femcore/fem_export.h"
#include <vector>

class FEMesh;

//-----------------------------------------------------------------------------
//! This class implements an algoritm that calculates a permutation of 
//! the node numbering in order to obtain a bandwidth reduced stiffness Matrix

//! The algorithm comes from "An algorithm for reducing the bandwidth and 
//! profile of a sparse Matrix", by N.E.Gibbs e.a. It applies the algorithm
//! on the node numberings in stead of the actual sparse Matrix since that
//! was easier to implement :). In the future I would like to extend it to
//! work with the actual sparse Matrix. 

class FEM_EXPORT FENodeReorder
{

public:
	//! default constructor
	FENodeReorder();

	//! destructor
	virtual ~FENodeReorder();

	//! calculates the permutation vector
	void Apply(FEMesh& m, std::vector<int>& P);
};
