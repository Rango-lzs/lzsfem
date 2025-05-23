#pragma once

#include "femcore/fem_export.h"
#include <vector>

class FEMesh;

//-----------------------------------------------------------------------------
//! This class implements an algoritm that calculates a permutation of 
//! the node numbering in order to obtain a bandwidth reduced stiffness matrix

//! The algorithm comes from "An algorithm for reducing the bandwidth and 
//! profile of a sparse matrix", by N.E.Gibbs e.a. It applies the algorithm
//! on the node numberings in stead of the actual sparse matrix since that
//! was easier to implement :). In the future I would like to extend it to
//! work with the actual sparse matrix. 

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
