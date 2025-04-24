#pragma once
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! This class describes a general purpose interaction between two surfaces.
class FEM_EXPORT FESurfacePairConstraintNL : public FENLConstraint
{
public:
    //! constructor
    FESurfacePairConstraintNL(FEModel* pfem);
    
public:
    //! return the primary surface
    virtual FESurface* GetPrimarySurface() = 0;
    
    //! return the secondary surface
    virtual FESurface* GetSecondarySurface () = 0;
    
    //! temporary construct to determine if contact interface uses nodal integration rule (or facet)
    virtual bool UseNodalIntegration() = 0;
    
    //! create a copy of this interface
    virtual void CopyFrom(FESurfacePairConstraintNL* pci) {}
    
    virtual int InitEquations(int neq);
    
    void Update(std::vector<double>& ui);
    
public:
    
    using FENLConstraint::Update;
};
