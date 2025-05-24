#pragma once
#include "FEElasticShellDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! shell element. Results indicate that using this class produces poorer convergence
//! with shells than the standard FEElasticShellDomain.  This class is included
//! only for development purposes.
class FEM_EXPORT FE3FieldElasticShellDomain : public FEElasticShellDomain
{
protected:
    struct ELEM_DATA
    {
        double    eJ;        // average element jacobian
        double    ep;        // average pressure
        double    Lk;        // Lagrangian multiplier

		void Serialize(DumpStream& ar);
    };
    
public:
    //! constructor
	FE3FieldElasticShellDomain(FEModel* pfem);
    
    //! \todo Is this really used?
	FE3FieldElasticShellDomain& operator = (FE3FieldElasticShellDomain& d);
    
    //! initialize class
	bool Init() override;
    
    //! Reset data
    void Reset() override;
    
    //! augmentation
    bool Augment(int naug) override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
public: // overridden from FEElasticDomain
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // calculate stiffness matrix
    void StiffnessMatrix(FELinearSystem& LS) override;
    
protected:
    //! Dilatational stiffness component for nearly-incompressible materials
    void ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke);
    
    //! material and geometrical stiffness components
    void ElementStiffness(int iel, matrix& ke);
    
    //! update the stress of an element
    void UpdateElementStress(int iel);

public:
	bool DoAugmentations() const;
    
protected:
    std::vector<ELEM_DATA>    m_Data;

	bool	m_blaugon;		//!< augmented lagrangian flag
	double	m_augtol;		//!< augmented lagrangian tolerance
	int		m_naugmin;		//!< minimum number of augmentations
	int		m_naugmax;		//!< max number of augmentations

	DECLARE_PARAM_LIST();
};
