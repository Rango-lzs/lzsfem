#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! The following domain implements the finite element formulation for a three-field
//! volume element. 
class FEM_EXPORT FE3FieldElasticSolidDomain : public FEElasticSolidDomain
{
protected:
	struct ELEM_DATA
	{
		double	eJ;		// average element jacobian at intermediate time
		double	ep;		// average pressure
		double	Lk;		// Lagrangian multiplier
        double  eJt;    // average element jacobian at current time
        double  eJp;    // average element jacobian at previous time

		void Serialize(DumpStream& ar);
	};

public:
	//! constructor
	FE3FieldElasticSolidDomain(FEModel* pfem);

	//! \todo Do I really use this?
	FE3FieldElasticSolidDomain& operator = (FE3FieldElasticSolidDomain& d);

	//! initialize class
	bool Init() override;

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
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

	//! material stiffness component
	void ElementMaterialStiffness(int iel, matrix& ke);

	//! geometrical stiffness (i.e. initial stress)
	void ElementGeometricalStiffness(int iel, matrix& ke);

	//! update the stress of an element
	void UpdateElementStress(int iel, const FETimeInfo& tp) override;

public:
	bool DoAugmentations() const;

protected:
	std::vector<ELEM_DATA>	m_Data;

	bool	m_blaugon;		//!< augmented lagrangian flag
	double	m_augtol;		//!< augmented lagrangian tolerance
	int		m_naugmin;		//!< minimum number of augmentations
	int		m_naugmax;		//!< max number of augmentations

	DECLARE_PARAM_LIST();
};
