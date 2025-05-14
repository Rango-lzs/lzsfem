
#pragma once
#include "femcore/Domain/FEElasticDomain.h"
#include "femcore/Domain/FESSIShellDomain.h"
#include "materials/FESolidMaterial.h"

class FEMaterial;
class FEElement;
class DumpStream;
class FEDofList;
class FEShellElementNew;
class FETimeInfo;

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FEElasticANSShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
    FEElasticANSShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticANSShellDomain& operator = (FEElasticANSShellDomain& d);
    
    //! Activate the domain
    void Activate() override;
    
    //! Unpack shell element data
    void UnpackLM(FEElement& el, std::vector<int>& lm) override;
    
    //! Set flag for update for dynamic quantities
    void SetDynamicUpdateFlag(bool b);
    
    //! serialization
    void Serialize(DumpStream& ar) override;
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat) override;

	//! get the total dof list
	const FEDofList& GetDOFList() const override;
    
public: // overrides from FEElasticDomain
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R) override;
    
    //! Calculates inertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;
    
    //! calculate body force
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    // update the element stress
    void UpdateElementStress(int iel, const FETimeInfo& tp);
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! calculates the global stiffness Matrix for this domain
    void StiffnessMatrix(FELinearSystem& LS) override;
    
    // inertial stiffness
    void MassMatrix(FELinearSystem& LS, double scale) override;
    
    // body force stiffness
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    
    // evaluate strain E and Matrix hu and hw
	void EvaluateEh(FEShellElementNew& el, const int n, const Vector3d* Gcnt, Matrix3ds& E, std::vector<Matrix>& hu, std::vector<Matrix>& hw, std::vector<Vector3d>& Nu, std::vector<Vector3d>& Nw);
    
public:
    
    // --- S T I F F N E S S ---
    
    //! calculates the shell element stiffness Matrix
    void ElementStiffness(int iel, Matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress std::vector for shell elements
	void ElementInternalForce(FEShellElementNew& el, std::vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElementNew& el, std::vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElementNew& el, std::vector<double>& fe);
    
    //! Calculates the inertial force for shell elements
    void ElementInertialForce(FEShellElementNew& el, std::vector<double>& fe);
    
    //! calculates the solid element mass Matrix
	void ElementMassMatrix(FEShellElementNew& el, Matrix& ke, double a);
    
    //! calculates the stiffness Matrix due to body forces
	void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElementNew& el, Matrix& ke);
    
public:
    
    // --- A N S  M E T H O D ---
    
    // Evaluate contravariant components of Matrix3ds tensor
    void mat3dsCntMat61(const Matrix3ds s, const Vector3d* Gcnt, Matrix& S);
    
    // Evaluate contravariant components of tens4ds tensor
    void tens4dsCntMat66(const tens4ds c, const Vector3d* Gcnt, Matrix& C);
    void tens4dmmCntMat66(const tens4dmm c, const Vector3d* Gcnt, Matrix& C);

    // Evaluate the strain using the ANS method
	void CollocationStrainsANS(FEShellElementNew& el, std::vector<double>& E, std::vector< std::vector<Vector3d>>& HU, std::vector< std::vector<Vector3d>>& HW, Matrix& NS, Matrix& NN);
    
	void EvaluateANS(FEShellElementNew& el, const int n, const Vector3d* Gcnt, Matrix3ds& Ec, std::vector<Matrix>& hu, std::vector<Matrix>& hw, std::vector<double>& E, std::vector< std::vector<Vector3d>>& HU, std::vector< std::vector<Vector3d>>& HW);
    
protected:
    FESolidMaterial*    m_pMat;
    int                 m_nEAS;
    bool                m_update_dynamic;    //!< flag for updating quantities only used in dynamic analysis
    
    bool    m_secant_stress;    //!< use secant approximation to stress
    bool    m_secant_tangent;   //!< flag for using secant tangent
    
    FEDofList   m_dofV;
    FEDofList   m_dofSV;
    FEDofList    m_dofSA;
    FEDofList    m_dofR;
    FEDofList    m_dof;
};
