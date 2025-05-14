#pragma once
#include "femcore/Domain/FEShellDomain.h"
#include "femcore/FEModelParam.h"
#include <functional>
#include "femcore/FEDofList.h"
class FEDataStream;

//-----------------------------------------------------------------------------
// This class extends the FEShellDomain and implements the solid-shell interface (SSI) logic.
// It is used by the new shell formulation.
class FEM_EXPORT FESSIShellDomain : public FEShellDomainNew
{
public:
	FESSIShellDomain(FEModel* pfem);

    //! initialize domain
    //! one-time initialization, called during model initialization
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;
    
	//! Update element data prior to solving time step
	void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! Initialize shell normals
    void InitShells() override;
    
protected:
	//! Find interfaces between solid element faces and shell elements
	void FindSSI();

public:
	//! calculates covariant basis vectors at an integration point
	void CoBaseVectors0(FEShellElement& el, int n, Vector3d g[3]);

    //! calculates covariant basis vectors at any point
    void CoBaseVectors0(FEShellElement& el, double r, double s, double t, Vector3d g[3]);

	//! calculates contravariant basis vectors at an integration point
	void ContraBaseVectors0(FEShellElement& el, int n, Vector3d g[3]);

    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors0(FEShellElement& el, double r, double s, double t, Vector3d g[3]);

	// inverse jacobian with respect to reference frame at an integration point
	double invjac0(FEShellElement& el, double J[3][3], int n);

    // inverse jacobian with respect to reference frame at any point
    double invjac0(FEShellElement& el, double J[3][3], double r, double s, double t);
    
	// jacobian with respect to reference frame
	double detJ0(FEShellElement& el, int n);

    // jacobian with respect to current frame at any point
    double detJ0(FEShellElement& el, double r, double s, double t);
    
public:
    //! calculates covariant basis vectors at an integration point
    void CoBaseVectors(FEShellElement& el, int n, Vector3d g[3]);
    
    //! calculates covariant basis vectors at any point
    void CoBaseVectors(FEShellElement& el, double r, double s, double t, Vector3d g[3]);
    
    //! calculates covariant basis vectors at any point
    void CoBaseVectors(FEShellElement& el, double r, double s, double t, Vector3d g[3], const double alpha);
    
    //! calculates covariant basis vectors at an integration point at previous time
    void CoBaseVectorsP(FEShellElement& el, int n, Vector3d g[3]);
    
    //! calculates covariant basis vectors at an integration point at intermediate time
    void CoBaseVectors(FEShellElement& el, int n, Vector3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors at an integration point
    void ContraBaseVectors(FEShellElement& el, int n, Vector3d g[3]);
    
    //! calculates contravariant basis vectors at an integration point at intermediate time
    void ContraBaseVectors(FEShellElement& el, int n, Vector3d g[3], const double alpha);
    
    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors(FEShellElement& el, double r, double s, double t, Vector3d g[3]);
    
    //! calculates contravariant basis vectors at any point
    void ContraBaseVectors(FEShellElement& el, double r, double s, double t, Vector3d g[3], const double alpha);
    
    // jacobian with respect to current frame at an integration point
    double detJ(FEShellElement& el, int n);
    
    // jacobian with respect to current frame at an integration point at intermediate time
    double detJ(FEShellElement& el, int n, const double alpha);
    
    // jacobian with respect to current frame at any point
    double detJ(FEShellElement& el, double r, double s, double t);
    
    // calculate deformation gradient at an integration point
    double defgrad(FEShellElement& el, Matrix3d& F, int n);
    
    // calculate deformation gradient at any point
    double defgrad(FEShellElement& el, Matrix3d& F, double r, double s, double t);
    
    // calculate deformation gradient at an integration point at previous time
    double defgradp(FEShellElement& el, Matrix3d& F, int n);
    
    // inverse jacobian with respect to current frame
    double invjact(FEShellElement& el, double J[3][3], int n);
    
    //! evaluate a std::vector function over the shell
    Vector3d evaluate(FEShellElement& el, Vector3d* vn, Vector3d* dvn, int n);
    
    //! evaluate a scalar function over the shell
    double evaluate(FEShellElement& el, double* pn, double* dpn, int n);
    
    //! calculate the gradient of a scalar function over the shell
    Vector3d gradient(FEShellElement& el, double* pn, double* dpn, int n);
    
    //! evaluate a scalar function over the shell
    double evaluate(FEShellElement& el, std::vector<double> pn, std::vector<double> dpn, int n);
    
    //! calculate the gradient of a scalar function over the shell
    Vector3d gradient(FEShellElement& el, std::vector<double> pn, std::vector<double> dpn, int n);
    
	//! Functions for element-DOF updates
	virtual void UpdateEAS(std::vector<double>& ui) {}
	virtual void UpdateIncrementsEAS(std::vector<double>& ui, const bool binc) {}

	void Update(const FETimeInfo& tp) override;

protected:
	FEDofList	m_dofU;		// displacement dofs
	FEDofList	m_dofSU;	// shell displacement dofs
	FEDofList	m_dofR;		// rigid rotation
    
public:
    bool        m_bnodalnormals; // flag for nodal (true) or element (false) normals

	DECLARE_PARAM_LIST();
};


//-----------------------------------------------------------------------------
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<double (const FEMaterialPoint& mp)> fnc);
void writeIntegratedElementValue(FESSIShellDomain& dom, FEDataStream& ar, std::function<Vector3d  (const FEMaterialPoint& mp)> fnc);
