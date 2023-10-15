#pragma once

#include "fem_export.h"

class ElementTraits;

//! Base class for all element classes
//! From this class the different element classes are derived.

/**
*@~English
* @brief brief - description - about - Element .
* @
*
*@~Chinese
* @brief brief - description - about - Element.
* Tasks:
*单元相关的数据，节点，材料等
* 计算单元刚度矩阵，载荷向量
* 计算单元应力，应变
* 结果输出
*
*/

class FEM_EXPORT FEElement
{
public:
	static const int MAX_NODES = 27;
	static const int MAX_INTPOINTS = 27;
	
public:
	//! default constructor
	FEElement();

	//! destructor
	virtual ~FEElement() {}

	//! get the element ID
	int GetID() const;

	//! set the element ID
	void SetID(int n);

	//! Get the element's material ID
	int GetMatID() const;

	//! Set the element's material ID
	void SetMatID(int id);

	//Get the mesh partition that contains this element
	FEMeshPartition * GetMeshPartition() const { return m_part; }

	//Set the mesh partition that contains this element
	void SetMeshPartition(FEMeshPartition* part){ m_part = part; }

	//! Set the Local ID
	void SetLocalID(int lid) { m_lid = lid; }

	//! Get the local ID
	int GetLocalID() const { return m_lid; }

	//! clear material point data
	void ClearData();

public:
	//! Set the type of the element
	void SetType(int ntype) { FEElementLibrary::SetElementTraits(*this, ntype); }

	//! Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits);

	//! Get the element traits
	ElementTraits* GetTraits() { return m_pT; }

	//! return number of nodes
	int Nodes() const { return m_pT->m_neln; }

	//! return the element class
	int Class() const { return m_pT->Class(); }

	//! return the element shape
	int Shape() const { return m_pT->Shape(); }

	//! return the type of element
	int Type() const { return m_pT->Type(); }

	//! return number of integration points
	int GaussPoints() const { return m_pT->m_nint; }

	//! shape function values
	double* H(int n) { return m_pT->m_H[n]; }
	const double* H(int n) const { return m_pT->m_H[n]; }

	//! return number of faces
	int Faces() const { return m_pT->Faces(); }

	//! return the nodes of the face
	int GetFace(int nface, int* nodeList) const;

public:
	//! Get the material point data
	FEMaterialPoint* GetMaterialPoint(int n) { return m_State[n]; }

	//! set the material point data
	void SetMaterialPointData(FEMaterialPoint* pmp, int n)
	{ 
		pmp->m_elem = this;
		pmp->m_index = n;
		m_State[n] = pmp; 
	}

	//! serialize
	//! NOTE: state data is not serialized by the element. This has to be done by the domains.
	virtual void Serialize(DumpStream& ar);

public:
	//! evaluate scalar field at integration point
	double Evaluate(double* fn, int n);
	double Evaluate(int order, double* fn, int n);

	//! evaluate scale field at integration point
	double Evaluate(std::vector<double>& fn, int n);

	//! evaluate vector field at integration point
	vec2d Evaluate(vec2d* vn, int n);

	//! evaluate vector field at integration point
	vec3d Evaluate(vec3d* vn, int n);

	// see if this element has the node n
    bool HasNode(int n) const;

    // see if this element has the list of nodes n
    int HasNodes(int* n, const int ns) const;
    
	// find local element index of node n
    int FindNode(int n) const;

	// project data to nodes, from gauss point to node 
	void project_to_nodes(double* ai, double* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(FloatArrayF<3>* ai, vec3d*  ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3ds* ai, mat3ds* ao) const { m_pT->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3d*  ai, mat3d*  ao) const { m_pT->project_to_nodes(ai, ao); }

	// evaluate scalar field at integration point using specific interpolation order
	double Evaluate(double* fn, int order, int n);

	int ShapeFunctions(int order);
	double* H(int order, int n);

public:
	void setStatus(unsigned int n) { m_status = n; }
	unsigned int status() const { return m_status; }
	bool isActive() const { return (m_status & ACTIVE); }
	void setActive() { m_status |= ACTIVE; }
	void setInactive() { m_status &= ~ACTIVE; }

private:
	int		m_nID;		//!< element ID
	int		m_lid;		//!< local ID
	int		m_mat;		//!< material index
	unsigned int	m_status;	//!< element status

	std::vector<int>		m_node;		//!< connectivity

	// This array stores the local node numbers, that is the node numbers
	// into the node list of a domain.
	std::vector<int>		m_lnode;	//!< local connectivity

	ElementTraits*	m_pT;		//!< pointer to element traits
};

//-----------------------------------------------------------------------------

class FECORE_API FETrussElement : public FEElement
{
public:
	FETrussElement();

	FETrussElement(const FETrussElement& el);

	FETrussElement& operator = (const FETrussElement& el);

	void Serialize(DumpStream& ar) override;

public:
	double	m_a0;	// cross-sectional area
	double	m_lam;	// current stretch ratio
	double	m_tau;	// Kirchoff stress
	double	m_L0;	// initial length
};

//-----------------------------------------------------------------------------
//! Discrete element class

class FECORE_API FEDiscreteElement : public FEElement
{
public:
	FEDiscreteElement(){}
	FEDiscreteElement(const FEDiscreteElement& e);
	FEDiscreteElement& operator = (const FEDiscreteElement& e);
};

//-----------------------------------------------------------------------------
//!  This class defines a 2D element
class FECORE_API FEElement2D : public FEElement
{
public:
	//! default constructor
	FEElement2D(){}

	//! copy constructor
	FEElement2D(const FEElement2D& el);

	//! assignment operator
	FEElement2D& operator = (const FEElement2D& el);

	double* GaussWeights() { return &((FE2DElementTraits*)(m_pT))->gw[0]; }			// weights of integration points

	double* Hr(int n) { return ((FE2DElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FE2DElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

    double* Hrr(int n) { return ((FE2DElementTraits*)(m_pT))->Grr[n]; }	// shape function 2nd derivative to rr
    double* Hsr(int n) { return ((FE2DElementTraits*)(m_pT))->Gsr[n]; }	// shape function 2nd derivative to sr
    
    double* Hrs(int n) { return ((FE2DElementTraits*)(m_pT))->Grs[n]; }	// shape function 2nd derivative to rs
    double* Hss(int n) { return ((FE2DElementTraits*)(m_pT))->Gss[n]; }	// shape function 2nd derivative to ss
    
	//! values of shape functions
	void shape_fnc(double* H, double r, double s) { ((FE2DElementTraits*)(m_pT))->shape(H, r, s); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) { ((FE2DElementTraits*)(m_pT))->shape_deriv(Hr, Hs, r, s); }
};

//-----------------------------------------------------------------------------
class FECORE_API FELineElement : public FEElement
{
public:
	FELineElement();

	FELineElement(const FELineElement& el);

	FELineElement& operator = (const FELineElement& el);

	void SetTraits(FEElementTraits* pt);
};
