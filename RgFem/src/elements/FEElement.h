#pragma once

#include "femcore/fem_export.h"
#include "elements/FEElementLibrary.h"
#include "elements/FEElementTraits.h"
#include <vector>

class FEMaterialPoint;

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
* 单元相关的数据，节点，材料等  ： 属性数据
* 计算单元刚度矩阵，载荷向量	： 物理特性
* 计算单元应力，应变等结果		： 结果数据
* 结果输出
*/

/* 单元如何分类
* 1、实体单元(连续体单元)，3D Solid，2D Plane, 
* 2、结构单元, Shell, Beam, Truss
* 3、连接单元, Spring, Cohesive
* 4、
*/


class FEM_EXPORT FEElementState
{
public:
    //! default constructor
    FEElementState()
    {
    }

    //! destructor
    ~FEElementState()
    {
        Clear();
    }

    //! copy constructor
    FEElementState(const FEElementState& s);

    //! assignment operator
    FEElementState& operator=(const FEElementState& s);

    //! clear state data
    void Clear()
    {
        for (size_t i = 0; i < m_data.size(); ++i)
            delete m_data[i];
        m_data.clear();
    }

    //! create
    void Create(int n)
    {
        m_data.assign(n, static_cast<FEMaterialPoint*>(0));
    }

    //! operator for easy access to element data
    FEMaterialPoint*& operator[](int n)
    {
        return m_data[n];
    }

private:
    std::vector<FEMaterialPoint*> m_data;
};

class FEM_EXPORT FEElement
{
public:
	static constexpr int MAX_NODES = 27;
    static constexpr int MAX_INTPOINTS = 27;

public:
	//! default constructor
	FEElement();

	//! destructor
	virtual ~FEElement() {}

	//! get the element ID
	int getID() const;

	//! set the element ID
	void setID(int n);

	//! Get the element's material ID
	int getMatID() const;

	//! Set the element's material ID
	void setMatID(int id);

	//! Set the Local ID, Local Id in the domain
	void SetLocalID(int lid) { m_loc_id = lid; }

	//! Get the local ID
	int GetLocalID() const { return m_loc_id; }

public:
	//! Set the type of the element and initialize the traits by type
	void SetType(int ntype) { FEElementLibrary::SetElementTraits(*this, ntype); }

	//! Set the traits of an element
	virtual void SetTraits(FEElementTraits* ptraits);

	//! Get the element traits
	FEElementTraits* GetTraits() { return m_pTraits; }

	//! return number of nodes
	int Nodes() const { return m_pTraits->m_neln; }

	//! return the element class
	int Class() const { return m_pTraits->Class(); }

	//! return the element shape
	int Shape() const { return m_pTraits->Shape(); }

	//! return the type of element
	int Type() const { return m_pTraits->Type(); }

	//! return number of integration points
	int GaussPoints() const { return m_pTraits->m_nint; }

	//! shape function values
	double* H(int n) { return m_pTraits->m_H[n]; }
	const double* H(int n) const { return m_pTraits->m_H[n]; }

	//! return number of faces
	int Faces() const { return m_pTraits->Faces(); }

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
	void project_to_nodes(double* ai, double* ao) const { m_pTraits->project_to_nodes(ai, ao); }
	void project_to_nodes(FloatArrayF<3>* ai, vec3d*  ao) const { m_pTraits->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3ds* ai, mat3ds* ao) const { m_pTraits->project_to_nodes(ai, ao); }
	void project_to_nodes(mat3d*  ai, mat3d*  ao) const { m_pTraits->project_to_nodes(ai, ao); }

	// evaluate scalar field at integration point using specific interpolation order
	double Evaluate(double* fn, int order, int n);

	int ShapeFunctions(int order);
	double* H(int order, int n);

protected:
	int		m_id;		//!< element ID
	int		m_loc_id;		//!< local ID
	int		m_mat_id;		//!< material index

	std::vector<int>		m_node;		//!< connectivity
	// This array stores the local node numbers, that is the node numbers
	// into the node list of a domain.
	std::vector<int>		m_loc_node;	//!< local connectivity

	FEElementState m_state;
	FEElementTraits*	m_pTraits;		//!< pointer to element traits

};

//-----------------------------------------------------------------------------

class FEM_EXPORT FETrussElement : public FEElement
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

class FEM_EXPORT FEDiscreteElement : public FEElement
{
public:
	FEDiscreteElement(){}
	FEDiscreteElement(const FEDiscreteElement& e);
	FEDiscreteElement& operator = (const FEDiscreteElement& e);
};

//-----------------------------------------------------------------------------
//!  This class defines a 2D element
class FEM_EXPORT FEElement2D : public FEElement
{
public:
	//! default constructor
	FEElement2D(){}

	//! copy constructor
	FEElement2D(const FEElement2D& el);

	//! assignment operator
	FEElement2D& operator = (const FEElement2D& el);

	double* GaussWeights() { return &((FE2DElementTraits*)(m_pTraits))->gw[0]; }			// weights of integration points

	double* Hr(int n) { return ((FE2DElementTraits*)(m_pTraits))->Gr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FE2DElementTraits*)(m_pTraits))->Gs[n]; }	// shape function derivative to s

    double* Hrr(int n) { return ((FE2DElementTraits*)(m_pTraits))->Grr[n]; }	// shape function 2nd derivative to rr
    double* Hsr(int n) { return ((FE2DElementTraits*)(m_pTraits))->Gsr[n]; }	// shape function 2nd derivative to sr
    
    double* Hrs(int n) { return ((FE2DElementTraits*)(m_pTraits))->Grs[n]; }	// shape function 2nd derivative to rs
    double* Hss(int n) { return ((FE2DElementTraits*)(m_pTraits))->Gss[n]; }	// shape function 2nd derivative to ss
    
	//! values of shape functions
	void shape_fnc(double* H, double r, double s) { ((FE2DElementTraits*)(m_pTraits))->shape(H, r, s); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) { ((FE2DElementTraits*)(m_pTraits))->shape_deriv(Hr, Hs, r, s); }
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELineElement : public FEElement
{
public:
	FELineElement();

	FELineElement(const FELineElement& el);

	FELineElement& operator = (const FELineElement& el);

	void SetTraits(FEElementTraits* pt);
};

