#pragma once
#include "femcore/fem_export.h"
#include "femcore/FEDofList.h"
#include "femcore/FEMesh.h"
#include <functional>
#include <assert.h>
#include "femcore/FEObjectBase.h"
#include "../Matrix/FEGlobalMatrix.h"

//-----------------------------------------------------------------------------
// forward declarations
class RgElement;
class FEMaterial;
class FENode;
class FEMesh;
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;
class FESolver;

struct FETimeInfo;

//-----------------------------------------------------------------------------
//! Base class for defining domains.
//!
class FEM_EXPORT RgDomain : public FEObjectBase
{
	DECLARE_META_CLASS(RgDomain, FEObjectBase);

public:
	//! constructor
	RgDomain(FEModel* pfem);

	//! destructor
	virtual ~RgDomain();

	//! copy constructor
	RgDomain(const RgDomain& d);

	//! assignment operator
	RgDomain& operator = (const RgDomain& d);

	//! Create domain data structures
	virtual bool Create(int nsize, FE_Element_Spec espec) = 0;

	virtual int domType() { return 1; }

	//! return number of elements
	virtual int Elements() const = 0;

	//! return a reference to an element
	virtual RgElement& ElementRef(int n) = 0;
	virtual const RgElement& ElementRef(int n) const = 0;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! Get the total number of degrees of freedom
	int GetTotalDofs();

public:
	// --- M E S H   I N T E R F A C E ---

	//! return number of nodes
	int Nodes() const { return (int)m_Node.size(); }

	//! return a node
	FENode& Node(int i) { return *m_Node[i]; }
	const FENode& Node(int i) const { return *m_Node[i]; }

	////! return a node index
	//int NodeIndex(int i) const { return m_Node[i] - &m_pMesh->Node(0); }

	//! return a node index from local node index
	int GetNodeIndex(int inode) const { return m_lnode[inode]; }

	//! return a node index from local element node index
	int GetElementNodeIndex(int iel, int inode) const;

	//! get the mesh
	FEMesh* GetMesh() { return m_pMesh; }

public:
	// --- M A T E R I A L   I N T E R F A C E ---

	//! Set the domain's material
	virtual void SetMaterial(FEMaterial* pmat);

	//! Get the domain's material
	FEMaterial* GetMaterial() { return m_pMat; }

	//! Get the domain's material
	const FEMaterial* GetMaterial() const { return m_pMat; }

public:
	// --- E L E M E N T   W A L K I N G ---

	//! loop over all elements
	virtual void ForEachElement(std::function<void(RgElement& el)> f);

public:
	// --- S U B C L A S S   I N T E R F A C E ---

	//! Initialize data
	virtual bool Init();

	//! Reset data
	virtual void Reset();

	//! Initialize elements
	virtual void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! Activate domain
	virtual void Activate();
	bool isActive() { return true; }

	//! Unpack element data
	virtual void UnpackLM(RgElement& el, std::vector<int>& lm) = 0;

	virtual void BuildMatrixProfile(FEGlobalMatrix& M);

protected:
	FEModel*		m_pfem;		//!< pointer to model
	FEMesh*		m_pMesh;		//!< pointer to mesh
	FEMaterial*	m_pMat;		//!< pointer to material
	std::vector<FENode*>	m_Node;		//!< list of nodes
	std::vector<int>		m_lnode;		//!< local node indices
	std::string		m_szname;		//!< domain name
};