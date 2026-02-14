#pragma once
#include "elements/RgElement/RgStructureElement.h"


//-----------------------------------------------------------------------------
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FEM_EXPORT RgShellElement : public RgStructureElement
{
public:
	RgShellElement();

	//! copy constructor
	RgShellElement(const RgShellElement& el);

	//! assignment operator
	RgShellElement& operator = (const RgShellElement& el);

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<double>	m_h0;	//!< initial shell thicknesses
	std::vector<double>	m_ht;	//!< current shell thickness
	std::vector<Vector3d>	m_d0;   //!< initial shell director

	std::vector<Vector3d>	m_g0[3];//!< reference covariant base vectors
	std::vector<Vector3d>	m_gt[3];//!< current covariant base vectors
	std::vector<Vector3d>	m_gp[3];//!< previous covariant base vectors

	std::vector<Vector3d>	m_G0[3];//!< reference contravariant base vectors
	std::vector<Vector3d>	m_Gt[3];//!< current contravariant base vectors

	// indices of solid elements this shell element is attached to.
	// the first element is attached to the back of the shell
	// and the second element is attached to the front.
	// the index is -1 if no solid is attached on that side.
	int        m_elem[2];
};

//-----------------------------------------------------------------------------
// Shell element used by old shell formulation
class FEM_EXPORT RgShellElementOld : public RgShellElement
{
public:
	RgShellElementOld();

	//! copy constructor
	RgShellElementOld(const RgShellElementOld& el);

	//! assignment operator
	RgShellElementOld& operator = (const RgShellElementOld& el);

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<Vector3d>	m_D0;	//!< initial shell directors
};

//-----------------------------------------------------------------------------
// Shell element used by new shell formulations
class FEM_EXPORT RgShellElementNew : public RgShellElement
{
public:
	RgShellElementNew();

	//! copy constructor
	RgShellElementNew(const RgShellElementNew& el);

	//! assignment operator
	RgShellElementNew& operator = (const RgShellElementNew& el);

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public: // EAS parameters

	Matrix          m_Kaai;
	Matrix          m_fa;
	Matrix          m_alpha;
	Matrix          m_alphat;
	Matrix          m_alphai;
	std::vector<Matrix>  m_Kua;
	std::vector<Matrix>  m_Kwa;
	std::vector<Matrix3ds>  m_E;
};