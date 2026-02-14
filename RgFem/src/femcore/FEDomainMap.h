#pragma once
#include "FEDataMap.h"
#include "femcore/Domain/RgDomain.h"
#include "elements/RgElementSet.h"
#include "materials/RgMaterialPoint.h"

class FEM_EXPORT FEDomainMap : public FEDataMap
{
public:
	//! default constructor
	FEDomainMap();
	FEDomainMap(FEDataType dataType, Storage_Fmt format = FMT_MULT);

	//! copy constructor
	FEDomainMap(const FEDomainMap& map);

	//! assignment operator
	FEDomainMap& operator = (const FEDomainMap& map);

	//! Create a surface data map for this surface
	bool Create(RgElementSet* ps, double val = 0.0);

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! get the value at a material point
	double value(const RgMaterialPoint& pt) override;
	Vector3d valueVec3d(const RgMaterialPoint& pt) override;
	Matrix3d valueMatrix3d(const RgMaterialPoint& pt) override;
	Matrix3ds valueMatrix3ds(const RgMaterialPoint& pt) override;

	//! get the value at a node (only works with FMT_NODE!!)
	double NodalValue(int nid);

	//! Get the element set
	const RgElementSet* GetElementSet() const { return m_elset; }

	//! return max nr of nodes
	int MaxNodes() const { return m_maxElemNodes; }

	//! return storage format
	int	StorageFormat() const { return m_fmt; }

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	// merge with another map
	bool Merge(FEDomainMap& map);

public:
	template <typename T> T value(int nelem, int node)
	{
		return get<T>(nelem*m_maxElemNodes + node);
	}

	template <typename T> void setValue(int nelem, int node, const T& v)
	{
		set<T>(nelem*m_maxElemNodes + node, v);
	}

	void setValue(int n, double v) override;
	void setValue(int n, const Vector2d& v) override;
	void setValue(int n, const Vector3d& v) override;
	void setValue(int n, const Matrix3d& v) override;
	void setValue(int n, const Matrix3ds& v) override;

	void fillValue(double v) override;
	void fillValue(const Vector2d& v) override;
	void fillValue(const Vector3d& v) override;
	void fillValue(const Matrix3d& v) override;
	void fillValue(const Matrix3ds& v) override;

private:
	void Realloc(int newElemSize, int newMaxElemNodes);

private:
	int					m_fmt;				//!< storage format
	int					m_maxElemNodes;		//!< max number of nodes for each element
	RgElementSet*		m_elset;			//!< the element set on which this map is defined

	std::vector<int>		m_NLT;		//!< node index lookup table for FMT_NODE
	int				m_imin;		//!< min index for lookup for FMT_NODE
};
