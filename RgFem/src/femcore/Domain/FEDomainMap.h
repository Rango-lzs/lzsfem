#pragma once
#include "femcore/FEDataMap.h"
#include "FEDomain.h"
#include "elements/FEElementSet.h"
#include "materials/FEMaterialPoint.h"

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
	bool Create(FEElementSet* ps, double val = 0.0);

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! get the value at a material point
	double value(const FEMaterialPoint& pt) override;
	Vector3d valueVec3d(const FEMaterialPoint& pt) override;
	Matrix3d valueMat3d(const FEMaterialPoint& pt) override;
	Matrix3ds valueMat3ds(const FEMaterialPoint& pt) override;

	//! Get the element set
	const FEElementSet* GetElementSet() const { return m_elset; }

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
	FEElementSet*		m_elset;			//!< the element set on which this map is defined

	std::vector<int>		m_NLT;		//!< node index lookup table for FMT_NODE
	int				m_imin;		//!< min index for lookup for FMT_NODE
};
