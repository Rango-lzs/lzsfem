#pragma once
#include <vector>
#include <string>
#include <assert.h>
#include "FEDataMap.h"

//-----------------------------------------------------------------------------
class FESurface;
class FEFacetSet;
class DumpStream;

//-----------------------------------------------------------------------------
typedef int FEFacetIndex;


//-----------------------------------------------------------------------------
// TODO: Perhaps I should rename this FESurfaceData
//       and then define FESurfaceMap as a tool for evaluating data across a surface (i.e. via shape functions)
class FEM_EXPORT FESurfaceMap : public FEDataMap
{
public:
	//! default constructor
	FESurfaceMap();
	FESurfaceMap(FEDataType dataType);

	//! copy constructor
	FESurfaceMap(const FESurfaceMap& map);

	//! assignment operator
	FESurfaceMap& operator = (const FESurfaceMap& map);

	//! Create a surface data map for this surface
	bool Create(const FEFacetSet* ps, double val = 0.0, Storage_Fmt fmt = FMT_MULT);

	//! serialization
	void Serialize(DumpStream& ar) override;

	const FEFacetSet* GetFacetSet() const { return m_surf; }

	int MaxNodes() const { return m_maxFaceNodes; }

	// return the item list associated with this map
	FEItemList* GetItemList() override;

	int StorageFormat() const;

public: // from FEDataMap
	double value(const FEMaterialPoint& pt) override;
	Vector3d valueVec3d(const FEMaterialPoint& mp) override;
	Matrix3d valueMat3d(const FEMaterialPoint& mp) override;
	Matrix3ds valueMat3ds(const FEMaterialPoint& mp) override;

public:
	template <typename T> T value(int nface, int node)
	{
		return get<T>(nface*m_maxFaceNodes + node);
	}

	template <typename T> void setValue(int nface, int node, const T& v)
	{
		set<T>(nface*m_maxFaceNodes + node, v);
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
	const FEFacetSet*	m_surf;		// the surface for which this data set is defined
	int					m_format;	// the storage format
	int	m_maxFaceNodes;				// number of nodes for each face
};
