
#pragma once
#include "datastructure/Vector2d.h"
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix3d.h"
#include "FEParam.h"
#include "femcore/fem_export.h"
#include "femcore/FEParamValidator.h"

#include <assert.h>
#include <list>
#include <memory>
#include <stdio.h>

//-----------------------------------------------------------------------------
class DumpStream;
class FEParamObject;
class FEParamDouble;
class FEParamVec3;
class FEParamMat3d;
class FEParamMat3ds;
class FEDataArray;
class tens3drs;
class FEMaterialPointProperty;

//-----------------------------------------------------------------------------
typedef std::list<FEParam>::iterator FEParamIterator;
typedef std::list<FEParam>::const_iterator FEParamIteratorConst;

//-----------------------------------------------------------------------------
//! A list of material parameters
class FEM_EXPORT FEParameterList
{
public:
	FEParameterList(FEParamObject* pc);
	virtual ~FEParameterList();

	//! assignment operator
	void operator = (FEParameterList& l);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType itype, int ndim, const char* sz, bool* watch = nullptr);

	//! Add a parameter to the list
	FEParam* AddParameter(void* pv, FEParamType type, int ndim, FEParamRange rng, double fmin, double fmax, const char* sz);

	//! find a parameter using the data pointer
	FEParam* FindFromData(void* pv);

	//! find a parameter using its name (the safe way)
	FEParam* FindFromName(const char* sz);

	//! get a parameter (the dangerous way)
	FEParam& operator [] (const char* sz) { return *FindFromName(sz); }

	//! returs the first parameter
	FEParamIterator first() { return m_pl.begin(); }

	//! returs the first parameter
	FEParamIteratorConst first() const { return m_pl.begin(); }

	//! number of parameters
	int Parameters() const { return (int) m_pl.size(); }

	//! return the parameter container
	FEParamObject* GetContainer() { return m_pc; }

public:
	int SetActiveGroup(const char* szgroup);
	int GetActiveGroup();
	int ParameterGroups() const;
	const char* GetParameterGroupName(int i);

protected:
	FEParamObject*	m_pc;	//!< parent container
	std::list<FEParam>		m_pl;	//!< the actual parameter list
	std::vector<const char*>	m_pg;	//!< parameter groups
	int	m_currentGroup;	//!< active parameter group (new parameters are assigned to the current group; can be -1)
};
