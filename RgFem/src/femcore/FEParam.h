/*****************************************************************//**
 * \file   FEParam.h
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once
#include "datastructure/vec3d.h"
#include "datastructure/mat3d.h"
#include "femcore/fem_export.h"

#include <assert.h>
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEParamValidator;
class DumpStream;
class FEParamContainer;

//FEParam一般是由FEObject所持有的，比较持久化的参数, 具有更丰富的信息，可能是数组类型
//FEParamList对其中的FEParam生命周期负责，FEParam对其中具体value的生命周期负责
//FEParamValue参数里的具体一个值，一般用于运行时取值
//做减法，而不是做加法

//-----------------------------------------------------------------------------
// Different supported parameter types
enum FEParamType {
	FE_PARAM_INVALID,
	FE_PARAM_INT,
	FE_PARAM_BOOL,
	FE_PARAM_DOUBLE,
	FE_PARAM_VEC2D,
	FE_PARAM_VEC3D,
	FE_PARAM_MAT3D,
	FE_PARAM_MAT3DS,
	FE_PARAM_STRING,
	FE_PARAM_DATA_ARRAY,
	FE_PARAM_TENS3DRS,
	FE_PARAM_STD_STRING,
	FE_PARAM_STD_VECTOR_INT,
	FE_PARAM_STD_VECTOR_DOUBLE,
	FE_PARAM_STD_VECTOR_VEC2D,
	FE_PARAM_STD_VECTOR_STRING,
	FE_PARAM_DOUBLE_MAPPED,
	FE_PARAM_VEC3D_MAPPED,
	FE_PARAM_MAT3D_MAPPED,
	FE_PARAM_MAT3DS_MAPPED,
	FE_PARAM_MATERIALPOINT
};

// class describing the value of parameter
class FEParamValue
{
private:
    void* m_pv;           // pointer to variable data
    FEParamType m_itype;  // type of variable (this is not the type of the param!)
    FEParam* m_param;     // the parameter (can be null if it is not a parameter)

public:
    FEParamValue()
    {
        m_pv = 0;
        m_itype = FE_PARAM_INVALID;
        m_param = 0;
    }

    explicit FEParamValue(FEParam* p, void* v, FEParamType itype)
    {
        m_pv = v;
        m_itype = itype;
        m_param = p;
    }

	FEParamValue(bool& v)
        : FEParamValue(0, &v, FE_PARAM_BOOL)
    {
    }

	FEParamValue(int& v)
        : FEParamValue(0, &v, FE_PARAM_INT)
    {
    }

    FEParamValue(double& v)
        : FEParamValue(0, &v, FE_PARAM_DOUBLE)
    {
    }

    FEParamValue(vec2d& v)
        : FEParamValue(0, &v, FE_PARAM_VEC2D)
    {
    }

    FEParamValue(vec3d& v)
        : FEParamValue(0, &v, FE_PARAM_VEC3D)
    {
    }

    FEParamValue(mat3ds& v)
        : FEParamValue(0, &v, FE_PARAM_MAT3DS)
    {
    }

    FEParamValue(mat3d& v)
        : FEParamValue(0, &v, FE_PARAM_MAT3D)
    {
    }

    bool isValid() const
    {
        return (m_pv != 0);
    }

    FEParamType type() const
    {
        return m_itype;
    }

    void* data_ptr() const
    {
        return m_pv;
    }

    FEParam* param()
    {
        return m_param;
    }

    template <typename T>
    T& value()
    {
        return *((T*)m_pv);
    }

    template <typename T>
    const T& value() const
    {
        return *((T*)m_pv);
    }

    FEM_EXPORT FEParamValue component(int n);
};

//-----------------------------------------------------------------------------
//! This class describes a user-defined parameter  ,是否能用Variant替代，
class FEM_EXPORT FEParam
{
private:
	void*			m_pv;		// pointer to variable data
	int				m_dim;		// dimension (in case data is array)
	FEParamType		m_type;		// type of variable
	unsigned int	m_flag;		// parameter flags
	bool*			m_watch;	// parameter watch (set to true if read in)
	int				m_group;	// index of parameter group (-1 by default)

	const char*	m_szname;	// name of the parameter
	const char*	m_szenum;	// enumerate values for ints
	const char* m_szunit;	// unit string
	const char* m_szlongname;	// a longer, more descriptive name (optional)

	// parameter validator
	FEParamValidator*	m_pvalid;

	FEParamContainer* m_parent;	// parent object of parameter

public:
	// constructor
	FEParam(void* pdata, FEParamType itype, int ndim, const char* szname, bool* watch = nullptr);
	FEParam(const FEParam& p);
	~FEParam();
	FEParam& operator = (const FEParam& p);

	// set the parameter's validator
	void SetValidator(FEParamValidator* pvalid);

	// see if the parameter's value is valid
	bool is_valid() const;

	// return the name of the parameter
	const char* name() const;

	// return the long name of the parameter
	const char* longName() const;

	// return the enum values
	const char* enums() const;

	// get the current enum value (or nullptr)
	const char* enumKey() const;

	// get the unit string
	const char* units() const;
	FEParam* setUnits(const char* szunit);

	// set the enum values (\0 separated. Make sure the end of the string has two \0's)
	FEParam* setEnums(const char* sz);

	// set the long name
	FEParam* setLongName(const char* sz);

	// parameter dimension
	int dim() const;

	// parameter type
	FEParamType type() const;

	// data pointer
	void* data_ptr() const;

	// get the param value
	FEParamValue paramValue(int i = -1);

	// Copy the state of one parameter to this parameter.
	// This requires that the parameters are compatible (i.e. same type, etc.)
	bool CopyState(const FEParam& p);

	void setParent(FEParamContainer* pc);
	FEParamContainer* parent();

	FEParam* SetFlags(unsigned int flags);
	unsigned int GetFlags() const;

	void SetWatchVariable(bool* watchVar);
	bool* GetWatchVariable();
	void SetWatchFlag(bool b);

	bool IsHidden() const;

	bool IsVolatile() const;

	FEParam* MakeVolatile(bool b);

	bool IsTopLevel() const;
	FEParam* MakeTopLevel(bool b);

public:
	int GetParamGroup() const;
	void SetParamGroup(int i);

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FEParam* p);
	static FEParam* LoadClass(DumpStream& ar, FEParam* p);

public:
	//! retrieves the value for a non-array item
	template <class T> T& value() { return *((T*) data_ptr()); }

	//! retrieves the value for a non-array item
	template <class T> const T& value() const { return *((T*) data_ptr()); }

	//! retrieves the value for an array item
	template <class T> T* pvalue() { return (T*) data_ptr(); }

	//! retrieves the value for an array item
	template <class T> T& value(int i) { return ((T*)data_ptr())[i]; }
	template <class T> T value(int i) const { return ((T*) data_ptr())[i]; }

	//! retrieves pointer to element in array
	template <class T> T* pvalue(int n);

	//! override the template for char pointers
	char* cvalue();
};

//-----------------------------------------------------------------------------
//! Retrieves a pointer to element in array
template<class T> inline T* FEParam::pvalue(int n)
{
	assert((n >= 0) && (n < m_dim));
	return &(pvalue<T>()[n]);
}

//-----------------------------------------------------------------------------
FEM_EXPORT FEParamValue GetParameterComponent(const std::string& paramName, FEParam* param);
