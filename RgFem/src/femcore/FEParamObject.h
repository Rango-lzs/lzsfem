#pragma once

#include "femcore/MetaObject.h"
#include "femcore/RTTI/RTTIMacroDefine.h"

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

class DumpStream;
class FEParameterList;

//-----------------------------------------------------------------------------
//!每一个FEM对象都需要由一些参数定义，这些参数可能是通过输入文件解析，也可以手动设置
class FEM_EXPORT FEParamObject : public MetaObject
{
    DECLARE_META_CLASS(FEParamObject, MetaObject);
public:

    FEParamObject();
    virtual ~FEParamObject();

    //return the material's parameter list
    FEParameterList& GetParameterList();
    const FEParameterList& GetParameterList() const;

    //find a parameter using it's name
    virtual FEParam* GetParameter(const std::string& strName);
    virtual FEParam* FindParameter(const std::string& strName);

    //find a parameter using a pointer to the variable
    virtual FEParam* FindParameterFromData(void* pv);

    //serialize parameter data
    virtual void Serialize(DumpStream& ar);

    //validate material parameters.
    virtual bool Validate();

public:
    void BeginParameterGroup(const char* szname);
    void EndParameterGroup();

public:
    //! This function will be overridden by each class that defines a parameter list
    virtual void BuildParamList() = 0;

    //! Add a parameter to the list
    FEParam* AddParameter(void* pv, FEParamType itype, int ndim, const char* sz, bool* watch = nullptr);

    //! Add a parameter to the list
    FEParam* AddParameter(void* pv, FEParamType type, int ndim, RANGE rng, const char* sz);

public:
    FEParam* AddParameter(int& v, const char* sz);
    FEParam* AddParameter(bool& v, const char* sz);
    FEParam* AddParameter(double& v, const char* sz);
    FEParam* AddParameter(vec2d& v, const char* sz);
    FEParam* AddParameter(vec3d& v, const char* sz);
    FEParam* AddParameter(mat3d& v, const char* sz);
    FEParam* AddParameter(mat3ds& v, const char* sz);
    FEParam* AddParameter(FEParamDouble& v, const char* sz);
    FEParam* AddParameter(FEParamVec3& v, const char* sz);
    FEParam* AddParameter(FEParamMat3d& v, const char* sz);
    FEParam* AddParameter(FEParamMat3ds& v, const char* sz);
    FEParam* AddParameter(FEDataArray& v, const char* sz);
    FEParam* AddParameter(tens3drs& v, const char* sz);
    FEParam* AddParameter(std::string& v, const char* sz);
    FEParam* AddParameter(std::vector<int>& v, const char* sz);
    FEParam* AddParameter(std::vector<double>& v, const char* sz);
    FEParam* AddParameter(std::vector<vec2d>& v, const char* sz);
    FEParam* AddParameter(std::vector<std::string>& v, const char* sz);
    FEParam* AddParameter(FEMaterialPointProperty& v, const char* sz);

    FEParam* AddParameter(int& v, RANGE rng, const char* sz);
    FEParam* AddParameter(double& v, RANGE rng, const char* sz);
    FEParam* AddParameter(FEParamDouble& v, RANGE rng, const char* sz);

    FEParam* AddParameter(double& v, const char* sz, bool& watch);

    FEParam* AddParameter(int* v, int ndim, const char* sz);
    FEParam* AddParameter(double* v, int ndim, const char* sz);
    FEParam* AddParameter(FEParamDouble* v, int ndim, const char* sz);

    FEParam* AddParameter(int* v, int ndim, RANGE rng, const char* sz);
    FEParam* AddParameter(double* v, int ndim, RANGE rng, const char* sz);
    FEParam* AddParameter(FEParamDouble* v, int ndim, RANGE rng, const char* sz);

    FEParam* AddParameter(int& v, const char* sz, unsigned int flags, const char* szenum);
    FEParam* AddParameter(std::vector<int>& v, const char* sz, unsigned int flags, const char* szenum);
    FEParam* AddParameter(std::string& s, const char* sz, unsigned int flags, const char* szenum = nullptr);

    template <typename T>
    void SetParameter(const char* sz, T v);

private:
    FEParameterList* m_pParam;  //!< parameter list
};

//-----------------------------------------------------------------------------
template <typename T>
void FEParamObject::SetParameter(const char* sz, T v)
{
    FEParam* p = m_pParam->FindFromName(sz);
    p->value<T>() = v;
}

//-----------------------------------------------------------------------------
// To add parameter list to a class, simply do the following two steps
// 1) add the DECLARE_FECORE_CLASS macro in the material class declaration
// 2) use the BEGIN_FECORE_CLASS, ADD_PARAM and END_FECORE_CLASS to
//    define a parameter list

// the following macro declares the parameter list for a material
#define DECLARE_PARAM_LIST()                                                                                         \
public:                                                                                                                \
    void BuildParamList() override;

#define FECORE_BASE_CLASS(theClass)                                                                                    \
public:                                                                                                                \
    static const char* BaseClassName()                                                                                 \
    {                                                                                                                  \
        return #theClass;                                                                                              \
    }

// the BEGIN_FECORE_CLASS defines the beginning of a parameter list
#define BEGIN_PARAM_DEFINE(theClass, baseClass)                                                                        \
    void theClass::BuildParamList()                                                                                    \
    {                                                                                                                  \
        baseClass::BuildParamList();

// the ADD_PARAMETER macro adds a parameter to the parameter list
#define ADD_PARAMETER(...)   AddParameter(__VA_ARGS__)

// the END_FECORE_CLASS defines the end of a parameter list
#define END_PARAM_DEFINE()   }

// macro for starting a parameter group
#define BEGIN_PARAM_GROUP(a) BeginParameterGroup(a)

// macro for ending a parameter group
#define END_PARAM_GROUP()    EndParameterGroup()
