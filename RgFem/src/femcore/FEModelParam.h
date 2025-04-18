#pragma once
#include "femcore/FEScalarValuator.h"
#include "femcore/FEVec3dValuator.h"
#include "femcore/FEMat3dValuator.h"
#include "femcore/FEMat3dsValuator.h"
#include "femcore/FEItemList.h"

//---------------------------------------------------------------------------------------
// 定义模型的一些参数，买这些参数可以通过Valuator进行计算，有点类似表达式.
class FEM_EXPORT FEModelParam
{
public:
	FEModelParam();
	virtual ~FEModelParam();

	// set the domain
	void SetItemList(FEItemList* itemList) { m_dom = itemList; }

	// get the domain list
	FEItemList* GetItemList() { return m_dom;  }

	// set the scale factor
	void SetScaleFactor(double s) { m_scl = s; }

	// return the scale factor
	double GetScaleFactor() const { return m_scl; }

	// serialization
	virtual void Serialize(DumpStream& ar);

protected:
	double			m_scl;	//!< scale factor. Used to store load curve value
	FEItemList*		m_dom;
};

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEParamDouble : public FEModelParam
{
public:
	FEParamDouble();
	~FEParamDouble();

	FEParamDouble(const FEParamDouble& p);

	// set the value
	void operator = (double v);
	void operator = (const FEParamDouble& p);

	// set the valuator
	void setValuator(FEScalarValuator* val);

	// get the valuator
	FEScalarValuator* valuator();

	// evaluate the parameter at a material point
	double operator () (const FEMaterialPoint& pt) { return m_scl*(*m_val)(pt); }

	// is this a const value
	bool isConst() const;

	// get the const value (return value undefined if param is not const)
	double& constValue();


//=======================================================================================

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEParamVec3 : public FEModelParam
{
public:
	FEParamVec3();
	~FEParamVec3();

	FEParamVec3(const FEParamVec3& p);

	bool Init();

	// set the value
	void operator = (const Vector3d& v);
	void operator = (const FEParamVec3& p);

	// set the valuator
	void setValuator(FEVec3dValuator* val);
	FEVec3dValuator* valuator();

	// evaluate the parameter at a material point
	Vector3d operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// return a unit vector
	Vector3d unitVector(const FEMaterialPoint& pt) { return (*this)(pt).normalized(); }

	// is this a const
	bool isConst() const { return m_val->isConst(); }

	// (return value undefined if param is not const)
	Vector3d& constValue() { assert(isConst()); return *m_val->constValue(); };

	void Serialize(DumpStream& ar) override;

private:
	FEVec3dValuator*	m_val;
};

//=======================================================================================

//---------------------------------------------------------------------------------------
class FEM_EXPORT FEParamMat3d : public FEModelParam
{
public:
	FEParamMat3d();
	~FEParamMat3d();

	FEParamMat3d(const FEParamMat3d& p);

	// set the value
	void operator = (const Matrix3d& v);
	void operator = (const FEParamMat3d& v);

	bool Init();

	// set the valuator
	void setValuator(FEMat3dValuator* val);

	// get the valuator
	FEMat3dValuator* valuator();

	// evaluate the parameter at a material point
	Matrix3d operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// is this a const
	bool isConst() const { return m_val->isConst(); }

	// (return value undefined if not constant)
	Matrix3d& constValue() { assert(isConst()); return *m_val->constValue(); };

	void Serialize(DumpStream& ar) override;

private:
	FEMat3dValuator*	m_val;
};


//---------------------------------------------------------------------------------------
class FEM_EXPORT FEParamMat3ds : public FEModelParam
{
public:
	FEParamMat3ds();
	~FEParamMat3ds();

	FEParamMat3ds(const FEParamMat3ds& p);

	// set the value
	void operator = (const mat3ds& v);
	void operator = (const FEParamMat3ds& v);

	// set the valuator
	void setValuator(FEMat3dsValuator* val);

	// get the valuator
	FEMat3dsValuator* valuator();

	// evaluate the parameter at a material point
	mat3ds operator () (const FEMaterialPoint& pt) { return (*m_val)(pt)*m_scl; }

	// is this a const
	bool isConst() const { return m_val->isConst(); }

	// (return value undefined if not constant)
	mat3ds& constValue() { assert(isConst()); return *m_val->constValue(); };

	void Serialize(DumpStream& ar) override;

private:
	FEMat3dsValuator*	m_val;
};
