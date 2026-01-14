/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#pragma once
#include "FEVec3dValuator.h"
#include "FEModelParam.h"

//---------------------------------------------------------------------------------------
// A constant valuator
class FEM_EXPORT FEConstValueVec3 : public FEVec3dValuator
{
public:
	FEConstValueVec3();

	FEVec3dValuator* copy() override;

	Vector3d operator()(const RgMaterialPoint& pt) override { return m_val; }

	// is this a const value
	bool isConst() override { return true; }

	// get the const value (returns 0 if param is not const)
	Vector3d* constValue() override { return &m_val; }

	Vector3d& value() { return m_val; }

	void setConstValue(const Vector3d& v) { m_val = v; }

private:
	Vector3d	m_val;

	DECLARE_PARAM_LIST();
};

//---------------------------------------------------------------------------------------
// The value is calculated using a mathematical expression
class FEM_EXPORT FEMathValueVec3 : public FEVec3dValuator
{
public:
	FEMathValueVec3();
	Vector3d operator()(const RgMaterialPoint& pt) override;

	bool Init() override;

	bool create(const std::string& sx, const std::string& sy, const std::string& sz);

	FEVec3dValuator* copy() override;

	bool UpdateParams() override;

private:
	std::string			m_expr;
	//FEMathExpression	m_math[3];

	DECLARE_PARAM_LIST();
};

//---------------------------------------------------------------------------------------
// The value is determined by a data map
class FEM_EXPORT FEMappedValueVec3 : public FEVec3dValuator
{
public:
	FEMappedValueVec3();

	void setDataMap(FEDataMap* val, Vector3d scl = Vector3d(1, 1, 1));

	Vector3d operator()(const RgMaterialPoint& pt) override;

	FEVec3dValuator* copy() override;

	void Serialize(DumpStream& ar) override;

	bool Init() override;

private:
	std::string		m_mapName;

private:
	FEDataMap*		m_val;

	DECLARE_PARAM_LIST();
};


//-----------------------------------------------------------------------------
// This class calculates a vector based on local element node numbering
class FEM_EXPORT FELocalVectorGenerator : public FEVec3dValuator
{
public:
	FELocalVectorGenerator();

	bool Init() override;

	Vector3d operator () (const RgMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	int	m_n[2];

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FECylindricalVectorGenerator : public FEVec3dValuator
{
public:
	FECylindricalVectorGenerator();

	bool Init() override;

	Vector3d operator () (const RgMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	Vector3d	m_center;		// center of map
	Vector3d	m_axis;			// cylinder axis
	Vector3d	m_vector;		// reference direction

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FESphericalVectorGenerator : public FEVec3dValuator
{
public:
	FESphericalVectorGenerator();

	bool Init() override;

	Vector3d operator () (const RgMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	Vector3d	m_center;
	Vector3d	m_vector;

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FESphericalAnglesVectorGenerator : public FEVec3dValuator
{
public:
	FESphericalAnglesVectorGenerator();

	Vector3d operator () (const RgMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

protected:
	FEParamDouble	m_theta;	// in-plane (x,y) angle from x-axis
	FEParamDouble	m_phi;		// angle from z-axis

	DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
// This class is mostly to support some older formats that used the "user" fiber
// generator option. It actually doesn't generate any vectors and should not be used. 
class FEM_EXPORT FEUserVectorGenerator : public FEVec3dValuator
{
public:
	FEUserVectorGenerator();

	Vector3d operator () (const RgMaterialPoint& mp) override;

	FEVec3dValuator* copy() override;

	DECLARE_PARAM_LIST();
};
