/*****************************************************************//**
 * \file   fe_interpol.cpp
 * \brief
 *
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#include "fe_interpol.h"
#include "element.h"
#include "gaussintegrationrule.h"

/*
插值模块需要知道高斯积分规则吗？
高斯积分规则，跟单元类型，单元阶数相关
*/

namespace fem
{
	double FEInterpolation::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		FloatMatrix jacobianMatrix;
		this->giveJacobianMatrixAt(jacobianMatrix, lcoords);
		return jacobianMatrix.giveDeterminant();
	}

	//根据Domain的类型，设置积分规则
	std::unique_ptr<IntegrationRule> FEInterpolation::giveIntegrationRule(int order) const
	{
		integrationDomain id = this->giveIntegrationDomain();
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

		int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
		iRule->SetUpPointsOnLine(points, _Unknown);
		return std::move(iRule);
	}
} // end namespace fem
