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

namespace fem
{
	double FEInterpolation::giveTransformationJacobian(const FloatArray& lcoords) const
	{
		FloatMatrix jacobianMatrix;
		this->giveJacobianMatrixAt(jacobianMatrix, lcoords);
		return jacobianMatrix.giveDeterminant();
	}

	//����Domain�����ͣ����û��ֹ���
	std::unique_ptr<IntegrationRule> FEInterpolation::giveIntegrationRule(int order) const
	{
		integrationDomain id = this->giveIntegrationDomain();
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

		int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
		iRule->SetUpPointsOnLine(points, _Unknown);
		return std::move(iRule);
	}

	std::unique_ptr<IntegrationRule> FEInterpolation::giveBoundaryIntegrationRule(int order, int boundary) const
	{
		integrationDomain id = this->giveBoundaryIntegrationDomain();
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

		int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
		iRule->setUpIntegrationPoints(id, points, _Unknown);
		return std::move(iRule);
	}

	std::unique_ptr<IntegrationRule>
		FEInterpolation::giveBoundaryEdgeIntegrationRule(int order, int boundary) const
	{
		integrationDomain id = this->giveBoundaryEdgeIntegrationDomain();
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

		int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
		iRule->setUpIntegrationPoints(id, points, _Unknown);
		return std::move(iRule);
	}

	std::unique_ptr<IntegrationRule>
		FEInterpolation::giveBoundarySurfaceIntegrationRule(int order, int boundary) const
	{
		integrationDomain id = this->giveBoundarySurfaceIntegrationDomain();
		auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);

		int points = iRule->getRequiredNumberOfIntegrationPoints(id, order + this->order);
		iRule->setUpIntegrationPoints(id, points, _Unknown);
		return std::move(iRule);
	}

} // end namespace fem
