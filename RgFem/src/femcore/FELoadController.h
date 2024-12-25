/*********************************************************************
 * \file   FELoadController.h
 * \brief  
 * 
 * \author Leizs
 * \date   December 2024
 *********************************************************************/

#pragma once
#include "FEModelComponent.h"

//-----------------------------------------------------------------------------
// Class that describes a load controller. A load controller can modify the value
// of model parameters during the analysis. 
class FECORE_API FELoadController : public FEModelComponent
{
	FECORE_SUPER_CLASS(FELOADCONTROLLER_ID)
	FECORE_BASE_CLASS(FELoadController);

public:
	FELoadController(FEModel* fem);

	//! evaluate the load controller 
	void Evaluate(double time);

	//! return the last calculated value
	double Value() const { return m_value; }

	//! serialization
	void Serialize(DumpStream& ar) override;

protected:
	// This must be implemented by derived classes
	virtual double GetValue(double time) = 0;

private:
	double	m_value;	//!< last calculated value
};
