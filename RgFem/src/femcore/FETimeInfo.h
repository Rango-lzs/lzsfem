#pragma once
#include "femcore/fem_export.h"
class DumpStream;

//-----------------------------------------------------------------------------
class FEM_EXPORT FETimeInfo
{
public:
	FETimeInfo();

	FETimeInfo(double time, double tinc);

	void Serialize(DumpStream& ar);

public:
	double	currentTime;		//!< current time value
	double	timeIncrement;		//!< current time step (difference between this time and previous one)
	int		timeStep;			//!< current time step
	int		currentIteration;	//!< iteration number
	int		augmentation;		//!< augmentation

	// HHT time integration parameters
	double	alpha;
	double	beta;
	double	gamma;
	double  alphaf;
	double  alpham;
};
