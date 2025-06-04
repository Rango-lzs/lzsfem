
#pragma once
#include "FEObjectBase.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
// The FECoreTask class is the base class for all tasks.
// A task is simply the highest level module which defines what the code will do.
class FEM_EXPORT FECoreTask : public FEObjectBase
{
    DECLARE_META_CLASS(FECoreTask, FEObjectBase);

public:
	FECoreTask();
	virtual ~FECoreTask(void);

	//! initialize the task
	//! make sure to call FEModel::Init at some point
	virtual bool Init(const char* szfile) = 0;

	//! Run the task.
	virtual bool Run() = 0;
};
