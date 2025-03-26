#pragma once

#include "femcore/FEObjectBase.h"
#include "FETimeInfo.h"

//-----------------------------------------------------------------------------
//! forward declaration of the FEModel class.
//! All classes inherited from FEModelComponent should take the model as a parameter
//! to the constructor.
class FENodeSet;
class FEMesh;

//-----------------------------------------------------------------------------
//! This class serves as a base class for many of the FECore classes. It defines
//! activation and deactivation functions which are used in multi-step analyses to determine which
//! components are active during an analysis.
//! A model component is basically anything that affects the state of a model.
//! For instance, boundary conditions, loads, contact definitions, etc.
class FEM_EXPORT FEModelComponent : public FEObjectBase
{
public:
	//! constructor
	FEModelComponent(FEModel* fem);

	//! destructor
	virtual ~FEModelComponent();

	//-----------------------------------------------------------------------------------
	//! Update the component
	//! This is called whenever the model is updated, i.e. the primary variables were updated.
	virtual void Update();

public: // some convenience functions (to pull data from FEModel without the need to include)
	double CurrentTime() const;
	double CurrentTimeIncrement() const;
	double GetGlobalConstant(const char* sz) const;
	int GetDOFIndex(const char* szvar, int n) const;
	int GetDOFIndex(const char* szdof) const;
	const FETimeInfo& GetTimeInfo() const;
	void AttachLoadController(const char* szparam, int lc);
	void AttachLoadController(void* pd, int lc);

	//! Get the model's mesh
	FEMesh& GetMesh();
};
