#pragma once
#include "femcore/fem_export.h"
#include <vector>

class FEModel;

class FEM_EXPORT FEModule
{
	class Impl;

public:
	enum Status {
		EXPERIMENTAL,
		RELEASED
	};

public:
	FEModule();
	FEModule(const char* szname, const char* szdescription = nullptr);

	virtual ~FEModule();

	// this function must be overridden by derived classes
	virtual void InitModel(FEModel* fem) = 0;

	void AddDependency(FEModule& mod);

	void ClearDependencies();

	std::vector<int> GetDependencies() const;

	int GetModuleID() const;

	const char* GetName() const;

	const char* GetDescription() const;

	int GetStatus() const;

	bool HasDependent(int modId) const;

protected:
	void SetStatus(FEModule::Status status);

private:
	void SetID(int newId);
	void SetName(const char* szname);
	void SetDescription(const char* szdesc);

public:
	Impl* im;
};
