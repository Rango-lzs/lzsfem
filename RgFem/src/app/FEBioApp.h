/*****************************************************************//**
 * \file   FEBioApp.h
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once


class FEModel;
class FEAppConfig;

class FEBioApp
{
public:
	FEBioApp();

	bool Init(int nargs, char* argv[]);

	bool Configure(const char* szconfig);

	int Run();

	void Finish();

	void ProcessCommands();

	// run an febio model
	int RunModel();

public:
	// get the current model
	FEBioModel* GetCurrentModel();

protected:
	// show FEBio prompt
	int prompt();

	// set the currently active model
	void SetCurrentModel(FEBioModel* fem);

	// apply configuration changes to model
	void ApplyConfig(FEBioModel& fem);

public:
	static FEBioApp* GetInstance();

private:

	FEAppConfig	m_config;			// configuration options

	FEModel*	m_fem;	
};
