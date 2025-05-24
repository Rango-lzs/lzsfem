#pragma once
#include "femcore/fem_export.h"
#include "input/XML/XMLReader.h"

class FEClassDescriptor;
class FEObjectBase;
class FEParameterList;

//This namespace defines some helper functions that facilitate processing the FEBio xml formatted files.
namespace fexml
{
//---------------------------------------------------------------------------------------
// Reads the value of a parameter.
// if paramName is zero, the tag's name will be used as the parameter name.
bool FEM_EXPORT readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName = 0);

//---------------------------------------------------------------------------------------
bool FEM_EXPORT readParameter(XMLTag& tag, FEObjectBase* pc);

//---------------------------------------------------------------------------------------
// reads the parameters and properties of a FECore class
bool FEM_EXPORT readParameterList(XMLTag& tag, FEObjectBase* pc);

//---------------------------------------------------------------------------------------
// read a list of integers
void FEM_EXPORT readList(XMLTag& tag, std::vector<int>& l);

//---------------------------------------------------------------------------------------
// create a class descriptor from the current tag
FEM_EXPORT FEClassDescriptor* readParameterList(XMLTag& tag);

}
