#pragma once
#include "MItem.h"
#include <vector>
#include "femcore/fem_export.h"

//-----------------------------------------------------------------------------
typedef std::vector<MVariable*>	MVarList;

//-----------------------------------------------------------------------------
// This class defines the base class for all math objects
// It also stores a list of all the variables
class FEM_EXPORT MathObject
{
public:
	MathObject();
	MathObject(const MathObject& mo);
	void operator = (const MathObject& mo);

	virtual ~MathObject();

	int Dim() { return (int)m_Var.size(); }

	MVariable* AddVariable(const std::string& var, double initVal = 0.0);
	void AddVariables(const std::vector<std::string>& varList);

	void AddVariable(MVariable* pv);
	MVariable* FindVariable(const std::string& s);
	int Variables() const { return (int)m_Var.size(); }

	MVariable* Variable(int i) { return m_Var[i]; }
	const MVariable* Variable(int i) const { return m_Var[i]; }

	virtual void Clear();

	virtual MathObject* copy() = 0;

protected:
	MVarList	m_Var;		// list of variables
};

//-----------------------------------------------------------------------------
// This class defines a simple epxression that can be evaluated by
// setting the values of the variables.
class FEM_EXPORT MSimpleExpression : public MathObject
{
public:
	MSimpleExpression() {}
	MSimpleExpression(const MSimpleExpression& mo);
	void operator = (const MSimpleExpression& mo);

	void SetExpression(MITEM& e) { m_item = e; }
	MITEM& GetExpression() { return m_item; }
	const MITEM& GetExpression() const { return m_item; }

	// Create a simple expression object from a string
	bool Create(const std::string& expr, bool autoVars = false);

	// copy the expression
	MathObject* copy() { return new MSimpleExpression(*this); }

	// These functions are not thread safe since variable values can be overridden by different threads
	// In multithreaded applications, use the thread safe functions below.
	double value() const { return value(m_item.ItemPtr());  }

	// combines Create and value. Not efficient usage! 
	double value(const std::string& s);

	// This is a thread safe function to evaluate the expression
	// The values of the variables are passed as an argument. This function
	// does not call MVariable->value, but uses these passed values insteads.
	// Make sure that the var array has the same size as the variable array of the expression
	double value_s(const std::vector<double>& var) const
	{ 
		assert(var.size() == m_Var.size());
		return value(m_item.ItemPtr(), var); 
	}

	int Items();

protected:
	double value(const MItem* pi) const;
	double value(const MItem* pi, const std::vector<double>& var) const;

protected:
	void fixVariableRefs(MItem* pi);

protected:
	MITEM	m_item;
};
