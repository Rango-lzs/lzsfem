/*****************************************************************//**
 * \file   FEObjectBase.h
 * \brief  
 * 
 * \author 11914
 * \date   December 2024
 *********************************************************************/

#pragma once

#include "femcore/meta_object.h"

#include "FEParameterList.h"
#include "fecore_enum.h"
#include "FEProperty.h"
#include <string>

//Use this clas to define the base impl class for the convenient of set and get facede class.
template <class Facade, class... Base>
class BaseImpl : public Base...
{
    static_assert(std::is_class<Facade>::value, "Facade should be a class");

protected:
    BaseImpl() = default;
    virtual ~BaseImpl() = default;
    template <class T>
    auto f_facade()
    {
        std::assert(m_pFacade);
        using Type = typename std::remove_cv<typename std::remove_pointer<T>::type>::type;
        return static_cast<Type*>(m_pFacade);
    }
    template <class T>
    auto f_facade() const
    {
        std::assert(m_pFacade);
        using Type = typename std::remove_cv<typename std::remove is_pointer<T>::type>::type;
        return static_cast<Type*>(m_pFacade);
    }

private:
    friend typename Facade;
    void setFacade(Facade* p)
    {
        m_pFacade = p;
    }
    Facade* m_pFacade = nullptr;
};


//-----------------------------------------------------------------------------
class FECoreFactory;

//-----------------------------------------------------------------------------
//! Base class for most classes in FECore library and the base class for all 
//! classes that can be registered with the framework.
class FEM_EXPORT FEObjectBase : public MetaObject
{
public:
	//! constructor
    explicit FEObjectBase(FEModel* pModel);

	//! destructor
	virtual ~FEObjectBase();

	//! Initialization
	virtual bool Init();
	
public:
	//! number of parameters
	int Parameters() const;
	//! return a parameter
	virtual FEParam* FindParameter(const ParamString& s) override;
	//! return the property (or this) that owns a parameter
	FEObjectBase* FindParameterOwner(void* pd);
	//! validates all properties and parameters
	bool Validate() override;
	//! call this after the parameters are changed
	virtual bool UpdateParams();

public:
	
	//! Get the FE model
	FEModel* GetFEModel() const;
	//! set the FEModel of this class (use with caution!)
	void SetFEModel(FEModel* fem);

	static void SaveClass(DumpStream& ar, FEObjectBase* p);
	static FEObjectBase* LoadClass(DumpStream& ar, FEObjectBase* p);
	//! data serialization
	void Serialize(DumpStream& ar) override;

private:
	std::string		m_name;			//!< user defined name of component
	FEModel*		m_fem;			//!< the model this class belongs to
	FEParamContainer* mp_param_cont;
	int		m_nID;			//!< component ID
};

