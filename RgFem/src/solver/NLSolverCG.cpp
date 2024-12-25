#include"NLSolverCG.h"

NLSolverCG::NLSolverCG(ImplAssembly& assembly, double tol/*=10e-10*/, 
					   int maxIterations/*=100*/, int stps/*=20*/):NLSolver(assembly, tol, maxIterations,stps)
{
	setLSolver();
}
void NLSolverCG::setLSolver()
{
	lSolver=&cg;
}


template <class Facade, class...Base>
class BmBaseImpl : public Base...
{
    static_assert(std::is_class<Facade>::value, "Facade should be a class");

protected:
    BmBaseImpl() = default;
    virtual ~BmBaseImpl() = default;
    template <class T>
    auto f_facade()
    {
        Zw_ASSERT(m _pFacade);
        using Type = typename std::remove_cv<typename std::remove_pointer<T>::type>::type;
        return static_cast<Type*>(m _pFacade);
    }
    template <class T>
    auto f_facade() const
    {
        ZW_ASSERT(m_pFacade);
        using Type = typename std::remove_cv<typename std::remove is_pointer<T>::type>::type;
        return static_cast<Type*>(m pFacade);
    }

private:
    friend typename Facade;
    void setFacade(Facade* p)
    {
        m pFacade = p;
    }
    Facade* m _pFacade = nullptr;
};
        
