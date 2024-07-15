/*****************************************************************//**
 * \file   fem_factory.cpp
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#include "fem_factory.h"
#include "engngm.h"

#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"

#include "sparsemtrx.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparselinsystemnm.h"

#include "function.h"
#include "material.h"
#include "crosssection.h"
#include "exportmodule.h"

#include "gaussintegrationrule.h"
#include "initialcondition.h"

#include <string>
#include <algorithm>
#include <cctype>

namespace fem 
{
	FemFactory& instance()
	{
		static FemFactory ans;
		return ans;
	}

	std::string conv2lower(std::string input)
	{
		std::transform(input.begin(), input.end(), input.begin(), ::tolower);
		return input;
	}

	// Non-string names (should be changed eventually):
	template< typename C, typename T, typename V, typename ... As> C* cf_create2(const T& list, V name, As ... args)
	{
		auto creator = list.find(name);
		return creator != list.end() ? creator->second(args...) : nullptr;
	}

	// Helper for storing creators
	template< typename T, typename V, typename C> bool cf_store2(T& list, V name, C& creator)
	{
		list[name] = creator;
		return true;
	}

	// Non string names (should be replaced with eventually)
	template< typename C, typename T, typename V, typename ... As> std::unique_ptr<C> cf_create4(const T& list, V name, As ... args)
	{
		auto creator = list.find(name);
		return creator != list.end() ? creator->second(args...) : nullptr;
	}

	// Helper for creating objects
	template< typename C, typename T, typename ... As> std::unique_ptr<C> cf_create(const T& list, const char* name, As ... args)
	{
		auto creator = list.find(conv2lower(name));
		return creator != list.end() ? creator->second(args...) : nullptr;
	}

	// Helper for storing creators
	template< typename T, typename C> bool cf_store(T& list, const char* name, C& creator)
	{
		list[conv2lower(name)] = creator;
		return true;
	}


	FemFactory::FemFactory()
	{
		// Fixed list for DOF types. No new components can register for these since these are part of the internal structure in OOFEM.
		dofList[DT_master] = dofCreator< MasterDof >;
		dofList[DT_simpleSlave] = dofCreator< SimpleSlaveDof >;
		dofList[DT_slave] = dofCreator< SlaveDof >;
		dofList[DT_active] = dofCreator< ActiveDof >;
	}

	std::unique_ptr<SparseMtrx> FemFactory::createSparseMtrx(SparseMtrxType name)
	{
		return cf_create4<SparseMtrx>(sparseMtrxList, name);
	}

	bool FemFactory::registerSparseMtrx(SparseMtrxType name, std::unique_ptr<SparseMtrx>(*creator)())
	{
		return cf_store2(sparseMtrxList, name, creator);
	}

	Dof* FemFactory::createDof(dofType name, DofIDItem dofid, DofManager* dman)
	{
		return cf_create2<Dof>(dofList, name, dofid, dman);
	}

	std::unique_ptr<SparseLinearSystemNM> FemFactory::createSparseLinSolver(LinSystSolverType name, Domain* domain, EngngModel* emodel)
	{
		return cf_create4<SparseLinearSystemNM>(sparseLinSolList, name, domain, emodel);
	}

	bool FemFactory::registerSparseLinSolver(LinSystSolverType name, std::unique_ptr<SparseLinearSystemNM>(*creator)(Domain*, EngngModel*))
	{
		return cf_store2(sparseLinSolList, name, creator);
	}

	std::unique_ptr<InitialCondition> FemFactory::createInitialCondition(const char* name, int number, Domain* domain)
	{
		if (conv2lower(name).compare("initialcondition") == 0) {
			return std::make_unique<InitialCondition>(number, domain);
		}
		return nullptr;
	}

	std::unique_ptr<Element> FemFactory::createElement(const char* name, int number, Domain* domain)
	{
		return cf_create<Element>(elemList, name, number, domain);
	}

	bool FemFactory::registerElement(const char* name, std::unique_ptr<Element>(*creator)(int, Domain*))
	{
		return cf_store(elemList, name, creator);
	}

	std::unique_ptr<DofManager> FemFactory::createDofManager(const char* name, int number, Domain* domain)
	{
		return cf_create<DofManager>(dofmanList, name, number, domain);
	}

	bool FemFactory::registerDofManager(const char* name, std::unique_ptr<DofManager>(*creator)(int, Domain*))
	{
		return cf_store(dofmanList, name, creator);
	}

	std::unique_ptr<GeneralBoundaryCondition> FemFactory::createBoundaryCondition(const char* name, int number, Domain* domain)
	{
		return cf_create<GeneralBoundaryCondition>(bcList, name, number, domain);
	}

	bool FemFactory::registerBoundaryCondition(const char* name, std::unique_ptr<GeneralBoundaryCondition>(*creator)(int, Domain*))
	{
		return cf_store(bcList, name, creator);
	}

	std::unique_ptr<CrossSection> FemFactory::createCrossSection(const char* name, int number, Domain* domain)
	{
		return cf_create<CrossSection>(csList, name, number, domain);
	}

	bool FemFactory::registerCrossSection(const char* name, std::unique_ptr<CrossSection>(*creator)(int, Domain*))
	{
		return cf_store(csList, name, creator);
	}

	std::unique_ptr<Material> FemFactory::createMaterial(const char* name, int number, Domain* domain)
	{
		return cf_create<Material>(matList, name, number, domain);
	}

	bool FemFactory::registerMaterial(const char* name, std::unique_ptr<Material>(*creator)(int, Domain*))
	{
		return cf_store(matList, name, creator);
	}

	std::unique_ptr<EngngModel> FemFactory::createEngngModel(const char* name, int number, EngngModel* master)
	{
		return cf_create<EngngModel>(engngList, name, number, master);
	}

	bool FemFactory::registerEngngModel(const char* name, std::unique_ptr<EngngModel>(*creator)(int, EngngModel*))
	{
		return cf_store(engngList, name, creator);
	}

	std::unique_ptr<Function> FemFactory::createFunction(const char* name, int number, Domain* domain)
	{
		return cf_create<Function>(funcList, name, number, domain);
	}

	bool FemFactory::registerFunction(const char* name, std::unique_ptr<Function>(*creator)(int, Domain*))
	{
		return cf_store(funcList, name, creator);
	}

	std::unique_ptr<NonlocalBarrier> FemFactory::createNonlocalBarrier(const char* name, int number, Domain* domain)
	{
		return cf_create<NonlocalBarrier>(nlbList, name, number, domain);
	}

	bool FemFactory::registerNonlocalBarrier(const char* name, std::unique_ptr<NonlocalBarrier>(*creator)(int, Domain*))
	{
		return cf_store(nlbList, name, creator);
	}

	std::unique_ptr<ExportModule> FemFactory::createExportModule(const char* name, int number, EngngModel* emodel)
	{
		return cf_create<ExportModule>(exportList, name, number, emodel);
	}

	bool FemFactory::registerExportModule(const char* name, std::unique_ptr<ExportModule>(*creator)(int, EngngModel*))
	{
		return cf_store(exportList, name, creator);
	}

	std::unique_ptr<Monitor> FemFactory::createMonitor(const char* name, int number)
	{
		return cf_create<Monitor>(monitorList, name, number);
	}

	bool FemFactory::registerMonitor(const char* name, std::unique_ptr<Monitor>(*creator)(int))
	{
		return cf_store(monitorList, name, creator);
	}

	std::unique_ptr<SparseNonLinearSystemNM>FemFactory::createNonLinearSolver(const char* name, Domain* domain, EngngModel* emodel)
	{
		return cf_create<SparseNonLinearSystemNM>(nonlinList, name, domain, emodel);
	}

	bool FemFactory::registerSparseNonLinearSystemNM(const char* name, std::unique_ptr<SparseNonLinearSystemNM>(*creator)(Domain*, EngngModel*))
	{
		return cf_store(nonlinList, name, creator);
	}

	std::unique_ptr<InitModule> FemFactory::createInitModule(const char* name, int number, EngngModel* emodel)
	{
		return cf_create<InitModule>(initList, name, number, emodel);
	}

	bool FemFactory::registerInitModule(const char* name, std::unique_ptr<InitModule>(*creator)(int, EngngModel*))
	{
		return cf_store(initList, name, creator);
	}

	std::unique_ptr<TopologyDescription> FemFactory::createTopology(const char* name, Domain* domain)
	{
		return cf_create<TopologyDescription>(topologyList, name, domain);
	}

	bool FemFactory::registerTopologyDescription(const char* name, std::unique_ptr<TopologyDescription>(*creator)(Domain*))
	{
		return cf_store(topologyList, name, creator);
	}


	// XFEM:
	std::unique_ptr<EnrichmentItem> FemFactory::createEnrichmentItem(const char* name, int number, XfemManager* xm, Domain* domain)
	{
		return cf_create<EnrichmentItem>(enrichItemList, name, number, xm, domain);
	}

	bool FemFactory::registerEnrichmentItem(const char* name, std::unique_ptr<EnrichmentItem>(*creator)(int, XfemManager*, Domain*))
	{
		return cf_store(enrichItemList, name, creator);
	}

	std::unique_ptr<NucleationCriterion> FemFactory::createNucleationCriterion(const char* name, Domain* domain)
	{
		return cf_create<NucleationCriterion>(nucleationCritList, name, domain);
	}

	bool FemFactory::registerNucleationCriterion(const char* name, std::unique_ptr<NucleationCriterion>(*creator)(Domain*))
	{
		return cf_store(nucleationCritList, name, creator);
	}

	std::unique_ptr<EnrichmentFunction> FemFactory::createEnrichmentFunction(const char* name, int number, Domain* domain)
	{
		return cf_create<EnrichmentFunction>(enrichFuncList, name, number, domain);
	}

	bool FemFactory::registerEnrichmentFunction(const char* name, std::unique_ptr<EnrichmentFunction>(*creator)(int, Domain*))
	{
		return cf_store(enrichFuncList, name, creator);
	}

#if 0
	std::unique_ptr<EnrichmentDomain> FemFactory::createEnrichmentDomain(const char* name)
	{
		return cf_create<EnrichmentDomain>(enrichmentDomainList, name);
	}

	bool FemFactory::registerEnrichmentDomain(const char* name, std::unique_ptr<EnrichmentDomain>(*creator)())
	{
		return cf_store(enrichmentDomainList, name, creator);
	}
#endif

	std::unique_ptr<EnrichmentFront> FemFactory::createEnrichmentFront(const char* name)
	{
		return cf_create<EnrichmentFront>(enrichmentFrontList, name);
	}

	bool FemFactory::registerEnrichmentFront(const char* name, std::unique_ptr<EnrichmentFront>(*creator)())
	{
		return cf_store(enrichmentFrontList, name, creator);
	}

	std::unique_ptr<PropagationLaw> FemFactory::createPropagationLaw(const char* name)
	{
		return cf_create<PropagationLaw>(propagationLawList, name);
	}

	bool FemFactory::registerPropagationLaw(const char* name, std::unique_ptr<PropagationLaw>(*creator)())
	{
		return cf_store(propagationLawList, name, creator);
	}

	std::unique_ptr<BasicGeometry> FemFactory::createGeometry(const char* name)
	{
		return cf_create<BasicGeometry>(geometryList, name);
	}

	bool FemFactory::registerGeometry(const char* name, std::unique_ptr<BasicGeometry>(*creator)())
	{
		return cf_store(geometryList, name, creator);
	}

	std::unique_ptr<XfemManager> FemFactory::createXfemManager(const char* name, Domain* domain)
	{
		return cf_create<XfemManager>(xManList, name, domain);
	}

	bool FemFactory::registerXfemManager(const char* name, std::unique_ptr<XfemManager>(*creator)(Domain*))
	{
		return cf_store(xManList, name, creator);
	}


	// Failure module:

	std::unique_ptr<FailureCriteria> FemFactory::createFailureCriteria(const char* name, int number, FractureManager* fracManager)
	{
		return cf_create<FailureCriteria>(failureCriteriaList, name, number, fracManager);
	}

	bool FemFactory::registerFailureCriteria(const char* name, std::unique_ptr<FailureCriteria>(*creator)(int, FractureManager*))
	{
		return cf_store(failureCriteriaList, name, creator);
	}

	std::unique_ptr<FailureCriteriaStatus> FemFactory::createFailureCriteriaStatus(const char* name, int number, FailureCriteria* fc)
	{
		return cf_create<FailureCriteriaStatus>(failureCriteriaStatusList, name, number, fc);
	}

	bool FemFactory::registerFailureCriteriaStatus(const char* name, std::unique_ptr<FailureCriteriaStatus>(*creator)(int, FailureCriteria*))
	{
		return cf_store(failureCriteriaStatusList, name, creator);
	}


	std::unique_ptr<ContactManager> FemFactory::createContactManager(const char* name, Domain* domain)
	{
		return cf_create<ContactManager>(contactManList, name, domain);
	}

	bool FemFactory::registerContactManager(const char* name, std::unique_ptr<ContactManager>(*creator)(Domain*))
	{
		return cf_store(contactManList, name, creator);
	}


	std::unique_ptr<ContactDefinition> FemFactory::createContactDefinition(const char* name, ContactManager* cMan)
	{
		return cf_create<ContactDefinition>(contactDefList, name, cMan);
	}

	bool FemFactory::registerContactDefinition(const char* name, std::unique_ptr<ContactDefinition>(*creator)(ContactManager*))
	{
		return cf_store(contactDefList, name, creator);
	}

	std::unique_ptr<SparseGeneralEigenValueSystemNM> FemFactory::createGeneralizedEigenValueSolver(GenEigvalSolverType name, Domain* domain, EngngModel* emodel)
	{
		return cf_create4<SparseGeneralEigenValueSystemNM>(generalizedEigenValueSolverList, name, domain, emodel);
	}

	bool FemFactory::registerGeneralizedEigenValueSolver(GenEigvalSolverType name, std::unique_ptr<SparseGeneralEigenValueSystemNM>(*creator)(Domain*, EngngModel*))
	{
		return cf_store2(generalizedEigenValueSolverList, name, creator);
	}

	std::unique_ptr<IntegrationRule> FemFactory::createIRule(IntegrationRuleType type, int number, Element* e)
	{
		if (type == IRT_Gauss) {
			return std::make_unique<GaussIntegrationRule>(number, e);
		}
		else if (type == IRT_Lobatto) {
			//return std::make_unique<LobattoIntegrationRule>(number, e);
		}
		return nullptr;
	}

	std::unique_ptr<LoadBalancerMonitor> FemFactory::createLoadBalancerMonitor(const char* name, EngngModel* emodel)
	{
		return cf_create<LoadBalancerMonitor>(loadMonitorList, name, emodel);
	}

	bool FemFactory::registerLoadBalancerMonitor(const char* name, std::unique_ptr<LoadBalancerMonitor>(*creator)(EngngModel*))
	{
		return cf_store(loadMonitorList, name, creator);
	}

	std::unique_ptr<LoadBalancer> FemFactory::createLoadBalancer(const char* name, Domain* domain)
	{
		return cf_create<LoadBalancer>(loadBalancerList, name, domain);
	}

	bool FemFactory::registerLoadBalancer(const char* name, std::unique_ptr<LoadBalancer>(*creator)(Domain*))
	{
		return cf_store(loadBalancerList, name, creator);
	}

} // End namespace fem
