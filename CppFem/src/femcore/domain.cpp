/*
*
*/

#include "domain.h"
#include "element.h"
#include "timestep.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "load.h"
#include "initialcondition.h"
#include "function.h"
#include "set.h"

#include "engngm.h"
#include "entityrenumberingscheme.h"
#include "datastream.h"
#include "contextioerr.h"
#include "verbose.h"
#include "connectivitytable.h"
#include "outputmanager.h"
#include "octreelocalizer.h"
#include "nodalrecoverymodel.h"
#include "nonlocalbarrier.h"
#include "classfactory.h"
#include "logger.h"

#include "topologydescription.h"
#include "errorestimator.h"
#include "range.h"
#include "fracturemanager.h"
#include "datareader.h"

#include "text_data_reader.h"
#include "initmodulemanager.h"
#include "exportmodulemanager.h"

#include "boundarycondition.h"
#include "activebc.h"
#include "simpleslavedof.h"
#include "masterdof.h"

#ifdef __PARALLEL_MODE
#include "parallel.h"
#include "processcomm.h"
#include "datastream.h"
#include "communicator.h"
#include "domaintransactionmanager.h"
#endif

#include <cstdarg>
#include <string>
#include <vector>
#include <set>

namespace fem 
{
	Domain::Domain(int n, int serNum, EngngModel* e)
	{
	}

	Domain :: ~Domain()
	{
	}

	int Domain::instanciateYourself(DataReader& dr)
		// Creates all objects mentioned in the data file.
	{
		int num;
		std::string name, topologytype;
		int nnode, nelem, nmat, nload, nic, nloadtimefunc, ncrossSections, nbarrier = 0, nset = 0;
		bool nxfemman = false;
		bool ncontactman = false;
		bool nfracman = false;

		// read type of Domain to be solved
		{
			auto& ir = dr.giveInputRecord   (DataReader::IR_domainRec, 1);
			IR_GIVE_FIELD(ir, name, _IFT_Domain_type); // This is inconsistent, "domain" isn't  exactly a field, but the actual record keyword.

			mDomainType = name;

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciating domain ", this->number);
#  endif

		resolveDomainDofsDefaults(name.c_str());

		// read output manager record
		{
			std::string tmp;
			auto& ir = dr.giveInputRecord(DataReader::IR_outManRec, 1);
			ir.giveRecordKeywordField(tmp);

			if (!giveEngngModel()->giveSuppressOutput()) {
				outputManager->initializeFrom(ir);
			}
			ir.finish();
		}

		// read domain description
		{
			auto& ir = dr.giveInputRecord(DataReader::IR_domainCompRec, 1);
			IR_GIVE_FIELD(ir, nnode, _IFT_Domain_ndofman);
			IR_GIVE_FIELD(ir, nelem, _IFT_Domain_nelem);
			IR_GIVE_FIELD(ir, ncrossSections, _IFT_Domain_ncrosssect);
			IR_GIVE_FIELD(ir, nmat, _IFT_Domain_nmat);
			IR_GIVE_FIELD(ir, nload, _IFT_Domain_nbc);
			IR_GIVE_FIELD(ir, nic, _IFT_Domain_nic);
			IR_GIVE_FIELD(ir, nloadtimefunc, _IFT_Domain_nfunct);
			IR_GIVE_OPTIONAL_FIELD(ir, nset, _IFT_Domain_nset);
			IR_GIVE_OPTIONAL_FIELD(ir, nxfemman, _IFT_Domain_nxfemman);
			IR_GIVE_OPTIONAL_FIELD(ir, ncontactman, _IFT_Domain_ncontactman);
			IR_GIVE_OPTIONAL_FIELD(ir, topologytype, _IFT_Domain_topology);
			this->nsd = -1; ///@todo Change this to default 0 when the domaintype record has been removed.
			IR_GIVE_OPTIONAL_FIELD(ir, this->nsd, _IFT_Domain_numberOfSpatialDimensions);
			this->axisymm = ir.hasField(_IFT_Domain_axisymmetric);
			IR_GIVE_OPTIONAL_FIELD(ir, nfracman, _IFT_Domain_nfracman);
			IR_GIVE_OPTIONAL_FIELD(ir, nbarrier, _IFT_Domain_nbarrier);
		}

		///@todo Eventually remove this backwards compatibility:
		//_HeatTransferMode _HeatMass1Mode // Are these deprecated?
		// set the number of spatial dimensions
		if (dType == _1dTrussMode) {
			nsd = 1;
		}
		else if (dType == _2dIncompressibleFlow || dType == _2dBeamMode || dType == _2dTrussMode || dType == _2dMindlinPlateMode || dType == _PlaneStrainMode || dType == _2dPlaneStressMode || dType == _2dPlaneStressRotMode || dType == _WarpingMode) {
			nsd = 2;
		}
		else if (dType == _3dIncompressibleFlow || dType == _3dShellMode || dType == _3dMode || dType == _3dDirShellMode) {
			nsd = 3;
		}
		else if (dType == _3dAxisymmMode) {
			nsd = 2;
			axisymm = true;
		}

		// read nodes
		dofManagerList.clear();
		dofManagerList.resize(nnode);
		for (int i = 1; i <= nnode; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_dofmanRec, i);
			// read type of dofManager
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			// assign component number according to record order
			// component number (as given in input record) becomes label
			std::unique_ptr< DofManager > dman(classFactory.createDofManager(name.c_str(), i, this));
			if (!dman) {
				FEM_ERROR("Couldn't create node of type: %s\n", name.c_str());
			}

			dman->initializeFrom(ir);
			dman->setGlobalNumber(num);    // set label
			dofManagerList[i - 1] = std::move(dman);

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated nodes & sides ", nnode)
#  endif

			BuildDofManPlaceInArrayMap();

		// read elements
		elementList.clear();
		elementList.resize(nelem);
		for (int i = 1; i <= nelem; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_elemRec, i);
			// read type of element
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< Element >elem(classFactory.createElement(name.c_str(), i, this));
			if (!elem) {
				FEM_ERROR("Couldn't create element: %s", name.c_str());
			}

			elem->initializeFrom(ir);
			elem->setGlobalNumber(num);
			elementList[i - 1] = std::move(elem);

			ir.finish();
		}

		BuildElementPlaceInArrayMap();

		// Support sets defined directly after the elements (special hack for backwards compatibility).
		setList.clear();
		if (dr.peakNext("set")) {
			setList.resize(nset);
			for (int i = 1; i <= nset; i++) {
				auto& ir = dr.giveInputRecord(DataReader::IR_setRec, i);
				// read type of set
				IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
				// Only one set for now (i don't see any need to ever introduce any other version)
				std::unique_ptr< Set > set = std::make_unique<Set>(num, this); //classFactory.createSet(name.c_str(), num, this)
				if (!set) {
					FEM_ERROR("Couldn't create set: %s", name.c_str());
				}

				set->initializeFrom(ir);

				// check number
				if (num < 1 || num > nset) {
					FEM_ERROR("Invalid set number (num=%d)", num);
				}

				if (!setList[num - 1]) {
					setList[num - 1] = std::move(set);
				}
				else {
					FEM_ERROR("Set entry already exist (num=%d)", num);
				}

				ir.finish();
			}
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated elements ", nelem);
#  endif

		// read cross sections
		crossSectionList.clear();
		crossSectionList.resize(ncrossSections);
		for (int i = 1; i <= ncrossSections; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_crosssectRec, i);
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< CrossSection >crossSection(classFactory.createCrossSection(name.c_str(), num, this));
			if (!crossSection) {
				FEM_ERROR("Couldn't create crosssection: %s", name.c_str());
			}

			crossSection->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > ncrossSections)) {
				FEM_ERROR("Invalid crossSection number (num=%d)", num);
			}

			if (!crossSectionList[num - 1]) {
				crossSectionList[num - 1] = std::move(crossSection);
			}
			else {
				FEM_ERROR("crossSection entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated cross sections ", ncrossSections)
#  endif

			// read materials
			materialList.clear();
		materialList.resize(nmat);
		for (int i = 1; i <= nmat; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_matRec, i);
			// read type of material
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< Material >mat(classFactory.createMaterial(name.c_str(), num, this));
			if (!mat) {
				FEM_ERROR("Couldn't create material: %s", name.c_str());
			}

			mat->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > nmat)) {
				FEM_ERROR("Invalid material number (num=%d)", num);
			}

			if (!materialList[num - 1]) {
				materialList[num - 1] = std::move(mat);
			}
			else {
				FEM_ERROR("material entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated materials ", nmat)
#  endif

			// read barriers
			nonlocalBarrierList.clear();
		nonlocalBarrierList.resize(nbarrier);
		for (int i = 1; i <= nbarrier; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_nlocBarRec, i);
			// read type of load
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< NonlocalBarrier >barrier(classFactory.createNonlocalBarrier(name.c_str(), num, this));
			if (!barrier) {
				FEM_ERROR("Couldn't create barrier: %s", name.c_str());
			}

			barrier->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > nbarrier)) {
				FEM_ERROR("Invalid barrier number (num=%d)", num);
			}

			if (!nonlocalBarrierList[num - 1]) {
				nonlocalBarrierList[num - 1] = std::move(barrier);
			}
			else {
				FEM_ERROR("barrier entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		if (nbarrier) {
			VERBOSE_PRINT0("Instanciated barriers ", nbarrier);
		}
#  endif

		// read boundary conditions
		bcList.clear();
		bcList.resize(nload);
		for (int i = 1; i <= nload; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_bcRec, i);
			// read type of bc
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< GeneralBoundaryCondition >bc(classFactory.createBoundaryCondition(name.c_str(), num, this));
			if (!bc) {
				FEM_ERROR("Couldn't create boundary condition: %s", name.c_str());
			}

			bc->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > nload)) {
				FEM_ERROR("Invalid boundary condition number (num=%d)", num);
			}

			if (!bcList[num - 1]) {
				bcList[num - 1] = std::move(bc);
			}
			else {
				FEM_ERROR("boundary condition entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated BCs ", nload)
#  endif

			// read initial conditions
			icList.clear();
		icList.resize(nic);
		for (int i = 1; i <= nic; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_icRec, i);
			// read type of load
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< InitialCondition >ic(new InitialCondition(num, this));
			if (!ic) {
				FEM_ERROR("Creation of IC no. %d failed", num);
			}

			ic->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > nic)) {
				FEM_ERROR("Invalid initial condition number (num=%d)", num);
			}

			if (!icList[num - 1]) {
				icList[num - 1] = std::move(ic);
			}
			else {
				FEM_ERROR("initial condition entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated ICs ", nic)
#  endif


			// read load time functions
			functionList.clear();
		functionList.resize(nloadtimefunc);
		for (int i = 1; i <= nloadtimefunc; i++) {
			auto& ir = dr.giveInputRecord(DataReader::IR_funcRec, i);
			// read type of func
			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

			std::unique_ptr< Function >func(classFactory.createFunction(name.c_str(), num, this));
			if (!func) {
				FEM_ERROR("Couldn't create time function: %s", name.c_str());
			}

			func->initializeFrom(ir);

			// check number
			if ((num < 1) || (num > nloadtimefunc)) {
				FEM_ERROR("Invalid Function number (num=%d)", num);
			}

			if (!functionList[num - 1]) {
				functionList[num - 1] = std::move(func);
			}
			else {
				FEM_ERROR("Function entry already exist (num=%d)", num);
			}

			ir.finish();
		}

#  ifdef VERBOSE
		VERBOSE_PRINT0("Instanciated load-time fncts ", nloadtimefunc)
#  endif

			// read sets
			if (setList.size() == 0) {
				setList.resize(nset);
				for (int i = 1; i <= nset; i++) {
					auto& ir = dr.giveInputRecord(DataReader::IR_setRec, i);
					// read type of set
					IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
					// Only one set for now (i don't see any need to ever introduce any other version)
					std::unique_ptr< Set > set = std::make_unique<Set>(num, this); //classFactory.createSet(name.c_str(), num, this)
					if (!set) {
						FEM_ERROR("Couldn't create set: %s", name.c_str());
					}

					set->initializeFrom(ir);

					// check number
					if ((num < 1) || (num > nset)) {
						FEM_ERROR("Invalid set number (num=%d)", num);
					}

					if (!setList[num - 1]) {
						setList[num - 1] = std::move(set);
					}
					else {
						FEM_ERROR("Set entry already exist (num=%d)", num);
					}

					ir.finish();
				}
			}

#  ifdef VERBOSE
		if (nset) {
			VERBOSE_PRINT0("Instanciated sets ", nset);
		}
#  endif

		if (nxfemman) {
			auto& ir = dr.giveInputRecord(DataReader::IR_xfemManRec, 1);

			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
			xfemManager = classFactory.createXfemManager(name.c_str(), this);
			if (!xfemManager) {
				FEM_ERROR("Couldn't create xfemmanager: %s", name.c_str());
			}

			xfemManager->initializeFrom(ir);
			xfemManager->instanciateYourself(dr);
#  ifdef VERBOSE
			VERBOSE_PRINT0("Instanciated xfem ", nxfemman);
#  endif
		}


		if (ncontactman) {
			// don't read any input yet
			auto& ir = dr.giveInputRecord(DataReader::IR_contactManRec, 1);

			IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
			contactManager = classFactory.createContactManager(name.c_str(), this);
			if (!contactManager) {
				FEM_ERROR("Couldn't create contact manager: %s", name.c_str());
			}

			contactManager->initializeFrom(ir);
			contactManager->instanciateYourself(dr);
		}
#  ifdef VERBOSE
		if (ncontactman) {
			VERBOSE_PRINT0("Instanciated contact manager ", ncontactman);
		}
#  endif


		this->topology = NULL;
		if (topologytype.length() > 0) {
			this->topology = classFactory.createTopology(topologytype.c_str(), this);
			if (!this->topology) {
				FEM_ERROR("Couldn't create topology of type '%s'", topologytype.c_str());
			}

			return this->topology->instanciateYourself(dr);
#  ifdef VERBOSE
			VERBOSE_PRINT0("Instanciated topologies ", topologytype.length());
#  endif
		}


		if (nfracman) {
			auto& ir = dr.giveInputRecord(DataReader::IR_fracManRec, 1);
			fracManager = std::make_unique<FractureManager>(this);
			fracManager->initializeFrom(ir);
			fracManager->instanciateYourself(dr);
#  ifdef VERBOSE
			VERBOSE_PRINT0("Instanciated fracture manager ", nxfemman);
#  endif
		}

		// change internal component references from labels to assigned local numbers
		MapBasedEntityRenumberingFunctor labelToLocNumFunctor(dofmanGlobal2LocalMap, elementGlobal2LocalMap);
		for (auto& dman : this->dofManagerList) {
			dman->updateLocalNumbering(labelToLocNumFunctor);
		}

		for (auto& element : this->elementList) {
			element->updateLocalNumbering(labelToLocNumFunctor);
		}

		for (auto& set : setList) {
			set->updateLocalNumbering(labelToLocNumFunctor);
		}


		BuildMaterialToElementMap();

		return 1;
	}


	void
		Domain::postInitialize()
	{
		// New  - in development /JB
		// set element cross sections based on element set definition and set the corresponding
		// material based on the cs
		for (int i = 1; i <= this->giveNumberOfCrossSectionModels(); i++) {
			if (int setNum = this->giveCrossSection(i)->giveSetNumber()) {
				Set* set = this->giveSet(setNum);
				for (int ielem : set->giveElementList()) {
					Element* element = this->giveElement(ielem);
					element->setCrossSection(i);
				}
			}
		}

		{
			spatialLocalizer = std::make_unique<OctreeSpatialLocalizer>(this);
			spatialLocalizer->init();
			connectivityTable = std::make_unique<ConnectivityTable>(this);
			FEM_LOG_INFO("Spatial localizer init done\n");
		}

		if (this->hasXfemManager()) {
			this->giveXfemManager()->postInitialize();
		}

		// Dofs must be created before dof managers due their post-initialization:
		this->createDofs();

		for (auto& dman : dofManagerList) {
			dman->postInitialize();
		}


		for (auto& el : elementList) {
			el->postInitialize();
		}

		for (auto& bc : bcList) {
			bc->postInitialize();
		}


	}


	std::string
		Domain::errorInfo(const char* func) const
	{
		return std::string("Domain::") + func + ", number: " + std::to_string(number);
	}


	const IntArray&
		Domain::giveDefaultNodeDofIDArry()
	{
		// returns default DofID array, defining physical meaning of particular DOFs
		// in Node Dof collection
		if (this->defaultNodeDofIDArry.giveSize()) {
			return defaultNodeDofIDArry;
		}

		if (dType == _2dPlaneStressRotMode) {
			defaultNodeDofIDArry = { D_u, D_v, R_w };
		}
		else if (dType == _2dPlaneStressMode) {
			defaultNodeDofIDArry = { D_u, D_v };
		}
		else if (dType == _PlaneStrainMode) {
			defaultNodeDofIDArry = { D_u, D_v };
		}
		else if (dType == _3dMode) {
			defaultNodeDofIDArry = { D_u, D_v, D_w };
		}
		else if (dType == _3dAxisymmMode) {
			defaultNodeDofIDArry = { D_u, D_v, R_w };
		}
		else if (dType == _2dMindlinPlateMode) {
			defaultNodeDofIDArry = { D_w, R_u, R_v };
		}
		else if (dType == _3dShellMode) {
			defaultNodeDofIDArry = { D_u, D_v, D_w, R_u, R_v, R_w };
		}
		else if (dType == _2dTrussMode) {
			defaultNodeDofIDArry = { D_u, D_w };
		}
		else if (dType == _1dTrussMode) {
			defaultNodeDofIDArry = { D_u };
		}
		else if (dType == _2dBeamMode) {
			defaultNodeDofIDArry = { D_u, D_w, R_v };
		}
		else if (dType == _2dLatticeMode) {
			defaultNodeDofIDArry = { D_u, D_v, R_w };
		}
		else if (dType == _HeatTransferMode) {
			defaultNodeDofIDArry = { T_f };
		}
		else if (dType == _Mass1TransferMode) {
			defaultNodeDofIDArry = { C_1 };
		}
		else if (dType == _HeatMass1Mode) {
			defaultNodeDofIDArry = { T_f, C_1 };
		}
		else if (dType == _2dIncompressibleFlow) {
			defaultNodeDofIDArry = { V_u, V_v, P_f };
		}
		else if (dType == _3dIncompressibleFlow) {
			defaultNodeDofIDArry = { V_u, V_v, V_w, P_f };
		}
		else if (dType == _3dDirShellMode) {
			defaultNodeDofIDArry = { D_u, D_v, D_w, W_u, W_v, W_w, Gamma };
		}
		else if (dType == _2dLatticeMassTransportMode) {
			defaultNodeDofIDArry = { P_f };
		}
		else if (dType == _3dLatticeMassTransportMode) {
			defaultNodeDofIDArry = { P_f };
		}
		else if (dType == _3dLatticeMode) {
			defaultNodeDofIDArry = { D_u, D_v, D_w, R_u, R_v, R_w };
		}
		else if (dType == _2dLatticeHeatTransferMode) {
			defaultNodeDofIDArry = { T_f };
		}
		else if (dType == _3dLatticeHeatTransferMode) {
			defaultNodeDofIDArry = { T_f };
		}
		else if (dType == _WarpingMode) {
			defaultNodeDofIDArry = { D_w };
		}
		else {
			FEM_ERROR("unknown domainType (%s)", __domainTypeToString(dType));
		}

		return defaultNodeDofIDArry;
	}


	int
		Domain::giveNumberOfSpatialDimensions()
	{
		return nsd;
	}


	bool
		Domain::isAxisymmetric()
	{
		return axisymm;
	}


	void
		Domain::resolveDomainDofsDefaults(const char* typeName)
		//
		// resolves default number of dofs per node according to domain type name.
		// and also resolves default dof mask according to domain type.
		//
	{
		if (!strncmp(typeName, "2dplanestressrot", 16)) {
			dType = _2dPlaneStressRotMode;
		}
		else if (!strncmp(typeName, "2dplanestress", 12)) {
			dType = _2dPlaneStressMode;
		}
		else if (!strncmp(typeName, "planestrain", 11)) {
			dType = _PlaneStrainMode;
		}
		else if (!strncmp(typeName, "3daxisymm", 9)) {
			dType = _3dAxisymmMode;
		}
		else if (!strncmp(typeName, "2dmindlinplate", 14)) {
			dType = _2dMindlinPlateMode;
		}
		else if (!strncmp(typeName, "3dshell", 7)) {
			dType = _3dShellMode;
		}
		else if (!strncmp(typeName, "2dtruss", 7)) {
			dType = _2dTrussMode;
		}
		else if (!strncmp(typeName, "1dtruss", 7)) {
			dType = _1dTrussMode;
		}
		else if (!strncmp(typeName, "2dbeam", 6)) {
			dType = _2dBeamMode;
		}
		else if (!strncmp(typeName, "2dlattice", 9)) {
			dType = _2dLatticeMode;
		}
		else if (!strncmp(typeName, "heattransfer", 12)) {
			dType = _HeatTransferMode;
		}
		else if (!strncmp(typeName, "mass1transfer", 13)) {
			dType = _Mass1TransferMode;
		}
		else if (!strncmp(typeName, "hema1", 5)) {
			dType = _HeatMass1Mode;
		}
		else if (!strncmp(typeName, "2dincompflow", 12)) {
			dType = _2dIncompressibleFlow;
		}
		else if (!strncmp(typeName, "3dincompflow", 12)) {
			dType = _3dIncompressibleFlow;
		}
		else if (!strncmp(typeName, "3ddirshell", 10)) {
			dType = _3dDirShellMode;
		}
		else if (!strncmp(typeName, "2dmasslatticetransport", 22)) {
			dType = _2dLatticeMassTransportMode;
		}
		else if (!strncmp(typeName, "3dlattice", 9)) {
			dType = _3dLatticeMode;
		}
		else if (!strncmp(typeName, "3dmasslatticetransport", 22)) {
			dType = _3dLatticeMassTransportMode;
		}
		else if (!strncmp(typeName, "2dheatlattice", 13)) {
			dType = _3dLatticeMassTransportMode;
		}
		else if (!strncmp(typeName, "3dheatlattice", 13)) {
			dType = _3dLatticeMassTransportMode;
		}
		else if (!strncmp(typeName, "3d", 2)) {
			dType = _3dMode;
		}
		else if (!strncmp(typeName, "warping", 7)) {
			dType = _WarpingMode;
		}
		else {
			FEM_ERROR("unknown domainType (%s)", typeName);
			return;
		}
	}


	NodalRecoveryModel*
		Domain::giveSmoother()
	{
		return this->smoother.get();
	}


	void
		Domain::setSmoother(NodalRecoveryModel* newSmoother, bool destroyOld)
	{
		if (!destroyOld) {
			this->smoother.release();
		}

		this->smoother.reset(newSmoother);
	}


	void
		Domain::setTopology(TopologyDescription* topo, bool destroyOld)
	{
		if (!destroyOld) {
			this->topology.release();
		}

		this->topology.reset(topo);
	}


	ConnectivityTable*
		Domain::giveConnectivityTable()
		//
		// return connectivity Table - if no defined - creates new one
		//
	{
		if (!connectivityTable) {
			//connectivityTable = std::make_unique<ConnectivityTable>(this);
			FEM_LOG_ERROR("Connectivity table init error");
		}

		return connectivityTable.get();
	}


	SpatialLocalizer*
		Domain::giveSpatialLocalizer()
		//
		// return connectivity Table - if no defined - creates new one
		//
	{
		//     if (spatialLocalizer == NULL) spatialLocalizer = new DummySpatialLocalizer(1, this);
		if (spatialLocalizer) {
			return spatialLocalizer.get();
		}
		else {
			FEM_LOG_ERROR("Spatial localizer init failure");
			return nullptr;
		}
	}

	void Domain::setSpatialLocalizer(std::unique_ptr<SpatialLocalizer> sl) {
		spatialLocalizer = std::move(sl);
	}


	void
		Domain::createDofs()
	{

		///////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// Step 1. Scan all required nodal dofs.
		std::vector< std::set< int > > node_dofs(this->giveNumberOfDofManagers());
		for (auto& element : this->elementList) {
			IntArray dofids;
			// Scan for all dofs needed by element.
			for (int j = 1; j <= element->giveNumberOfNodes(); ++j) {
				element->giveDofManDofIDMask(j, dofids);
				for (int k = 1; k <= dofids.giveSize(); k++) {
					node_dofs[element->giveNode(j)->giveNumber() - 1].insert(dofids.at(k));
				}
			}
		}
		for (auto& dman : this->dofManagerList) {
			// Nodes can also contain their own list of dofs (typical usecase: RigidArmNode )
			const IntArray* dofids = dman->giveForcedDofIDs();
			if (dofids) {
				for (int k = 1; k <= dofids->giveSize(); ++k) {
					node_dofs[dman->giveNumber() - 1].insert(dofids->at(k));
				}
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		// Step 2. Scan all Dirichlet b.c.s (or active dofs). For every node we store a map from the dofid to it's b.c. number.
		// This loop won't check for slave dofs or so, and will give a bc id for every single relevant dof.
		// This must be a separate step since we store the inverse mapping (bc->dof instead of dof->bc) so we want to loop over all b.c.s to invert this.
		std::vector< std::map< int, int > > dof_bc(this->giveNumberOfDofManagers());
		for (int i = 1; i <= this->giveNumberOfBoundaryConditions(); ++i) {
			GeneralBoundaryCondition* gbc = this->giveBc(i);
			if (gbc->giveSetNumber() > 0) { ///@todo This will eventually not be optional.
				// Loop over nodes in set and store the bc number in each dof.
				Set* set = this->giveSet(gbc->giveSetNumber());
				ActiveBoundaryCondition* active_bc = dynamic_cast<ActiveBoundaryCondition*>(gbc);
				BoundaryCondition* bc = dynamic_cast<BoundaryCondition*>(gbc);
				if (bc || (active_bc && active_bc->requiresActiveDofs())) {
					const IntArray& appliedDofs = gbc->giveDofIDs();
					const IntArray& nodes = set->giveNodeList();
					for (int inode = 1; inode <= nodes.giveSize(); ++inode) {
						for (int idof = 1; idof <= appliedDofs.giveSize(); ++idof) {
							dof_bc[nodes.at(inode) - 1][appliedDofs.at(idof)] = i;
						}
					}
				}
			}
		}
		// Step 2b. This step asks nodes for their bc-vector, which is the old approach to dirichlet b.c.s (i.e. this is for backwards compatibility)
		///@todo Remove this input method whenever we decide on deprecating the old approach.
		for (int i = 1; i <= this->giveNumberOfDofManagers(); ++i) {
			DofManager* dman = this->giveDofManager(i);
			const std::map< int, int >* dmanBcs = dman->giveBcMap();
			if (dmanBcs) {
				dof_bc[i - 1].insert(dmanBcs->begin(), dmanBcs->end());     // This will ignore duplicated dofiditems.
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		// Step 3. Same for initial conditions as for boundary conditions in step 2.
		std::vector< std::map< int, int > > dof_ic(this->giveNumberOfDofManagers());
		for (int i = 1; i <= this->giveNumberOfInitialConditions(); ++i) {
			InitialCondition* ic = this->giveIc(i);
			if (ic->giveSetNumber() > 0) { ///@todo This will eventually not be optional.
				// Loop over nodes in set and store the bc number in each dof.
				Set* set = this->giveSet(ic->giveSetNumber());
				const IntArray& appliedDofs = ic->giveDofIDs();
				const IntArray& nodes = set->giveNodeList();
				for (int inode = 1; inode <= nodes.giveSize(); ++inode) {
					for (int idof = 1; idof <= appliedDofs.giveSize(); ++idof) {
						dof_ic[nodes.at(inode) - 1][appliedDofs.at(idof)] = i;
					}
				}
			}
		}
		// Step 3b. This step asks nodes for their bc-vector, which is the old approach to dirichlet b.c.s (i.e. this is for backwards compatibility)
		///@todo Remove this input method whenever we decide on deprecating the old approach.
		for (int i = 1; i <= this->giveNumberOfDofManagers(); ++i) {
			DofManager* dman = this->giveDofManager(i);
			const std::map< int, int >* dmanIcs = dman->giveIcMap();
			if (dmanIcs) {
				dof_ic[i - 1].insert(dmanIcs->begin(), dmanIcs->end());     // This will ignore duplicated dofiditems.
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////
		// Step 3. Create the dofs. This involves obtaining the correct
		for (int i = 1; i <= this->giveNumberOfDofManagers(); ++i) {
			DofManager* dman = this->giveDofManager(i);
			//printf("Dofs in node %d (of %d) = %d\n", i, this->giveNumberOfDofManagers(), node_dofs[i-1].size());

			/* do not delete existing DOFs; that may be created during adaptive solution scheme (mesh generator applies DOFs) */
			if (0) dman->setNumberOfDofs(0);

			for (int id : node_dofs[i - 1]) {
				// Find bc and ic if there are any, otherwise zero.
				int bcid = dof_bc[i - 1].find(id) != dof_bc[i - 1].end() ? dof_bc[i - 1][id] : 0;
				int icid = dof_ic[i - 1].find(id) != dof_ic[i - 1].end() ? dof_ic[i - 1][id] : 0;

				// Determine the doftype:
				dofType dtype = DT_master;
				const std::map< int, int >* dmanTypes = dman->giveDofTypeMap();
				if (dmanTypes) {
					std::map< int, int > ::const_iterator it = dmanTypes->find(id);
					if (it != dmanTypes->end()) {
						dtype = (dofType)it->second;
					}
				}
				// Check if active dofs are needed:
				if (bcid > 0) {
					// What should take precedence here if there is a slave node?
					// Right now the active b.c. overrides anything set prior, if necessary.
					// This seems like the most suitable choice, but it could possibly be changed.
					ActiveBoundaryCondition* active_bc = dynamic_cast<ActiveBoundaryCondition*>(this->giveBc(bcid));
					if (active_bc && active_bc->requiresActiveDofs()) {
						dtype = DT_active;
					}
				}

				if (!dman->isDofTypeCompatible(dtype)) {
					FEM_ERROR("Incompatible dof type (%d) in node %d", dtype, i);
				}

				// Finally create the new DOF: 
				//printf("Creating: node %d, id = %d, dofType = %d, bc = %d, ic = %d\n", i, id, dtype, bcid, icid);
				if (!dman->hasDofID((DofIDItem)id)) {

					Dof* dof = classFactory.createDof(dtype, (DofIDItem)id, dman);
					dof->setBcId(bcid); // Note: slave dofs and such will simple ignore this.
					dof->setIcId(icid);
					// Slave dofs obtain their weights post-initialization, simple slave dofs must have their master node specified.
					if (dtype == DT_simpleSlave) {
						static_cast<SimpleSlaveDof*>(dof)->setMasterDofManagerNum((*dman->giveMasterMap())[id]);
					}
					dman->appendDof(dof);
				}
			}
		}

		// XFEM manager create additional dofs themselves:
		if (this->hasXfemManager()) {
			xfemManager->createEnrichedDofs();
		}

		if (this->hasContactManager()) {
			contactManager->createContactDofs();
		}
	}


	int
		Domain::checkConsistency()
		// this function transverse tree of all objects and invokes
		// checkConsistency on this objects
		// currently this function checks noly consistency
		// of internal object structures, mainly whether referenced other objects
		// are having required support
		//
	{
		int result = 1;

		for (auto& dman : this->dofManagerList) {
			result &= dman->checkConsistency();
		}

		for (auto& element : this->elementList) {
			result &= element->checkConsistency();
		}

		for (auto& material : this->materialList) {
			result &= material->checkConsistency();
		}

		return result;
	}

	double
		Domain::giveArea()
	{
		double area = 0.0;
		for (auto& element : this->elementList) {
			area += element->computeArea();
		}

		return area;
	}

	double
		Domain::giveVolume()
	{
		double volume = 0.0;
		for (auto& element : this->elementList) {
			volume += element->computeVolume();
		}

		return volume;
	}

	double
		Domain::giveSize()
	{
		double volume = 0.0;
		for (auto& element : this->elementList) {
			volume += element->computeVolumeAreaOrLength();
		}

		return volume;
	}

	int
		Domain::giveNextFreeDofID(int increment)
	{
		if (this->engineeringModel->isParallel()) {
			FEM_ERROR("Additional dof id's not implemented/tested for parallel problems");
		}

		int freeID = this->freeDofID;
		this->freeDofID += increment;
		return freeID;
	}

	void
		Domain::resetFreeDofID()
	{
		this->freeDofID = MaxDofID;
	}

	ErrorEstimator*
		Domain::giveErrorEstimator()
	{
		return engineeringModel->giveDomainErrorEstimator(this->number);
	}


	OutputManager*
		Domain::giveOutputManager()
	{
		return outputManager.get();
	}


	TopologyDescription*
		Domain::giveTopology()
	{
		return topology.get();
	}


	template< typename T >
	void save_components(T& list, DataStream& stream, ContextMode mode)
	{
		if (!stream.write((int)list.size())) {
			THROW_CIOERR(CIO_IOERR);
		}
		for (const auto& object : list) {
			if ((mode & CM_Definition) != 0) {
				if (stream.write(std::string(object->giveInputRecordName())) == 0) {
					THROW_CIOERR(CIO_IOERR);
				}
			}
			object->saveContext(stream, mode);
		}
	}


	template< typename T, typename C >
	void restore_components(T& list, DataStream& stream, ContextMode mode, const C& creator)
	{
		int size = 0;
		if (!stream.read(size)) {
			THROW_CIOERR(CIO_IOERR);
		}
		if (mode & CM_Definition) {
			list.clear();
			list.resize(size);
		}
		for (int i = 1; i <= size; i++) {
			if (mode & CM_Definition) {
				std::string name;
				if (!stream.read(name)) {
					THROW_CIOERR(CIO_IOERR);
				}
				list[i - 1] = creator(name, i);
			}
			list[i - 1]->restoreContext(stream, mode);
		}
	}


	void
		Domain::saveContext(DataStream& stream, ContextMode mode)
	{
		if (!stream.write(this->giveSerialNumber())) {
			THROW_CIOERR(CIO_IOERR);
		}

		if ((mode & CM_Definition)) {
			save_components(this->setList, stream, mode);
			save_components(this->materialList, stream, mode);
			save_components(this->crossSectionList, stream, mode);
			save_components(this->icList, stream, mode);
			save_components(this->functionList, stream, mode);
			save_components(this->nonlocalBarrierList, stream, mode);
		}

		save_components(this->dofManagerList, stream, mode);
		save_components(this->elementList, stream, mode);
		save_components(this->bcList, stream, mode);

		auto ee = this->giveErrorEstimator();
		if (ee) {
			ee->saveContext(stream, mode);
		}
	}


	void
		Domain::restoreContext(DataStream& stream, ContextMode mode)
	{
		bool domainUpdated = false;
		int serNum = this->giveSerialNumber();
		// restore domain serial number
		if (!stream.read(this->serialNumber)) {
			THROW_CIOERR(CIO_IOERR);
		}

		if ((mode & CM_Definition)) {
			// clear cached data:
			elementGlobal2LocalMap.clear();
			dofmanGlobal2LocalMap.clear();

			///@todo Saving and restoring xfemmanagers.
			xfemManager = nullptr;

			restore_components(this->setList, stream, mode,
				[this](std::string& x, int i) { return std::make_unique<Set>(i, this); });
			restore_components(this->materialList, stream, mode,
				[this](std::string& x, int i) { return classFactory.createMaterial(x.c_str(), i, this); });
			restore_components(this->crossSectionList, stream, mode,
				[this](std::string& x, int i) { return classFactory.createCrossSection(x.c_str(), i, this); });
			restore_components(this->icList, stream, mode,
				[this](std::string& x, int i) { return classFactory.createInitialCondition(x.c_str(), i, this); });
			restore_components(this->functionList, stream, mode,
				[this](std::string& x, int i) { return classFactory.createFunction(x.c_str(), i, this); });
			restore_components(this->nonlocalBarrierList, stream, mode,
				[this](std::string& x, int i) { return classFactory.createNonlocalBarrier(x.c_str(), i, this); });

			domainUpdated = true;
		}
		else {
			if (serNum != this->giveSerialNumber()) {
				// read corresponding domain
				FEM_LOG_INFO("restoring domain %d.%d\n", this->number, this->giveSerialNumber());
				OOFEMTXTDataReader domainDr(this->engineeringModel->giveDomainFileName(1, this->giveSerialNumber()));
				this->clear();

				if (!this->instanciateYourself(domainDr)) {
					FEM_ERROR("domain Instanciation failed");
				}

				domainUpdated = true;
			}
		}

		restore_components(this->dofManagerList, stream, mode,
			[this](std::string& x, int i) { return classFactory.createDofManager(x.c_str(), i, this); });
		restore_components(this->elementList, stream, mode,
			[this](std::string& x, int i) { return classFactory.createElement(x.c_str(), i, this); });
		restore_components(this->bcList, stream, mode,
			[this](std::string& x, int i) { return classFactory.createBoundaryCondition(x.c_str(), i, this); });

		auto ee = this->giveErrorEstimator();
		if (ee) {
			if (domainUpdated) {
				ee->setDomain(this);
			}
			ee->restoreContext(stream, mode);
		}

		if (domainUpdated) {
			if (this->smoother) {
				this->smoother->clear();
			}
		}
	}

	void Domain::BuildElementPlaceInArrayMap()
	{
		elementGlobal2LocalMap.clear();

		int nelem = giveNumberOfElements();

		for (int i = 1; i <= nelem; i++) {
			Element* elem = this->giveElement(i);
			int globnum = elem->giveGlobalNumber();
			if (elementGlobal2LocalMap.find(globnum) == elementGlobal2LocalMap.end()) {
				// label does not exist yet
				elementGlobal2LocalMap[globnum] = i;
			}
			else {
				FEM_ERROR("Element with same global number already exist (label=%d)", globnum);
			}
		}
	}

	void Domain::BuildDofManPlaceInArrayMap()
	{
		dofmanGlobal2LocalMap.clear();

		int ndman = giveNumberOfDofManagers();

		for (int i = 1; i <= ndman; i++) {
			DofManager* dMan = this->giveDofManager(i);
			int globnum = dMan->giveGlobalNumber();
			if (dofmanGlobal2LocalMap.find(globnum) == dofmanGlobal2LocalMap.end()) {
				// label does not exist yet
				dofmanGlobal2LocalMap[globnum] = i;
			}
			else {
				FEM_ERROR("DofManagaer with same global number already exist (label=%d)", globnum);
			}
		}
	}

	void Domain::BuildMaterialToElementMap()
	{
		materialNum2ElMap.clear();

		int nelem = giveNumberOfElements();

		for (int i = 1; i <= nelem; i++) {
			Element* elem = this->giveElement(i);
			int matNum = elem->giveMaterialNumber();
			materialNum2ElMap[matNum].followedBy(i);
		}
	}

} // end namespace fem
