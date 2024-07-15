#ifndef modulemanager_h
#define modulemanager_h

#include "datareader.h"
#include "error.h"
#include "exportmodule.h"


#include <memory>
#include <vector>

///@name Input fields for module managers
//@{
#define _IFT_ModuleManager_nmodules "nmodules"
//@}

namespace fem
{
	class EngngModel;

	/**
	 * Class representing and implementing ModuleManager. It is attribute of EngngModel.
	 * It manages the modules of given type.
	 */
	template< class M >
	class FEM_EXPORT ModuleManager
	{
	protected:
		/// Module list.
		std::vector< std::unique_ptr< M > > moduleList;
		/// Number of modules.
		int numberOfModules;
		/// Associated Engineering model.
		EngngModel* emodel;

	public:
		ModuleManager(EngngModel* emodel) {
			this->emodel = emodel;
			numberOfModules = 0;
		}
		virtual ~ModuleManager() {
		}
		ModuleManager(const ModuleManager&) = delete;
		ModuleManager& operator=(const ModuleManager&) = delete;
		/**
		 * Creates new instance of module.
		 * @param name Name of module.
		 * @param n Number associated with module.
		 * @param emodel Engineering model which receiver belongs to.
		 */
		virtual std::unique_ptr<M> CreateModule(const char* name, int n, EngngModel* emodel) = 0;
		/**
		 *
		 * Reads receiver description from input stream and creates corresponding modules components accordingly.
		 * It scans input file, each line is assumed to be single record describing particular module.
		 * The record line is converted to lowercase letters.
		 * After new output module object is created, its initializeForm member function is
		 * called with its record as parameter.
		 * @param dr Data reader for input records.
		 * @param ir Record for receiver.
		 * @return Nonzero if o.k.
		 */
		virtual int instanciateYourself(DataReader& dr, InputRecord& ir)
		{
			std::string name;

			// read modules
			moduleList.reserve(numberOfModules);
			for (int i = 0; i < numberOfModules; i++) {
				auto& mir = dr.giveInputRecord(DataReader::IR_expModuleRec, i + 1);
				mir.giveRecordKeywordField(name);

				// read type of module
				std::unique_ptr< M > module = this->CreateModule(name.c_str(), i, emodel);
				if (!module) {
					FEM_ERROR("unknown module (%s)", name.c_str());
				}

				module->initializeFrom(mir);
				registerModule(module);
				//             moduleList.push_back(std :: move(module));
			}

#  ifdef VERBOSE
			VERBOSE_PRINT0("Instanciated output modules ", numberOfModules)
#  endif
				return 1;
		}

		/**
		 * Stores a module in moduleList. Useful when
		 * adding modules externally, e.g. from Python
		 */

		virtual void registerModule(std::unique_ptr< M >& module) {
			moduleList.push_back(std::move(module));
		}


		/**
		 * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
		 * from corresponding record. Should be called before instanciateYourself.
		 */
		virtual void initializeFrom(InputRecord& ir) = 0;
		/// Returns class name of the receiver.
		virtual const char* giveClassName() const = 0;

		/**
		 * Returns the required module.
		 * @param num Module number.
		 */
		M* giveModule(int num) {
			M* elem = NULL;

			if (num >= 1 && num <= (int)moduleList.size()) {
				elem = moduleList[num - 1].get();
			}
			else {
				FEM_ERROR("No module no. %d defined", num);
			}

			return elem;
		}

		int giveNumberOfModules() const { return (int)moduleList.size(); }
	};

	template class FEM_EXPORT ModuleManager<ExportModule>;

} // end namespace fem
#endif // modulemanager_h
