

#ifndef exportmodulemanager_h
#define exportmodulemanager_h

#include "modulemanager.h"
#include "exportmodule.h"

namespace fem
{
	class EngngModel;

	/**
	 * Class representing and implementing ExportModuleManager. It is attribute of EngngModel.
	 * It manages the export output modules, which perform module - specific output operations.
	 */
	class FEM_EXPORT ExportModuleManager : public ModuleManager< ExportModule >
	{
	public:
		ExportModuleManager(EngngModel* emodel);
		virtual ~ExportModuleManager();

		void initializeFrom(InputRecord& ir) override;
		std::unique_ptr<ExportModule> CreateModule(const char* name, int num, EngngModel* emodel) override;

		/**
		 * Writes the output. Loops over all modules and calls corresponding doOutput module service.
		 * @param tStep Time step.
		 * @param substepFlag is set to true, only the modules with substepFlag set to true will be processed.
		 */
		void doOutput(TimeStep* tStep, bool substepFlag = false);
		/**
		 * Initializes output manager. The corresponding initialize module services are called.
		 */
		void initialize();
		/**
		 * Terminates the receiver, the corresponding terminate module services are called.
		 */
		void terminate();
		const char* giveClassName() const override { return "ExportModuleManager"; }
	};
} // end namespace oofem
#endif // exportmodulemanager_h
