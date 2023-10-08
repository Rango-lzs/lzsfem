
/*****************************************************************//**
 * \file   timestep.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef timestep_h
#define timestep_h

#include "contextioresulttype.h"
#include "contextmode.h"
#include "statecountertype.h"
#include "timediscretizationtype.h"
#include "inputrecord.h"

namespace fem 
{
	class EngngModel;
	class DataStream;

	typedef long StateCounterType;
	/**
	 * Class representing solution step. The timeStep instance may represent either
	 * time step, load increment, or load case depending on used Engineering model.
	 * See corresponding EngngModel reference for details. TimeStep maintains a reference to corresponding
	 * Engineering model class instance.
	 *
	 * The class hold target time, which may represent the end of time interval. In addition, there is
	 * intrinsic time, which is normally set the same as target time. Intrinsic (assembly) time is used
	 * especially in constitutive laws, where the material law is not meant to be evaluated at target time. Also,
	 * imposing boundary conditions occurs at intrinsic time by default. This reflects changes of static system and proper
	 * equation numbering during each timeStep.
	 *
	 * Some components (typically integration points real stresses or integration points nonlocal values)
	 * are computationally very demanding. Because in typical code, there are number of requests for same value
	 * during the computation process, it may be efficient to store these values and compute them only once.
	 * The principal problem is to recognize, when is necessary to re-compute these stored values to reflect
	 * newly reached state. This cannot be determined form solution step "time", because solution step may
	 * represent for example load increment, inside which typically many iterations are needed to reach
	 * convergence. For this purpose, a concept of solution state counters is provided.
	 * Whenever the solution state changes, the engineering model updates the solution state counter.
	 * The solution state counter is guaranteed to grow up smoothly (it newer decreases) during solution process.
	 * Other components of program (integration points) can then store their computationally expensive values
	 * but have to store also corresponding solution state counter value valid when these were computed.
	 * Then, easy check is done for finding differences between frozen solution state counter and their value with
	 * current solution state requested from solution step and recompute the values if necessary.
	 */
	class FEM_EXPORT TimeStep
	{
	protected:
		/// Engineering model reference.
		EngngModel* eModel;
		/// Current target time, which represents time at the end of a time step.
		double targetTime;
		/// Current intrinsic time, which may represents imposing time of boundary condition or time entering constitutive laws.
		double intrinsicTime;
		/// Current intrinsic time increment.
		double deltaT;
		/// Solution state counter.
		StateCounterType solutionStateCounter;
		/// Receiver's number.
		int number;
		/**
		 * Receiver's version, used for special applications; default set to 0.
		 * Typically, new version of same step is generated after adaptive restart, when
		 * the restarted step is equilibrated on new domain.
		 */
		int version;
		/**
		 * Receiver's substep (iteration) number.
		 */
		int subStepNumber;
		/// Corresponding meta step number.
		int mStepNumber;
		/// Time discretization.
		TimeDiscretizationType timeDiscretization;

	public:
		/**
		 * Creates a new solution step.
		 * @param n Solution step number.
		 * @param e Reference to corresponding engng model.
		 * @param mn Meta step number.
		 * @param tt Intrinsic time.
		 * @param dt Intrinsic time increment.
		 * @param counter Solution state counter.
		 * @param td Time discretization.
		 */
		TimeStep(int n, EngngModel* e, int mn, double tt, double dt, StateCounterType counter, TimeDiscretizationType td = TD_Unspecified);
		TimeStep(const TimeStep&);
		/// Convenience ctor for constructing the next timestep based on the previous one.
		TimeStep(const TimeStep& previous, double dt);
		TimeStep(EngngModel* e);
		TimeStep& operator = (const TimeStep&);

		/// Returns receiver's number.
		int giveNumber() { return number; }
		/// Set receiver's number.
		void setNumber(int i) { number = i; }
		/// Returns receiver's version.
		int giveVersion() { return version; }
		/// Returns receiver's meta step number.
		int giveMetaStepNumber() { return mStepNumber; }
		/// Returns receiver's substep number.
		int giveSubStepNumber() { return subStepNumber; }
		/**
		 * Returns class name of receiver.
		 * @return Pointer to s parameter filled with name.
		 */
		const char* giveClassName() const { return "TimeStep"; }
		/// Returns pointer to previous solution step.
		TimeStep* givePreviousStep();
		/// Returns target time.
		double giveTargetTime() { return targetTime; }
		/// Returns intrinsic time, e.g. time in which constitutive model is evaluated.
		double giveIntrinsicTime() { return intrinsicTime; }
		/// Returns solution step associated time increment.
		double giveTimeIncrement() { return deltaT; }
		/// Sets solution step time increment.
		void setTimeIncrement(double newDt) { deltaT = newDt; }
		/// Sets target and intrinsic time to be equal.
		void setTime(double newt) {
			targetTime = newt;
			intrinsicTime = newt;
		}
		/// Sets only target time.
		void setTargetTime(double newt) { targetTime = newt; }
		/// Sets only intrinsic time.
		void setIntrinsicTime(double newt) { intrinsicTime = newt; }
		/// Sets time discretization.
		void setTimeDiscretization(TimeDiscretizationType td) { timeDiscretization = td; }

		/**
		 * Check if solution step is not the last step.
		 * @return True if not last step, false otherwise.
		 */
		bool isNotTheLastStep();
		/**
		 * Check if receiver is first step.
		 * @return True if receiver is the first step, false otherwise.
		 */
		bool isTheFirstStep();
		/**
		 * Check if receiver is current solution step.
		 * @returns True if receiver is current step, false otherwise.
		 */
		bool isTheCurrentTimeStep();
		/**
		 * Check if receiver is solution step when initial conditions should apply.
		 * @return True if ic apply, false otherwise.
		 */
		bool isIcApply();
		/**
		 * Returns current solution state counter.
		 */
		StateCounterType giveSolutionStateCounter() { return solutionStateCounter; }
		/// Updates solution state counter.
		void incrementStateCounter() { solutionStateCounter++; }
		/// Increments receiver's version.
		void incrementVersion() { version++; }
		/// Increments receiver's substep number.
		void incrementSubStepNumber() { subStepNumber++; }
		/// Returns time discretization.
		TimeDiscretizationType giveTimeDiscretization() { return timeDiscretization; }

		void initializeFrom(InputRecord& ir) { }
		/**
		 * Stores receiver state to output stream.
		 * @param stream Output stream.
		 * @exception ContextIOERR If error encountered.
		 */
		void saveContext(DataStream& stream);
		/**
		 * Restores the receiver state previously written in stream.
		 * @see saveContext member function.
		 */
		void restoreContext(DataStream& stream);

		std::string errorInfo(const char* func) { return std::string("TimeStep::") + func; }
	};
} // end namespace fem
#endif // timestep_h
