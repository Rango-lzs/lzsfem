
#include "femcmpnn.h"
#include "contextioerr.h"
#include "dynamic_input_record.h"

#include <cstdarg>


namespace fem
{
	void FEMComponent::giveInputRecord(DynamicInputRecord& input)
	{
		input.setRecordKeywordField(this->giveInputRecordName(), this->giveNumber());
	}


	std::string
		FEMComponent::errorInfo(const char* func) const
	{
		return std::string(this->giveClassName()) + "::" + func + ", number: " + std::to_string(this->giveNumber());
	}

	void FEMComponent::initializeFrom(InputRecord& ir)
	{
	}

	int FEMComponent::checkConsistency()
	{
		return 1;
	}

} // end namespace oofem
