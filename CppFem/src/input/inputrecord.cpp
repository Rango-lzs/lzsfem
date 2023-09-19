/*****************************************************************//**
 * \file   inputrecord.cpp
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#include "inputrecord.h"

namespace fem {

	void InputRecord::giveOptionalField(int& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(double& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(bool& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(std::string& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(FloatArray& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(IntArray& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(FloatMatrix& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(std::vector< std::string >& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(Dictionary& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(std::list< Range >& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}

	void InputRecord::giveOptionalField(ScalarFunction& answer, InputFieldType id)
	{
		if (this->hasField(id)) {
			try {
				this->giveField(answer, id);
			}
			catch (MissingKeywordInputException&) {}
		}
	}


	InputException::InputException(const InputRecord& ir, std::string keyword, int number) :
		record(ir.giveRecordAsString()), keyword(std::move(keyword)), number(number)
	{ }


	MissingKeywordInputException::MissingKeywordInputException(const InputRecord& ir, std::string kw, int n) :
		InputException(ir, std::move(kw), n)
	{
		msg = "Missing keyword \"" + keyword + "\" on input " + std::to_string(number) + \
			"\nRecord: \"" + record + "\"";
	}


	BadFormatInputException::BadFormatInputException(const InputRecord& ir, std::string kw, int n) :
		InputException(ir, std::move(kw), n)
	{
		msg = "Bad format for keyword \"" + keyword + "\" on input " + std::to_string(number) + \
			"\nRecord: \"" + record + "\"";
	}


	ValueInputException::ValueInputException(const InputRecord& ir, std::string kw, const std::string& reason) :
		InputException(ir, std::move(kw), -1)
	{
		msg = "Value input error for keyword \"" + keyword + "\"" + \
			"\nReason: \"" + reason + "\"\n" + \
			"\nRecord: \"" + record + "\"";
	}


	const char* MissingKeywordInputException::what() const noexcept
	{
		return msg.c_str();
	}


	const char* BadFormatInputException::what() const noexcept
	{
		return msg.c_str();
	}

	const char* ValueInputException::what() const noexcept
	{
		return msg.c_str();
	}

} // end namespace fem
