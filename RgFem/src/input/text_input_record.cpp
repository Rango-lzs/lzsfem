/*****************************************************************//**
 * \file   text_input_record.cpp
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#include "text_input_record.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "dictionary.h"
#include "range.h"
#include "scalarfunction.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <ostream>
#include <sstream>

namespace fem {
	FEMTXTInputRecord::FEMTXTInputRecord() : tokenizer(), record()
	{ }

	FEMTXTInputRecord::FEMTXTInputRecord(const FEMTXTInputRecord& src) : tokenizer(),
		record(src.record), lineNumber(src.lineNumber)
	{
		tokenizer.tokenizeLine(this->record);
		int ntok = tokenizer.giveNumberOfTokens();
		readFlag.resize(ntok);
		for (int i = 0; i < ntok; i++) {
			readFlag[i] = src.readFlag[i];
		}
	}

	FEMTXTInputRecord::FEMTXTInputRecord(int linenumber, std::string source) : tokenizer(),
		record(std::move(source)), lineNumber(linenumber)
	{
		tokenizer.tokenizeLine(this->record);
		int ntok = tokenizer.giveNumberOfTokens();
		readFlag.resize(ntok);
		for (int i = 0; i < ntok; i++) {
			readFlag[i] = false;
		}
	}

	FEMTXTInputRecord&
		FEMTXTInputRecord :: operator = (const FEMTXTInputRecord& src)
	{
		this->record = src.record;
		tokenizer.tokenizeLine(this->record);
		int ntok = tokenizer.giveNumberOfTokens();
		readFlag.resize(ntok);
		for (int i = 0; i < ntok; i++) {
			readFlag[i] = src.readFlag[i];
		}

		return *this;
	}

	void
		FEMTXTInputRecord::setRecordString(std::string newRec)
	{
		this->record = std::move(newRec);
		tokenizer.tokenizeLine(this->record);
		int ntok = tokenizer.giveNumberOfTokens();
		readFlag.resize(ntok);
		for (int i = 0; i < ntok; i++) {
			readFlag[i] = false;
		}
	}

	void
		FEMTXTInputRecord::giveRecordKeywordField(std::string& answer, int& value)
	{
		if (tokenizer.giveNumberOfTokens() > 0) {
			answer = std::string(tokenizer.giveToken(1));
			setReadFlag(1);
			auto ptr = scanInteger(tokenizer.giveToken(2), value);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, "RecordID", lineNumber);
			}
			setReadFlag(2);
		}
		else {
			throw BadFormatInputException(*this, "RecordID", lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveRecordKeywordField(std::string& answer)
	{
		if (tokenizer.giveNumberOfTokens() > 0) {
			answer = std::string(tokenizer.giveToken(1));
			setReadFlag(1);
		}
		else {
			throw BadFormatInputException(*this, "RecordID", lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(int& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			auto ptr = scanInteger(tokenizer.giveToken(indx + 1), answer);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			setReadFlag(indx + 1);
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(double& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			auto ptr = scanDouble(tokenizer.giveToken(indx + 1), answer);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			setReadFlag(indx + 1);
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(bool& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int val;
			auto ptr = scanInteger(tokenizer.giveToken(indx + 1), val);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			setReadFlag(indx + 1);
			answer = val != 0;
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(std::string& answer, InputFieldType id)
	{
		int indx = 0;
		if (id) {
			if ((indx = this->giveKeywordIndx(id)) == 0) {
				throw MissingKeywordInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			indx++;
		}
		else {
			indx = 1;
		}

		const char* _token = tokenizer.giveToken(indx);
		if (_token) {
			answer = std::string(tokenizer.giveToken(indx));
			setReadFlag(indx);
		}
		else {
			answer = "";
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(IntArray& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int size;
			setReadFlag(indx);
			auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			answer.resize(size);
			setReadFlag(indx);

			for (int i = 1; i <= size; i++) {
				int value;
				ptr = scanInteger(tokenizer.giveToken(indx + i), value);
				if (ptr == nullptr || *ptr != 0) {
					throw BadFormatInputException(*this, id, lineNumber);
				}

				answer.at(i) = value;
				setReadFlag(indx + i);
			}

		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(FloatArray& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int size;
			setReadFlag(indx);
			auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			answer.resize(size);
			setReadFlag(indx);

			for (int i = 1; i <= size; i++) {
				double value;
				auto ptr = scanDouble(tokenizer.giveToken(indx + i), value);
				if (ptr == nullptr || *ptr != 0) {
					throw BadFormatInputException(*this, id, lineNumber);
				}

				answer.at(i) = value;
				setReadFlag(indx + i);
			}

		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(FloatMatrix& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int nrows, ncols;
			setReadFlag(indx);

			auto ptr = scanInteger(tokenizer.giveToken(++indx), nrows);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			ptr = scanInteger(tokenizer.giveToken(++indx), ncols);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);

			if (readMatrix(tokenizer.giveToken(++indx), nrows, ncols, answer) == 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(std::vector< std::string >& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int size;
			setReadFlag(indx);
			auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}
			answer.reserve(size);
			setReadFlag(indx);
			for (int i = 1; i <= size; i++) {
				answer.push_back(tokenizer.giveToken(indx + i));
				setReadFlag(indx + i);
			}

		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(Dictionary& answer, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			setReadFlag(indx);
			int size;
			auto ptr = scanInteger(tokenizer.giveToken(++indx), size);
			if (ptr == nullptr || *ptr != 0) {
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);

			answer.clear();
			for (int i = 1; i <= size; i++) {
				char key = (tokenizer.giveToken(++indx))[0];
				double value;
				setReadFlag(indx);
				auto ptr = scanDouble(tokenizer.giveToken(++indx), value);
				if (ptr == nullptr || *ptr != 0) {
					throw BadFormatInputException(*this, id, lineNumber);
				}

				setReadFlag(indx);
				answer.add(key, value);
			}
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void
		FEMTXTInputRecord::giveField(std::list< Range >& list, InputFieldType id)
	{
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			int li, hi;
			setReadFlag(indx);
			const char* rec = tokenizer.giveToken(++indx);
			if (*rec != '{') {
				FEM_WARNING("missing left '{'");
				list.clear();
				throw BadFormatInputException(*this, id, lineNumber);
			}

			setReadFlag(indx);
			rec++;
			// read ranges
			while (readRange(&rec, li, hi)) {
				Range range(li, hi);
				list.push_back(range);
			}

			// skip whitespaces after last range
			while (isspace(*rec)) {
				rec++;
			}

			// test for enclosing bracket
			if (*rec != '}') {
				FEM_WARNING("missing end '}'");
				list.clear();
				throw BadFormatInputException(*this, id, lineNumber);
			}
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	void FEMTXTInputRecord::giveField(ScalarFunction& answer, InputFieldType id)
	{
		const char* rec;
		int indx = this->giveKeywordIndx(id);

		if (indx) {
			setReadFlag(indx);
			rec = tokenizer.giveToken(++indx);
			if (*rec == '@') {
				// reference to function
				int refVal;
				auto ptr = scanInteger(rec + 1, refVal);
				if (ptr == nullptr || *ptr != 0) {
					throw BadFormatInputException(*this, id, lineNumber);
				}
				setReadFlag(indx);
				answer.setReference(refVal);
			}
			else if (*rec == '$') {
				// simple expression
				std::string expr;

				expr = std::string(tokenizer.giveToken(indx));
				setReadFlag(indx);
				std::string _v = expr.substr(1, expr.size() - 2);

				answer.setSimpleExpression(_v); // get rid of enclosing '$'
			}
			else {
				double val;
				auto ptr = scanDouble(tokenizer.giveToken(indx), val);
				if (ptr == nullptr || *ptr != 0) {
					throw BadFormatInputException(*this, id, lineNumber);
				}

				setReadFlag(indx);
				answer.setValue(val);
			}
		}
		else {
			throw MissingKeywordInputException(*this, id, lineNumber);
		}
	}

	bool
		FEMTXTInputRecord::hasField(InputFieldType id)
	{
		//returns nonzero if id is present in source
		int indx = this->giveKeywordIndx(id);
		if (indx) {
			setReadFlag(indx);
		}

		return (indx > 0) ? true : false;
	}

	void
		FEMTXTInputRecord::printYourself()
	{
		printf("%s", this->record.c_str());
	}

	const char*
		FEMTXTInputRecord::scanInteger(const char* source, int& value)
	{
		//
		// reads integer value from source, returns pointer to char after this number
		//
		char* endptr;

		if (source == nullptr) {
			value = 0;
			return nullptr;
		}

		value = strtol(source, &endptr, 10);
		return endptr;
	}

	const char*
		FEMTXTInputRecord::scanDouble(const char* source, double& value)
	{
		//
		// reads double value from source, returns pointer to char after this number
		//
		char* endptr;

		if (source == nullptr) {
			value = 0;
			return nullptr;
		}

		value = strtod(source, &endptr);
		return endptr;
	}

	int
		FEMTXTInputRecord::giveKeywordIndx(const char* kwd)
	{
		int ntokens = tokenizer.giveNumberOfTokens();
		for (int i = 1; i <= ntokens; i++) {
			if (strcmp(kwd, tokenizer.giveToken(i)) == 0) {
				return i;
			}
		}

		return 0;
	}

	void
		FEMTXTInputRecord::finish(bool wrn)
	{
		if (!wrn) {
			return;
		}

		std::ostringstream buff;
		bool pf = true, wf = false;
		int ntokens = tokenizer.giveNumberOfTokens();
		for (int i = 0; i < ntokens; i++) {
			//fprintf (stderr, "[%s] ", tokenizer.giveToken(i+1));
			if (!readFlag[i]) {
				if (pf) {
					buff << "Unread token(s) detected in the following record\n\"";
					for (int j = 0; j < 40; j++) {
						if (this->record[j] == '\n' || this->record[j] == '\0') {
							break;
						}
						else {
							buff << this->record[j];
						}
					}
					if (this->record.size() > 41) {
						buff << "...";
					}
					buff << "\":\n";

					pf = false;
					wf = true;
				}

				buff << "[" << tokenizer.giveToken(i + 1) << "]";
			}
		}

		if (wf) {
			FEM_WARNING(buff.str().c_str());
		}
	}

	int FEMTXTInputRecord::readRange(const char** helpSource, int& li, int& hi)
	{
		char* endptr;
		// skip whitespaces
		while (isspace(**helpSource)) {
			(*helpSource)++;
		}

		// test if character is digit
		if (isdigit(**helpSource)) {
			// digit character - read one value range
			li = hi = strtol(*helpSource, &endptr, 10);
			*helpSource = endptr;
			return 1;
		}
		else if (**helpSource == '(') {
			// range left parenthesis found
			(*helpSource)++;
			// read lower index
			li = strtol(*helpSource, &endptr, 10);
			*helpSource = endptr;
			// test whitespaces
			if (**helpSource != ' ' && **helpSource != '\t') {
				FEM_WARNING("unexpected token while reading range value");
				return 0;
			}

			// read end index
			hi = strtol(*helpSource, &endptr, 10);
			*helpSource = endptr;
			// skip whitespaces
			while (isspace(**helpSource)) {
				(*helpSource)++;
			}

			// test for enclosing bracket
			if (**helpSource == ')') {
				(*helpSource)++;
				return 1;
			}
			else {
				FEM_WARNING("end ')' missing while parsing range value");
				return 0;
			}
		}

		return 0;
	}

	int
		FEMTXTInputRecord::readMatrix(const char* helpSource, int r, int c, FloatMatrix& ans)
	{
		const char* endptr = helpSource;

		if (helpSource == NULL) {
			ans.clear();
			return 0;
		}

		ans.resize(r, c);
		// skip whitespaces
		while (isspace(*endptr)) {
			(endptr)++;
		}

		if (*endptr == '{') {
			// range left parenthesis found
			(endptr)++;
			// read row by row separated by semicolon
			for (int i = 1; i <= r; i++) {
				for (int j = 1; j <= c; j++) {
					endptr = scanDouble(endptr, ans.at(i, j));
				}

				if (i < r) {
					// skip whitespaces
					while (isspace(*endptr)) {
						(endptr)++;
					}

					// test for row terminating semicolon
					if (*endptr == ';') {
						(endptr)++;
					}
					else {
						FEM_WARNING("missing row terminating semicolon");
						return 0;
					}
				}
			}

			// skip whitespaces
			while (isspace(*endptr)) {
				(endptr)++;
			}

			// test for enclosing bracket
			if (*endptr == '}') {
				return 1;
			}
			else {
				FEM_WARNING("end '}' missing while parsing matrix value");
				return 0;
			}
		}
		else {
			return 0;
		}

	}

} // end namespace fem
