/*****************************************************************//**
 * \file   dynamic_data_reader.cpp
 * \brief  
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#include "dynamic_data_reader.h"
#include "inputrecord.h"
#include "error.h"

#include <memory>
#include <fstream>

namespace fem {
DynamicDataReader::DynamicDataReader(std :: string name) : DataReader(), name(std :: move(name))
{
    this->it = recordList.end();
}

DynamicDataReader :: ~DynamicDataReader()
{
}

void
DynamicDataReader :: insertInputRecord(InputRecordType type, std::unique_ptr<InputRecord> record)
{
    // Should care about the record type, but its a hassle.
    this->recordList.push_back(std::move(record));
    this->it = recordList.end();
}

InputRecord &
DynamicDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    // Ignores recordId in favor of having a dynamic list (just incremental access). typeId could be supported, but its a hassle.
    // The txt data reader makes the same assumptions.
    if ( this->it == this->recordList.end() ) {
        this->it = this->recordList.begin();
    } else {
        ++this->it;
    }
    return **(this->it);
}

bool DynamicDataReader :: peakNext(const std :: string &keyword)
{
    std :: string nextKey;
    auto temp = this->it;
    temp++;
    (*temp)->giveRecordKeywordField(nextKey);
    return keyword.compare( nextKey ) == 0;
}

void DynamicDataReader :: finish()
{
    this->recordList.clear();
}

void
DynamicDataReader :: writeToFile(const char *fileName)
{
    std :: ofstream fout(fileName);

    fout << get << '\n';
    fout << this->giveOutputFileName() << '\n';
    for ( auto &rec: this->recordList ) {
        fout << rec->giveRecordAsString() << "\n";
    }
    fout.close();
}
} // end namespace fem
