/*****************************************************************//**
 * \file   dynamic_data_reader.h
 * \brief  
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef dynamicdatareader_h
#define dynamicdatareader_h

#include "datareader.h"
#include <vector>
#include <memory>

namespace fem {
class InputRecord;

/**
 * Class representing the implementation of a dynamic data reader for in-code use.
 */
class FEM_EXPORT DynamicDataReader: public DataReader
{
protected:
    std :: string name;
    /// Keeps track of the current position in the list
    std :: vector< std :: unique_ptr< InputRecord > > :: iterator it;
    /// All record types will be appended to this list, no split in terms of InputRecordType is implemented yet.
    std :: vector< std :: unique_ptr< InputRecord > >recordList;

public:
    /// Constructor.
    DynamicDataReader(std :: string name);
    virtual ~DynamicDataReader();

    DynamicDataReader(const DynamicDataReader& src) = delete;
    DynamicDataReader &operator = (const DynamicDataReader &src) = delete;

    /**
     * Main purpose of this class it the possibility to add new input records in code.
     * The input records can be any implementation, but the intended use would be together with the DynamicInputRecord class.
     * @param type Currently ignored, but left here for consistency with giveInputRecord
     * @param record New record to be added at the end. New input records have to be added in the same order as the text input files do it.
     */
    void insertInputRecord(InputRecordType type, std::unique_ptr<InputRecord> record);

    
    InputRecord &giveInputRecord(InputRecordType, int recordId) override;

    bool peakNext(const std :: string &keyword) override;
    void finish() override;
    std :: string giveReferenceName() const override { return name; }

    /**
     * Writes all input records to file.
     * @param fileName Name of file to dump data to.
     */
    void writeToFile(const char *fileName);
};
} // end namespace fem
#endif // dynamicdatareader_h
