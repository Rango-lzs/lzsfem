
#include "text_data_reader.h"
#include "error.h"

#include <string>

namespace fem {
FEMTXTDataReader :: FEMTXTDataReader(std :: string inputfilename) : DataReader(),
    dataSourceName(std :: move(inputfilename)), recordList()
{
    std :: list< std :: pair< int, std :: string > >lines;
    // Read all the lines in the main input file:
    {
        std :: ifstream inputStream(dataSourceName);
        if ( !inputStream.is_open() ) {
            FEM_ERROR("Can't open input stream (%s)", dataSourceName.c_str());
        }

        int lineNumber = 0;
        std :: string line;

        this->giveRawLineFromInput(inputStream, lineNumber, outputFileName);
        this->giveRawLineFromInput(inputStream, lineNumber, description);

        while (this->giveLineFromInput(inputStream, lineNumber, line)) {
            lines.emplace_back(make_pair(lineNumber, line));
        }
    }
    // Check for included files: @include "somefile"
    for ( auto it = lines.begin(); it != lines.end(); ++it ) {
        if ( it->second.compare(0, 8, "@include") == 0 ) {
            std :: string fname = it->second.substr(10, it->second.length()-11);
            FEM_LOG_INFO("Reading included file: %s\n", fname.c_str());

            // Remove the include line
            lines.erase(it++);
            // Add all the included lines:
            int includedLine = 0;
            std :: string line;
            std :: ifstream includedStream(fname);
            if ( !includedStream.is_open() ) {
                FEM_ERROR("Can't open input stream (%s)", fname.c_str());
            }
            while (this->giveLineFromInput(includedStream, includedLine, line)) {
                lines.emplace(it, make_pair(includedLine, line));
            }
        }
    }
    ///@todo This could be parallelized, but I'm not sure it is worth it 
    /// (might make debugging faulty input files harder for users as well)
    for ( auto &line: lines ) {
        //printf("line: %s\n", line.second.c_str());
        this->recordList.emplace_back(line.first, line.second);
    }
    this->it = this->recordList.begin();
}

FEMTXTDataReader :: FEMTXTDataReader(const FEMTXTDataReader &x) : FEMTXTDataReader(x.dataSourceName) {}

FEMTXTDataReader :: ~FEMTXTDataReader()
{
}

InputRecord &
FEMTXTDataReader :: giveInputRecord(InputRecordType typeId, int recordId)
{
    if ( this->it == this->recordList.end() ) {
        FEM_ERROR("Out of input records, file contents must be missing");
    }
    return *this->it++;
}

bool
FEMTXTDataReader :: peakNext(const std :: string &keyword)
{
    std :: string nextKey;
    this->it->giveRecordKeywordField(nextKey);
    return keyword.compare( nextKey ) == 0;
}

void
FEMTXTDataReader :: finish()
{
    if ( this->it != this->recordList.end() ) {
        FEM_WARNING("There are unread lines in the input file\n"
            "The most common cause are missing entries in the domain record, e.g. 'nset'");
    }
    this->recordList.clear();
}

bool
FEMTXTDataReader :: giveLineFromInput(std :: ifstream &stream, int &lineNum, std :: string &line)
{
    bool flag = false; //0-tolower, 1-remain with capitals

    bool read = this->giveRawLineFromInput(stream, lineNum, line);
    if ( !read ) {
        return false;
    }

    for ( auto &c: line ) {
        if ( c == '"' ) { //do not change to lowercase inside quotation marks
            flag = !flag; // switch flag
        }

        if ( !flag ) {
            c = (char)tolower(c); // convert line to lowercase
        }
    }
    return true;
}

bool
FEMTXTDataReader :: giveRawLineFromInput(std :: ifstream &stream, int &lineNum, std :: string &line)
{
    do {
        lineNum++;
        std :: getline(stream, line);
        if ( !stream ) {
            return false;
        } if ( line.length() > 0 ) {
            if ( line.back() == '\\' ) {
                std :: string continuedLine;
                do {
                    lineNum++;
                    std :: getline(stream, continuedLine);
                    if ( !stream ) {
                        return false;
                    }
                    line.pop_back();
                    line += continuedLine;
                } while ( continuedLine.back() == '\\' );
            }
        }
    } while ( line.length() == 0 || line [ 0 ] == '#' ); // skip comments
    return true;
}
} // end namespace fem