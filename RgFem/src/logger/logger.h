/*****************************************************************//**
 * \file   logger.h
 * \brief  
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef logger_h
#define logger_h

#include "femcore/fem_export.h"

#include <cstdio>
#include <string>
#ifdef __PARALLEL_MODE
#include <mpi.h>
#endif

// MSVC doesn't properly implement C99. (might need to wrap __func__ behind a macro to support all platforms correctly(?))
#ifdef _MSC_VER
# define __func__ __FUNCTION__
#endif

namespace fem {
/**
 * Logger class used by OOFEM to print information during analysis.
 * Prints warnings and errors into a separate stream from normal output.
 */
class FEM_EXPORT Logger
{
public:
    /// Type defining basic log levels.
    enum logLevelType {
        LOG_LEVEL_FORCED=-1,
        LOG_LEVEL_FATAL=0, LOG_LEVEL_ERROR=0,
        LOG_LEVEL_WARNING = 1,
        LOG_LEVEL_RELEVANT = 2,
        LOG_LEVEL_INFO = 3,
        LOG_LEVEL_ALL = 4, LOG_LEVEL_DEBUG = 4
    };
protected:
    /// Stream used for logging.
    FILE *logStream, *errStream;
    /// flag indicating whether to close mylogStream.
    bool closeFlag, errCloseFlag;
    /// Current log level, messages with higher level are not reported.
    logLevelType logLevel;
    /// Counter of all warning and error messages.
    int numberOfWrn, numberOfErr;
#ifdef __PARALLEL_MODE
    /// Parallell comm
    MPI_Comm comm;
#endif
public:
    Logger(logLevelType level);
    ~Logger();
    /// Redirects log output to given file name (with path).
    void appendLogTo(const std :: string &fname);
    /// Redirects error output to given file name (with path).
    void appendErrorTo(const std :: string &fname);
    /// Redirects log output to given stream.
    void appendLogTo(FILE* stream);
    /// Redirects error output to given stream.
    void appendErrorTo(FILE* stream);
#ifdef __PARALLEL_MODE
    /// Parallell comm
    void setComm(MPI_Comm comm);
#endif

    /// Writes the normal log message.
    void writeLogMsg(logLevelType level, const char *format, ...);
    /// Writes extended log message with file and line info.
    void writeELogMsg(logLevelType level, const char *_func, const char *_file, int _line, const char *format, ...);
    /// Flushes the log stream.
    void flush() { fflush(logStream); fflush(errStream); }

    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(logLevelType level) { logLevel = level; }
    /// Sets log level to given one. Only log messages with level less or equal given threshold will be printed.
    void setLogLevel(int level);
    /// Increment error count by one
    void incrementErrorCounter() {numberOfErr++;}
    /// Prints number of errors and warning logged.
    void printStatistics();

protected:
    const char *giveLevelName(logLevelType l) const;
};

extern FEM_EXPORT Logger FEM_logger;

/**
 * Log reporting macros
 */
//@{
#define FEM_LOG_FATAL(...) FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_FATAL, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define FEM_LOG_ERROR(...) FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __func__, __FILE__, __LINE__, __VA_ARGS__)
#define FEM_LOG_WARNING(...) FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING,  __func__, __FILE__, __LINE__, __VA_ARGS__)

#define FEM_LOG_FORCED(...) FEM_logger.writeLogMsg(Logger :: LOG_LEVEL_FORCED, __VA_ARGS__)
#define FEM_LOG_RELEVANT(...) FEM_logger.writeLogMsg(Logger :: LOG_LEVEL_RELEVANT, __VA_ARGS__)
#define FEM_LOG_INFO(...) FEM_logger.writeLogMsg(Logger :: LOG_LEVEL_INFO, __VA_ARGS__)
#define FEM_LOG_DEBUG(...) FEM_logger.writeLogMsg(Logger :: LOG_LEVEL_DEBUG, __VA_ARGS__)
//@}
} // end namespace fem
#endif // logger_h
