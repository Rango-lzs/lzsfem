/*****************************************************************//**
 * \file   Element.h
 * \brief  The abstract base class of element
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef error_h
#define error_h

#include "fem_export.h"
#include "logger.h"

#include <string>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>


namespace fem {
/** Cause oofem program termination by calling exit. */
#define FEM_EXIT(code) \
    FEM_logger.printStatistics(); \
    fprintf(stderr, "oofem exit code %d\n", code); \
    exit(code);


class RuntimeException : public std::exception
{
public:
    std::string msg;

    RuntimeException(const char* _func, const char* _file, int _line, const char *format, ...);
    const char* what() const noexcept override;
};


/**
 * Macros for printing errors.
 * This macro can be used only within classes that implement errorInfo function.
 */
//@{
#define FEM_FATAL(...) { throw RuntimeException(__func__, __FILE__, __LINE__, __VA_ARGS__);}
#define FEM_ERROR(...) { throw RuntimeException(__func__, __FILE__, __LINE__, __VA_ARGS__);}
#define FEM_WARNING(...) FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, errorInfo(__func__).c_str(), __FILE__, __LINE__, __VA_ARGS__)
#define FEM_SERROR(...) { FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_ERROR, __func__, __FILE__, __LINE__, __VA_ARGS__); FEM_EXIT(1); }
#define FEM_SWARNING(...) FEM_logger.writeELogMsg(Logger :: LOG_LEVEL_WARNING, __func__, __FILE__, __LINE__, __VA_ARGS__)
//@}

FEM_EXPORT std::string errorInfo(const char *func);

} // end namespace fem
#endif // error_h
