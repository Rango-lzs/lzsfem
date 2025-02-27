#pragma once
#include "femcore/fem_export.h"

class FEModel;

FEM_EXPORT void write_log(FEModel* fem, int ntag, const char* szmsg, ...);

#define feLog(...) write_log(GetFEModel(), 0, __VA_ARGS__)
#define feLogWarning(...) write_log(GetFEModel(), 1, __VA_ARGS__)
#define feLogError(...) write_log(GetFEModel(), 2, __VA_ARGS__)
#define feLogInfo(...) write_log(GetFEModel(), 3, __VA_ARGS__)
#define feLogDebug(...) write_log(GetFEModel(), 4, __VA_ARGS__)

#define feLogEx(fem, ...) write_log(fem, 0, __VA_ARGS__)
#define feLogWarningEx(fem, ...) write_log(fem, 1, __VA_ARGS__)
#define feLogErrorEx(fem, ...) write_log(fem, 2, __VA_ARGS__)
#define feLogInfoEx(fem, ...) write_log(fem, 3, __VA_ARGS__)
#define feLogDebugEx(fem, ...) write_log(fem, 4, __VA_ARGS__)
