#include "logger/log.h"
#include "femcore/FEModel.h"
#include <stdarg.h>

void write_log(FEModel* fem, int ntag, const char* szmsg, ...)
{
	assert(fem);
	//if (fem->LogBlocked()) return;

	// get a pointer to the argument list
	va_list	args;

	// make the message
	char sztxt[2048] = { 0 };
	va_start(args, szmsg);
	vsprintf(sztxt, szmsg, args);
	va_end(args);
	//Rango TODO
	//fem->Log(ntag, sztxt);
}
