
#include "field.h"

#include <cstdarg>

namespace fem 
{
std :: string Field :: errorInfo(const char *func) const
{
    return std :: string(this->giveClassName()) + "::" + func;
}

} // end namespace oofem
