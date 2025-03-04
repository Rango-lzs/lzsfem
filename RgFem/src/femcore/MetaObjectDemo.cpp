#include "MetaObjectDemo.h"
#include "femcore/RTTI/MetaClass.h"

const MetaClass* DemoObject::mpMeta = ConcretMeta<DemoObject>::instance();
