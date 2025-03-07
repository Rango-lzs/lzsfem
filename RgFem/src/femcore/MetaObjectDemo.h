#include "femcore/MetaObject.h"
#include "femcore/RTTI/RTTIMacroDefine.h"

class DemoObject : public MetaObject
{
    DECLARE_META_CLASS(DemoObject, MetaObject);

public:
    DemoObject() = default;

public:
    void member();
};

// MetaClass* DemoObject::mpMeta = new MetaClass("DemoObject", MetaObject::meta(), []() { return new DemoObject(); });
