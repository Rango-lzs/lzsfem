#include "femcore/MetaObject.h"
#include "femcore/RTTI/MetaClass.h"

class DemoObject : public MetaObject
{
public:
	DemoObject() {}

	//定义如下函数，然后将其提取为宏
	virtual const const MetaClass* meta()
	{
		return mpMeta;
	}

	static const MetaClass* staic_meta()
	{
		return mpMeta;
	}

	static std::string class_name()
    {
        return "DemoObject";
    }          

	static MetaObject* meta_cast(MetaObject* pOther)
	{
        if (!pOther)
            return nullptr;
        return pOther->isKindOf(MetaObject::staic_meta()) ? static_cast<MetaObject*>(pOther) : nullptr;
	}

public:
	void member()
	{

	}
private:
	static const MetaClass* mpMeta;
};

//MetaClass* DemoObject::mpMeta = new MetaClass("DemoObject", MetaObject::meta(), []() { return new DemoObject(); });
