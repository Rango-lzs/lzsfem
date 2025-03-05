
class MetaClass;
// �궨�帨��ע������
#define DECLARE_RTTI_CLASS(ClassName, BaseClass)                                                                       \
public:                                                                                                                \
    static const MetaClass* static_meta()                                                                             \
    {                                                                                                                  \
        static const MetaClass MetaClass(#ClassName, BaseClass::static_meta(),                                         \
                                       []() -> std::unique_ptr<MetaObject> { return std::make_unique<ClassName>(); }); \
        return &MetaClass;                                                                                              \
    }                                                                                                                  \
    const MetaClass* meta() const override                                                                           \
    {                                                                                                                  \
        return static_meta();                                                                                        \
    }                                                                                                                  \
    void RegisterType() const                                                                                          \
    {                                                                                                                  \
        RTTIFactory::GetInstance().RegisterType(GetStaticType());                                                      \
    }                                                                                                                  \
    ClassName()                                                                                                        \
    {                                                                                                                  \
        RegisterType();                                                                                                \
    }


#define DEFINE_RTTI_CLASS(ClassName)                                                                                   \
    static const MetaClass* GetStaticType()                                                                             \
    {                                                                                                                  \
        static const MetaClass MetaClass(#ClassName, nullptr,                                                            \
                                       []() -> std::unique_ptr<RTTIObject> { return std::make_unique<ClassName>(); }); \
        return &MetaClass;                                                                                              \
    }                                                                                                                  \
    const MetaClass* GetType() const override                                                                           \
    {                                                                                                                  \
        return GetStaticType();                                                                                        \
    }                                                                                                                  \
    void RegisterType() const                                                                                          \
    {                                                                                                                  \
        RTTIFactory::GetInstance().RegisterType(GetStaticType());                                                      \
    }                                                                                                                  \
    ClassName()                                                                                                        \
    {                                                                                                                  \
        RegisterType();                                                                                                \
    }



// ʾ�����νṹ
class Base : public MetaObject
{
    DEFINE_RTTI_CLASS(Base)
};

class Derived : public Base
{
    DECLARE_RTTI_CLASS(Derived, Base)
};

class FurtherDerived : public Derived
{
    DECLARE_RTTI_CLASS(FurtherDerived, Derived)
};

int main()
{
    // ͨ����������ʵ��
    auto base = RTTIFactory::GetInstance().Create("Base");
    auto derived = RTTIFactory::GetInstance().Create("Derived");

    std::cout << "Base type: " << base->GetType()->GetName() << std::endl;
    std::cout << "Derived type: " << derived->GetType()->GetName() << std::endl;

    // ��̬ת������
    if (auto* derivedPtr = derived->DynamicCast<Base>())
    {
        std::cout << "Successfully cast to Base" << std::endl;
    }

    return 0;
}