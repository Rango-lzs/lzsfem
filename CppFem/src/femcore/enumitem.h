#ifndef enumitem_h
#define enumitem_h

namespace fem
{
#ifndef ENUM_ITEM
#define ENUM_ITEM(element) element,
#define ENUM_ITEM_WITH_VALUE(element, val) element = val,
#endif
} // end namespace fem
#endif // enumitem_h
