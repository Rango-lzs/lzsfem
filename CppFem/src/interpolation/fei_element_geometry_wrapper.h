
#ifndef fei_element_geometry_wrapper_hh
#define fei_element_geometry_wrapper_hh

#include "fem_export.h"
#include "error.h"
#include "floatarray.h"
#include "element.h"
#include "node.h"

namespace fem
{
    /**
     * Class representing a general abstraction for cell geometry.
     * The motivation for this class is that the interpolation classes require to pass underlying cell geometry.
     * The aim here is to hide and encapsulate as much as possible from actual cell geometry specification,
     * elements describe its geometry using nodes, which are independent objects, some cells may be
     * directly specified using vertices, etc.
     */
    class FEM_EXPORT FEIElementGeometry
    {
    public:
        FEIElementGeometry() { }
        virtual ~FEIElementGeometry() { }
        virtual int giveNumberOfVertices() const = 0;
        virtual const FloatArray& giveVertexCoordinates(int i) const = 0;
    };


    /**
     * Void cell geometry wrapper.
     * Allows to use some interpolation services not needing the reference to cell geometry.
     */
    class FEM_EXPORT FEIVoidElementGeometry : public FEIElementGeometry
    {
        FloatArray tmp;
    public:
        FEIVoidElementGeometry() : FEIElementGeometry() { }
        virtual ~FEIVoidElementGeometry() { }
        int giveNumberOfVertices() const override
        {
            FEM_ERROR("no reference geometry");
            return 0;
        }
        const FloatArray& giveVertexCoordinates(int i) const override
        {
            FEM_ERROR("no reference geometry");
            return tmp;
        }
        std::string errorInfo(const char* func) const { return func; } ///@todo Class name?
    };

    /**
     * Wrapper around element definition to provide FEIElementGeometry interface.
     * 插值函数进行插值计算的时候需要单元的几何信息
     */
    class FEM_EXPORT FEIElementGeometryWrapper : public FEIElementGeometry
    {
    protected:
        const Element* elem;
    public:
        FEIElementGeometryWrapper(const Element* elem) :
            FEIElementGeometry(), elem(elem) { }
        virtual ~FEIElementGeometryWrapper() { }
        int giveNumberOfVertices() const override;
        const FloatArray& giveVertexCoordinates(int i) const override
        {
            return elem->giveNode(i)->giveCoordinates();
        }
    };

} //end namespace fem

#endif
