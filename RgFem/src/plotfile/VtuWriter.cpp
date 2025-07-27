#include "plotfile/VtuWriter.h"

#include "datastructure/Vector3d.h"
#include "elements/RgElemTypeDefine.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"


VTUWriter::VTUWriter(FEModel* pModel, std::ostream& out)
    : mp_model(pModel)
    , m_out(out)
{

}

void VTUWriter::write(const std::string& filename, const std::vector<NodeData>& nodes,
                      const std::vector<ElementData>& elements)
{
    // 0. write header
    writeHeader();

    // 1. 写入节点坐标
    writePoints(nodes);

    // 2. 写入单元信息
    writeCells(elements);

    // 3. 写入节点数据（位移、速度、加速度等）
    m_out << "      <PointData>\n";
    writeVectorData("Displacement", nodes, [](const NodeData& n) { return n.displacement; });
    writeVectorData("Velocity", nodes, [](const NodeData& n) { return n.velocity; });
    writeVectorData("Acceleration", nodes, [](const NodeData& n) { return n.acceleration; });
    writeScalarData("Temperature", nodes, [](const NodeData& n) { return n.temperature; });
    m_out << "      </PointData>\n";

    // 4. 写入单元数据（应力、应变等）
    m_out << "      <CellData>\n";
    writeScalarData("VonMisesStress", elements, [](const ElementData& e) { return e.vonMisesStress; });
    writeScalarData("PlasticStrain", elements, [](const ElementData& e) { return e.plasticStrain; });
    m_out << "      </CellData>\n";

    // 5. 写入节点相关的单元数据（需要特殊处理）
    writeNodalCellData(elements);

    // 结束文件
    m_out << "    </Piece>\n";
    m_out << "  </UnstructuredGrid>\n";
    m_out << "</VTKFile>\n";
}

void VTUWriter::writeHeader()
{
    // XML 头
    int node = mp_model->GetMesh().Nodes();
    int elem = mp_model->GetMesh().Elements();
    m_out << "<?xml version=\"1.0\"?>\n";
    m_out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    m_out << "  <UnstructuredGrid>\n";
    m_out << "    <Piece NumberOfPoints=\"" << node << "\" NumberOfCells=\"" << elem << "\">\n";
}
// 写入节点坐标
void VTUWriter::writePoints(const std::vector<NodeData>& nodes)
{
    m_out << "      <Points>\n";
    m_out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : nodes)
    {
        m_out << "          " << node.coordinates[0] << " " << node.coordinates[1] << " " << node.coordinates[2] << "\n";
    }
    m_out << "        </DataArray>\n";
    m_out << "      </Points>\n";
}

// 写入单元信息
void VTUWriter::writeCells(const std::vector<ElementData>& elements)
{
    m_out << "      <Cells>\n";

    // 连接关系
    m_out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        for (int nodeId : elem.nodeIds)
        {
            m_out << (nodeId - 1) << " ";  // VTK 使用0-based索引
        }
        m_out << "\n";
    }
    m_out << "        </DataArray>\n";

    // 偏移量
    m_out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& elem : elements)
    {
        offset += elem.nodeIds.size();
        m_out << offset << " ";
    }
    m_out << "\n        </DataArray>\n";

    // 单元类型
    m_out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        m_out << getVTKCellType(elem.type) << " ";
    }
    m_out << "\n        </DataArray>\n";

    m_out << "      </Cells>\n";
}


// 写入节点相关的单元数据（特殊处理）
void VTUWriter::writeNodalCellData(const std::vector<ElementData>& elements)
{
    // 应力张量 (6个分量: xx, yy, zz, xy, yz, xz)
    m_out << "      <PointData>\n";
    m_out << "        <DataArray type=\"Float64\" Name=\"Stress\" NumberOfComponents=\"6\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        for (const auto& stress : elem.nodalStresses)
        {
            m_out << "          ";
            for (int i = 0; i < 6; i++)
            {
                m_out << stress[i] << " ";
            }
            m_out << "\n";
        }
    }
    m_out << "        </DataArray>\n";

    // 应变张量
    m_out << "        <DataArray type=\"Float64\" Name=\"Strain\" NumberOfComponents=\"6\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        for (const auto& strain : elem.nodalStrains)
        {
            m_out << "          ";
            for (int i = 0; i < 6; i++)
            {
                m_out << strain[i] << " ";
            }
            m_out << "\n";
        }
    }
    m_out << "        </DataArray>\n";
    m_out << "      </PointData>\n";
}

// 获取VTK单元类型代码
int VTUWriter::getVTKCellType(ElementType type) const
{
    static std::map<ElementType, int> typeMap = {
        {ElementType::FE_HEX8G8, 10},
        {ElementType::FE_HEX8G8, 12},
        {ElementType::FE_HEX8G8,  9},
        {ElementType::FE_HEX8G8,  5}
    };
    return typeMap[type];
}
