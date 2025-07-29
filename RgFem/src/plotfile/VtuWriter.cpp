#include "plotfile/VtuWriter.h"

#include "datastructure/Vector3d.h"
#include "elements/RgElemTypeDefine.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"

#define SPACE2 "  "
#define SPACE4 "    "
#define SPACE6 "      "
#define SPACE8 "        "
#define SPACE10 "          "

VTUWriter::VTUWriter(FEModel* pModel, std::ostream& out)
    : mp_model(pModel)
    , m_out(out)
{

}

void VTUWriter::write(const std::vector<NodeData>& nodes,
                      const std::vector<ElementData>& elements)
{
    // 0. write header
    writeHeader();

    // 1. 写入节点坐标
    writePoints();

    // 2. 写入单元信息
    writeCells();

    // 3. 写入节点数据（位移、速度、加速度等）
    const auto& allNode = mp_model->GetMesh().AllNode();
    m_out << "      <PointData>\n";
    writeVectorData("Displacement", allNode, [](const FENode& n) { return n.get_vec3d(0, 1, 2); });
    //writeVectorData("Velocity", nodes, [](const NodeData& n) { return n.velocity; });
    //writeVectorData("Acceleration", nodes, [](const NodeData& n) { return n.acceleration; });
    //writeScalarData("Temperature", nodes, [](const NodeData& n) { return n.temperature; });
    m_out << "      </PointData>\n";

    // 4. 写入单元数据（应力、应变等）
    //m_out << "      <CellData>\n";
    //writeScalarData("VonMisesStress", elements, [](const ElementData& e) { return e.vonMisesStress; });
    //writeScalarData("PlasticStrain", elements, [](const ElementData& e) { return e.plasticStrain; });
    //m_out << "      </CellData>\n";

    // 5. 写入节点相关的单元数据（需要特殊处理）
    //writeNodalCellData(elements);

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
    m_out <<SPACE2 <<"<UnstructuredGrid>\n";
    m_out <<SPACE4<<"<Piece NumberOfPoints=\"" << node << "\" NumberOfCells=\"" << elem << "\">\n";
}
// 写入节点坐标
void VTUWriter::writePoints()
{
    m_out << SPACE6<< "<Points>\n";
    m_out << SPACE8<< "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    
    for (int i = 0; i < mp_model->GetMesh().Nodes(); i++)
    {
        const FENode& node = mp_model->GetMesh().Node(i);
        m_out << "          " << node.m_r0[0] << " " << node.m_r0[1] << " " << node.m_r0[2] << "\n";
    }
    m_out << SPACE8<< "</DataArray>\n";
    m_out << SPACE6<< "</Points>\n";
}

// 写入单元信息
void VTUWriter::writeCells()
{
    m_out << SPACE6<< "<Cells>\n";

    // 连接关系
    m_out << SPACE8<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";


    for (int i =0; i< mp_model->GetMesh().Elements();i++)
    {
        m_out << SPACE10; 
        const FEElement& elem = *mp_model->GetMesh().Element(i);
        for (int nodeId : elem.getNodeIds())
        {
            m_out << nodeId << " ";  // VTK 使用0-based索引
        }
        m_out << "\n";
    }
    m_out << SPACE8<< "</DataArray>\n";

    // 偏移量
    m_out << SPACE8<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    m_out << SPACE10; 
    int offset = 0;
    for (int i = 0; i < mp_model->GetMesh().Elements(); i++)
    {
        const FEElement& elem = *mp_model->GetMesh().Element(i);
        offset += elem.NodeSize();
        m_out << offset << " ";
    }
    m_out << "\n" << SPACE8<< "</DataArray>\n";

    // 单元类型
    m_out << SPACE8<< "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    m_out << SPACE10; 
    for (int i = 0; i < mp_model->GetMesh().Elements(); i++)
    {
        const FEElement& elem = *mp_model->GetMesh().Element(i);
        m_out << 12 << " ";
    }
    m_out << "\n" << SPACE8 <<"</ DataArray >\n ";

    m_out << SPACE6<<"</Cells>\n";
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
        {ElementType::FE_HEX8G8, 12},
        {ElementType::FE_HEX8G8, 12},
        {ElementType::FE_HEX8G8,  9},
        {ElementType::FE_HEX8G8,  5}
    };
    return typeMap[type];
}
