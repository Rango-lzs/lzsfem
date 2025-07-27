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

    // 1. д��ڵ�����
    writePoints(nodes);

    // 2. д�뵥Ԫ��Ϣ
    writeCells(elements);

    // 3. д��ڵ����ݣ�λ�ơ��ٶȡ����ٶȵȣ�
    m_out << "      <PointData>\n";
    writeVectorData("Displacement", nodes, [](const NodeData& n) { return n.displacement; });
    writeVectorData("Velocity", nodes, [](const NodeData& n) { return n.velocity; });
    writeVectorData("Acceleration", nodes, [](const NodeData& n) { return n.acceleration; });
    writeScalarData("Temperature", nodes, [](const NodeData& n) { return n.temperature; });
    m_out << "      </PointData>\n";

    // 4. д�뵥Ԫ���ݣ�Ӧ����Ӧ��ȣ�
    m_out << "      <CellData>\n";
    writeScalarData("VonMisesStress", elements, [](const ElementData& e) { return e.vonMisesStress; });
    writeScalarData("PlasticStrain", elements, [](const ElementData& e) { return e.plasticStrain; });
    m_out << "      </CellData>\n";

    // 5. д��ڵ���صĵ�Ԫ���ݣ���Ҫ���⴦��
    writeNodalCellData(elements);

    // �����ļ�
    m_out << "    </Piece>\n";
    m_out << "  </UnstructuredGrid>\n";
    m_out << "</VTKFile>\n";
}

void VTUWriter::writeHeader()
{
    // XML ͷ
    int node = mp_model->GetMesh().Nodes();
    int elem = mp_model->GetMesh().Elements();
    m_out << "<?xml version=\"1.0\"?>\n";
    m_out << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    m_out << "  <UnstructuredGrid>\n";
    m_out << "    <Piece NumberOfPoints=\"" << node << "\" NumberOfCells=\"" << elem << "\">\n";
}
// д��ڵ�����
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

// д�뵥Ԫ��Ϣ
void VTUWriter::writeCells(const std::vector<ElementData>& elements)
{
    m_out << "      <Cells>\n";

    // ���ӹ�ϵ
    m_out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        for (int nodeId : elem.nodeIds)
        {
            m_out << (nodeId - 1) << " ";  // VTK ʹ��0-based����
        }
        m_out << "\n";
    }
    m_out << "        </DataArray>\n";

    // ƫ����
    m_out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& elem : elements)
    {
        offset += elem.nodeIds.size();
        m_out << offset << " ";
    }
    m_out << "\n        </DataArray>\n";

    // ��Ԫ����
    m_out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& elem : elements)
    {
        m_out << getVTKCellType(elem.type) << " ";
    }
    m_out << "\n        </DataArray>\n";

    m_out << "      </Cells>\n";
}


// д��ڵ���صĵ�Ԫ���ݣ����⴦��
void VTUWriter::writeNodalCellData(const std::vector<ElementData>& elements)
{
    // Ӧ������ (6������: xx, yy, zz, xy, yz, xz)
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

    // Ӧ������
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

// ��ȡVTK��Ԫ���ʹ���
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
