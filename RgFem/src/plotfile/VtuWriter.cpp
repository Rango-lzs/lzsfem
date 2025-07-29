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

    // 1. д��ڵ�����
    writePoints();

    // 2. д�뵥Ԫ��Ϣ
    writeCells();

    // 3. д��ڵ����ݣ�λ�ơ��ٶȡ����ٶȵȣ�
    const auto& allNode = mp_model->GetMesh().AllNode();
    m_out << "      <PointData>\n";
    writeVectorData("Displacement", allNode, [](const FENode& n) { return n.get_vec3d(0, 1, 2); });
    //writeVectorData("Velocity", nodes, [](const NodeData& n) { return n.velocity; });
    //writeVectorData("Acceleration", nodes, [](const NodeData& n) { return n.acceleration; });
    //writeScalarData("Temperature", nodes, [](const NodeData& n) { return n.temperature; });
    m_out << "      </PointData>\n";

    // 4. д�뵥Ԫ���ݣ�Ӧ����Ӧ��ȣ�
    //m_out << "      <CellData>\n";
    //writeScalarData("VonMisesStress", elements, [](const ElementData& e) { return e.vonMisesStress; });
    //writeScalarData("PlasticStrain", elements, [](const ElementData& e) { return e.plasticStrain; });
    //m_out << "      </CellData>\n";

    // 5. д��ڵ���صĵ�Ԫ���ݣ���Ҫ���⴦��
    //writeNodalCellData(elements);

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
    m_out <<SPACE2 <<"<UnstructuredGrid>\n";
    m_out <<SPACE4<<"<Piece NumberOfPoints=\"" << node << "\" NumberOfCells=\"" << elem << "\">\n";
}
// д��ڵ�����
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

// д�뵥Ԫ��Ϣ
void VTUWriter::writeCells()
{
    m_out << SPACE6<< "<Cells>\n";

    // ���ӹ�ϵ
    m_out << SPACE8<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";


    for (int i =0; i< mp_model->GetMesh().Elements();i++)
    {
        m_out << SPACE10; 
        const FEElement& elem = *mp_model->GetMesh().Element(i);
        for (int nodeId : elem.getNodeIds())
        {
            m_out << nodeId << " ";  // VTK ʹ��0-based����
        }
        m_out << "\n";
    }
    m_out << SPACE8<< "</DataArray>\n";

    // ƫ����
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

    // ��Ԫ����
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
        {ElementType::FE_HEX8G8, 12},
        {ElementType::FE_HEX8G8, 12},
        {ElementType::FE_HEX8G8,  9},
        {ElementType::FE_HEX8G8,  5}
    };
    return typeMap[type];
}
