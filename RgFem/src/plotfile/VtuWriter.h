#include "datastructure/Vector3d.h"
#include "elements/RgElemTypeDefine.h"

#include <fstream>
#include <vector>
#include <map>

class FEModel;

// �ڵ����ݽṹ
struct NodeData
{
    Vector3d coordinates;   // ����
    Vector3d displacement;  // λ��
    Vector3d velocity;      // �ٶ�
    Vector3d acceleration;  // ���ٶ�
    double temperature;     // �¶ȣ���ѡ��
};

// ��Ԫ���ݽṹ
struct ElementData
{
    ElementType type;                     // ��Ԫ����
    std::vector<int> nodeIds;             // �ڵ�ID
    std::vector<Vector3d> nodalStresses;  // �ڵ�Ӧ����ÿ���ڵ�һ����
    std::vector<Vector3d> nodalStrains;   // �ڵ�Ӧ�䣨ÿ���ڵ�һ����
    double vonMisesStress;                // ��ЧӦ������Ԫƽ����
    double plasticStrain;                 // ����Ӧ�䣨��ѡ��
};

class VTUWriter
{
public:
    VTUWriter(FEModel* pModel, std::ostream& out);

    void write(const std::vector<NodeData>& nodes,
               const std::vector<ElementData>& elements);
  
private:

    void writeHeader();
   
    // д��ڵ�����
    void writePoints();
   
    // д�뵥Ԫ��Ϣ
    void writeCells();
    
    // ͨ�ñ�������д��
    template <typename T, typename Func>
    void writeScalarData(const std::string& name, const std::vector<T>& items, Func getValue)
    {
        m_out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (const auto& item : items)
        {
            m_out << "          " << getValue(item) << "\n";
        }
        m_out << "        </DataArray>\n";
    }

    // ͨ��ʸ������д��
    template <typename T, typename Func>
    void writeVectorData(const std::string& name, const std::vector<T>& items, Func getVector)
    {
        m_out << "        <DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& item : items)
        {
            const auto& vec = getVector(item);
            m_out << "          " << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
        }
        m_out << "        </DataArray>\n";
    }
   
    // д��ڵ���صĵ�Ԫ���ݣ����⴦��
    void writeNodalCellData(const std::vector<ElementData>& elements);
    
    // ��ȡVTK��Ԫ���ʹ���
    int getVTKCellType(ElementType type) const;
    
private:
    std::ostream& m_out;
    FEModel* mp_model = nullptr;
};
