#include "datastructure/Vector3d.h"
#include "elements/RgElemTypeDefine.h"

#include <fstream>
#include <vector>
#include <map>

class FEModel;

// 节点数据结构
struct NodeData
{
    Vector3d coordinates;   // 坐标
    Vector3d displacement;  // 位移
    Vector3d velocity;      // 速度
    Vector3d acceleration;  // 加速度
    double temperature;     // 温度（可选）
};

// 单元数据结构
struct ElementData
{
    ElementType type;                     // 单元类型
    std::vector<int> nodeIds;             // 节点ID
    std::vector<Vector3d> nodalStresses;  // 节点应力（每个节点一个）
    std::vector<Vector3d> nodalStrains;   // 节点应变（每个节点一个）
    double vonMisesStress;                // 等效应力（单元平均）
    double plasticStrain;                 // 塑性应变（可选）
};

class VTUWriter
{
public:
    VTUWriter(FEModel* pModel, std::ostream& out);

    void write(const std::vector<NodeData>& nodes,
               const std::vector<ElementData>& elements);
  
private:

    void writeHeader();
   
    // 写入节点坐标
    void writePoints();
   
    // 写入单元信息
    void writeCells();
    
    // 通用标量数据写入
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

    // 通用矢量数据写入
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
   
    // 写入节点相关的单元数据（特殊处理）
    void writeNodalCellData(const std::vector<ElementData>& elements);
    
    // 获取VTK单元类型代码
    int getVTKCellType(ElementType type) const;
    
private:
    std::ostream& m_out;
    FEModel* mp_model = nullptr;
};
