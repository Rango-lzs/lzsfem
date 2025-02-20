#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <sstream>
#include <set>

// 定义节点、单元、材料等结构
struct Node {
    int id;
    double x, y, z;
};

struct Element {
    int id;
    std::vector<int> connectivity;
};

struct Material {
    std::string name;
    double density;
    double elasticModulus;
    double poissonRatio;
};

struct ElementSet {
    std::string name;
    std::set<int> elements;
};

struct NodeSet {
    std::string name;
    std::set<int> nodes;
};

struct BoundaryCondition {
    std::string type;
    int nodeId;
    std::string dof;
    double value;
};

struct Load {
    int nodeId;
    std::string type;
    double magnitude;
};

// 数据存储结构
struct ModelData {
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::unordered_map<std::string, Material> materials;
    std::unordered_map<std::string, ElementSet> elementSets;
    std::unordered_map<std::string, NodeSet> nodeSets;
    std::vector<BoundaryCondition> boundaries;
    std::vector<Load> loads;
};

// 解析节点数据
void parseNode(std::ifstream &file, ModelData &model) {
    std::string line;
    while (std::getline(file, line) && line[0] != '*') {
        Node node;
        std::istringstream iss(line);
        iss >> node.id >> node.x >> node.y >> node.z;
        model.nodes.push_back(node);
    }
}

// 解析单元数据
void parseElement(std::ifstream &file, ModelData &model) {
    std::string line;
    while (std::getline(file, line) && line[0] != '*') {
        Element element;
        int nodeId;
        std::istringstream iss(line);
        iss >> element.id;
        while (iss >> nodeId) {
            element.connectivity.push_back(nodeId);
        }
        model.elements.push_back(element);
    }
}

// 解析材料数据
void parseMaterial(std::ifstream &file, ModelData &model) {
    std::string line;
    Material material;
    while (std::getline(file, line)) {
        if (line.find("*DENSITY") != std::string::npos) {
            std::getline(file, line);
            material.density = std::stod(line);
        } else if (line.find("*ELASTIC") != std::string::npos) {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> material.elasticModulus >> material.poissonRatio;
        } else if (line[0] == '*') {
            break;
        } else {
            material.name = line;
        }
    }
    model.materials[material.name] = material;
}

// 解析单元集合
void parseElementSet(std::ifstream &file, ModelData &model, const std::string& setName) {
    std::string line;
    ElementSet elset;
    elset.name = setName;
    while (std::getline(file, line) && line[0] != '*') {
        int elementId;
        std::istringstream iss(line);
        while (iss >> elementId) {
            elset.elements.insert(elementId);
        }
    }
    model.elementSets[setName] = elset;
}

// 解析节点集合
void parseNodeSet(std::ifstream &file, ModelData &model, const std::string& setName) {
    std::string line;
    NodeSet nset;
    nset.name = setName;
    while (std::getline(file, line) && line[0] != '*') {
        int nodeId;
        std::istringstream iss(line);
        while (iss >> nodeId) {
            nset.nodes.insert(nodeId);
        }
    }
    model.nodeSets[setName] = nset;
}

// 解析边界条件
void parseBoundary(std::ifstream &file, ModelData &model) {
    std::string line;
    while (std::getline(file, line) && line[0] != '*') {
        BoundaryCondition boundary;
        std::istringstream iss(line);
        iss >> boundary.nodeId >> boundary.dof >> boundary.value;
        model.boundaries.push_back(boundary);
    }
}

// 解析载荷
void parseLoad(std::ifstream &file, ModelData &model) {
    std::string line;
    while (std::getline(file, line) && line[0] != '*') {
        Load load;
        std::istringstream iss(line);
        iss >> load.nodeId >> load.type >> load.magnitude;
        model.loads.push_back(load);
    }
}

// 主解析函数
void parseAbaqusInp(const std::string &filename, ModelData &model) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    // 关键字到解析函数的映射
    std::unordered_map<std::string, std::function<void(std::ifstream&, ModelData&)>> parseFunctionMap = {
        {"*NODE", parseNode},
        {"*ELEMENT", parseElement},
        {"*MATERIAL", parseMaterial},
        {"*BOUNDARY", parseBoundary},
        {"*CLOAD", parseLoad}
    };

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("*ELSET") != std::string::npos) {
            auto pos = line.find("ELSET=");
            if (pos != std::string::npos) {
                std::string setName = line.substr(pos + 6);
                parseElementSet(file, model, setName);
            }
        } else if (line.find("*NSET") != std::string::npos) {
            auto pos = line.find("NSET=");
            if (pos != std::string::npos) {
                std::string setName = line.substr(pos + 5);
                parseNodeSet(file, model, setName);
            }
        } else {
            // 查找标准关键字
            for (const auto& [keyword, parseFunc] : parseFunctionMap) {
                if (line.find(keyword) == 0) {
                    parseFunc(file, model);
                    break;
                }
            }
        }
    }

    file.close();
}

int main() {
    ModelData model;
    parseAbaqusInp("model.inp", model);

    // 输出解析到的数据，检查结果
    std::cout << "节点数: " << model.nodes.size() << std::endl;
    std::cout << "单元数: " << model.elements.size() << std::endl;
    std::cout << "材料数: " << model.materials.size() << std::endl;
    std::cout << "边界条件数: " << model.boundaries.size() << std::endl;
    std::cout << "载荷数: " << model.loads.size() << std::endl;

    return 0;
}
