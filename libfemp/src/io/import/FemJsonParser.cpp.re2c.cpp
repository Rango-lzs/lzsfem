#include "FemJsonParser.hpp"

#include <iostream>
#include <string.h> // for memcpy

#include <boost/lexical_cast.hpp>

#include <libfemp/Element.hpp>
#include <libfemp/FemException.hpp>
#include <libfemp/Material.hpp>

FemJsonParser::FemJsonParser()
{
    limit = buffer;
    tok = limit;
    pos = limit;
    marker = limit;
    *limit = '\000';
    line_number = 1;

    setParserTable();

    // resets the error message
    this->error.code = Error::Type::ERR_OK;
    this->error.message.clear();
}

void FemJsonParser::operator()(std::istream& file, fem::Model& model)
{
    parse(file, model);
}

enum FemJsonParser::Error::Type FemJsonParser::parse(std::istream& file, fem::Model& model)
{
    // temp variables used by the parser
    double temp_material_E = 0;
    double temp_material_nu = 0;
    double temp_material_fy = 0;

    ss.push(NT_START); // initializes the symbol stack with  the start symbol

    enum TerminalTokens terminal_token = lexer(file);
    int nonterminal_token;

    // resets the error message
    this->error.code = Error::Type::ERR_OK;
    this->error.message.clear();

    // temporary object
    fem::Element element;

    while (!ss.empty()) {
        if (terminal_token == ss.top()) {
            ss.pop();
            //std::cout << names(terminal_token) << " << ";
            terminal_token = lexer(file);
            //std::cout << names(terminal_token) << std::endl;
        } else {
            nonterminal_token = ss.top();
            ss.pop();
            //std::cout << "table[\t" << names(nonterminal_token) << "\t]\t[\t" << names(terminal_token) << "\t]" << std::endl;
            switch (table[(enum NonTerminalTokens)nonterminal_token][terminal_token]) {
            case PR_LABEL:
                ss.push(NT_LABEL_SET);
                ss.push(TT_TEXT_STRING);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_LABEL);
                break;

            case PR_LABEL_SET:
                label = temp_token;
                break;

            case PR_VECTOR3D:
                ss.push(NT_VECTOR3D_SET);
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_NUMBER);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_NUMBER);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_NUMBER);
                ss.push(TT_ARRAY_OPEN);
                break;

            case PR_VECTOR3D_SET: {
                fem::Point3D p;
                p.data[2] = temp_float.back();
                temp_float.pop_back();
                p.data[1] = temp_float.back();
                temp_float.pop_back();
                p.data[0] = temp_float.back();
                temp_float.pop_back();

                temp_point.push_back(p);
            } break;

            case PR_REFERENCE_1:
                ss.push(NT_REFERENCE_SET);
                ss.push(TT_NATURAL_NUMBER);
                break;

            case PR_REFERENCE_ERROR:
                this->error.code = Error::Type::ERR_UNKNOWN; //TODO return a dedicated error
                this->error.message = "Expected a reference number at line " + boost::lexical_cast<std::string>(line_number);
                return this->error.code;
                break;

            case PR_REFERENCE_SET:
                try {
                    reference.push_back(boost::lexical_cast<size_t>(temp_token));
                    temp_token.clear();
                } catch (boost::bad_lexical_cast& e) {
                    std::cerr << "PR_REFERENCE_SET" << e.what() << std::endl;
                    std::cerr << temp_token << std::endl;
                }
                break;

            case PR_REFERENCE_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_REFERENCE_FOLLOW);
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                break;

            case PR_REFERENCE_FOLLOW_1:
                ss.push(NT_REFERENCE_FOLLOW);
                ss.push(NT_REFERENCE);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_REFERENCE_FOLLOW_2:
                break;

            case PR_MULTI_REFERENCE_1:
                ss.push(NT_REFERENCE);
                break;

            case PR_MULTI_REFERENCE_2:
                ss.push(NT_MULTI_REFERENCE_SET);
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_REFERENCE);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                break;

            case PR_MULTI_REFERENCE_SET: {
                size_t last = reference.back();
                reference.pop_back();
                size_t first = reference.back();
                reference.pop_back();
                for (size_t i = first; i <= last; i++) {
                    reference.push_back(i);
                }
            } break;

            case PR_BOOLEAN_TRUE:
                //TODO push boolean value
                ss.push(TT_TRUE);
                break;

            case PR_BOOLEAN_FALSE:
                //TODO push boolean value
                ss.push(TT_FALSE);
                break;

            case PR_NUMBER_1:
                ss.push(NT_FLOAT_SET);
                ss.push(TT_FLOAT_NUMBER);
                break;

            case PR_NUMBER_2:
                ss.push(NT_INT_TO_FLOAT_SET);
                ss.push(TT_NATURAL_NUMBER);
                break;

            case PR_FLOAT_SET:
                try {
                    temp_float.push_back(boost::lexical_cast<double>(temp_token));
                    temp_token.clear();
                } catch (boost::bad_lexical_cast& e) {
                    std::cerr << "PR_FLOAT_SET" << e.what() << std::endl;
                    std::cerr << temp_token << std::endl;
                }
                break;

            case PR_INT_TO_FLOAT_SET:
                try {
                    temp_float.push_back(boost::lexical_cast<double>(temp_token));
                    temp_token.clear();
                } catch (boost::bad_lexical_cast& e) {
                    std::cerr << "PR_INT_TO_FLOAT_SET" << e.what() << std::endl;
                    std::cerr << temp_token << std::endl;
                }
                break;

            case PR_START_1:
                ss.push(TT_EOF);
                ss.push(NT_DOCUMENT);
                break;

            case PR_START_2:
                ss.push(TT_EOF);
                break;

            case PR_DOCUMENT_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_ELEMENT_LIST_FOLLOW);
                ss.push(NT_ELEMENT_LIST);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_NODE_LIST);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_MATERIAL_LIST);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_HEADER);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_HEADER_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_MODEL_TYPE);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_FORMAT_VERSION);
                ss.push(TT_OBJECT_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FEM);
                break;

            case PR_FORMAT_VERSION_1:
                ss.push(TT_LABEL_V1_0);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_VERSION);
                break;

            case PR_MODEL_TYPE_1:
                ss.push(NT_MODEL_TYPE_SET); // sets the model type
                ss.push(TT_LABEL_3D_SOLID);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_TYPE);
                break;

            case PR_MODEL_TYPE_SET:
                // due to poor design, it isn't possible to set the model type
                break;

            case PR_MATERIAL_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_MATERIAL_FOLLOW);
                ss.push(NT_MATERIAL_SET);
                ss.push(NT_MATERIAL);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_MATERIALS);
                break;

            case PR_MATERIAL_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_ELASTIC_LIMIT);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_POISSON_COEFFICIENT);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_YOUNG_MODULUS);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_LABEL);
                ss.push(TT_VALUE_SEPARATOR);

                // material is linear elastic

                ss.push(NT_MATERIAL_TYPE);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_MATERIAL_TYPE_1:
                ss.push(TT_LABEL_LINEAR_ELASTIC);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_TYPE);
                break;

            case PR_YOUNG_MODULUS_1:
                ss.push(NT_YOUNG_MODULUS_SET);
                ss.push(NT_NUMBER);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_E);
                break;

            case PR_YOUNG_MODULUS_SET:
                temp_material_E = temp_float.back();
                temp_float.pop_back();
                break;

            case PR_POISSON_COEFFICIENT_1:
                ss.push(NT_POISSON_COEFFICIENT_SET);
                ss.push(NT_NUMBER);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NU);
                break;

            case PR_POISSON_COEFFICIENT_SET:
                temp_material_nu = temp_float.back();
                temp_float.pop_back();
                break;

            case PR_ELASTIC_LIMIT_1:
                ss.push(NT_ELASTIC_LIMIT_SET);
                ss.push(NT_NUMBER);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FY);
                break;

            case PR_ELASTIC_LIMIT_SET:
                temp_material_fy = temp_float.back();
                temp_float.pop_back();
                break;

            case PR_MATERIAL_SET: {
                fem::Material material(label, temp_material_E, temp_material_nu);
                model.pushMaterial(material);
            } break;

            case PR_MATERIAL_FOLLOW_1:
                ss.push(NT_MATERIAL_FOLLOW);
                ss.push(NT_MATERIAL_SET);
                ss.push(NT_MATERIAL);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_MATERIAL_FOLLOW_2:
                break;

            case PR_NODE_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_NODE_FOLLOW);
                ss.push(NT_NODE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_NODE_1:
                ss.push(NT_NODE_SET);
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_VECTOR3D);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                break;

            case PR_NODE_SET:
                model.setNode(reference.back(), temp_point.back());
                reference.pop_back();
                temp_point.pop_back();
                break;

            case PR_NODE_FOLLOW_1:
                ss.push(NT_NODE_FOLLOW);
                ss.push(NT_NODE);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_FOLLOW_2:
                break;

            case PR_ELEMENT_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_ELEMENT_FOLLOW);
                ss.push(NT_ELEMENT);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_ELEMENTS);
                break;

            case PR_ELEMENT_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_ELEMENT_SET);
                ss.push(NT_ELEMENT_TYPE_FOLLOW);
                ss.push(NT_ELEMENT_TYPE);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_TYPE);
                ss.push(TT_OBJECT_OPEN);

                reference.clear();
                break;

            case PR_ELEMENT_SET_1:
                model.pushElement(element);
                break;

            case PR_ELEMENT_TYPE_1:
                ss.push(NT_ELEMENT_HEXAHEDRON8);
                break;

            case PR_ELEMENT_TYPE_2:
                ss.push(NT_ELEMENT_HEXAHEDRON20);
                break;

            case PR_ELEMENT_TYPE_3:
                ss.push(NT_ELEMENT_HEXAHEDRON27);
                break;

            case PR_ELEMENT_TYPE_4:
                ss.push(NT_ELEMENT_TETRAHEDRON4);
                break;

            case PR_ELEMENT_TYPE_5:
                ss.push(NT_ELEMENT_TETRAHEDRON10);
                break;
            case PR_ELEMENT_TYPE_6:
                ss.push(NT_ELEMENT_PRISM6);
                break;

            case PR_ELEMENT_TYPE_7:
                ss.push(NT_ELEMENT_PRISM15);
                break;

            case PR_ELEMENT_TYPE_8:
                ss.push(NT_ELEMENT_PRISM18);
                break;

            case PR_ELEMENT_HEXAHEDRON8_1:
                ss.push(NT_ELEMENT_HEXAHEDRON8_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 7; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_HEXAHEDRON8);
                break;

            case PR_ELEMENT_HEXAHEDRON20_1:
                ss.push(NT_ELEMENT_HEXAHEDRON20_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 19; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_HEXAHEDRON20);
                break;

            case PR_ELEMENT_HEXAHEDRON27_1:
                ss.push(NT_ELEMENT_HEXAHEDRON27_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 26; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_HEXAHEDRON27);
                break;

            case PR_ELEMENT_TETRAHEDRON4_1:
                ss.push(NT_ELEMENT_TETRAHEDRON4_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 3; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_TETRAHEDRON4);
                break;

            case PR_ELEMENT_TETRAHEDRON10_1:
                ss.push(NT_ELEMENT_TETRAHEDRON10_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 9; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_TETRAHEDRON10);
                break;

            case PR_ELEMENT_PRISM6_1:
                ss.push(NT_ELEMENT_PRISM6_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 5; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_PRISM6);
                break;

            case PR_ELEMENT_PRISM15_1:
                ss.push(NT_ELEMENT_PRISM15_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 14; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_PRISM15);
                break;

            case PR_ELEMENT_PRISM18_1:
                ss.push(NT_ELEMENT_PRISM18_SET);
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 17; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_PRISM18);
                break;

            case PR_ELEMENT_HEXAHEDRON8_SET:
                element.clear();
                element.type = fem::Element::FE_HEXAHEDRON8;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_HEXAHEDRON20_SET:
                element.clear();
                element.type = fem::Element::FE_HEXAHEDRON20;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_HEXAHEDRON27_SET:
                element.clear();
                element.type = fem::Element::FE_HEXAHEDRON27;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_TETRAHEDRON4_SET:
                element.clear();
                element.type = fem::Element::FE_TETRAHEDRON4;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_TETRAHEDRON10_SET:
                element.clear();
                element.type = fem::Element::FE_TETRAHEDRON10;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_PRISM6_SET:
                element.clear();
                element.type = fem::Element::FE_PRISM6;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_PRISM15_SET:
                element.clear();
                element.type = fem::Element::FE_PRISM15;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_PRISM18_SET:
                element.clear();
                element.type = fem::Element::FE_PRISM18;
                element.nodes = reference;
                reference.clear();
                break;

            case PR_ELEMENT_TYPE_FOLLOW_1:
                ss.push(NT_ELEMENT_MATERIAL);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_ELEMENT_TYPE_FOLLOW_2:
                break;

            case PR_ELEMENT_MATERIAL_1:
                ss.push(NT_ELEMENT_MATERIAL_SET);
                ss.push(NT_REFERENCE);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_MATERIAL);
                break;

            case PR_ELEMENT_MATERIAL_SET_1:
                model.setDefaultMaterial(reference.back());
                reference.pop_back();
                break;

            case PR_ELEMENT_FOLLOW_1:
                ss.push(NT_ELEMENT_FOLLOW);
                ss.push(NT_ELEMENT);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_ELEMENT_FOLLOW_2:
                break;

            case PR_ELEMENT_LIST_FOLLOW_1:
                ss.push(NT_EXTRA_DEFINITIONS);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_ELEMENT_LIST_FOLLOW_2:
                break;

            case PR_EXTRA_DEFINITIONS_1:
                ss.push(NT_GROUPS_LIST_FOLLOW);
                ss.push(NT_GROUPS_LIST);
                break;

            case PR_EXTRA_DEFINITIONS_2:
                ss.push(NT_NODE_RESTRICTIONS_LIST_FOLLOW);
                ss.push(NT_NODE_RESTRICTIONS_LIST);
                break;

            case PR_EXTRA_DEFINITIONS_3:
                ss.push(NT_LOAD_PATTERNS_LIST_FOLLOW);
                ss.push(NT_LOAD_PATTERNS_LIST);
                break;

            case PR_GROUPS_LIST_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_GROUP_TYPE);
                ss.push(TT_OBJECT_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_GROUPS);
                break;

            case PR_GROUP_TYPE_1:
                ss.push(NT_NODE_GROUP_LIST_FOLLOW);
                ss.push(NT_NODE_GROUP_LIST);
                break;

            case PR_GROUP_TYPE_2:
                ss.push(NT_ELEMENT_GROUP_LIST_FOLLOW);
                ss.push(NT_ELEMENT_GROUP_LIST);
                break;

            case PR_GROUPS_LIST_FOLLOW_1:
                ss.push(NT_GROUPS_LIST_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_GROUPS_LIST_FOLLOW_2:
                break;

            case PR_GROUPS_LIST_FOLLOW_FOLLOW_1:
                ss.push(NT_NODE_RESTRICTIONS_LIST_FOLLOW);
                ss.push(NT_NODE_RESTRICTIONS_LIST);
                break;

            case PR_GROUPS_LIST_FOLLOW_FOLLOW_2:
                ss.push(NT_LOAD_PATTERNS_LIST_FOLLOW);
                ss.push(NT_LOAD_PATTERNS_LIST);
                break;

            case PR_NODE_GROUP_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_NODE_GROUP_FOLLOW);
                ss.push(NT_NODE_GROUP_SET);
                ss.push(NT_NODE_GROUP);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_NODE_GROUP_LIST_FOLLOW_1:
                ss.push(NT_ELEMENT_GROUP_LIST);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_GROUP_LIST_FOLLOW_2:
                break;

            case PR_NODE_GROUP_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_REFERENCE_LIST);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_LABEL);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_NODE_GROUP_SET: {
                fem::NodeGroup n;
                n.setLabel(label);
                for (std::vector<size_t>::const_iterator i = reference.begin(); i != reference.end(); i++) {
                    if (*i < model.getNodeMap().size())
                        n.pushNode(*i);
                    else {
                        //TODO warning message
                        std::cerr << "FemJsonParser: node group " << label << " refers to invalid node: " << *i << std::endl;
                    }
                }
                reference.clear();
                model.pushNodeGroup(n);
            } break;

            case PR_NODE_GROUP_FOLLOW_1:
                ss.push(NT_NODE_GROUP_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_GROUP_FOLLOW_2:
                break;

            case PR_NODE_GROUP_FOLLOW_FOLLOW_1:
                ss.push(NT_NODE_GROUP_FOLLOW);
                ss.push(NT_NODE_GROUP_SET);
                ss.push(NT_NODE_GROUP);
                break;

            case PR_ELEMENT_GROUP_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_ELEMENT_GROUP_FOLLOW);
                ss.push(NT_ELEMENT_GROUP_SET);
                ss.push(NT_ELEMENT_GROUP);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_ELEMENTS);
                break;

            case PR_ELEMENT_GROUP_LIST_FOLLOW_1:
                break;

            case PR_ELEMENT_GROUP_1:
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_REFERENCE_LIST);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_ELEMENTS);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_LABEL);
                ss.push(TT_OBJECT_OPEN);
                reference.clear();
                break;

            case PR_ELEMENT_GROUP_SET: {
                fem::ElementGroup e;
                e.setLabel(label);
                for (std::vector<size_t>::const_iterator i = reference.begin(); i != reference.end(); i++) {
                    e.pushElement(*i);
                }
                reference.clear();
                model.pushElementGroup(e);
            } break;

            case PR_ELEMENT_GROUP_FOLLOW_1:
                ss.push(NT_ELEMENT_GROUP_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_ELEMENT_GROUP_FOLLOW_2:
                break;

            case PR_ELEMENT_GROUP_FOLLOW_FOLLOW_1:
                ss.push(NT_ELEMENT_GROUP_FOLLOW);
                ss.push(NT_ELEMENT_GROUP_SET);
                ss.push(NT_ELEMENT_GROUP);
                break;

            case PR_NODE_RESTRICTIONS_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_NODE_RESTRICTION_FOLLOW);
                ss.push(NT_NODE_RESTRICTION);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODE_RESTRICTIONS);
                break;

            case PR_NODE_RESTRICTIONS_LIST_FOLLOW_1:
                ss.push(NT_NODE_RESTRICTIONS_LIST_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_RESTRICTIONS_LIST_FOLLOW_2:
                break;

            case PR_NODE_RESTRICTIONS_LIST_FOLLOW_FOLLOW_1:
                ss.push(NT_LOAD_PATTERNS_LIST_FOLLOW);
                ss.push(NT_LOAD_PATTERNS_LIST);
                break;

            case PR_NODE_RESTRICTION_1:
                ss.push(NT_NODE_RESTRICTION_SET);
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_NODE_DOF);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_NODE_REFERENCE);
                ss.push(TT_OBJECT_OPEN);

                temp_node_restrictions.reset();
                break;

            case PR_NODE_REFERENCE_1:
                ss.push(NT_NODE_REFERENCE_SET);
                ss.push(NT_MULTI_REFERENCE);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODE);
                reference.clear();
                break;

            case PR_NODE_REFERENCE_SET:
                break;

            case PR_NODE_DOF_1:
                ss.push(NT_NODE_DOF_DX_FOLLOW);
                ss.push(NT_NODE_DOF_DX);
                break;

            case PR_NODE_DOF_2:
                ss.push(NT_NODE_DOF_DY_FOLLOW);
                ss.push(NT_NODE_DOF_DY);
                break;

            case PR_NODE_DOF_3:
                ss.push(NT_NODE_DOF_DZ_FOLLOW);
                ss.push(NT_NODE_DOF_DZ);
                break;

            case PR_NODE_DOF_DX_1:
                ss.push(NT_NODE_DOF_DX_SET);
                ss.push(NT_BOOLEAN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_DX);
                break;

            case PR_NODE_DOF_DX_SET:
                temp_node_restrictions.setdx();
                break;

            case PR_NODE_DOF_DX_FOLLOW_1:
                ss.push(NT_NODE_DOF_DX_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_DOF_DX_FOLLOW_2:
                break;

            case PR_NODE_DOF_DX_FOLLOW_FOLLOW_1:
                ss.push(NT_NODE_DOF_DY_FOLLOW);
                ss.push(NT_NODE_DOF_DY);
                break;

            case PR_NODE_DOF_DX_FOLLOW_FOLLOW_2:
                ss.push(NT_NODE_DOF_DZ_FOLLOW);
                ss.push(NT_NODE_DOF_DZ);
                break;

            case PR_NODE_DOF_DX_FOLLOW_FOLLOW_3:
                break;

            case PR_NODE_DOF_DY_1:
                ss.push(NT_NODE_DOF_DY_SET);
                ss.push(NT_BOOLEAN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_DY);
                break;

            case PR_NODE_DOF_DY_SET:
                temp_node_restrictions.setdy();
                break;

            case PR_NODE_DOF_DY_FOLLOW_1:
                ss.push(NT_NODE_DOF_DY_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_DOF_DY_FOLLOW_2:
                break;

            case PR_NODE_DOF_DY_FOLLOW_FOLLOW_1:
                ss.push(NT_NODE_DOF_DZ_FOLLOW);
                ss.push(NT_NODE_DOF_DZ);
                break;

            case PR_NODE_DOF_DY_FOLLOW_FOLLOW_2:
                break;

            case PR_NODE_DOF_DZ_1:
                ss.push(NT_NODE_DOF_DZ_SET);
                ss.push(NT_BOOLEAN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_DZ);
                break;

            case PR_NODE_DOF_DZ_SET:
                temp_node_restrictions.setdz();
                break;

            case PR_NODE_DOF_DZ_FOLLOW_1:
                break;

            case PR_NODE_RESTRICTION_SET:
                //TODO finish this
                for (std::vector<size_t>::iterator i = reference.begin(); i != reference.end(); i++) {
                    model.pushNodeRestrictions(*i, temp_node_restrictions);
                }
                reference.clear();
                temp_node_restrictions.reset();
                break;

            case PR_NODE_RESTRICTION_FOLLOW_1:
                ss.push(NT_NODE_RESTRICTION_FOLLOW);
                ss.push(NT_NODE_RESTRICTION);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODE_RESTRICTION_FOLLOW_2:
                break;

            case PR_LOAD_PATTERNS_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_LOAD_PATTERN_FOLLOW);
                ss.push(NT_LOAD_PATTERN);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_LOAD_PATTERNS);
                break;

            case PR_LOAD_PATTERN_1:
                ss.push(NT_LOAD_PATTERN_SET);
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_LOAD_PATTERN_LABEL_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_LABEL);
                ss.push(TT_OBJECT_OPEN);

                temp_load_pattern.clear();
                break;

            case PR_LOAD_PATTERN_SET:
                temp_load_pattern.setLabel(label);
                model.pushLoadPattern(temp_load_pattern);
                break;

            case PR_LOAD_PATTERN_LABEL_FOLLOW_1:
                ss.push(NT_NODAL_LOADS_LIST_FOLLOW);
                ss.push(NT_NODAL_LOADS_LIST);
                break;

            case PR_LOAD_PATTERN_LABEL_FOLLOW_2:
                ss.push(NT_SURFACE_LOADS_LIST_FOLLOW);
                ss.push(NT_SURFACE_LOADS_LIST);
                break;

            case PR_LOAD_PATTERN_LABEL_FOLLOW_3:
                ss.push(NT_DOMAIN_LOADS_LIST_FOLLOW);
                ss.push(NT_DOMAIN_LOADS_LIST);
                break;

            case PR_NODAL_LOADS_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_NODAL_LOAD_FOLLOW);
                ss.push(NT_NODAL_LOAD);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODAL_LOADS);
                break;

            case PR_NODAL_LOAD_1:
                ss.push(NT_NODAL_LOAD_SET);
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_NODAL_LOAD_FORCE);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_NODE_REFERENCE);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_NODAL_LOAD_SET:
                for (std::vector<size_t>::iterator i = reference.begin(); i != reference.end(); i++) {
                    temp_load_pattern.addNodalLoad(*i, temp_point.back());
                }
                reference.clear();
                temp_point.pop_back();
                break;

            case PR_NODAL_LOAD_FOLLOW_1:
                ss.push(NT_NODAL_LOAD_FOLLOW);
                ss.push(NT_NODAL_LOAD);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODAL_LOAD_FOLLOW_2:
                break;

            case PR_NODAL_LOAD_FORCE_1:
                ss.push(NT_VECTOR3D);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_NODAL_LOADS_LIST_FOLLOW_1:
                ss.push(NT_NODAL_LOADS_LIST_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_NODAL_LOADS_LIST_FOLLOW_2:
                break;

            case PR_SURFACE_LOADS_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_SURFACE_LOAD_FOLLOW);
                ss.push(NT_SURFACE_LOAD);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_SURFACE_LOADS);
                break;

            case PR_SURFACE_LOAD_1:
                ss.push(NT_SURFACE_LOAD_SET);
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_SURFACE_LOAD_TYPE);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_TYPE);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_SURFACE_LOAD_TYPE_1:
                ss.push(NT_SURFACE_TRIANGLE3_FORCES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_SURFACE_TRIANGLE3_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_TRIANGLE3);

                element.type = fem::Element::FE_TRIANGLE3;
                break;

            case PR_SURFACE_LOAD_TYPE_2:
                ss.push(NT_SURFACE_TRIANGLE6_FORCES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_SURFACE_TRIANGLE6_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_TRIANGLE6);
                element.type = fem::Element::FE_TRIANGLE6;
                break;

            case PR_SURFACE_LOAD_TYPE_3:
                ss.push(NT_SURFACE_QUADRANGLE4_FORCES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_SURFACE_QUADRANGLE4_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_QUADRANGLE4);
                element.type = fem::Element::FE_QUADRANGLE4;
                break;

            case PR_SURFACE_LOAD_TYPE_4:
                ss.push(NT_SURFACE_QUADRANGLE8_FORCES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_SURFACE_QUADRANGLE8_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_QUADRANGLE8);
                element.type = fem::Element::FE_QUADRANGLE8;
                break;

            case PR_SURFACE_LOAD_TYPE_5:
                ss.push(NT_SURFACE_QUADRANGLE9_FORCES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_SURFACE_QUADRANGLE9_NODES);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(TT_LABEL_QUADRANGLE9);
                element.type = fem::Element::FE_QUADRANGLE9;
                break;

            case PR_SURFACE_LOAD_SET: {
                element.nodes = reference;

                fem::SurfaceLoad* surface_load = fem::SurfaceLoad::makeSurfaceLoad(&element);
                surface_load->surface_forces = temp_point;

                temp_load_pattern.addSurfaceLoad(surface_load);
                reference.clear();
                temp_point.clear();
            } break;

            case PR_SURFACE_TRIANGLE3_NODES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 2; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_SURFACE_TRIANGLE6_NODES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 5; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_SURFACE_QUADRANGLE4_NODES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 3; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_SURFACE_QUADRANGLE8_NODES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 7; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_SURFACE_QUADRANGLE9_NODES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 8; i++) {
                    ss.push(NT_REFERENCE);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_REFERENCE);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_NODES);
                break;

            case PR_SURFACE_TRIANGLE3_FORCES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 2; i++) {
                    ss.push(NT_VECTOR3D);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_VECTOR3D);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_SURFACE_TRIANGLE6_FORCES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 5; i++) {
                    ss.push(NT_VECTOR3D);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_VECTOR3D);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_SURFACE_QUADRANGLE4_FORCES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 3; i++) {
                    ss.push(NT_VECTOR3D);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_VECTOR3D);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_SURFACE_QUADRANGLE8_FORCES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 7; i++) {
                    ss.push(NT_VECTOR3D);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_VECTOR3D);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_SURFACE_QUADRANGLE9_FORCES_1:
                ss.push(TT_ARRAY_CLOSE);
                for (int i = 0; i < 8; i++) {
                    ss.push(NT_VECTOR3D);
                    ss.push(TT_VALUE_SEPARATOR);
                }
                ss.push(NT_VECTOR3D);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_SURFACE_LOADS_LIST_FOLLOW_1:
                ss.push(NT_SURFACE_LOADS_LIST_FOLLOW_FOLLOW);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_SURFACE_LOADS_LIST_FOLLOW_2:
                break;

            case PR_SURFACE_LOADS_LIST_FOLLOW_FOLLOW_1:
                ss.push(NT_DOMAIN_LOADS_LIST_FOLLOW);
                ss.push(NT_DOMAIN_LOADS_LIST);
                break;

            case PR_SURFACE_LOAD_FOLLOW_1:
                ss.push(NT_SURFACE_LOAD_FOLLOW);
                ss.push(NT_SURFACE_LOAD);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_SURFACE_LOAD_FOLLOW_2:
                break;

            case PR_DOMAIN_LOADS_LIST_1:
                ss.push(TT_ARRAY_CLOSE);
                ss.push(NT_DOMAIN_LOAD_FOLLOW);
                ss.push(NT_DOMAIN_LOAD);
                ss.push(TT_ARRAY_OPEN);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_DOMAIN_LOADS);
                break;

            case PR_DOMAIN_LOADS_LIST_FOLLOW_1:
                break;

            case PR_DOMAIN_LOAD_1:
                ss.push(NT_DOMAIN_LOAD_SET);
                ss.push(TT_OBJECT_CLOSE);
                ss.push(NT_DOMAIN_FORCE);
                ss.push(TT_VALUE_SEPARATOR);
                ss.push(NT_ELEMENT_REFERENCE);
                ss.push(TT_OBJECT_OPEN);
                break;

            case PR_DOMAIN_LOAD_SET:
                for (std::vector<size_t>::iterator i = reference.begin(); i != reference.end(); i++) {
                    if (*i < model.numberOfElements()) {
                        temp_load_pattern.addDomainLoad(*i, temp_point.back());
                    } else {
                        /* this model tries to set a domain load on an element which doesn't exist */
                        std::cerr << "Error: domain load applied to non-existing element #" << *i << std::endl;
                    }
                }
                reference.clear();
                temp_point.pop_back();
                break;

            case PR_ELEMENT_REFERENCE_1:
                ss.push(NT_MULTI_REFERENCE);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_ELEMENT);
                break;

            case PR_DOMAIN_FORCE_1:
                ss.push(NT_VECTOR3D);
                ss.push(TT_NAME_SEPARATOR);
                ss.push(TT_LABEL_FORCE);
                break;

            case PR_DOMAIN_LOAD_FOLLOW_1:
                ss.push(NT_DOMAIN_LOAD_FOLLOW);
                ss.push(NT_DOMAIN_LOAD);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_DOMAIN_LOAD_FOLLOW_2:
                break;

            case PR_NODAL_LOADS_LIST_FOLLOW_FOLLOW_1:
                ss.push(NT_SURFACE_LOADS_LIST_FOLLOW);
                ss.push(NT_SURFACE_LOADS_LIST);
                break;

            case PR_NODAL_LOADS_LIST_FOLLOW_FOLLOW_2:
                ss.push(NT_DOMAIN_LOADS_LIST_FOLLOW);
                ss.push(NT_DOMAIN_LOADS_LIST);
                break;

            case PR_LOAD_PATTERN_FOLLOW_1:
                ss.push(NT_LOAD_PATTERN_FOLLOW);
                ss.push(NT_LOAD_PATTERN);
                ss.push(TT_VALUE_SEPARATOR);
                break;

            case PR_LOAD_PATTERN_FOLLOW_2:
                break;

            case PR_LOAD_PATTERNS_LIST_FOLLOW_1:
                break;

            default:
                this->error.code = Error::Type::ERR_UNKNOWN;
                this->error.message = "Unknown error while parsing line " + boost::lexical_cast<std::string>(line_number);
                return this->error.code;
                break;
            }
        }
    }

    // all is good
    this->error.code = Error::Type::ERR_OK;
    this->error.message.clear();
    return this->error.code;
}

enum FemJsonParser::TerminalTokens
FemJsonParser::lexer(std::istream& file)
{
#define YYFILL fill(file);
starting_point:
    tok = pos;
    /*!re2c
        re2c:define:YYCTYPE  	= "unsigned char";
        re2c:define:YYCURSOR 	= pos;
        re2c:define:YYLIMIT 	= limit;
        re2c:define:YYMARKER 	= marker;
        re2c:yyfill:enable	= 1;
        re2c:yych:conversion	= 0;
        re2c:indent:top		= 1;

        natural = "0"|[1-9][0-9]*;

        "\"label\""		{ return lexerReturnProcedures(TT_LABEL_LABEL);	}
        "\"fem\""		{ return lexerReturnProcedures(TT_LABEL_FEM);	}
        "\"version\""		{ return lexerReturnProcedures(TT_LABEL_VERSION);	}
        "\"1.0\""		{ return lexerReturnProcedures(TT_LABEL_V1_0);	}
        "\"type\""		{ return lexerReturnProcedures(TT_LABEL_TYPE);	}
        "\"3D solid\""		{ return lexerReturnProcedures(TT_LABEL_3D_SOLID);	}
        "\"materials\""		{ return lexerReturnProcedures(TT_LABEL_MATERIALS);	}
        "\"E\""			{ return lexerReturnProcedures(TT_LABEL_E);	}
        "\"nu\""		{ return lexerReturnProcedures(TT_LABEL_NU);	}
        "\"fy\""		{ return lexerReturnProcedures(TT_LABEL_FY);	}
        "\"linear elastic\""	{ return lexerReturnProcedures(TT_LABEL_LINEAR_ELASTIC);	}
        "\"nodes\""		{ return lexerReturnProcedures(TT_LABEL_NODES);	}
        "\"elements\""		{ return lexerReturnProcedures(TT_LABEL_ELEMENTS);	}
        "\"tetrahedron4\""	{ return lexerReturnProcedures(TT_LABEL_TETRAHEDRON4);	}
        "\"tetrahedron10\""	{ return lexerReturnProcedures(TT_LABEL_TETRAHEDRON10);	}
        "\"hexahedron8\""	{ return lexerReturnProcedures(TT_LABEL_HEXAHEDRON8);	}
        "\"hexahedron20\""	{ return lexerReturnProcedures(TT_LABEL_HEXAHEDRON20);	}
        "\"hexahedron27\""	{ return lexerReturnProcedures(TT_LABEL_HEXAHEDRON27);	}
        "\"prism6\""		{ return lexerReturnProcedures(TT_LABEL_PRISM6);	}
        "\"prism15\""		{ return lexerReturnProcedures(TT_LABEL_PRISM15);	}
        "\"prism18\""		{ return lexerReturnProcedures(TT_LABEL_PRISM18);	}
        "\"material\""		{ return lexerReturnProcedures(TT_LABEL_MATERIAL);	}
        "\"groups\""		{ return lexerReturnProcedures(TT_LABEL_GROUPS);	}
        "\"material\""		{ return lexerReturnProcedures(TT_LABEL_MATERIAL);	}
        "\"node restrictions\""	{ return lexerReturnProcedures(TT_LABEL_NODE_RESTRICTIONS);	}
        "\"load patterns\""	{ return lexerReturnProcedures(TT_LABEL_LOAD_PATTERNS);	}
        "\"node\""		{ return lexerReturnProcedures(TT_LABEL_NODE);	}
        "\"dx\""		{ return lexerReturnProcedures(TT_LABEL_DX);	}
        "\"dy\""		{ return lexerReturnProcedures(TT_LABEL_DY);	}
        "\"dz\""		{ return lexerReturnProcedures(TT_LABEL_DZ);	}
        "\"nodal loads\""	{ return lexerReturnProcedures(TT_LABEL_NODAL_LOADS);	}
        "\"surface loads\""	{ return lexerReturnProcedures(TT_LABEL_SURFACE_LOADS);	}
        "\"domain loads\""	{ return lexerReturnProcedures(TT_LABEL_DOMAIN_LOADS);	}
        "\"triangle3\""		{ return lexerReturnProcedures(TT_LABEL_TRIANGLE3);	}
        "\"triangle6\""		{ return lexerReturnProcedures(TT_LABEL_TRIANGLE6);	}
        "\"quadrangle4\""	{ return lexerReturnProcedures(TT_LABEL_QUADRANGLE4);	}
        "\"quadrangle8\""	{ return lexerReturnProcedures(TT_LABEL_QUADRANGLE8);	}
        "\"quadrangle9\""	{ return lexerReturnProcedures(TT_LABEL_QUADRANGLE9);	}
        "\"element\""		{ return lexerReturnProcedures(TT_LABEL_ELEMENT);	}
        "\"force\""		{ return lexerReturnProcedures(TT_LABEL_FORCE);		}
        "\"forces\""		{ return lexerReturnProcedures(TT_LABEL_FORCE);		}

        natural			{ return lexerReturnProcedures(TT_NATURAL_NUMBER);	}
        [-+]?natural ("." [0-9]+)? ([Ee] [+-]?[0-9]+)?	{ return lexerReturnProcedures(TT_FLOAT_NUMBER);	}
        "-"natural		{ return lexerReturnProcedures(TT_FLOAT_NUMBER);	}
        "\"" [^\"\n]* "\""	{ return lexerReturnProcedures(TT_TEXT_STRING);	}
        "{"			{ return lexerReturnProcedures(TT_OBJECT_OPEN); 	}
        "}"			{ return lexerReturnProcedures(TT_OBJECT_CLOSE); 	}
        "["			{ return lexerReturnProcedures(TT_ARRAY_OPEN); 	}
        "]"			{ return lexerReturnProcedures(TT_ARRAY_CLOSE); 	}
        ":"			{ return lexerReturnProcedures(TT_NAME_SEPARATOR);	}
        ","			{ return lexerReturnProcedures(TT_VALUE_SEPARATOR);	}
        "true"			{ return lexerReturnProcedures(TT_TRUE);	}
        "false"			{ return lexerReturnProcedures(TT_FALSE);	}
        "null"			{ return lexerReturnProcedures(TT_NULL);	}

        "\000"			{ return lexerReturnProcedures(TT_EOF);		}

        [ \t]+			{ goto starting_point;				}
        "\n"			{ line_number++; goto starting_point;		}
        [^]			{ return lexerReturnProcedures(TT_UNKNOWN); 		}
*/
#undef YYFILL
}

void FemJsonParser::setParserTable()
{
    table[NT_LABEL][TT_LABEL_LABEL] = PR_LABEL;
    table[NT_LABEL_SET][TT_VALUE_SEPARATOR] = PR_LABEL_SET;
    table[NT_LABEL_SET][TT_OBJECT_CLOSE] = PR_LABEL_SET;

    table[NT_VECTOR3D][TT_ARRAY_OPEN] = PR_VECTOR3D;
    table[NT_VECTOR3D_SET][TT_VALUE_SEPARATOR] = PR_VECTOR3D_SET;
    table[NT_VECTOR3D_SET][TT_ARRAY_CLOSE] = PR_VECTOR3D_SET;
    table[NT_VECTOR3D_SET][TT_OBJECT_CLOSE] = PR_VECTOR3D_SET;

    table[NT_REFERENCE][TT_NATURAL_NUMBER] = PR_REFERENCE_1;
    //	table	[NT_REFERENCE]		[TT_VALUE_SEPARATOR]	=	PR_REFERENCE_ERROR;	// error state

    table[NT_REFERENCE_SET][TT_VALUE_SEPARATOR] = PR_REFERENCE_SET;
    table[NT_REFERENCE_SET][TT_ARRAY_CLOSE] = PR_REFERENCE_SET;
    table[NT_REFERENCE_SET][TT_OBJECT_CLOSE] = PR_REFERENCE_SET;

    table[NT_REFERENCE_LIST][TT_ARRAY_OPEN] = PR_REFERENCE_LIST_1;

    table[NT_REFERENCE_FOLLOW][TT_VALUE_SEPARATOR] = PR_REFERENCE_FOLLOW_1;
    table[NT_REFERENCE_FOLLOW][TT_ARRAY_CLOSE] = PR_REFERENCE_FOLLOW_2;

    table[NT_MULTI_REFERENCE][TT_NATURAL_NUMBER] = PR_MULTI_REFERENCE_1;
    table[NT_MULTI_REFERENCE][TT_ARRAY_OPEN] = PR_MULTI_REFERENCE_2;

    table[NT_MULTI_REFERENCE_SET][TT_VALUE_SEPARATOR] = PR_MULTI_REFERENCE_SET;
    table[NT_MULTI_REFERENCE_SET][TT_ARRAY_CLOSE] = PR_MULTI_REFERENCE_SET;
    table[NT_MULTI_REFERENCE_SET][TT_OBJECT_CLOSE] = PR_MULTI_REFERENCE_SET;

    table[NT_BOOLEAN][TT_TRUE] = PR_BOOLEAN_TRUE;
    table[NT_BOOLEAN][TT_FALSE] = PR_BOOLEAN_FALSE;

    table[NT_NUMBER][TT_FLOAT_NUMBER] = PR_NUMBER_1;
    table[NT_NUMBER][TT_NATURAL_NUMBER] = PR_NUMBER_2;

    table[NT_FLOAT_SET][TT_VALUE_SEPARATOR] = PR_FLOAT_SET;
    table[NT_FLOAT_SET][TT_ARRAY_CLOSE] = PR_FLOAT_SET;
    table[NT_FLOAT_SET][TT_OBJECT_CLOSE] = PR_FLOAT_SET;

    table[NT_INT_TO_FLOAT_SET][TT_VALUE_SEPARATOR] = PR_INT_TO_FLOAT_SET;
    table[NT_INT_TO_FLOAT_SET][TT_ARRAY_CLOSE] = PR_INT_TO_FLOAT_SET;
    table[NT_INT_TO_FLOAT_SET][TT_OBJECT_CLOSE] = PR_INT_TO_FLOAT_SET;

    table[NT_START][TT_OBJECT_OPEN] = PR_START_1;
    table[NT_START][TT_EOF] = PR_START_2;

    table[NT_DOCUMENT][TT_OBJECT_OPEN] = PR_DOCUMENT_1;

    table[NT_HEADER][TT_LABEL_FEM] = PR_HEADER_1;
    table[NT_FORMAT_VERSION][TT_LABEL_VERSION] = PR_FORMAT_VERSION_1;

    table[NT_MODEL_TYPE][TT_LABEL_TYPE] = PR_MODEL_TYPE_1;
    table[NT_MODEL_TYPE_SET][TT_OBJECT_CLOSE] = PR_MODEL_TYPE_SET;

    table[NT_MATERIAL_LIST][TT_LABEL_MATERIALS] = PR_MATERIAL_LIST_1;
    table[NT_MATERIAL][TT_OBJECT_OPEN] = PR_MATERIAL_1;
    table[NT_MATERIAL_TYPE][TT_LABEL_TYPE] = PR_MATERIAL_TYPE_1;
    table[NT_YOUNG_MODULUS][TT_LABEL_E] = PR_YOUNG_MODULUS_1;
    table[NT_YOUNG_MODULUS_SET][TT_VALUE_SEPARATOR] = PR_YOUNG_MODULUS_SET;
    table[NT_POISSON_COEFFICIENT][TT_LABEL_NU] = PR_POISSON_COEFFICIENT_1;
    table[NT_POISSON_COEFFICIENT_SET][TT_VALUE_SEPARATOR] = PR_POISSON_COEFFICIENT_SET;
    table[NT_ELASTIC_LIMIT][TT_LABEL_FY] = PR_ELASTIC_LIMIT_1;
    table[NT_ELASTIC_LIMIT_SET][TT_OBJECT_CLOSE] = PR_ELASTIC_LIMIT_SET;
    table[NT_MATERIAL_SET][TT_VALUE_SEPARATOR] = PR_MATERIAL_SET;
    table[NT_MATERIAL_SET][TT_ARRAY_CLOSE] = PR_MATERIAL_SET;
    table[NT_MATERIAL_FOLLOW][TT_VALUE_SEPARATOR] = PR_MATERIAL_FOLLOW_1;
    table[NT_MATERIAL_FOLLOW][TT_ARRAY_CLOSE] = PR_MATERIAL_FOLLOW_2;

    table[NT_NODE_LIST][TT_LABEL_NODES] = PR_NODE_LIST_1;
    table[NT_NODE][TT_ARRAY_OPEN] = PR_NODE_1;
    table[NT_NODE_SET][TT_VALUE_SEPARATOR] = PR_NODE_SET;
    table[NT_NODE_SET][TT_ARRAY_CLOSE] = PR_NODE_SET;
    table[NT_NODE_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_FOLLOW_1;
    table[NT_NODE_FOLLOW][TT_ARRAY_CLOSE] = PR_NODE_FOLLOW_2;

    table[NT_ELEMENT_LIST][TT_LABEL_ELEMENTS] = PR_ELEMENT_LIST_1;
    table[NT_ELEMENT][TT_OBJECT_OPEN] = PR_ELEMENT_1;
    table[NT_ELEMENT_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_SET_1;

    table[NT_ELEMENT_TYPE][TT_LABEL_HEXAHEDRON8] = PR_ELEMENT_HEXAHEDRON8_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_HEXAHEDRON20] = PR_ELEMENT_HEXAHEDRON20_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_HEXAHEDRON27] = PR_ELEMENT_HEXAHEDRON27_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_TETRAHEDRON4] = PR_ELEMENT_TETRAHEDRON4_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_TETRAHEDRON10] = PR_ELEMENT_TETRAHEDRON10_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_PRISM6] = PR_ELEMENT_PRISM6_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_PRISM15] = PR_ELEMENT_PRISM15_1;
    table[NT_ELEMENT_TYPE][TT_LABEL_PRISM18] = PR_ELEMENT_PRISM18_1;

    table[NT_ELEMENT_HEXAHEDRON8][TT_LABEL_HEXAHEDRON8] = PR_ELEMENT_HEXAHEDRON8_1;
    table[NT_ELEMENT_HEXAHEDRON20][TT_LABEL_HEXAHEDRON20] = PR_ELEMENT_HEXAHEDRON20_1;
    table[NT_ELEMENT_HEXAHEDRON27][TT_LABEL_HEXAHEDRON27] = PR_ELEMENT_HEXAHEDRON27_1;
    table[NT_ELEMENT_TETRAHEDRON4][TT_LABEL_TETRAHEDRON4] = PR_ELEMENT_TETRAHEDRON4_1;
    table[NT_ELEMENT_TETRAHEDRON10][TT_LABEL_TETRAHEDRON10] = PR_ELEMENT_TETRAHEDRON10_1;
    table[NT_ELEMENT_PRISM6][TT_LABEL_PRISM6] = PR_ELEMENT_PRISM6_1;
    table[NT_ELEMENT_PRISM15][TT_LABEL_PRISM15] = PR_ELEMENT_PRISM15_1;
    table[NT_ELEMENT_PRISM18][TT_LABEL_PRISM18] = PR_ELEMENT_PRISM18_1;

    table[NT_ELEMENT_HEXAHEDRON8_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_HEXAHEDRON8_SET;
    table[NT_ELEMENT_HEXAHEDRON8_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_HEXAHEDRON8_SET;
    table[NT_ELEMENT_HEXAHEDRON20_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_HEXAHEDRON20_SET;
    table[NT_ELEMENT_HEXAHEDRON20_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_HEXAHEDRON20_SET;
    table[NT_ELEMENT_HEXAHEDRON27_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_HEXAHEDRON27_SET;
    table[NT_ELEMENT_HEXAHEDRON27_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_HEXAHEDRON27_SET;
    table[NT_ELEMENT_TETRAHEDRON4_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_TETRAHEDRON4_SET;
    table[NT_ELEMENT_TETRAHEDRON4_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_TETRAHEDRON4_SET;
    table[NT_ELEMENT_TETRAHEDRON10_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_TETRAHEDRON10_SET;
    table[NT_ELEMENT_TETRAHEDRON10_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_TETRAHEDRON10_SET;
    table[NT_ELEMENT_PRISM6_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_PRISM6_SET;
    table[NT_ELEMENT_PRISM6_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_PRISM6_SET;
    table[NT_ELEMENT_PRISM15_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_PRISM15_SET;
    table[NT_ELEMENT_PRISM15_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_PRISM15_SET;
    table[NT_ELEMENT_PRISM18_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_PRISM18_SET;
    table[NT_ELEMENT_PRISM18_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_PRISM18_SET;

    table[NT_ELEMENT_TYPE_FOLLOW][TT_VALUE_SEPARATOR] = PR_ELEMENT_TYPE_FOLLOW_1;
    table[NT_ELEMENT_TYPE_FOLLOW][TT_OBJECT_CLOSE] = PR_ELEMENT_TYPE_FOLLOW_2;
    table[NT_ELEMENT_MATERIAL][TT_LABEL_MATERIAL] = PR_ELEMENT_MATERIAL_1;
    table[NT_ELEMENT_MATERIAL_SET][TT_OBJECT_CLOSE] = PR_ELEMENT_MATERIAL_SET_1;
    table[NT_ELEMENT_FOLLOW][TT_VALUE_SEPARATOR] = PR_ELEMENT_FOLLOW_1;
    table[NT_ELEMENT_FOLLOW][TT_ARRAY_CLOSE] = PR_ELEMENT_FOLLOW_2;

    table[NT_ELEMENT_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_ELEMENT_LIST_FOLLOW_1;
    table[NT_ELEMENT_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_ELEMENT_LIST_FOLLOW_2;

    table[NT_EXTRA_DEFINITIONS][TT_LABEL_GROUPS] = PR_EXTRA_DEFINITIONS_1;
    table[NT_EXTRA_DEFINITIONS][TT_LABEL_NODE_RESTRICTIONS] = PR_EXTRA_DEFINITIONS_2;
    table[NT_EXTRA_DEFINITIONS][TT_LABEL_LOAD_PATTERNS] = PR_EXTRA_DEFINITIONS_3;

    table[NT_GROUPS_LIST][TT_LABEL_GROUPS] = PR_GROUPS_LIST_1;
    table[NT_GROUP_TYPE][TT_LABEL_NODES] = PR_GROUP_TYPE_1;
    table[NT_GROUP_TYPE][TT_LABEL_ELEMENTS] = PR_GROUP_TYPE_2;

    table[NT_GROUPS_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_GROUPS_LIST_FOLLOW_1;
    table[NT_GROUPS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_GROUPS_LIST_FOLLOW_2;

    table[NT_GROUPS_LIST_FOLLOW_FOLLOW][TT_LABEL_NODE_RESTRICTIONS] = PR_GROUPS_LIST_FOLLOW_FOLLOW_1;
    table[NT_GROUPS_LIST_FOLLOW_FOLLOW][TT_LABEL_LOAD_PATTERNS] = PR_GROUPS_LIST_FOLLOW_FOLLOW_2;

    table[NT_NODE_GROUP_LIST][TT_LABEL_NODES] = PR_NODE_GROUP_LIST_1;
    table[NT_NODE_GROUP_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_GROUP_LIST_FOLLOW_1;
    table[NT_NODE_GROUP_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_GROUP_LIST_FOLLOW_2;

    table[NT_NODE_GROUP][TT_OBJECT_OPEN] = PR_NODE_GROUP_1;

    table[NT_NODE_GROUP_SET][TT_VALUE_SEPARATOR] = PR_NODE_GROUP_SET;
    table[NT_NODE_GROUP_SET][TT_ARRAY_CLOSE] = PR_NODE_GROUP_SET;

    table[NT_NODE_GROUP_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_GROUP_FOLLOW_1;
    table[NT_NODE_GROUP_FOLLOW][TT_ARRAY_CLOSE] = PR_NODE_GROUP_FOLLOW_2;

    table[NT_NODE_GROUP_FOLLOW_FOLLOW][TT_OBJECT_OPEN] = PR_NODE_GROUP_FOLLOW_FOLLOW_1;

    table[NT_ELEMENT_GROUP_LIST][TT_LABEL_ELEMENTS] = PR_ELEMENT_GROUP_LIST_1;
    table[NT_ELEMENT_GROUP_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_ELEMENT_GROUP_LIST_FOLLOW_1;
    table[NT_ELEMENT_GROUP][TT_OBJECT_OPEN] = PR_ELEMENT_GROUP_1;
    table[NT_ELEMENT_GROUP_SET][TT_VALUE_SEPARATOR] = PR_ELEMENT_GROUP_SET;
    table[NT_ELEMENT_GROUP_SET][TT_ARRAY_CLOSE] = PR_ELEMENT_GROUP_SET;
    table[NT_ELEMENT_GROUP_FOLLOW][TT_VALUE_SEPARATOR] = PR_ELEMENT_GROUP_FOLLOW_1;
    table[NT_ELEMENT_GROUP_FOLLOW][TT_ARRAY_CLOSE] = PR_ELEMENT_GROUP_FOLLOW_2;

    table[NT_ELEMENT_GROUP_FOLLOW_FOLLOW][TT_OBJECT_OPEN] = PR_ELEMENT_GROUP_FOLLOW_FOLLOW_1;

    table[NT_NODE_RESTRICTIONS_LIST][TT_LABEL_NODE_RESTRICTIONS] = PR_NODE_RESTRICTIONS_LIST_1;

    table[NT_NODE_RESTRICTIONS_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_RESTRICTIONS_LIST_FOLLOW_1;
    table[NT_NODE_RESTRICTIONS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_RESTRICTIONS_LIST_FOLLOW_2;
    table[NT_NODE_RESTRICTIONS_LIST_FOLLOW_FOLLOW][TT_LABEL_LOAD_PATTERNS] = PR_NODE_RESTRICTIONS_LIST_FOLLOW_FOLLOW_1;

    table[NT_NODE_RESTRICTION][TT_OBJECT_OPEN] = PR_NODE_RESTRICTION_1;

    table[NT_NODE_REFERENCE][TT_LABEL_NODE] = PR_NODE_REFERENCE_1;
    table[NT_NODE_REFERENCE_SET][TT_VALUE_SEPARATOR] = PR_NODE_REFERENCE_SET;

    table[NT_NODE_DOF][TT_LABEL_DX] = PR_NODE_DOF_1;
    table[NT_NODE_DOF][TT_LABEL_DY] = PR_NODE_DOF_2;
    table[NT_NODE_DOF][TT_LABEL_DZ] = PR_NODE_DOF_3;
    table[NT_NODE_DOF_DX][TT_LABEL_DX] = PR_NODE_DOF_DX_1;

    table[NT_NODE_DOF_DX_SET][TT_VALUE_SEPARATOR] = PR_NODE_DOF_DX_SET;
    table[NT_NODE_DOF_DX_SET][TT_OBJECT_CLOSE] = PR_NODE_DOF_DX_SET;
    table[NT_NODE_DOF_DX_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_DOF_DX_FOLLOW_1;
    table[NT_NODE_DOF_DX_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DX_FOLLOW_2;

    table[NT_NODE_DOF_DX_FOLLOW_FOLLOW][TT_LABEL_DY] = PR_NODE_DOF_DX_FOLLOW_FOLLOW_1;
    table[NT_NODE_DOF_DX_FOLLOW_FOLLOW][TT_LABEL_DZ] = PR_NODE_DOF_DX_FOLLOW_FOLLOW_2;
    table[NT_NODE_DOF_DX_FOLLOW_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DX_FOLLOW_FOLLOW_3;

    table[NT_NODE_DOF_DY][TT_LABEL_DY] = PR_NODE_DOF_DY_1;
    table[NT_NODE_DOF_DY_SET][TT_VALUE_SEPARATOR] = PR_NODE_DOF_DY_SET;
    table[NT_NODE_DOF_DY_SET][TT_OBJECT_CLOSE] = PR_NODE_DOF_DY_SET;
    table[NT_NODE_DOF_DY_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_DOF_DY_FOLLOW_1;
    table[NT_NODE_DOF_DY_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DY_FOLLOW_2;
    table[NT_NODE_DOF_DY_FOLLOW_FOLLOW][TT_LABEL_DZ] = PR_NODE_DOF_DY_FOLLOW_FOLLOW_1;
    table[NT_NODE_DOF_DY_FOLLOW_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DY_FOLLOW_FOLLOW_2;

    table[NT_NODE_DOF_DZ][TT_LABEL_DZ] = PR_NODE_DOF_DZ_1;
    table[NT_NODE_DOF_DZ_SET][TT_OBJECT_CLOSE] = PR_NODE_DOF_DZ_SET;
    table[NT_NODE_DOF_DZ_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DZ_FOLLOW_1;
    table[NT_NODE_DOF_DZ_FOLLOW_FOLLOW][TT_OBJECT_CLOSE] = PR_NODE_DOF_DZ_FOLLOW_FOLLOW_1;

    table[NT_NODE_RESTRICTION_SET][TT_VALUE_SEPARATOR] = PR_NODE_RESTRICTION_SET;
    table[NT_NODE_RESTRICTION_SET][TT_ARRAY_CLOSE] = PR_NODE_RESTRICTION_SET;
    table[NT_NODE_RESTRICTION_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODE_RESTRICTION_FOLLOW_1;
    table[NT_NODE_RESTRICTION_FOLLOW][TT_ARRAY_CLOSE] = PR_NODE_RESTRICTION_FOLLOW_2;

    table[NT_LOAD_PATTERNS_LIST][TT_LABEL_LOAD_PATTERNS] = PR_LOAD_PATTERNS_LIST_1;
    table[NT_LOAD_PATTERN][TT_OBJECT_OPEN] = PR_LOAD_PATTERN_1;

    table[NT_LOAD_PATTERN_SET][TT_VALUE_SEPARATOR] = PR_LOAD_PATTERN_SET;
    table[NT_LOAD_PATTERN_SET][TT_ARRAY_CLOSE] = PR_LOAD_PATTERN_SET;

    table[NT_LOAD_PATTERN_LABEL_FOLLOW][TT_LABEL_NODAL_LOADS] = PR_LOAD_PATTERN_LABEL_FOLLOW_1;
    table[NT_LOAD_PATTERN_LABEL_FOLLOW][TT_LABEL_SURFACE_LOADS] = PR_LOAD_PATTERN_LABEL_FOLLOW_2;
    table[NT_LOAD_PATTERN_LABEL_FOLLOW][TT_LABEL_DOMAIN_LOADS] = PR_LOAD_PATTERN_LABEL_FOLLOW_3;
    table[NT_NODAL_LOADS_LIST][TT_LABEL_NODAL_LOADS] = PR_NODAL_LOADS_LIST_1;
    table[NT_NODAL_LOAD][TT_OBJECT_OPEN] = PR_NODAL_LOAD_1;
    table[NT_NODAL_LOAD_SET][TT_VALUE_SEPARATOR] = PR_NODAL_LOAD_SET;
    table[NT_NODAL_LOAD_SET][TT_ARRAY_CLOSE] = PR_NODAL_LOAD_SET;

    table[NT_NODAL_LOAD_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODAL_LOAD_FOLLOW_1;
    table[NT_NODAL_LOAD_FOLLOW][TT_ARRAY_CLOSE] = PR_NODAL_LOAD_FOLLOW_2;
    table[NT_NODAL_LOAD_FORCE][TT_LABEL_FORCE] = PR_NODAL_LOAD_FORCE_1;
    table[NT_NODAL_LOADS_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_NODAL_LOADS_LIST_FOLLOW_1;
    table[NT_NODAL_LOADS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_NODAL_LOADS_LIST_FOLLOW_2;
    table[NT_SURFACE_LOADS_LIST][TT_LABEL_SURFACE_LOADS] = PR_SURFACE_LOADS_LIST_1;
    table[NT_SURFACE_LOAD][TT_OBJECT_OPEN] = PR_SURFACE_LOAD_1;
    table[NT_SURFACE_LOAD_TYPE][TT_LABEL_TRIANGLE3] = PR_SURFACE_LOAD_TYPE_1;
    table[NT_SURFACE_LOAD_TYPE][TT_LABEL_TRIANGLE6] = PR_SURFACE_LOAD_TYPE_2;
    table[NT_SURFACE_LOAD_TYPE][TT_LABEL_QUADRANGLE4] = PR_SURFACE_LOAD_TYPE_3;
    table[NT_SURFACE_LOAD_TYPE][TT_LABEL_QUADRANGLE8] = PR_SURFACE_LOAD_TYPE_4;
    table[NT_SURFACE_LOAD_TYPE][TT_LABEL_QUADRANGLE9] = PR_SURFACE_LOAD_TYPE_5;
    table[NT_SURFACE_LOAD_SET][TT_VALUE_SEPARATOR] = PR_SURFACE_LOAD_SET;
    table[NT_SURFACE_LOAD_SET][TT_ARRAY_CLOSE] = PR_SURFACE_LOAD_SET;

    table[NT_SURFACE_TRIANGLE3_NODES][TT_LABEL_NODES] = PR_SURFACE_TRIANGLE3_NODES_1;
    table[NT_SURFACE_TRIANGLE6_NODES][TT_LABEL_NODES] = PR_SURFACE_TRIANGLE6_NODES_1;
    table[NT_SURFACE_QUADRANGLE4_NODES][TT_LABEL_NODES] = PR_SURFACE_QUADRANGLE4_NODES_1;
    table[NT_SURFACE_QUADRANGLE8_NODES][TT_LABEL_NODES] = PR_SURFACE_QUADRANGLE8_NODES_1;
    table[NT_SURFACE_QUADRANGLE9_NODES][TT_LABEL_NODES] = PR_SURFACE_QUADRANGLE9_NODES_1;

    table[NT_SURFACE_TRIANGLE3_FORCES][TT_LABEL_FORCE] = PR_SURFACE_TRIANGLE3_FORCES_1;
    table[NT_SURFACE_TRIANGLE6_FORCES][TT_LABEL_FORCE] = PR_SURFACE_TRIANGLE6_FORCES_1;
    table[NT_SURFACE_QUADRANGLE4_FORCES][TT_LABEL_FORCE] = PR_SURFACE_QUADRANGLE4_FORCES_1;
    table[NT_SURFACE_QUADRANGLE8_FORCES][TT_LABEL_FORCE] = PR_SURFACE_QUADRANGLE8_FORCES_1;
    table[NT_SURFACE_QUADRANGLE9_FORCES][TT_LABEL_FORCE] = PR_SURFACE_QUADRANGLE9_FORCES_1;

    table[NT_SURFACE_LOADS_LIST_FOLLOW][TT_VALUE_SEPARATOR] = PR_SURFACE_LOADS_LIST_FOLLOW_1;
    table[NT_SURFACE_LOADS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_SURFACE_LOADS_LIST_FOLLOW_2;

    table[NT_SURFACE_LOADS_LIST_FOLLOW_FOLLOW][TT_LABEL_DOMAIN_LOADS] = PR_SURFACE_LOADS_LIST_FOLLOW_FOLLOW_1;

    table[NT_SURFACE_LOAD_FOLLOW][TT_VALUE_SEPARATOR] = PR_SURFACE_LOAD_FOLLOW_1;
    table[NT_SURFACE_LOAD_FOLLOW][TT_ARRAY_CLOSE] = PR_SURFACE_LOAD_FOLLOW_2;
    table[NT_DOMAIN_LOADS_LIST][TT_LABEL_DOMAIN_LOADS] = PR_DOMAIN_LOADS_LIST_1;
    table[NT_DOMAIN_LOADS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_DOMAIN_LOADS_LIST_FOLLOW_1;
    table[NT_DOMAIN_LOAD][TT_OBJECT_OPEN] = PR_DOMAIN_LOAD_1;
    table[NT_DOMAIN_LOAD_SET][TT_VALUE_SEPARATOR] = PR_DOMAIN_LOAD_SET;
    table[NT_DOMAIN_LOAD_SET][TT_ARRAY_CLOSE] = PR_DOMAIN_LOAD_SET;
    table[NT_ELEMENT_REFERENCE][TT_LABEL_ELEMENT] = PR_ELEMENT_REFERENCE_1;
    table[NT_DOMAIN_FORCE][TT_LABEL_FORCE] = PR_DOMAIN_FORCE_1;
    table[NT_DOMAIN_LOAD_FOLLOW][TT_VALUE_SEPARATOR] = PR_DOMAIN_LOAD_FOLLOW_1;
    table[NT_DOMAIN_LOAD_FOLLOW][TT_ARRAY_CLOSE] = PR_DOMAIN_LOAD_FOLLOW_2;

    table[NT_NODAL_LOADS_LIST_FOLLOW_FOLLOW][TT_LABEL_SURFACE_LOADS] = PR_NODAL_LOADS_LIST_FOLLOW_FOLLOW_1;
    table[NT_NODAL_LOADS_LIST_FOLLOW_FOLLOW][TT_LABEL_DOMAIN_LOADS] = PR_NODAL_LOADS_LIST_FOLLOW_FOLLOW_2;
    table[NT_LOAD_PATTERN_FOLLOW][TT_VALUE_SEPARATOR] = PR_LOAD_PATTERN_FOLLOW_1;
    table[NT_LOAD_PATTERN_FOLLOW][TT_ARRAY_CLOSE] = PR_LOAD_PATTERN_FOLLOW_2;
    table[NT_LOAD_PATTERNS_LIST_FOLLOW][TT_OBJECT_CLOSE] = PR_LOAD_PATTERNS_LIST_FOLLOW_1;
}

enum FemJsonParser::TerminalTokens FemJsonParser::lexerReturnProcedures(TerminalTokens tt)
{
    switch (tt) {
    case TT_TEXT_STRING:
        temp_token.assign(tok + 1, pos - tok - 2);
        return TT_TEXT_STRING;
        break;

    case TT_NATURAL_NUMBER:
        temp_token.assign(tok, pos - tok);
        return TT_NATURAL_NUMBER;
        break;

    case TT_FLOAT_NUMBER:
        temp_token.assign(tok, pos - tok);
        return TT_FLOAT_NUMBER;
        break;

    case TT_OBJECT_OPEN:
        return TT_OBJECT_OPEN;
        break;
    case TT_OBJECT_CLOSE:
        return TT_OBJECT_CLOSE;
        break;
    case TT_ARRAY_OPEN:
        return TT_ARRAY_OPEN;
        break;
    case TT_ARRAY_CLOSE:
        return TT_ARRAY_CLOSE;
        break;
    case TT_NAME_SEPARATOR:
        return TT_NAME_SEPARATOR;
        break;
    case TT_VALUE_SEPARATOR:
        return TT_VALUE_SEPARATOR;
        break;
    case TT_TRUE:
        temp_bool = true;
        return TT_TRUE;
        break;
    case TT_FALSE:
        temp_bool = false;
        return TT_FALSE;
        break;
    case TT_NULL:
        return TT_NULL;
        break;

    case TT_LABEL_LABEL:
        return TT_LABEL_LABEL;
        break;

    case TT_LABEL_FEM:
        return TT_LABEL_FEM;
        break;

    case TT_LABEL_VERSION:
        return TT_LABEL_VERSION;
        break;
    case TT_LABEL_V1_0:
        return TT_LABEL_V1_0;
        break;
    case TT_LABEL_TYPE:
        return TT_LABEL_TYPE;
        break;
    case TT_LABEL_3D_SOLID:
        return TT_LABEL_3D_SOLID;
        break;
    case TT_LABEL_MATERIALS:
        return TT_LABEL_MATERIALS;
        break;
    case TT_LABEL_E:
        return TT_LABEL_E;
        break;
    case TT_LABEL_NU:
        return TT_LABEL_NU;
        break;
    case TT_LABEL_FY:
        return TT_LABEL_FY;
        break;
    case TT_LABEL_LINEAR_ELASTIC:
        return TT_LABEL_LINEAR_ELASTIC;
        break;
    case TT_LABEL_NODES:
        return TT_LABEL_NODES;
        break;
    case TT_LABEL_ELEMENTS:
        return TT_LABEL_ELEMENTS;
        break;
    case TT_LABEL_TETRAHEDRON4:
        return TT_LABEL_TETRAHEDRON4;
        break;
    case TT_LABEL_TETRAHEDRON10:
        return TT_LABEL_TETRAHEDRON10;
        break;
    case TT_LABEL_HEXAHEDRON8:
        return TT_LABEL_HEXAHEDRON8;
        break;
    case TT_LABEL_HEXAHEDRON20:
        return TT_LABEL_HEXAHEDRON20;
        break;
    case TT_LABEL_HEXAHEDRON27:
        return TT_LABEL_HEXAHEDRON27;
        break;
    case TT_LABEL_PRISM6:
        return TT_LABEL_PRISM6;
        break;
    case TT_LABEL_PRISM15:
        return TT_LABEL_PRISM15;
        break;
    case TT_LABEL_PRISM18:
        return TT_LABEL_PRISM18;
        break;
    case TT_LABEL_MATERIAL:
        return TT_LABEL_MATERIAL;
        break;
    case TT_LABEL_GROUPS:
        return TT_LABEL_GROUPS;
        break;
    case TT_LABEL_NODE_RESTRICTIONS:
        return TT_LABEL_NODE_RESTRICTIONS;
        break;
    case TT_LABEL_LOAD_PATTERNS:
        return TT_LABEL_LOAD_PATTERNS;
        break;
    case TT_LABEL_NODE:
        return TT_LABEL_NODE;
        break;
    case TT_LABEL_DX:
        return TT_LABEL_DX;
        break;
    case TT_LABEL_DY:
        return TT_LABEL_DY;
        break;
    case TT_LABEL_DZ:
        return TT_LABEL_DZ;
        break;
    case TT_LABEL_NODAL_LOADS:
        return TT_LABEL_NODAL_LOADS;
        break;
    case TT_LABEL_SURFACE_LOADS:
        return TT_LABEL_SURFACE_LOADS;
        break;
    case TT_LABEL_DOMAIN_LOADS:
        return TT_LABEL_DOMAIN_LOADS;
        break;
    case TT_LABEL_TRIANGLE3:
        return TT_LABEL_TRIANGLE3;
        break;
    case TT_LABEL_TRIANGLE6:
        return TT_LABEL_TRIANGLE6;
        break;
    case TT_LABEL_QUADRANGLE4:
        return TT_LABEL_QUADRANGLE4;
        break;
    case TT_LABEL_QUADRANGLE8:
        return TT_LABEL_QUADRANGLE8;
        break;
    case TT_LABEL_QUADRANGLE9:
        return TT_LABEL_QUADRANGLE9;
        break;
    case TT_LABEL_ELEMENT:
        return TT_LABEL_ELEMENT;
        break;
    case TT_LABEL_FORCE:
        return TT_LABEL_FORCE;
        break;

    case TT_EOF:
        return TT_EOF;
        break;
    default:
        return TT_UNKNOWN;
        break;
    }
}

void FemJsonParser::fill(std::istream& file)
{
    if (file.good()) {
        // move the remaining unprocessed buffer to the start
        char count = limit - tok;
        memcpy(buffer, tok, count);
        pos = &buffer[pos - tok];
        marker = &buffer[marker - tok];
        tok = buffer;

        // fill the vacant space
        file.read(buffer + count, 1024 - count - 1);
        limit = buffer + count + file.gcount();
        *limit = '\000';
    }
}