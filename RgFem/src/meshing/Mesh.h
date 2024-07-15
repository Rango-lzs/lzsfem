
#ifndef MESH_H
#define MESH_H

#include "Element.h"
#include "Node.h"
#include "Vector.h"

//定义简单的mesh数据
class Mesh
{
public:
	/*!\brief Print MESH information
	 */
	virtual void print();
	//*****SETTERS AND GETTERS
	/*!\brief Pointer to a node depending on the node local number starts at 1
	 */
	Node* getNode(int key);
	Element* getElement(int key);
	virtual void setNumElements(int n, int m) = 0;
	int getNumElements();
	virtual void setNumNodes(int n, int m) = 0;
	int getNumNodes();
	Vector<Node*>& getNodes();
	Vector<Element*>& getElements();
	virtual void setElementMaterial(Material* mat) = 0;
protected:
	//Structures to store our data
	Vector<Node*> nodes;/**<Stores nodes in MESH*/
	Vector<Element*> elements;/**<Stores elements in MESH*/
	int numNodes;
	int numElements;
	std::string typeOfElement;

};
#endif 