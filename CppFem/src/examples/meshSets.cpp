#include<cstdlib>
#include"../meshing/NaiveMesh.h"
#include"../physics/MechanicalBoundaryConditions.h"
#include"../materials/VonMisesPlaneStress.h"
#include"../materials/VonMisesPlaneStrain.h"
#include"../materials/PlaneStrain.h"
#include"../plotter/gnuplot_i.hpp"
#include"../fem/ImplAssembly.h"
#include"../solver/descent/ConjugateGradientDescent.h"
#include"../solver/NLSolverCG.h"
#include"meshing/CustomMesh.h"
#include"tensors/Tensor.h"
#include"tensors/Stress.h"

using std::cout;

#define mesh_
#ifdef mesh_

int main(int argc, char* argv[])
{
	
	PlaneStrain mat1(210,0.25);
	
	NaiveMesh mesh(2,2);
	Vector<int> leftInclusive=mesh.setLeftInclusive();
	Vector<int> rightInclusive=mesh.setRightInclusive();
	//mesh.print();
	mesh.setElementMaterial(&mat1);
	
	Force f(10,0);
	for(int i=0;i<rightInclusive.size();++i)
	{
		mesh.getNode(rightInclusive[i])->setPointForce(&f);
	}
	
	Discretization disc(2,mesh);
	
	MechanicalBoundaryConditions encastre(true,true,false);
	encastre.setNodes(leftInclusive);
	
	//mesh.getNode(2)->setPointForce(&f);
	//mesh.getNode(4)->setPointForce(&f);
	//mesh.print();
	disc.addBoundaryCondition(encastre);
	disc.DOFenum();
	//mesh.print();
	ImplAssembly ass(&disc);
	
	ass.printImplAssembly();
	
	ass.matrixAssemblyRoutine();
	ass.vectorAssemblyRoutine();
	
	
	/* 
	ConjugateGradientDescent solver(ass.getGlobalMatrix(), ass.getGlobalVector(), 10e-10, 200, true);
	
	solver.solve(); */
	
	NLSolverCG solver(ass,10e-10,8,10);
	solver.printNLSolver();
	solver.solve();
	
	
	
 	/* ass.localSolutionVectorAssemblyRoutine(solver.getU());
	ass.globalInternalForceAssembly();
	cout<<"CALCULATED VIA Ku\n";
	cout<<ass.getGlobalVector()-ass.getGlobalMatrix()*solver.getU();
	cout<<(ass.getGlobalVector()-ass.getGlobalMatrix()*solver.getU()).norm();
	cout<<"\nASSEMBLED\n";
	cout<<ass.getGlobalVector()-ass.getGlobalInternalForce();
	cout<<(ass.getGlobalVector()-ass.getGlobalInternalForce()).norm(); */
	
	mesh.print();
}

#endif