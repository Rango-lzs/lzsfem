#include "libfemp/io/import/ModelImporterFactory.hpp"
#include "libfemp/io/import/Parser.hpp"
#include "libfemp/LinearAnalysis.hpp"
#include "libfemp/LoadPattern.hpp"
#include "la/ProgressIndicatorStrategy.hpp"

#include <Eigen/Sparse>
#include "libfemp/solvers/CholeskySolver.hpp"
//#include <iostream>

int main() {
	std::string file_name = "D:/Lzs/PhysicEngine/FEM_Project/femp/data/models/surfaces/surface hexa8x2.fem.json";
	std::fstream file;
	file.open(file_name, std::fstream::in);

	// TODO react to failed file open


	/*try {
		Eigen::Matrix3d JI = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d invJ = JI.inverse();
		std::cout << "Identity Matrix:\n" << JI << std::endl;
		std::cout << "Inverse of Identity Matrix:\n" << invJ << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Exception caught: " << e.what() << std::endl;
	}*/

	fem::Model model;
	model.clear();
	auto parser = fem::ModelImporterFactory::makeFemJsonParser();
	parser->parse(file, model);

	fem::Analysis<double>* p_analysis = new fem::LinearAnalysis<double>();
	fem::AnalysisResult result;
	fem::Solver<double>* pSolver = new fem::CholeskySolver<double>();
	((fem::LinearAnalysis<double>*)(p_analysis))->set(model, model.getLoadPatternList().front(), result, ProgressIndicatorStrategy(), pSolver);

	p_analysis->run(model, model.getLoadPatternList().front(), result, ProgressIndicatorStrategy());

	return 0;
}
