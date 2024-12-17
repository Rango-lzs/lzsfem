#include "libfemp/io/import/ModelImporterFactory.hpp"
#include "libfemp/io/import/Parser.hpp"
#include "libfemp/LinearAnalysis.hpp"
#include "libfemp/LoadPattern.hpp"
#include "la/ProgressIndicatorStrategy.hpp"

#include <Eigen/Sparse>
#include "libfemp/solvers/CholeskySolver.hpp"
//#include <iostream>

int main() {
	std::string file_name = "../../data/models/surfaces/surface hexa8x2.fem.json";
	std::fstream file;
	file.open(file_name, std::fstream::in);

	// TODO react to failed file open

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
