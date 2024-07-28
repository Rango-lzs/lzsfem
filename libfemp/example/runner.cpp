#include "libfemp/io/import/ModelImporterFactory.hpp"
#include "libfemp/io/import/Parser.hpp"
#include "libfemp/LinearAnalysis.hpp"

int main() {
	std::string file_name = "D:/Lzs/PhysicEngine/FEM_Project/libfemp/libfemp/example/minimal.fem.json";
	std::fstream file;
	file.open(file_name, std::fstream::in);

	// TODO react to failed file open

	fem::Model model;
	model.clear();
	auto parser = fem::ModelImporterFactory::makeFemJsonParser();
	parser->parse(file, model);

	fem::Analysis<double>* p_analysis = new fem::LinearAnalysis<double>();
  return 0;
}
