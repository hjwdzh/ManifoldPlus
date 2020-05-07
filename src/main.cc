#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include "Manifold.h"
#include "Parser.h"
#include "types.h"

int main(int argc, char**argv) {
	Parser parser;
	parser.AddArgument("input", "../examples/input.obj");
	parser.AddArgument("output", "../examples/output.obj");
	parser.AddArgument("resolution", "20000");
	parser.ParseArgument(argc, argv);
	parser.Log();

	MatrixX V, out_V;
	MatrixXi F, out_F;
	igl::readOBJ(parser["input"], V, F);

	int resolution = 0;
	sscanf(parser["resolution"].c_str(), "%d", &resolution);

	Manifold manifold;
	manifold.ProcessManifold(V, F, resolution, &out_V, &out_F);

	igl::writeOBJ(parser["output"], out_V, out_F);

	return 0;
}