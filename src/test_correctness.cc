#include <igl/readOBJ.h>
#include <Eigen/Core>
#include <map>
int main(int argc, char** argv) {
	const char* input = argv[1];
	const char* output = argv[2];

	Eigen::MatrixXd V, VN, FN;
	Eigen::MatrixXi F;
	igl::readOBJ(input, V, F);

	std::map<std::pair<int, int>, int> edge_counts;
	for (int i = 0; i < F.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v0 = F(i, j);
			int v1 = F(i, (j + 1) % 3);
			auto k = std::make_pair(v0, v1);
			if (edge_counts.count(k) == 0)
				edge_counts[k] = 1;
			else
				edge_counts[k] += 1;
		}
	}

	int non_manifold = 0;
	int boundary = 0;
	for (auto& p : edge_counts) {
		if (p.second > 1)
			non_manifold += 1;
		auto rk = std::make_pair(p.first.second, p.first.first);
		if (edge_counts.count(rk) == 0)
			boundary += 1;
	}

	double area = 0;
	std::vector<double> areas;
	areas.resize(F.rows());
	for (int i = 0; i < F.rows(); ++i) {
		Eigen::Vector3d d1 = V.row(F(i, 1)) - V.row(F(i, 0));
		Eigen::Vector3d d2 = V.row(F(i, 2)) - V.row(F(i, 0));
		Eigen::Vector3d n = d1.cross(d2);
		areas[i] = n.norm();
		area += areas[i];
	}
	igl::per_vertex_normals(V, F, VN);
	igl::per_face_normals(V, F, FN);
	int inversion = 0;
	for (int i = 0; i < F.rows(); ++i) {
		if (areas[i] < area * 1e-2)
			continue;
		for (int j = 0; j < 3; ++j) {
			if (VN.row(F(i, j)).dot(FN.row(i)) < -1e-2) {
				inversion += 1;
			}
		}
	}
	printf("%d %d %d\n", boundary, non_manifold, inversion);
	std::ofstream os(output);
	os << boundary << " " << non_manifold << " " << inversion << "\n";
	os.close();
	return 0;
}