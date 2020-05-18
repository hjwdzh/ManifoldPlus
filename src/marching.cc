#include <igl/readOBJ.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/copyleft/marching_cubes.h>
#include <fstream>

class UniformGrid
{
public:
	UniformGrid()
	: N(0)
	{}
	UniformGrid(int _N) {
		N = _N;
		distances.resize(N);
		for (auto& d : distances) {
			d.resize(N);
			for (auto& v : d)
				v.resize(N, 1e30);
		}
	}
	template <class T>
	T distance(const T* const p) const {
		int px = *(double*)&p[0] * N;
		int py = *(double*)&p[1] * N;
		int pz = *(double*)&p[2] * N;
		if (px < 0 || py < 0 || pz < 0 || px >= N - 1 || py >= N - 1 || pz >= N - 1) {
			T l = (T)0;
			if (px < 0)
				l = l + -p[0] * (T)N;
			else if (px >= N)
				l = l + (p[0] * (T)N - (T)(N - 1 - 1e-3));

			if (py < 0)
				l = l + -p[1] * (T)N;
			else if (py >= N)
				l = l + (p[1] * (T)N - (T)(N - 1 - 1e-3));

			if (pz < 0)
				l = l + -p[2] * (T)N;
			else if (pz >= N)
				l = l + (p[2] * (T)N - (T)(N - 1 - 1e-3));

			return l;
		}
		T wx = p[0] * (T)N - (T)px;
		T wy = p[1] * (T)N - (T)py;
		T wz = p[2] * (T)N - (T)pz;
		T w0 = ((T)1 - wx) * ((T)1 - wy) * ((T)1 - wz) * distances[pz    ][py    ][px    ];
		T w1 = wx 		   * ((T)1 - wy) * ((T)1 - wz) * distances[pz    ][py    ][px + 1];
		T w2 = ((T)1 - wx) * wy 		 * ((T)1 - wz) * distances[pz    ][py + 1][px    ];
		T w3 = wx 		   * wy 		 * ((T)1 - wz) * distances[pz    ][py + 1][px + 1];
		T w4 = ((T)1 - wx) * ((T)1 - wy) * wz 		   * distances[pz + 1][py    ][px    ];
		T w5 = wx 		   * ((T)1 - wy) * wz 		   * distances[pz + 1][py    ][px + 1];
		T w6 = ((T)1 - wx) * wy 		 * wz		   * distances[pz + 1][py + 1][px    ];
		T w7 = wx 		   * wy 		 * wz 		   * distances[pz + 1][py + 1][px + 1];
		T res = w0 + w1 + w2 + w3 + w4 + w5 + w6 + w7;
		//T thres = (T)0.02;

		//commented out for deform_tune
		//if (res > thres)
		//	res = thres;
		//
		
		return res;
	}
	int N;
	std::vector<std::vector<std::vector<double> > > distances;
};
class Mesh
{
public:
	Mesh() : scale(1.0) {}
	std::vector<Eigen::Vector3d> V;
	std::vector<Eigen::Vector3i> F;
	void ReadOBJ(const char* filename) {
		Eigen::MatrixXd v;
		Eigen::MatrixXi f;
		igl::readOBJ(filename, v, f);
		V.resize(v.rows());
		F.resize(f.rows());
		for (int i = 0; i < V.size(); ++i)
			V[i] = v.row(i);
		for (int i = 0; i < F.size(); ++i)
			F[i] = f.row(i);
	}

	void WriteOBJ(const char* filename) {
		std::ofstream os(filename);
		for (int i = 0; i < V.size(); ++i) {
			auto v = V[i];
			v = v * scale + pos;
			os << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
		}
		for (int i = 0; i < F.size(); ++i) {
			os << "f " << F[i][0] + 1 << " " << F[i][1] + 1 << " " << F[i][2] + 1 << "\n";
		}
		os.close();
	}
	double scale;
	Eigen::Vector3d pos;
	void Normalize() {
		double min_p[3], max_p[3];
		for (int j = 0; j < 3; ++j) {
			min_p[j] = 1e30;
			max_p[j] = -1e30;
			for (int i = 0; i < V.size(); ++i) {
				if (V[i][j] < min_p[j])
					min_p[j] = V[i][j];
				if (V[i][j] > max_p[j])
					max_p[j] = V[i][j];
			}
		}
		scale = std::max(max_p[0] - min_p[0], std::max(max_p[1] - min_p[1], max_p[2] - min_p[2])) * 1.1;
		for (int j = 0; j < 3; ++j)
			pos[j] = min_p[j] - 0.05 * scale;
		for (auto& v : V) {
			v = (v - pos) / scale;
		}
		for (int j = 0; j < 3; ++j) {
			min_p[j] = 1e30;
			max_p[j] = -1e30;
			for (int i = 0; i < V.size(); ++i) {
				if (V[i][j] < min_p[j])
					min_p[j] = V[i][j];
				if (V[i][j] > max_p[j])
					max_p[j] = V[i][j];
			}
		}
	}
	void ApplyTransform(Mesh& m) {
		pos = m.pos;
		scale = m.scale;
		for (auto& v : V) {
			v = (v - pos) / scale;
		}
	}
	void ConstructDistanceField(UniformGrid& grid, int thres) {
		Eigen::MatrixXd P(grid.N * grid.N * grid.N, 3);
		int offset = 0;
		for (int i = 0; i < grid.N; ++i) {
			for (int j = 0; j < grid.N; ++j) {
				for (int k = 0; k < grid.N; ++k) {
					P.row(offset) = Eigen::Vector3d(double(k) / grid.N, double(j) / grid.N, double(i) / grid.N);
					offset += 1;
				}
			}
		}

		Eigen::MatrixXd V2(V.size(), 3);
		for (int i = 0; i < V.size(); ++i)
			V2.row(i) = V[i];

		Eigen::MatrixXi F2(F.size(), 3);
		for (int i = 0; i < F.size(); ++i)
			F2.row(i) = F[i];

		Eigen::MatrixXd N(F.size(), 3);
		for (int i = 0; i < F.size(); ++i) {
			Eigen::Vector3d x = V[F[i][1]] - V[F[i][0]];
			Eigen::Vector3d y = V[F[i][2]] - V[F[i][0]];
			N.row(i) = x.cross(y).normalized();
		}

		Eigen::VectorXd sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		igl::point_mesh_squared_distance(P,V2,F2,sqrD,I,C);

		offset = 0;

		for (int i = 0; i < grid.N; ++i) {
			for (int j = 0; j < grid.N; ++j) {
				for (int k = 0; k < grid.N; ++k) {
					if (i < 3 || i + 3 >= grid.N || j < 3 || j + 3 >= grid.N
						|| k < 3 || k + 3 >= grid.N)
						grid.distances[i][j][k] = 1;
					else {
						Eigen::Vector3d n = N.row(I[offset]);
						Eigen::Vector3d off = P.row(offset);
						off -= V[F[I[offset]][0]];
						double d = n.dot(off);
						if (d > 0)
							grid.distances[i][j][k] = sqrt(sqrD[offset]);
						else
							grid.distances[i][j][k] = -sqrt(sqrD[offset]);
					}
					offset += 1;
				}
			}
		}
		if (thres) {
			std::queue<Eigen::Vector3i> q;
			for (int i = 0; i < grid.N; ++i) {
				for (int j = 0; j < grid.N; ++j) {
					for (int k = 0; k < grid.N; ++k) {
						grid.distances[i][j][k] = -std::abs(grid.distances[i][j][k]);
					}
				}
			}
			int counter = 0;
			for (int i = 0; i < grid.N; ++i) {
				if (i > 0 && i < grid.N - 1)
					continue;
				for (int j = 0; j < grid.N; ++j) {
					if (j > 0 && j < grid.N - 1)
						continue;
					for (int k = 0; k < grid.N; ++k) {
						if (k > 0 && k < grid.N - 1)
							continue;
						q.push(Eigen::Vector3i(i, j, k));
						grid.distances[i][j][k] = std::abs(grid.distances[i][j][k]);
						counter += 1;
					}
				}
			}
			while (!q.empty()) {
				auto info = q.front();
				q.pop();
				int i = info[0];
				int j = info[1];
				int k = info[2];
				int di[] = {0, 0, 0, 0, -1, 1};
				int dj[] = {0, 0, -1, 1, 0, 0};
				int dk[] = {-1, 1, 0, 0, 0, 0};
				for (int x = 0; x < 6; ++x) {
					int ii = i + di[x];
					int jj = j + dj[x];
					int kk = k + dk[x];
					if (ii < 0 || ii >= grid.N || jj < 0 || jj >= grid.N || kk < 0 || kk >= grid.N)
						continue;
					if (grid.distances[ii][jj][kk] > 0)
						continue;
					if (std::abs(grid.distances[ii][jj][kk]) > 8e-3) {
						q.push(Eigen::Vector3i(ii, jj, kk));
						grid.distances[ii][jj][kk] = std::abs(grid.distances[ii][jj][kk]);
						counter += 1;
					} else {
						grid.distances[ii][jj][kk] = std::abs(grid.distances[ii][jj][kk]);
					}
				}
			}

		}
	}

	void FromDistanceField(UniformGrid& grid) {
		Eigen::VectorXd S(grid.N * grid.N * grid.N);
		Eigen::MatrixXd GV(grid.N * grid.N * grid.N, 3);
		int offset = 0;
		for (int i = 0; i < grid.N; ++i) {
			for (int j = 0; j < grid.N; ++j) {
				for (int k = 0; k < grid.N; ++k) {
					S[offset] = grid.distances[i][j][k];
					GV.row(offset) = Eigen::Vector3d(k, j, i);
					offset += 1;
				}
			}
		}
		Eigen::MatrixXd SV;
		Eigen::MatrixXi SF;
		igl::copyleft::marching_cubes(S,GV,grid.N,grid.N,grid.N,SV,SF);
		printf("%d %d\n", SV.rows(), SF.rows());
		V.resize(SV.rows());
		F.resize(SF.rows());
		for (int i = 0; i < SV.rows(); ++i)
			V[i] = SV.row(i) / (double)grid.N;
		for (int i = 0; i < SF.rows(); ++i)
			F[i] = SF.row(i);
	}
};

int main(int argc, char** argv) {
	Mesh m;
	m.ReadOBJ(argv[1]);
	m.Normalize();
	UniformGrid grid(128);

	if (argc > 3)
		m.ConstructDistanceField(grid, 1);
	else
		m.ConstructDistanceField(grid, 0);

	Mesh n;
	n.FromDistanceField(grid);
	n.scale = m.scale;
	n.pos = m.pos;
	n.WriteOBJ(argv[2]);
	return 0;
}
