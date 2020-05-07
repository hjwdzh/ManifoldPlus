#include "MeshProjector.h"

#include <map>
#include <unordered_set>

#include <Eigen/Dense>

#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/point_mesh_squared_distance.h>

#include "Intersection.h"
#include "Loss.h"

#define ZERO_THRES 1e-9
MeshProjector::MeshProjector()
{}

void MeshProjector::ComputeHalfEdge()
{
	V2E_.resize(out_V_.rows());
	E2E_.resize(out_F_.rows() * 3);

	for (int i = 0; i < V2E_.size(); ++i)
		V2E_[i] = -1;
	for (int i = 0; i < E2E_.size(); ++i)
		E2E_[i] = -1;

	std::map<std::pair<int, int>, int> dedges;
	for (int i = 0; i < out_F_.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v0 = out_F_(i, j);
			int v1 = out_F_(i, (j + 1) % 3);
			V2E_[v0] = i * 3 + j;
			auto k = std::make_pair(v1, v0);
			auto it = dedges.find(k);
			if (it == dedges.end()) {
				dedges[std::make_pair(v0, v1)] = i * 3 + j;
			} else {
				int rid = it->second;
				E2E_[i * 3 + j] = rid;
				E2E_[rid] = i * 3 + j;
			}
		}
	}
	for (int i = 0; i < V2E_.size(); ++i) {
		if (V2E_[i] == -1) {
			printf("independent vertex! %d\n", i);
			exit(0);
		}
	}
	for (int i = 0; i < E2E_.size(); ++i) {
		if (E2E_[i] == -1) {
			printf("Wrong edge!\n");
			exit(0);
		}
		if (E2E_[E2E_[i]] != i) {
			printf("Wrong edge 2!\n");
			exit(0);
		}
	}
}

void MeshProjector::SplitVertices() {
	std::vector<std::unordered_set<int> > vlinks(out_V_.rows());
	for (int i = 0; i < out_F_.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v0 = out_F_(i, j);
			vlinks[v0].insert(i * 3 + j);
		}
	}
	int invalid_vertex = 0;
	std::vector<std::pair<int, int> > insert_vertex_info;
	int num_vertices = out_V_.rows();
	for (int i = 0; i < out_V_.rows(); ++i) {
		int deid = V2E_[i];
		int deid0 = deid;
		int vertex_count = 0;
		do {
			vertex_count += 1;
			deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
		} while (deid0 != deid);
		if (vertex_count != vlinks[i].size()) {
			int group_id = 0;
			invalid_vertex += 1;
			while (!vlinks[i].empty()) {
				int deid = *vlinks[i].begin();
				int deid0 = deid;
				std::vector<int> dedges;
				do {
					dedges.push_back(deid);
					deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
				} while (deid0 != deid);
				for (auto& p : dedges)
					vlinks[i].erase(p);

				if (group_id != 0) {
					insert_vertex_info.push_back(
						std::make_pair(num_vertices, i));
					for (auto& p : dedges) {
						out_F_(p / 3, p % 3) = num_vertices;
					}
					num_vertices += 1;
				}
				group_id += 1;
			}
		}
	}
	out_V_.conservativeResize(num_vertices, 3);
	for (auto& p : insert_vertex_info) {
		out_V_.row(p.first) = out_V_.row(p.second);
	}
}

void MeshProjector::ComputeIndependentSet() {
	int marked_vertices = 0;
	int group_id = 0;
	std::vector<int> vertex_colors(out_V_.rows(), -1);
	while (marked_vertices < out_V_.rows()) {
		for (int i = 0; i < vertex_colors.size(); ++i) {
			vertex_groups_.push_back(std::vector<int>());
			auto& group = vertex_groups_.back();
			if (vertex_colors[i] != -1)
				continue;
			int deid = V2E_[i];
			int deid0 = deid;
			bool conflict = false;
			do {
				int next_v = out_F_(deid / 3, (deid + 1) % 3);
				if (vertex_colors[next_v] == group_id) {
					conflict = true;
					break;
				}
				deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
			} while (deid0 != deid);
			if (!conflict) {
				vertex_colors[i] = group_id;
				group.push_back(i);
				marked_vertices += 1;
			}
			std::random_shuffle(group.begin(), group.end());
		}
		group_id += 1;
	}	
}

void MeshProjector::Project(const MatrixX& V, const MatrixXi& F,
	MatrixX* out_V, MatrixXi* out_F)
{
	V_ = V;
	F_ = F;
	out_V_ = *out_V;
	out_F_ = *out_F;

	FT len = (V.row(out_F_(0,0)) - V.row(out_F_(0,1))).norm();

	ComputeHalfEdge();
	SplitVertices();
	ComputeHalfEdge();
	ComputeIndependentSet();
	printf("Optimize...\n");
	IterativeOptimize(len);

	out_V->conservativeResize(out_V_.rows(), 3);
	out_F->conservativeResize(out_F_.rows(), 3);
	for (int i = 0; i < out_V_.rows(); ++i) {
		out_V->row(i) = out_V_.row(i);
	}
	for (int i = 0; i < out_F_.rows(); ++i) {
		out_F->row(i) = out_F_.row(i);
	}
}

void MeshProjector::UpdateNearestDistance()
{
	igl::point_mesh_squared_distance(out_V_, V_, F_, sqrD_, I_, target_V_);
}

void MeshProjector::UpdateVertexNormal(int conservative)
{
	if (out_N_.rows() != out_V_.rows())
		out_N_.resize(out_V_.rows(), 3);
	for (int i = 0; i < out_V_.rows(); ++i) {
		int deid = V2E_[i];
		int deid0 = deid;
		Vector3 n(0,0,0);
		do {
			int f = deid / 3;
			int v0 = out_F_(f, deid % 3);
			int v1 = out_F_(f, (deid + 1) % 3);
			int v2 = out_F_(f, (deid + 2) % 3);
			Vector3 d0 = out_V_.row(v1) - out_V_.row(v0);
			Vector3 d1 = out_V_.row(v2) - out_V_.row(v0);
			d0.normalize();
			d1.normalize();
			auto vn = d0.cross(d1);
			double l = vn.norm();
			vn = vn * (asin(l) / l);
			n += vn;
			deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
		} while (deid0 != deid);

		if (conservative) {
			do {
				int f = deid / 3;
				int v0 = out_F_(f, deid % 3);
				int v1 = out_F_(f, (deid + 1) % 3);
				int v2 = out_F_(f, (deid + 2) % 3);
				Vector3 d0 = out_V_.row(v1) - out_V_.row(v0);
				Vector3 d1 = out_V_.row(v2) - out_V_.row(v0);
				auto vn = d0.cross(d1).normalized();
				if (n.dot(vn) < 0) {
					n -= n.dot(vn) * vn;
				}
				deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
			} while (deid0 != deid);
		}

		out_N_.row(i) = n.normalized();
	}
}

int MeshProjector::BoundaryCheck() {
	igl::per_face_normals(out_V_, out_F_, out_FN_);
	int consistent = 0;
	int inconsistent = 0;
	for (int i = 0; i < out_V_.rows(); ++i) {
		Vector3 n = out_N_.row(i);
		int deid = V2E_[i];
		int deid0 = deid;
		do {
			Vector3 fn = out_FN_.row(deid / 3);
			if (n.dot(fn) < -ZERO_THRES) {
				inconsistent += 1;
				printf("%d %d %f: <%f %f %f> <%f %f %f>\n",
					i, deid / 3, n.dot(fn),
					n[0], n[1], n[2],
					fn[0], fn[1], fn[2]);
			}
			else
				consistent += 1;
			deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
		} while (deid0 != deid);
	}
	return inconsistent;
}

void MeshProjector::IterativeOptimize(FT len) {
	UpdateVertexNormal(1);
	std::vector<std::pair<FT, int> > indices(out_V_.size());
	UpdateNearestDistance();
	for (int iter = 0; iter < 4; ++iter) {
		for (int i = 0; i < sharp_vertices_.size(); ++i) {
			if (sharp_vertices_[i] > 0) {
				sqrD_[i] = 10000 * sharp_vertices_[i] + (sharp_positions_[i]
					- Vector3(out_V_.row(i))).squaredNorm();
				target_V_.row(i) = sharp_positions_[i];
				out_V_.row(i) = sharp_positions_[i];
			}
		}
		for (int i = 0; i < out_V_.rows(); ++i) {
			indices[i] = std::make_pair(sqrD_[i], i);
		}
		std::sort(indices.rbegin(), indices.rend());
		for (int i = 0; i < out_V_.rows(); ++i) {
			OptimizePosition(indices[i].second,
				target_V_.row(indices[i].second), len);
		}
		OptimizeNormals();
		UpdateNearestDistance();
		if (iter == 2)
			PreserveSharpFeatures(len);
	}
}

void MeshProjector::OptimizePosition(int v, const Vector3& p, FT len) {
	std::vector<Vector3> A;
	std::vector<FT> B;
	int deid = V2E_[v];
	int deid0 = deid;
	do {
		int v0 = out_F_(deid / 3, deid % 3);
		int v1 = out_F_(deid / 3, (deid + 1) % 3);
		int v2 = out_F_(deid / 3, (deid + 2) % 3);

		Vector3 vn[] = {out_N_.row(v0), out_N_.row(v1), out_N_.row(v2)};
		for (int i = 0; i < 3; ++i) {
			Vector3 d = (out_V_.row(v2) - out_V_.row(v1));
			d = d.cross(vn[i]).normalized();
			FT b = d.dot(out_V_.row(v1) - out_V_.row(v0));
			
			b -= 1e-5;

			A.push_back(d);
			B.push_back(b);

		}
		deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
	} while (deid0 != deid);

	std::vector<int> attached_dimensions(A.size(), 0);
	std::vector<Vector3> constraints;
	constraints.reserve(3);
	
	for (int i = 0; i < A.size(); ++i) {
		Vector3 offset = p - Vector3(out_V_.row(v));

		FT tar_step = offset.norm();
		if (tar_step < ZERO_THRES)
			return;

		Vector3 tar_dir = offset / tar_step;

#ifdef PROJ_THREE_TIMES
		if (constraints.size() == 1) {
			tar_dir = tar_dir - tar_dir.dot(constraints[0]) * constraints[0];
			FT n = tar_dir.norm();
			if (n < ZERO_THRES)
				return;
			tar_step *= n;
			tar_dir /= n;
		} 
		else if (constraints.size() == 2) {
			Vector3 dir = constraints[0].cross(constraints[1]).normalized();
			tar_dir = tar_dir.dot(dir) * dir;
			FT n = tar_dir.norm();
			if (n < ZERO_THRES)
				return;
			tar_step *= n;
			tar_dir /= n;			
		}
		else if (constraints.size() == 3) {
			return;
		}
#else
		if (constraints.size() > 0) {
			Vector3 c = constraints.back();
			Vector3 temp_dir = tar_dir - tar_dir.dot(c) * c;
			FT n = temp_dir.norm();
			if (n < ZERO_THRES)
				return;
			temp_dir /= n;
			int boundary_constraint = 0;
			Vector3 temp_boundary[3];
			for (int j = 0; j < constraints.size(); ++j) {
				FT denominator = constraints[j].dot(temp_dir);
				if (denominator > -1e-3) {
					temp_boundary[boundary_constraint] = constraints[j];
					boundary_constraint += 1;
				}
			}
			if (boundary_constraint == 3) {
				return;
			}
			if (boundary_constraint == 2) {
				temp_dir = temp_boundary[0].cross(temp_boundary[1]);
				if (temp_dir.dot(tar_dir) < 0)
					temp_dir = -temp_dir;
				FT n = temp_dir.norm();
				if (n < ZERO_THRES)
					return;
				temp_dir /= n;
				boundary_constraint = 0;
				for (int j = 0; j < constraints.size(); ++j) {
					FT denominator = constraints[j].dot(temp_dir);
					if (denominator > -1e-3)
						boundary_constraint += 1;
				}
				if (boundary_constraint == 3)
					return;
			}

			int top = 0;
			for (int j = 0; j < constraints.size(); ++j) {
				FT denominator = constraints[j].dot(temp_dir);
				if (denominator > -1e-3) {
					constraints[top++] = constraints[j];
				}
			}
			constraints.resize(top);
			if (top == 3) {
				return;
			}
			if (top == 2) {
				Vector3 c1 = constraints[0];
				Vector3 c2 = constraints[1];
				Vector3 dir = c1.cross(c2).normalized();
				temp_dir = tar_dir.dot(dir) * dir;
				n = temp_dir.norm();
				if (n < ZERO_THRES)
					return;
				temp_dir /= n;
			}
			tar_step *= n;
			tar_dir = temp_dir;
		}
#endif

		FT max_step = tar_step;
		int constraint_id = -1;
		for (int j = 0; j < A.size(); ++j) {
			if (attached_dimensions[j])
				continue;
			FT denominator = A[j].dot(tar_dir);
			if (denominator < ZERO_THRES)
				continue;
			FT step = B[j] / denominator;
			if (step < max_step) {
				constraint_id = j;
				max_step = step;
			}
		}

		if (max_step < 1e-6)
			max_step = 0;

		out_V_.row(v) += max_step * tar_dir;

		if (max_step == tar_step)
			return;

		int constraint_size = constraints.size();
		int new_element = 0;
		for (int j = 0; j < A.size(); ++j) {
			if (attached_dimensions[j])
				continue;
			FT denominator = A[j].dot(tar_dir);
			B[j] -= denominator * max_step;

			if (B[j] < ZERO_THRES && denominator >= ZERO_THRES) {
				bool linear_dependent = false;
				if (constraint_size == 1
					&& constraints[0].cross(A[j]).norm() < ZERO_THRES)
					linear_dependent = true;
				if (constraint_size == 2) {
					Vector3 n = constraints[0].cross(constraints[1]);
					if (std::abs(n.normalized().dot(A[j])) < ZERO_THRES) {
						linear_dependent = true;
					}
				} 
				if (!linear_dependent) {
					if (new_element == 0) {
						constraints.push_back(A[j]);
						new_element = 1;
						attached_dimensions[j] = 1;
					}
				} else {
					attached_dimensions[j] = 1;
				}
			}
		}
	}
}

void MeshProjector::OptimizeNormals() {
	MatrixX prev_norm = out_N_;
	UpdateVertexNormal(0);
	igl::per_face_normals(out_V_, out_F_, out_FN_);
	for (int i = 0; i < out_N_.rows(); ++i) {
		Vector3 vn = prev_norm.row(i);
		Vector3 target_vn = out_N_.row(i);
		Vector3 d = target_vn - vn;
		FT max_step = 1.0;
		int deid = V2E_[i];
		int deid0 = deid;
		do {
			Vector3 fn = out_FN_.row(deid / 3);
			FT denominator = d.dot(fn);
			if (denominator < -ZERO_THRES) {
				FT step = -fn.dot(vn) / denominator;
				if (step < max_step) {
					max_step = step;
				}
			}
			deid = E2E_[deid / 3 * 3 + (deid + 2) % 3];
		} while (deid0 != deid);
		if (max_step < 0) {
			max_step = 0;
		}
		out_N_.row(i) = vn + max_step * d;
	}
}

void MeshProjector::PreserveSharpFeatures(FT len_thres) {
	/*
	UpdateNearestDistance();
	MatrixX origin_FN;
	igl::per_face_normals(V_, F_, origin_FN);
	igl::per_face_normals(out_V_, out_F_, out_FN_);
	auto consistent = [&](int src_f0, int src_f1) {
		if (src_f0 == src_f1)
			return true;
		Vector3 n1 = origin_FN.row(src_f0);
		Vector3 n2 = origin_FN.row(src_f1);		
		FT norm_angle = std::abs(n1.dot(n2));
		if (norm_angle < std::cos(30 / 180.0 * 3.141592654)) {
			return false;
		}
		return true;
	};

	std::vector<int> sharp_edges;
	std::vector<Vector3> vertex_positions;
	for (int i = 0; i < out_F_.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v0 = out_F_(i, j);
			int v1 = out_F_(i, (j + 1) % 3);
			if (v0 > v1)
				continue;
			int src_f0 = I_[v0];
			int src_f1 = I_[v1];
			if (consistent(src_f0, src_f1))
				continue;
			Vector3 p0 = V_.row(F_(src_f0, 0));
			Vector3 n0 = origin_FN.row(src_f0);

			Vector3 p1 = V_.row(F_(src_f1, 0));
			Vector3 n1 = origin_FN.row(src_f1);

			Vector3 o, t;
			if (!PlaneIntersect(p0, n0, p1, n1, &o, &t)) {
				continue;
			}

			Vector3 u1 = out_V_.row(v0);
			Vector3 u2 = out_V_.row(v1);
			Vector3 ntarget1 = (u1 - o).dot(t) * t + o;
			Vector3 ntarget2 = (u2 - o).dot(t) * t + o;

			if ((u1 - ntarget1).norm() > sqrt(3.0) * len_thres)
				continue;
			if ((u2 - ntarget2).norm() > sqrt(3.0) * len_thres)
				continue;
			sharp_edges.push_back(i * 3 + j);
			vertex_positions.push_back(ntarget1);
			vertex_positions.push_back(ntarget2);
		}
	}


	MatrixX P(vertex_positions.size(), 3);
	memcpy(P.data(), vertex_positions.data(),
		sizeof(Vector3) * vertex_positions.size());
	VectorX sqrD;
	VectorXi I;
	MatrixX tarP;
	*/
	UpdateNearestDistance();
	MatrixX origin_FN;
	igl::per_face_normals(V_, F_, origin_FN);
	auto consistent = [&](int src_f0, int src_f1) {
		if (src_f0 == src_f1)
			return true;
		Vector3 n1 = origin_FN.row(src_f0);
		Vector3 n2 = origin_FN.row(src_f1);		
		FT norm_angle = std::abs(n1.dot(n2));
		if (norm_angle < std::cos(60 / 180.0 * 3.141592654)) {
			return false;
		}
		return true;
	};

	std::vector<std::set<std::pair<int, int> > > vfeatures(out_V_.rows());

	for (int i = 0; i < out_F_.rows(); ++i) {
		for (int j = 0; j < 3; ++j) {
			int v0 = out_F_(i, j);
			int v1 = out_F_(i, (j + 1) % 3);
			if (v1 < v0)
				continue;
			int src_f0 = I_[v0];
			int src_f1 = I_[v1];
			if (consistent(src_f0, src_f1))
				continue;
			if (src_f0 > src_f1)
				std::swap(src_f0, src_f1);
			auto key = std::make_pair(src_f0, src_f1);
			vfeatures[v0].insert(key);
			vfeatures[v1].insert(key);
		}
	}

	std::vector<std::pair<int, int> > face_elements;
	int vid = -1;
	std::vector<int> vertex_to_update;
	std::vector<Vector3> vertex_update_position;

	for (auto& v : vfeatures) {
		vid += 1;
		if (v.size() > 0) {
			face_elements.clear();
			for (auto& e : v)
				face_elements.push_back(e);

			Vector3 p0 = V_.row(F_(face_elements[0].first, 0));
			Vector3 n0 = origin_FN.row(face_elements[0].first);

			Vector3 p1 = V_.row(F_(face_elements[0].second, 0));
			Vector3 n1 = origin_FN.row(face_elements[0].second);

			Vector3 o, t;
			if (!PlaneIntersect(p0, n0, p1, n1, &o, &t)) {
				continue;
			}
			Vector3 target = out_V_.row(vid), ntarget;
			bool solved = false;

			if (face_elements.size() > 1) {
				std::set<int> fset;
				fset.insert(face_elements[0].first);
				fset.insert(face_elements[0].second);
				fset.insert(face_elements[1].first);
				fset.insert(face_elements[1].second);
				fset.erase(face_elements[0].first);
				fset.erase(face_elements[1].second);
				if (fset.size() > 1) {
					printf("...\n");
				}
				FT max_len = 1e30;
				for (auto& p : fset) {
					Vector3 p2 = V_.row(F_(p, 0));
					Vector3 n2 = origin_FN.row(p);
					if (std::abs(t.dot(n2)) > 0.1) {
						FT lambda = (p2 - o).dot(n2) / t.dot(n2);
						Vector3 nt = o + lambda * t;
						FT len = (nt - target).norm();
						if (len < max_len) {
							max_len = len;
							ntarget = nt;
						}
						solved = true;
					}
				}
			}
			if (!solved) {
				ntarget = (Vector3(target) - o).dot(t) * t + o;
			}
			if ((ntarget - target).norm() < len_thres || true) {
				vertex_to_update.push_back((solved) ? -(vid + 1) : vid + 1);
				vertex_update_position.push_back(ntarget);
			}
		}
	}
	MatrixX P(vertex_update_position.size(), 3);
	memcpy(P.data(), vertex_update_position.data(),
		sizeof(Vector3) * vertex_update_position.size());
	VectorX sqrD;
	VectorXi I;
	MatrixX tarP;
	sharp_vertices_.resize(out_V_.rows(), 0);
	sharp_positions_.resize(out_V_.rows());
	igl::point_mesh_squared_distance(P, V_, F_, sqrD, I, tarP);

	for (int i = 0; i < sqrD.size(); i++) {
		if (sqrt(sqrD[i]) < 3e-2 * len_thres) {
			int v = vertex_to_update[i];
			int id = 1;
			if (v < 0) {
				id = 2;
				v = -v;
			}
			v -= 1;
			sharp_vertices_[v] = id;
			sharp_positions_[v] = vertex_update_position[i];
		}
	}

	/*
	std::ofstream os("../examples/sharps.obj");
	sharp_vertices_.resize(out_V_.rows(), 0);
	sharp_positions_.resize(out_V_.rows());
	igl::point_mesh_squared_distance(P, V_, F_, sqrD, I, tarP);
	for (int i = 0; i < sqrD.size(); i += 2) {
		int dedge = sharp_edges[i / 2];
		int v0 = out_F_(dedge / 3, dedge % 3);
		int v1 = out_F_(dedge / 3, (dedge + 1) % 3);
		if (sqrt(sqrD[i]) < len_thres * 1e-1) {
			sharp_vertices_[v0] = 1;
			sharp_positions_[v0] = vertex_positions[i];
			os << "v " << vertex_positions[i][0] << " " << vertex_positions[i][1] << " " << vertex_positions[i][2] << "\n";
		}
		if (sqrt(sqrD[i + 1]) < len_thres * 1e-1) {
			sharp_vertices_[v1] = 1;
			sharp_positions_[v1] = vertex_positions[i + 1];
			os << "v " << vertex_positions[i+1][0] << " " << vertex_positions[i+1][1] << " " << vertex_positions[i+1][2] << "\n";
		}
	}
	os.close();
	*/
}