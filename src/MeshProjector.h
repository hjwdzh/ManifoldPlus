#ifndef MANIFOLD2_MESH_PROJECTOR_H_
#define MANIFOLD2_MESH_PROJECTOR_H_

#include <vector>

#include "types.h"

class MeshProjector
{
public:
	MeshProjector();
	void ComputeHalfEdge();
	void ComputeIndependentSet();
	void UpdateVertexNormal(int conservative);
	void IterativeOptimize(FT len);
	void Project(const MatrixX& V, const MatrixXi& F,
		MatrixX* out_V, MatrixXi* out_F);
	void UpdateNearestDistance();
	int BoundaryCheck();
	void SplitVertices();
	void OptimizePosition(int v, const Vector3& target_p, FT len);
	void OptimizeNormals();
	void PreserveSharpFeatures(FT len_thres);
private:
	std::vector<std::vector<int> > vertex_groups_;
	MatrixX V_, out_V_, target_V_, out_N_, out_FN_;
	MatrixXi F_, out_F_;
	VectorXi V2E_, E2E_;

	VectorX sqrD_;
	VectorXi I_;

	std::vector<int> sharp_vertices_;
	std::vector<Vector3> sharp_positions_;
};
#endif