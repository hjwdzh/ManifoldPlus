#ifndef MANIFOLD2_LOSS_H_
#define MANIFOLD2_LOSS_H_

#include <ceres/ceres.h>

#include "types.h"

template <typename T>
void Minus(const T* v0, const T* v1, T* v2) {
	v2[0] = v0[0] - v1[0];
	v2[1] = v0[1] - v1[1];
	v2[2] = v0[2] - v1[2];
}

template <typename T>
void Cross(const T* n1, const T* n2, T* n3) {
	n3[0] = n1[1] * n2[2] - n1[2] * n2[1];
	n3[1] = n1[2] * n2[0] - n1[0] * n2[2];
	n3[2] = n1[0] * n2[1] - n1[1] * n2[0];
}

template <typename T>
void Normalize(T* n) {
	T l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	if (l > 1e-6) {
		n[0] /= l;
		n[1] /= l;
		n[2] /= l;
	}
}

template <typename T>
void ComputeNormal(const T* const v0,
	const T* const v1,
	const T* const v2,
	T* n) {
	T d1[3], d2[3];
	Minus(v1, v0, d1);
	Minus(v2, v0, d2);
	Cross(d1, d2, n);
	Normalize(n);
}

template <typename T>
T dot(const T* const v0, const T* const v1) {
	return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

struct FaceError {
	FaceError(double lambda)
	: lambda_(lambda)
	{}

	template <typename T>
	bool operator()(const T* const v0,
	          const T* const v1,
	          const T* const v2,
	          const T* const v3,
	          T* residuals) const {
		T n1[3], n2[3], n3[3], d[3];
		ComputeNormal(v0, v1, v2, n1);
		ComputeNormal(v1, v0, v3, n2);
		Cross(n1, n2, n3);
		Minus(v1, v0, d);
		Normalize(d);
		T sin_theta = dot(d, n3);
		T cos_theta = -dot(n1, n2);
		T theta = atan2(sin_theta, cos_theta) * (T)(180.0/3.141592654);
		residuals[0] = (T)0;
		if (theta < (T)0. && theta > (T)-85.) {
			residuals[0] = (theta + (T)90.) * lambda_;
		}

		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(double lambda) {
	return (new ceres::AutoDiffCostFunction<FaceError, 1, 3, 3, 3, 3>(
	         new FaceError(lambda)));
	}

	double lambda_;
};

struct VertexError {
	VertexError(const Vector3& v)
	: v_(v)
	{}

	template <typename T>
	bool operator()(const T* const v0,
	          T* residuals) const {
		residuals[0] = v0[0] - (T)v_[0];
		residuals[1] = v0[1] - (T)v_[1];
		residuals[2] = v0[2] - (T)v_[2];
		return true;
	}

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const Vector3& v) {
	return (new ceres::AutoDiffCostFunction<VertexError, 3, 3>(
	         new VertexError(v)));
	}

	Vector3 v_;
};

struct EdgeError {
	EdgeError(const double l)
	: l_(l)
	{}

	template <typename T>
	bool operator()(const T* const v0,
		const T* const v1,
		T* residuals) const {
		T d[3];
		Minus(v0, v1, d);
		T l = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		residuals[0] = (T)0;
		if (l < l_ * 0.5)
			residuals[0] = (l - l_ * 0.5);
		else if (l > l_)
			residuals[0] = (l - l_);
		return true;
	}

	static ceres::CostFunction* Create(const double l) {
	return (new ceres::AutoDiffCostFunction<EdgeError, 1, 3, 3>(
	         new EdgeError(l)));
	}

	double l_;

};

#endif