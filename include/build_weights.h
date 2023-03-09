#ifndef BUILD_WEIGHTS_H
#define BUILD_WEIGHTS_H
#include <Eigen/Core>
#include <vector>
// Build two weight matrices as described in the paper. It will be used in "build_A_bar" to build
// a matrix A_bar which is used for computing Sobolev inner products.
//
// W corresponds to the discrete high-order fractional Sobolev inner product
// W0 corresponds to the discrete low-order term (which will build a matrix B0 which will act
// like a regularizer for the more important matrix B). B and B0 will be used to construct A_bar.
// A_bar is described in "build_A_bar.h"
//
// Inputs:
//   alpha  double which changes the nature of the tangent point energy
//   beta  double which changes the nature of the tangent point energy like alpha
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   V  #V by 3 list of curve vertex positions
//	 T  #E by 3 list of unit tangent vectors for each edge of the curve. 
//	   Namely, T.row(e) = (V.row(E(e, 1)) - v.row(E(e, 0)))/L(e) where L is the next input
//	 L  #E by 1 list of the length of each edge in the curve
// Outputs:
//   W  #E by #E list of weights used to construct the B matrix in build_A_bar.cpp
//   W0  #E by #E list of weights used to construct the B0 matrix in build_A_bar.cpp
void build_weights(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::MatrixXd& W,
	Eigen::MatrixXd& W0);
#endif
