#ifndef BUILD_A_BAR_H
#define BUILD_A_BAR_H
#include <Eigen/Core>
#include <vector>
// Build the Sobolev inner product matrix A_bar
//
// A_bar satisfies the equation:
// A_bar*g = (the derivative of the tangent point energy with respect to vertex positions)
// where g is the discrete fractional Sobolev gradient
//
// Inputs:
//   pt_num  the number of vertices in the curve (equal to #V)
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//	 T  #E by 3 list of unit tangent vectors for each edge of the curve. 
//	   Namely, T.row(e) = (V.row(E(e, 1)) - v.row(E(e, 0)))/L(e) where L is the next input
//	 L  #E by 1 list of the length of each edge in the curve
//   W  #E by #E list of weights used to construct the B matrix
//   W0  #E by #E list of weights used to construct the B0 matrix
// Outputs:
//   A_bar  (3*#V) by (3*#V) matrix corresponding to the fractional Sobolev inner product
//	  for a given curve, and given energy.
//   
void build_A_bar(
	const int& pt_num,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	const Eigen::MatrixXd& W,
	const Eigen::MatrixXd& W0,
	Eigen::MatrixXd& A_bar);
#endif

