#ifndef REPULSION_ITERATION_H
#define REPULSION_ITERATION_H

#include <Eigen/Core>
#include <vector>

// Conduct a single iteration of fractional Sobolev descent on a curve
// to minimize the tangent-point energy
//
// The first 6 inputs to this function are explained in much more
// detail in the main.cpp file, where they are meant to be tuned
//
// Inputs:
//   alpha  double which changes the nature of the tangent point energy
//   beta  double which changes the nature of the tangent point energy like alpha
//   a_const  parameter for my implementation of backtracking line search (must be between 0 and 0.5)
//	   used as the constant for the Armijo condition 
//   b_const  parameter for my implementation of backtracking line search (must be between 0 and 1)
//	   where in each step of line search, the step size gets multiplied by b_const
//   threshold  enforces that after the descent step, the norm of the constraints at the new vertex
//	   positions is less than this threshold
//   max_iters  max number of iterations to use when iteratively projecting 
//   E  #E by 2 list of edge indices for the curve
//   Ac  a vector of vectors containing a vector for each edge "I". This vector consists of
//	   the edge indices for all edges which do not intersect edge "I"
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   lengths  #E by 1 list of the length of each edge in the curve before any iterations were done
// Outputs:
//   V  #V by 3 list of curve vertex positions after the descent step 
void repulsion_iteration(
	const double& alpha,
	const double& beta,
	const double& a_const,
	const double& b_const,
	const double& threshold,
	const int& max_iters,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXd& lengths,
	Eigen::MatrixX3d& V);

#endif
