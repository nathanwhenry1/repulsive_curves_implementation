#ifndef CONSTRAINT_DERIVATIVE_H
#define CONSTRAINT_DERIVATIVE_H
#include <Eigen/Core>
#include <vector>

// Compute the derivative of each constraint with respect to vert positions.
// The constraints used are the total length constraint and the 
// individual edge length constraint. The formulas for these constraints are on
// page 13 of the paper, but I computed the derivative formulas by hand and then
// confirmed the result with Mathematica).
//
// Inputs:
//   E  #E by 2 list of edge indices for the curve
//   V  #V by 3 list of curve vertex positions
//   L  #E by 1 list of the length of each edge in the curve
//   E_adj  a vector of vectors where there is a vector for each point "p" consisting of 
//	   the edge indices for all edges containing the point "p"
//   constr_edges  list of the indices of the edges you want to constrain the lengths of.
//	   Note that this is different from constraining the total length of all edges. In the
//	   code, it is automatically set so that the edges which are in this list are the ones
//	   connected to a verticy of index 3 or more.
// Outputs:
//   C  (1 + #constr_edges) by 3 * #V list of the derivative of two constraints: the total length
//	   constraint, and the individual edge length constraints. The first row is the derivative
//	   of the total length constraint, and each following row corresponds to the derivative of
//	   the fixed edge length constraint for a different edge in "constr_edges"
void constraint_derivative(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	const Eigen::VectorXd& L,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXi& constr_edges,
	Eigen::MatrixXd& C);
#endif
