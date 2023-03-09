#include "constraint_derivative.h"
#include <Eigen/Dense>

#include <vector>
#include <iostream>

void constraint_derivative(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	const Eigen::VectorXd& L,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::VectorXi& constr_edges,
	Eigen::MatrixXd& C)
{
	int pt_num = V.rows();

	// First we will find the derivative of the "total length constraint"
	// Phi_l = l0 - Sum over all edges I of {L(I)} 

	Eigen::RowVectorXd C_l;
	C_l.resize(3 * pt_num);

	// loop through all points
	for (size_t p = 0; p < pt_num; p++)
	{
		Eigen::RowVector3d deriv_p(0, 0, 0); // will store the derivative wrt phi_p

		// an edge will contribute to the derivative with respect to phi_p
		// if the edge is adjacent to p

		// loop through all edges adjacent to p
		for (size_t I_ind = 0; I_ind < E_adj[p].size(); I_ind++)
		{
			int I = E_adj[p][I_ind];

			// add or subtract from deriv_p based on the derivative formula

			if (E(I, 0) == p)
			{
				deriv_p -= (V.row(E(I, 0)) - V.row(E(I, 1))) / L(I);
			}
			else if (E(I, 1) == p)
			{
				deriv_p -= (V.row(E(I, 1)) - V.row(E(I, 0))) / L(I);
			}
		}

		// update the derivative matrix
		C_l(3 * p) = deriv_p(0);
		C_l(3 * p + 1) = deriv_p(1);
		C_l(3 * p + 2) = deriv_p(2);
	}


	// We do the same thing as above but now for constraining
	// the length of a specific set of edges "constr_edges"

	Eigen::MatrixXd C_e;
	C_e.setZero(constr_edges.rows(), 3 * pt_num);

	for (size_t i = 0; i < C_e.rows(); i++)
	{
		int I = constr_edges(i);

		int i1 = E(I, 0);
		int i2 = E(I, 1);

		Eigen::Vector3d diff = V.row(i1) - V.row(i2);

		// the formula for the derivative was computed by hand
		// the following code implements it

		C_e(i, i1 * 3 + 0) -= diff(0) / L(I);
		C_e(i, i1 * 3 + 1) -= diff(1) / L(I);
		C_e(i, i1 * 3 + 2) -= diff(2) / L(I);

		C_e(i, i2 * 3 + 0) += diff(0) / L(I);
		C_e(i, i2 * 3 + 1) += diff(1) / L(I);
		C_e(i, i2 * 3 + 2) += diff(2) / L(I);
	}

	// put each constraint in a row of C (the overall constraint derivative matrix which we are returning)
	C.resize(C_l.rows() + C_e.rows(), pt_num * 3);
	C << C_l, C_e;

	
	// the paper recommends that you add in a constraint for the barycenter
	// but I found it was faster to just adjust the barycenter manually after each iteration
	// for the applications I am using repulsive curves for

}
