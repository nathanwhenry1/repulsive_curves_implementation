#include "build_weights.h"
#include <Eigen/Dense>

#include <iostream>
#include <vector>

void build_weights(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::MatrixXd& W,
	Eigen::MatrixXd& W0)
{
	int edge_num = E.rows();

	// sigma is defined in the paper
	// it is used to compute the value of the repulsion energy
	double sigma = (beta - 1) / alpha;

	// Now, when constructing W (which will happen on each iteration), we iterate through Ac
	W.setZero(edge_num, edge_num);
	// We will construct W0 at the same time (the weights used to make B0)
	W0.setZero(edge_num, edge_num);
	
	// loop through all edges
	for (size_t I = 0; I < Ac.size(); I++)
	{
		// loop through all edges which don't intersect edge I
		for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
		{
			int J = Ac[I][J_ind]; // the index of the edge which doesn't intersect I

			double elt1 = 0;
			double elt2 = 0;

			// iterate through all combinations of endpoints of these 2 edges
			// use the coordinates of each combination of endpoints for formula from the paper for W and W0

			for (size_t a = 0; a < 2; a++)
			{
				for (size_t b = 0; b < 2; b++)
				{
					int i = E(I, a);
					int j = E(J, b);
					Eigen::Vector3d p = V.row(i);
					Eigen::Vector3d q = V.row(j);
					double diff_norm = (p - q).norm();

					// sometimes when diff_norm gets super small (by 2 vertices overlapping by chance), we divide by 0 and things mess up
					// this is a quick fix
					if (abs(diff_norm) < 0.00000001) {
						elt1 += 1000;
						elt2 += 1000;
					}

					elt1 += 1 / std::pow(diff_norm, 2 * sigma + 1);


					// use alpha = 2 and beta = 4, as specified for B0
					double alph = 2;
					double bet = 4;

					double k = std::pow(((p - q).cross(T.row(I))).norm(), alph) / std::pow(diff_norm, bet);
					elt2 += k / std::pow(diff_norm, 2 * sigma + 1);
				}
			}

			W(I, J) = 0.25 * L(I) * L(J) * elt1; // update the weight matrix W0
			W0(I, J) = 0.25 * L(I) * L(J) * elt2; // update the weight matrix W
		}
	}

}
