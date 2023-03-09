#include "loss_derivative.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


void loss_derivative(
	const double& alpha,
	const double& beta,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const std::vector<std::vector<int>>& E_adj,
	const Eigen::MatrixX3d& V,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	Eigen::RowVectorXd& Deriv)
{
	// this file only works if we assume that there are no edges (u,v) such that u=v

	// I derived the derivative by hand so that the code would be faster
	// and we wouldn't have to waste computation power on computing the derivative
	// It is very messy because it had to be broken up into 4 cases
	// But, I checked with mathematica afterwards and it is correct.


	// we want to build the derivative of the repulsion energy with respect to the point gamma_i
	// the energy is a sum of many terms, and the derivative will be 0 if gamma_i isn't in a term
	// Thus, we loop through all disjoint pairs of edges where one of the edges contains gamma_i


	int pt_num = V.rows();


	// first loop through all the points we will be differentiating with respect to
	for (size_t p = 0; p < pt_num; p++)
	{
		Eigen::RowVector3d deriv_p(0, 0, 0);

		// now, loop through all the edges "I" containing vertex "i"
		for (size_t I_ind = 0; I_ind < E_adj[p].size(); I_ind++)
		{
			int I = E_adj[p][I_ind];

			// loop through all edges "J" not intersecting "I"
			for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
			{
				int J = Ac[I][J_ind];

				// iterate through combinations of verts from the 2 edges
				for (size_t i = 0; i < 2; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{

						// note that by construction, only I will have point p in it
						// the following if statement is just so that I can know which index of i has p

						if (E(I, i) == p) {

							// relabel the indices so "i1" is the vertex index equal to p
							// "i2" is the one that isn't p
							// this is just so it matches the math I wrote on paper
							int i1 = E(I, i);
							int i2 = E(I, (i + 1) % 2);

							// to fit with notation, we will also set
							int j1 = E(J, j);

							// what I'm about to write won't depend on j, so we don't need to break it up into cases for j

							// Case 1,1

							// the messy cross product in the numerator
							Eigen::Vector3d cross_term = (V.row(i2) - V.row(j1)).cross(V.row(i1)) - V.row(i2).cross(V.row(j1));

							// the difference on the denominator
							Eigen::Vector3d denom_diff = V.row(i1) - V.row(j1);

							Eigen::Vector3d e1(1, 0, 0);
							Eigen::Vector3d e2(0, 1, 0);
							Eigen::Vector3d e3(0, 0, 1);

							Eigen::Vector3d matrix_cross = (Eigen::Vector3d)(V.row(i2)-V.row(j1));

							Eigen::Matrix3d TMat = (matrix_cross.cross(e1)) * e1.transpose() + (matrix_cross.cross(e2)) * e2.transpose() + (matrix_cross.cross(e3)) * e3.transpose();

							Eigen::RowVector3d term1 = (1 - alpha)* std::pow(L(I), -1 * alpha - 1)* (V.row(i1) - V.row(i2))* std::pow(cross_term.norm(), alpha)* std::pow(denom_diff.norm(), -1 * beta);

							Eigen::RowVector3d term2 = alpha * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2) * cross_term.transpose() * TMat;

							Eigen::RowVector3d term3 = -1 * beta * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta - 2) * std::pow(cross_term.norm(), alpha) * denom_diff.transpose();

							deriv_p += 0.25 * L(J) * (term1 + term2 + term3);

							
							// Case 2,1

							denom_diff = V.row(i2) - V.row(j1);

							term1 = (1 - alpha) * std::pow(L(I), -1 * alpha - 1) * (V.row(i1) - V.row(i2)) * std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta);

							term2 = alpha * std::pow(L(I), 1 - alpha) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2) * cross_term.transpose() * TMat;

							deriv_p += 0.25 * L(J) * (term1 + term2);


							// Now cases with J
							// Case J : 1, 1

							Eigen::Vector3d TJ = T.row(J);
							TMat = (TJ.cross(e1)) * e1.transpose() + (TJ.cross(e2)) * e2.transpose() + (TJ.cross(e3)) * e3.transpose();

							denom_diff = V.row(i1) - V.row(j1);

							cross_term = TJ.cross(denom_diff);

							term1 = std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta) * denom_diff.transpose() / L(I);

							term2 = cross_term.transpose() * TMat;
							term2 = term2 * alpha * L(I) * std::pow(denom_diff.norm(), -1 * beta) * std::pow(cross_term.norm(), alpha - 2);

							term3 = -1 * beta * L(I) * std::pow(denom_diff.norm(), -1 * beta - 2) * std::pow(cross_term.norm(), alpha) * denom_diff.transpose();

							deriv_p += 0.25 * L(J) * (term1 + term2 + term3);


							// Case J : 1, 2

							denom_diff = V.row(i2) - V.row(j1);

							cross_term = TJ.cross(denom_diff);

							deriv_p += 0.25 * (L(J) / L(I)) * std::pow(cross_term.norm(), alpha) * std::pow(denom_diff.norm(), -1 * beta) * (V.row(i1) - V.row(i2));

						}
						
					}
				}

			}

		}

		// put values in the derivative vector for the derivative with respect to p (which is a vector) so it corresponds to 3 elements
		Deriv(3 * p) = deriv_p(0);
		Deriv(3 * p + 1) = deriv_p(1);
		Deriv(3 * p + 2) = deriv_p(2);
	}

	
}
