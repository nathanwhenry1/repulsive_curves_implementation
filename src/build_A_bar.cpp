#include "build_A_bar.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

void build_A_bar(
	const int& pt_num,
	const Eigen::MatrixXi& E,
	const std::vector<std::vector<int>>& Ac,
	const Eigen::MatrixX3d& T,
	const Eigen::VectorXd& L,
	const Eigen::MatrixXd& W,
	const Eigen::MatrixXd& W0,
	Eigen::MatrixXd& A_bar)
{

	Eigen::MatrixXd B;
	Eigen::MatrixXd B0;
	// Now we will build B and B0

	// Start with B and B0 as zero matrices
	B.setZero(pt_num, pt_num);
	B0.setZero(pt_num, pt_num);

	// B is the inner product matrix for the sobolev inner product
	// B0 acts like a regularizer (it's a lower order derivative) which makes the sum A = B + B0 nicer to work with

	// to build B and B0, we first iterate through all edges I, and then all edges J which don't intersect I
	for (size_t I = 0; I < Ac.size(); I++)
	{
		for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
		{
			int J = Ac[I][J_ind];

			// then we just apply a formula from the paper to construct B and B0 based on W and W0
			for (int a = 0; a < 2; a++)
			{
				for (int b = 0; b < 2; b++)
				{
					B(E(I, a), E(I, b)) += std::pow(-1, a + b) * W(I, J) / std::pow(L(I), 2);
					B(E(J, a), E(J, b)) += std::pow(-1, a + b) * W(I, J) / std::pow(L(J), 2);

					B(E(I, a), E(J, b)) -= std::pow(-1, a + b) * W(I, J) * T.row(I).dot(T.row(J)) / (L(I) * L(J));
					B(E(J, a), E(I, b)) -= std::pow(-1, a + b) * W(I, J) * T.row(J).dot(T.row(I)) / (L(J) * L(I));


					B0(E(I, a), E(I, b)) += 0.25 * W0(I, J);
					B0(E(J, a), E(J, b)) += 0.25 * W0(I, J);

					B0(E(I, a), E(J, b)) -= 0.25 * W0(I, J);
					B0(E(J, a), E(I, b)) -= 0.25 * W0(I, J);
				}
			}
		}
	}

	// Now we construct matrix A by adding B and B0
	Eigen::MatrixXd A = B + B0;

	// Construct matrix A_bar
	A_bar.resize(A.cols() * 3, A.rows() * 3);
	A_bar << A, 0 * A, 0 * A,
		0 * A, A, 0 * A,
		0 * A, 0 * A, A;


}

