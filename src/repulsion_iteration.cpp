#include "repulsion_iteration.h"  

#include "loss_derivative.h"  
#include "build_A_bar.h"  
#include "read_curve_file.h"  

#include "build_weights.h"

#include "constraint_derivative.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


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
	Eigen::MatrixX3d& V)
{
	
	std::cout << "Iteration in Progress" << std::endl;


	// number of points and number of edges
	int pt_num = V.rows();
	int edge_num = E.rows();

	// We want to constrain the length of edges connected to vertices
	// with 3 or more edges connected to the vertex.
	// Otherwise the tangent point energy will casue the graph will explode

	// number of verts with 3 or more edges
	int three_way = 0;
	// Count the number of constrain edges for graph embeddings
	for (size_t i = 0; i < pt_num; i++)
	{
		if (E_adj[i].size() > 2) {
			three_way += E_adj[i].size();
		}
	}

	// vector containing the indices of the edges whose lengths we want to fix
	Eigen::VectorXi constr_edges;
	constr_edges.resize(three_way);

	int n = 0;
	for (size_t i = 0; i < pt_num; i++)
	{
		// if there are more than 2 edges adjacent to this edge
		if (E_adj[i].size() > 2) {
			for (int I_ind = 0; I_ind < E_adj[i].size(); I_ind++) {

				int I = E_adj[i][I_ind];
				constr_edges(n) = I;
				n++;

			}			
		}
	}

	// total length of the original mesh before any of the iterations
	// it will be used later for the "total length constraint"
	double l0 = lengths.sum();

	// Get all the edge lengths
	Eigen::VectorXd L;
	L.resize(edge_num);
	for (size_t i = 0; i < edge_num; i++)
	{
		L(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
	}

	Eigen::MatrixX3d T;
	T.resize(edge_num, 3);
	// for each edge I, we have a corresponding T_I vector (representing the tangent vector)
	for (size_t i = 0; i < edge_num; i++)
	{
		T.row(i) = (V.row(E(i, 1)) - V.row(E(i, 0))) / L(i);
	}

	// Build a weight matrix as described in the paper. It will be used to build
	// a matrix A_bar which is used for computing Sobolev inner products
	Eigen::MatrixXd W;
	Eigen::MatrixXd W0;
	build_weights(alpha, beta, E, Ac, V, T, L, W, W0);

	// Construct matrix A_bar
	Eigen::MatrixXd A_bar;
	A_bar.resize(pt_num * 3, pt_num * 3);
	build_A_bar(pt_num, E, Ac, T, L, W, W0, A_bar);

	
	// build derivative of the loss function
	Eigen::RowVectorXd Deriv;
	Deriv.setZero(3 * pt_num);
	loss_derivative(alpha, beta, E, Ac, E_adj, V, T, L, Deriv);

	// Make constraint derivative matrix
	Eigen::MatrixXd C;
	constraint_derivative(E, V, L, E_adj, constr_edges, C);

	int k = C.rows(); // number of constraints 

	Eigen::MatrixXd Z;
	Z.setZero(k, k);


	// SOLVE FOR UNKNOWN in equation
	// Left*Unknown = Right

	// the first 3*pt_num elements of "Unknown" give the descent direction

	// make the matrix on the left in the equation
	Eigen::MatrixXd Left;
	Left.resize(k + pt_num * 3, k + pt_num * 3);
	Left << A_bar, C.transpose(), C, Z;

	// make the vector on the right in the equation we are solving
	Eigen::VectorXd Right;
	Right.setZero(k + pt_num * 3);
	for (size_t i = 0; i < pt_num * 3; i++)
	{
		Right(i) = Deriv(i);
	}

	Eigen::VectorXd Unknown;
	
	Unknown = Left.fullPivLu().solve(Right);
	// we use full PivLu because partial pivoting sometimes results in "unknown" having 
	// undefined values

	// "g_tilde" from the paper is the first pt_num * 3 elements of Unknown
	
	
	// The step size is t. We start it at 1 and use backtracking line search to reduce it
	// Note that starting at 1 works because we will first normalize the descent direction and the derivative
	double t = 1;

	// Matrix to store the updated vertex position so we can compare it to the original
	Eigen::MatrixX3d V_new;
	V_new.resize(V.rows(), 3);

	// the descent direction as computed from the equation where we solved for "Unknown"
	Eigen::VectorXd descent_dir;
	descent_dir.resize(pt_num * 3);
	for (size_t i = 0; i < pt_num; i++)
	{
		descent_dir(3 * i + 0) = -1 * Unknown(3 * i + 0);
		descent_dir(3 * i + 1) = -1 * Unknown(3 * i + 1);
		descent_dir(3 * i + 2) = -1 * Unknown(3 * i + 2);
	}
	
	// normalize the gradient and descent direction before doing backtracking line search
	descent_dir = descent_dir/descent_dir.norm();
	Deriv = Deriv / Deriv.norm(); 

	// backtracking line search
	while (1) {

		// first take a step with step size t and store new vertex positions in V_new
		for (size_t i = 0; i < pt_num; i++)
		{
			V_new(i, 0) = V(i, 0) - t * Unknown(3 * i + 0);
			V_new(i, 1) = V(i, 1) - t * Unknown(3 * i + 1);
			V_new(i, 2) = V(i, 2) - t * Unknown(3 * i + 2);
		}

		// this vector will store the value of each constraint function at the new vertices	
		Eigen::VectorXd constraint_vals;
		constraint_vals.resize(k);

		// compute new edge lengths
		Eigen::VectorXd L_new;
		L_new.resize(L.rows());
		for (size_t i = 0; i < edge_num; i++)
		{
			L_new(i) = (V_new.row(E(i, 0)) - V_new.row(E(i, 1))).norm();
		}

		// compute new tangent vectors
		Eigen::MatrixX3d T_new;
		T_new.resize(edge_num, 3);
		for (size_t i = 0; i < edge_num; i++)
		{
			T_new.row(i) = (V_new.row(E(i, 1)) - V_new.row(E(i, 0))) / L_new(i);
		}


		// We project the step we have taken onto the constraint space in the following loop
		// This is done by solving an iterative equation
		// We terminate either when the constriants are below "threshold" or when we have done a max number of iterations
		// The max number of iterations is here so that we don't waste time making the constraint error small when 
		// just decreasing the step size will work better
		int iter = 0;
		while (iter < max_iters) {
			// the first constraint is that the total length doesn't change
			// note that l0 is the length of the mesh when it was initialized
			constraint_vals(0) = l0 - L_new.sum();

			// the remaining constraints are that the
			for (size_t i = 1; i < k; i++)
			{
				constraint_vals(i) = lengths(constr_edges(i - 1)) - L_new(constr_edges(i - 1));
			}

			if (constraint_vals.norm() <= threshold) {
				break;
			}

			// to do the projection, we solve a matrix equation
			// Left * Unknonw2 = Right2

			// "Left" is the same as the one from the previous equation

			// Now, construct the "Right2" vector
			Eigen::VectorXd Right2;
			Right2.setZero(k + 3 * pt_num);
			for (size_t i = 0; i < k; i++)
			{
				Right2(i + 3 * pt_num) = -1 * constraint_vals(i);
			}

			Left.resize(k + pt_num * 3, k + pt_num * 3);
			Left << A_bar, C.transpose(), C, Z;

			// Solve the matrix equation and update the vertex positions accordingly
			Eigen::VectorXd Unknown2;
			//Unknown2 = Left.lu().solve(Right2);
			Unknown2 = Left.fullPivLu().solve(Right2);
			for (size_t i = 0; i < pt_num; i++)
			{
				V_new(i, 0) += Unknown2(3 * i + 0);
				V_new(i, 1) += Unknown2(3 * i + 1);
				V_new(i, 2) += Unknown2(3 * i + 2);
			}

			iter++;

			// recompute lengths and tangents based on the new vertex positions before the next iter
			for (size_t i = 0; i < edge_num; i++)
			{
				L_new(i) = (V_new.row(E(i, 0)) - V_new.row(E(i, 1))).norm();
			}
			for (size_t i = 0; i < edge_num; i++)
			{
				T_new.row(i) = (V_new.row(E(i, 1)) - V_new.row(E(i, 0))) / L_new(i);
			}
		}


		// for backtracking line search, we want f_delta < "right side"
		// f_delta = f(x + tdelta(x)) where delta(x) is the direction of descent
		// right side = f(x) + a_const * t * the dot product between the gradient and the direction of descent
		// where f is the energy we want to minimize

		double f_delta = 0;
		double right_side = descent_dir.dot((Eigen::VectorXd)Deriv);

		right_side *= t * a_const;

		// we still need to add f(x) to "right_side"
		// and add f(x + t*delta(x)) to f_delta

		// we do this in the following nested for loop
		for (size_t I = 0; I < edge_num; I++)
		{
			for (size_t J_ind = 0; J_ind < Ac[I].size(); J_ind++)
			{
				int J = Ac[I][J_ind];

				for (size_t i = 0; i < 2; i++)
				{
					for (size_t j = 0; j < 2; j++)
					{
						Eigen::Vector3d d = V.row(E(I, i)) - V.row(E(J, j));
						double d_norm = d.norm();

						double cross_norm = (T.row(I)).cross(d).norm();

						right_side += L(I) * L(J) * 0.25 * std::pow(cross_norm, alpha) * std::pow(d_norm, -1*beta);



						Eigen::Vector3d d_new = V_new.row(E(I, i)) - V_new.row(E(J, j));
						double d_norm_new = d_new.norm();

						double cross_norm_new = (T_new.row(I)).cross(d_new).norm();

						f_delta += L_new(I) * L_new(J) * 0.25 * std::pow(cross_norm_new, alpha) * std::pow(d_norm_new, -1 * beta);
					}
				}
			}
		}

		if (f_delta <= right_side && constraint_vals.norm() < threshold) {
			V = V_new;
			break;
		}

		// decrease the step size
		t = t * b_const;
	}


	// in updating our vertex positions, we might have moved the entire mesh
	// we fix this by translating the mesh so the barycenter is 0

	// compute barycenter
	Eigen::Vector3d x0_new(0, 0, 0); 
	for (size_t i = 0; i < edge_num; i++)
	{
		x0_new += L(i) * 0.5 * (V.row(E(i, 0)) + V.row(E(i, 1)));
	}
	x0_new = x0_new / L.sum();

	// subtract off barycenter
	for (size_t i = 0; i < pt_num; i++)
	{
		V.row(i) = V.row(i) - (Eigen::RowVector3d)x0_new;
	}
	
	std::cout << "Iteration Complete" << std::endl;

}