#include "init_curve.h"

#include <iostream>
#include <cmath>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>


void init_curve(
	const Eigen::MatrixXi& E,
	const Eigen::MatrixX3d& V,
	std::vector<std::vector<int>>& Ac,
	std::vector<std::vector<int>>& E_adj,
	Eigen::VectorXd& lengths)
{
	int pt_num = V.rows();
	int edge_num = E.rows();

	// Later, in each iteration of repulsion we will have to construct a weight matrix w_{I,J}
	// To do this we will need to iterate through all pairs of edges which are disjoing

	// In order to save time later, we will just store all edges J which are disjoint from
	// an edge I in a vector Ac[I]
	// So Ac is a vector of vectors of the form Ac[I]

	// Ac stands for "adjacency complement" because it's all the edges not adjacent to other edges
	// This will only depend on the topology of the mesh
	// So we won't need to recompute it

	// We use the "vector" class for this since different rows will have different lengths

	// build Ac
	for (size_t i = 0; i < edge_num; i++)
	{
		std::vector<int> row;
		for (size_t j = 0; j < edge_num; j++)
		{
			// if edge i and j don't share any vertices, add edge j to the row for edge i
			if (!(E(i, 0) == E(j, 0) || E(i, 0) == E(j, 1) || E(i, 1) == E(j, 0) || E(i, 1) == E(j, 1))) {
				row.push_back(j);
			}
		}
		Ac.push_back(row);
	}

	// Note that some of the rows of Ac will be empty, and this is fine


	// Define a vector of vectors called E_Adj
	// One row for each vertex
	// row i contains an element for each edge containing vertex i

	// Similar to "Ac", we only need to build this once since it only depends on the topology
	// of the curve, and it will save us time during each iteration when we need to find the
	// set of edges adjacent to a point.

	// iterate thorugh all points
	for (size_t i = 0; i < pt_num; i++)
	{
		std::vector<int> row;
		// iterate through all edges
		for (size_t j = 0; j < edge_num; j++)
		{
			// if edge i and j don't share any vertices, add edge j to the row for edge i
			if (E(j, 0) == i || E(j, 1) == i) {
				row.push_back(j);
			}

		}

		E_adj.push_back(row);
	}


	// We will get all the edge lengths now
	// This is important because during each iteartion we will want to compare
	// new edge lenghts to the edge lengths of the mesh when it was initialized

	Eigen::VectorXd L;
	L.resize(edge_num);
	for (size_t i = 0; i < edge_num; i++)
	{
		L(i) = (V.row(E(i, 0)) - V.row(E(i, 1))).norm();
	}

	// "lengths" is what gets outputted by the function
	lengths = L;
}
