#ifndef READ_CURVE_FILE_H
#define READ_CURVE_FILE_H
#include <Eigen/Core>
#include <string>

// Given a file path to an OBJ file containing a curve,
// find the vertex position and edges of this curve.
//
// THE FILE MUST BE AN OBJ FILE --- not another file format
//
// Inputs:
//   input  the file path to an OBJ file containing the curve to be used
// Outputs:
//   V  #V by 3 list of 3D positions of each vertex
//   E  #E by 2 list of indices for the virtices of each edge
//   
void read_curve_file(
    const std::string& input,
    Eigen::MatrixX3d& V,
    Eigen::MatrixXi& E);
#endif


