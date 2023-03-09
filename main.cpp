#include "init_curve.h"
#include "repulsion_iteration.h"
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

#include <vector>

#include "read_curve_file.h"

// This project was made by Nathan Henry
// The main file was adapted from the registration assignment in CSC419
// PLEASE READ THE README.pdf

int main(int argc, char* argv[])
{

    // ----------------------------------------------------------------------------------------------------------------
    // THESE ARE THE PARAMETERS THAT CAN BE TUNED -- more details in the README pdf

    // copy the file path of the curve you want to use the program on
    // since I wrote a custom function to parse the file, IT ONLY ACCEPTS OBJ FILES!!!!!!!!
    std::string curve_file_path = "C:/Users/Nathan Henry/Desktop/Geometry_processing/nathan-henry-repulsive-curves/data/petersen.obj";

    // -------------

    // alpha and beta are parameters that change the nature of the tangent-point energy
    // (See Figure 4 in the paper for intuition on what they do)
    // values of 3 and 6 are recommended
    double alpha = 3;
    double beta = 6;

    // -------------

    double a_const = 0.01;
    double b_const = 0.9;
    // These are 2 parameters used in my implementation of line search
    // a_const must be between 0 and 0.5
    // b_const must be between 0 and 1
   
    // decreasing a_const gives preference to larger step sizes
    // though, these steps will have a tendency to overshoot

    // increasing b_const can increase the step sizes in a way where it won't overshoot as much if it does overshoot
    // however, increasing b_const results in each step taking a longer time to compute

    // empirically, I have found that for the examples I tried a_const = 0.01 and b_const = 0.9 work well

    
    // -------------

    int max_iters = 20;
    // THIS IS NOT THE MAXIMUM NUMBER OF DESCENT STEPS
    // for each descent step, we will have to project the result of the step back onto constrained configuration space
    // this projection is iterative, and the maximum number is something we can set here
    // usually it should only take 3 iterations, but this is set to 20 to deal with edge cases
    // on much bigger meshes, it should be set to be higher (ex. 50)
    // on very small meshes (less than 10 vertices), I have found that it's faster if you set max_iters to 10
    // though, this is not an absolute rule

    // -------------

    double threshold = 0.001;
    // This is the maximum absolute value of the constraint function
    // In this case, it means that the sum of squares of the total length deviation with the length
    // deviation of each constrained edge must be less than 0.001
    
    // ----------------------------------------------------------------------------------------------------------------


    Eigen::MatrixXd OVX, VX;
    Eigen::MatrixXi FX, E;
    Eigen::MatrixX3d OV, V;

    std::vector<std::vector<int>> Ac;
    std::vector<std::vector<int>> E_adj;
    Eigen::VectorXd lengths;


    // Load mesh which shows x,y,z axis
    igl::read_triangle_mesh(
        ("C:/Users/Nathan Henry/Desktop/Geometry_processing/nathan-henry-repulsive-curves/data/x_y_z_axis.obj"), OVX, FX);

    // I made a cpp file to read OBJ files which contain curves (i.e. no faces)
    read_curve_file(argc > 1 ? argv[1] : curve_file_path, OV, E);

    V.resize(OV.rows(), 3);

    double factor = OV.maxCoeff() / OVX.maxCoeff();
    OVX = OVX * factor / 1.25; // scale axis so it's a bit smaller than the curve



    bool show_samples = true;
    igl::opengl::glfw::Viewer viewer;
    const int xid = viewer.selected_data_index;
    viewer.append_mesh();

    std::cout << R"(
  [space]  toggle animation
  h        to a single step with fractional Sobolev descent
  r        reset the curve to its initial state
  p        hide points and edges of curve
)";

    // predefined color
    const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
    const auto& set_points = [&]()
    {
        if (show_samples)
        {
            // show verts and edges
            viewer.data_list[xid].set_points(V, (1. - (1. - orange.array()) * .8));
            viewer.data_list[xid].set_edges(V, E, Eigen::RowVector3d(1, 1, 1));
        }
        else
        {
            viewer.data_list[xid].clear_points();
            viewer.data_list[xid].clear_edges();
        }
    };
    const auto& reset = [&]()
    {
        V = OV;
        Ac.clear();
        E_adj.clear();

        // based on E and V, compute Ac, E_adj, and lengths (of edges in the curve's initial state)
        // Ac and E_adj as well as their motivation are described in the init_curve.cpp file
        init_curve(E, V, Ac, E_adj, lengths);

        set_points();
    };
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        if (viewer.core().is_animating)
        {
            // do a sobolev descent step and update V
            repulsion_iteration(alpha, beta, a_const, b_const, threshold, max_iters, E, Ac, E_adj, lengths, V);

            set_points();
        }
        return false;
    };
    viewer.callback_key_pressed =
        [&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
    {
        switch (key)
        {
        case ' ':
            viewer.core().is_animating ^= 1;
            break;
        case 'H':
        case 'h':
            // do a sobolev descent step and update V
            repulsion_iteration(alpha, beta, a_const, b_const, threshold, max_iters, E, Ac, E_adj, lengths, V);
            set_points();
            break;
        case 'M':
        case 'm':
        case 'P':
        case 'p':
            show_samples ^= 1;
            set_points();
            break;
        case 'R':
        case 'r':
            reset();
            break;
        case 'S':
        case 's':
        default:
            return false;
        }
        return true;
    };

    viewer.data_list[xid].set_mesh(OVX, FX);
    viewer.data_list[xid].set_colors(orange);

    reset();
    viewer.core().is_animating = false;
    viewer.data().point_size = 10;
    viewer.launch();

    return EXIT_SUCCESS;
}
