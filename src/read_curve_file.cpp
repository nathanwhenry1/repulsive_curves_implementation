#include "read_curve_file.h"  
#include <cmath>  
#include <Eigen/Dense>  
#include <iostream>  

#include <vector>
#include <fstream>
#include <string>
#include <sstream>

#include <iterator>

void read_curve_file(
    const std::string& input,
    Eigen::MatrixX3d& V,
    Eigen::MatrixXi& E)
{

    // First, we find the number of vertices and number of edges by looping through lines in the obj file
    std::ifstream file(input);

    int v_ind = 0;
    int e_ind = 0;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            // if the line starts with 'v' (corresponding to a vert)
            if (line[0] == 'v') {
                v_ind++;
            }

            // if the line starts with 'l' (corresponding to an edge)
            if (line[0] == 'l') {
                e_ind++;
            }


        }
    }

    file.close();


    // set the sizes of V and E
    V.resize(v_ind, 3);
    E.resize(e_ind, 2);


    // now, we will go through the file again and actually set the elements of V and E to the values in the file
    std::ifstream file2(input);

    v_ind = 0;
    e_ind = 0;

    if (file2.is_open()) {
        std::string line;
        while (std::getline(file2, line)) {
            // if the line starts with 'v' (corresponding to a vert)
            if (line[0] == 'v') {
                line.erase(0, 1);
                std::istringstream iss(line);

                std::vector<std::string> tokens;

                // split up string by "space" delimeter
                std::copy(std::istream_iterator<std::string>(iss),
                    std::istream_iterator<std::string>(),
                    std::back_inserter(tokens));


                V(v_ind, 0) = std::stod(tokens[0]);
                V(v_ind, 1) = std::stod(tokens[1]);
                V(v_ind, 2) = std::stod(tokens[2]);

                v_ind++;
            }

            // if the line starts with 'l' (corresponding to an edge)
            if (line[0] == 'l') {
                line.erase(0, 1);
                std::istringstream iss(line);

                std::vector<std::string> tokens;

                // split up string by "space" delimeter
                std::copy(std::istream_iterator<std::string>(iss),
                    std::istream_iterator<std::string>(),
                    std::back_inserter(tokens));

                E(e_ind, 0) = std::stoi(tokens[0]) - 1;
                E(e_ind, 1) = std::stoi(tokens[1]) - 1;

                e_ind++;
            }


        }
    }

    file.close();


}