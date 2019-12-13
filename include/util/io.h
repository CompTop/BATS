#pragma once
/*
Utilities for I/O
*/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <topology/data.h>

// read point cloud from csv file.
// header=True indicates first line is a header
Matrix<double> read_point_cloud(std::string &fname, bool header=false) {
	
}
