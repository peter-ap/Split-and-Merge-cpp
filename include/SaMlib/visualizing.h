#ifndef VISUALIZING_H
#define VISUALIZING_H

#include <vector>

#include <Eigen/Dense>

#include <SaMlib/structs.h>

void view_cluster_set(std::vector<Eigen::MatrixXf> data);

void view_cluster_and_corners(ClusteredData dataSet);

#endif
