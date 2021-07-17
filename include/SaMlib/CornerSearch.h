#ifndef CORNERSEARCH_H_
#define CORNERSEARCH_H_

#include <vector>

#include <Eigen/Dense>

#include <SaMlib/structs.h>

//set rays from polar coordinates to carhesian coordinates (all measurements must be present)
Eigen::MatrixXf polar_to_carthesian(std::vector<double> rays, double view_angle);
//remove all zero points from carthesian measurments and return carthesian matrix 
Eigen::MatrixXf remove_zeros(Eigen::MatrixXf &matrix);
//Determine first clusters based on data jumps
std::vector<Eigen::MatrixXf> determine_rough_clusters(Eigen::MatrixXf scan);
// Remove outlier clusters
std::vector<Eigen::MatrixXf> remove_outlier_clusters(std::vector<Eigen::MatrixXf> clusterSet);
//Splits the given cluster once!
std::vector<Eigen::MatrixXf> split_and_merge_single_set(Eigen::MatrixXf objectCluster);
//Further split clusters
std::vector<Eigen::MatrixXf> split_and_merge(std::vector<Eigen::MatrixXf> objectClusters);
//Further split clusters
std::vector<Eigen::MatrixXf> split_and_merge_total(std::vector<Eigen::MatrixXf> objectClusters);
//calculate line based on linear regression + calculate error
LineValues get_line(Eigen::MatrixXf data);
//calculate regression error of data i.f.o. ransac line
LineValues calulate_ransac_full_regression(Eigen::MatrixXf data, LineValues line);
//calculate line using ransac based on regression line
LineValues get_line_ransac(Eigen::MatrixXf data, float setSize, int iterations);
//calculate intersection/corner point between two lines
cornerPoint get_intersection(LineValues line_1, LineValues line_2);
//check if corners are actual corners and not the intersection of two lines which do not create a physical corner
bool check_virtual_corners(std::vector<Eigen::MatrixXf> clusters, cornerPoint corner);


//Calulate clusters, lines and cornerpoints
//THE ONE FUNCTION TO CALL THEM ALL!!
ClusteredData calculate_clusters_lines_corners(std::vector<double> rays, double view_angle);


//find file size
FileData get_file_size(std::string filename);
//read file data
Laser get_laser_data(std::string laser_file, int i);

//visualize all clusters
void view_cluster_set(std::vector<Eigen::MatrixXf> data);
//visualize clusters and Corners  
void view_cluster_and_corners(ClusteredData dataSet);

#endif
