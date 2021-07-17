#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <cmath>

#include <matplotlib/matplotlibcpp.h>

#include <SaMlib/structs.h>
#include <SaMlib/CornerSearch.h>


int main()
{
    std::string laser_file = "/data/laserscanner_data.txt";
 
    for (int i = 0; i < 92; i++)
    {

        std::cout << "laserscan: " << i << std::endl;
        Laser curScan = get_laser_data(laser_file, i);
        ClusteredData data = calculate_clusters_lines_corners(curScan.rays, curScan.angle);
        view_cluster_set(data.clusters);
        view_cluster_and_corners(data);
    }
    return 0;
}
