#include <SaMlib/visualizing.h>

#include <Eigen/Dense>
#include <vector>

#include <matplotlib/matplotlibcpp.h>
#include <SaMlib/structs.h>

namespace plt = matplotlibcpp;

void view_cluster_set(std::vector<Eigen::MatrixXf> data)
{
    plt::figure_size(1200, 780);

    for (int i = 0; i < data.size(); i++)
    {
        int n = data[i].cols();
        std::vector<double> x_data(n), y_data(n);
        for (int k = 0; k < n; k++)
        {
            x_data.at(k) = data[i](0, k);
            y_data.at(k) = data[i](1, k);
        }
        plt::plot(x_data, y_data, ".");
        x_data.clear();
        y_data.clear();
    }

    // Plot line from given x and y data. Color is selected automatically.
    plt::title("clustered set");

    // Set x-axis to interval [0,1000000]
    plt::xlim(-10, 10);
    plt::ylim(-10, 10);

    // Add graph title

    plt::show();
    plt::close();
}

void view_cluster_and_corners(ClusteredData dataSet)
{
    std::vector<Eigen::MatrixXf> data = dataSet.clusters;
    plt::figure_size(1200, 780);

    for (int i = 0; i < data.size(); i++)
    {
        int n = data[i].cols();
        std::vector<double> x_data(n), y_data(n);
        for (int k = 0; k < n; k++)
        {
            x_data.at(k) = data[i](0, k);
            y_data.at(k) = data[i](1, k);
        }
        plt::plot(x_data, y_data, ".");
        x_data.clear();
        y_data.clear();
    }

    int n = dataSet.corners.size();
    std::vector<double> x_data(n), y_data(n);
    for (int k = 0; k < n; k++)
    {
        x_data.at(k) = dataSet.corners[k].x;
        y_data.at(k) = dataSet.corners[k].y;
    }
    plt::plot(x_data, y_data, ".");
    x_data.clear();
    y_data.clear();
    // Plot line from given x and y data. Color is selected automatically.
    plt::title("clustered set and corners");

    // Set x-axis to interval [0,1000000]
    plt::xlim(-10, 10);
    plt::ylim(-10, 10);

    // Add graph title

    plt::show();
    plt::close();
}
