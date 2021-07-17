
#ifndef MYSTRUCT_H
#define MYSTRUCT_H

#include <vector>

#include <Eigen/Dense>

struct LineValues
{
    float slope = -1;
    float intercept = -1;
    float regressionError = -1;
};

struct lineLink
{
    std::vector<int> lineNumber;
};

struct FileData
{
    int lines;
    int linewords;
};

struct Laser
{
    double angle;
    std::vector<double> rays;

};

struct cornerPoint
{
    float x;
    float y;
    int lineID1;
    int lineID2;
};

struct ClusteredData
{
     std::vector<Eigen::MatrixXf> clusters;
     std::vector<cornerPoint> corners;
     std::vector<LineValues> lines;
};

#endif
