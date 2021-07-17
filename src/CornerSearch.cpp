#include <SaMlib/CornerSearch.h>

#include <math.h>
#include <vector>
#include <iostream>
#include <random>
#include <thread>

#include <Eigen/Dense>

#include <SaMlib/structs.h>

//calculate line based on start and end point
LineValues get_line(float x1, float x2, float y1, float y2)
{
    LineValues lineEquation;

    float dx = x2 - x1;
    float dy = y2 - y1;

    float slope = dy / dx;
    float intercept = y1 - slope * x1;

    if (isinf(slope))
    {
        lineEquation.slope = slope;
        lineEquation.intercept = x1;
    }
    else if (isnan(abs(slope)) or isnan(slope))
    {
        lineEquation.slope = slope;
        lineEquation.intercept = x1;
    }
    else
    {
        lineEquation.slope = slope;
        lineEquation.intercept = intercept;
    }

    return lineEquation;
}

//Polar to carthesian coordinates
Eigen::MatrixXf polar_to_carthesian(std::vector<double> rays, double view_angle)
{
    int num_rays = rays.size();
    Eigen::MatrixXf array(2, num_rays);
    float start_angle = -(view_angle / 2);
    float angle = view_angle / (num_rays - 1);

    for (int i = 0; i < num_rays; i++)
    {
        float omega = start_angle + (angle * i);
        float x = rays[i] * cos(omega);
        float y = rays[i] * sin(omega);
        array(0, i) = x;
        array(1, i) = y;
    }
    return array;
}

//Remove non measurments
Eigen::MatrixXf remove_zeros(Eigen::MatrixXf &matrix)
{
    Eigen::MatrixXf new_matrix;
    std::vector<int> zero_points;

    for (int i = 0; i < matrix.cols(); i++)
    {
        if (matrix(0, i) == 0 or matrix(1, i) == 0)
        {
            zero_points.push_back(i);
        }
        else
        {
            new_matrix.conservativeResize(2, new_matrix.cols() + 1);
            new_matrix.col(new_matrix.cols() - 1) = matrix.col(i);
        }
    }

    return new_matrix;
}

/***
*ABD: Deelt scan op in ruwe clusters wanneer grote afstandsverschillen zijn tussen opeenvolgende punten adhv een drempelwaarde.
***/
std::vector<Eigen::MatrixXf> determine_rough_clusters(Eigen::MatrixXf scan)
{

    float jumpSize = 0.35; // size in which we determine if point belongs to the same cluster or not. (max_length*sin(angle/num_beams+1))

    std::vector<Eigen::MatrixXf> clusterSet; // list with all break numbers
    Eigen::MatrixXf c_scan;
    Eigen::MatrixXf cluster(2, 0);
    //Delete all 0 points!!
    for (int j = 0; j < scan.cols(); j++)
    {
        if (scan(0, j) != 0. && scan(1, j) != 0)
        {
            c_scan.conservativeResize(2, c_scan.cols() + 1);
            c_scan.col(c_scan.cols() - 1) = scan.col(j);
        }
    }
    //Split where there are large jumps in the location of the data
    //First point added into cluster:
    cluster.conservativeResize(2, cluster.cols() + 1);
    cluster(0, cluster.cols() - 1) = c_scan(0, 0);
    cluster(1, cluster.cols() - 1) = c_scan(1, 0);
    for (int i = 1; i < c_scan.cols(); i++)
    {

        float distPrevPoint = sqrt(pow((c_scan(0, i) - c_scan(0, i - 1)), 2) + pow((c_scan(1, i) - c_scan(1, i - 1)), 2)); // Stelling van pythagoras
        if (distPrevPoint < jumpSize)
        {
            cluster.conservativeResize(2, cluster.cols() + 1);
            cluster(0, cluster.cols() - 1) = c_scan(0, i);
            cluster(1, cluster.cols() - 1) = c_scan(1, i);
        }
        else
        {
            clusterSet.push_back(cluster);
            cluster.resize(2, 1);
            cluster(0, 0) = c_scan(0, i);
            cluster(1, 0) = c_scan(1, i);
        }
    }

    // Add last cluster to clusterSet and check if the last cluster does not belong together with the first cluster!
    float distFirstPoint = sqrt(pow((c_scan(0, 0) - c_scan(0, c_scan.cols() - 1)), 2) + pow((c_scan(1, 0) - c_scan(1, c_scan.cols() - 1)), 2));
    if (distFirstPoint < jumpSize)
    {
        // view_custerPoint_set(clusterSet);
        Eigen::MatrixXf joined(2, (clusterSet[0].cols() + cluster.cols()));
        // joined << clusterSet[0].rowwise().reverse(), cluster;
        // joined << cluster, clusterSet[0].rowwise().reverse();
        joined << cluster, clusterSet[0];

        clusterSet.erase(clusterSet.begin());
        clusterSet.insert(clusterSet.begin(), joined);
    }
    else
    {
        clusterSet.push_back(cluster);
    }

    // view_custer_set(clusterSet);
    return clusterSet;
}

// Remove outlier clusters
std::vector<Eigen::MatrixXf> remove_outlier_clusters(std::vector<Eigen::MatrixXf> clusterSet)
{
    float jumpSize = 0.35; // size in which we determine if point belongs to the same cluster or not. (max_length*sin(angle/num_beams+1))
    int clusterAmount = 3; // size of cluster which is considered outliers

    //check for outliers
    for (int i = 0; i < clusterSet.size(); i++)
    {
        if (clusterSet[i].cols() <= clusterAmount)
        {
            clusterSet.erase(clusterSet.begin() + i);

            if (i != 0)
            {

                float x1 = clusterSet[i - 1](0, clusterSet[i - 1].cols() - 1);
                float x2 = clusterSet[i](0, 0);
                float y1 = clusterSet[i - 1](1, clusterSet[i - 1].cols() - 1);
                float y2 = clusterSet[i](1, 0);

                float distance = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2));
                if (distance < jumpSize)
                {
                    Eigen::MatrixXf joined(2, (clusterSet[i - 1].cols() + clusterSet[i].cols()));
                    joined << clusterSet[i - 1], clusterSet[i];
                    clusterSet.insert(clusterSet.begin() + i, joined);
                    clusterSet.erase(clusterSet.begin() + (i + 1));
                    clusterSet.erase(clusterSet.begin() + (i - 1));
                }
            }
            i = i - 1;
        }
    }
    return clusterSet;
}

//Further split clusters
std::vector<Eigen::MatrixXf> split_and_merge(std::vector<Eigen::MatrixXf> objectClusters)
{
    Eigen::MatrixXf lineEquation(2, 0);
    // std::vector<Eigen::MatrixXf> objectClusters;
    std::vector<Eigen::MatrixXf> totalObjectCLusters;
    std::vector<Eigen::MatrixXf> temp_objectClusters;
    std::vector<bool> split;

    /***
         * Split-and-Merge (SaM): Deelt clusters verder op zodanig dat punten binnen één cluster tot dezelfde rechte lijn behoren.
        ***/

    //calculate line equation
    float slope, intercept;
    float x1, x2, y1, y2, dx, dy;
    float startPoint, endPoint;
    bool executeSaM = true;
    // objectClusters.push_back(points[points.size() - 1]);

    while (executeSaM)
    {
        temp_objectClusters.resize(0);
        lineEquation.resize(2, 0);
        for (int i = 0; i < objectClusters.size(); i++)
        {
            objectClusters[i];
            startPoint = 0;
            endPoint = objectClusters[i].cols() - 1;

            x1 = objectClusters[i](0, startPoint);
            x2 = objectClusters[i](0, endPoint);
            y1 = objectClusters[i](1, startPoint);
            y2 = objectClusters[i](1, endPoint);
            LineValues lVal = get_line(x1, x2, y1, y2);

            lineEquation.conservativeResize(2, lineEquation.cols() + 1);
            lineEquation(0, lineEquation.cols() - 1) = lVal.slope;
            lineEquation(1, lineEquation.cols() - 1) = lVal.intercept;
            // view_custer_set_and_start_end_points(objectClusters[i], x1, y1, x2, y2);
        }
        // calculate distance to line
        float breakValue = 0.21;
        float xPoint;
        float yPoint;
        std::vector<float> distance;
        int numAppended = 0;
        for (int i = 0; i < lineEquation.cols(); i++)
        {
            int index = i;
            float a = lineEquation(0, i);
            float b = -1;
            float c = lineEquation(1, i);

            for (int k = 0; k < objectClusters[index].cols(); k++)
            {
                xPoint = objectClusters[index](0, k);
                yPoint = objectClusters[index](1, k);
                distance.push_back(abs(a * xPoint + b * yPoint + c) / sqrt(pow(a, 2) + pow(b, 2)));
            }
            //get index of max element
            int maxElementIndex = std::max_element(distance.begin(), distance.end()) - distance.begin();
            float maxElement = *std::max_element(distance.begin(), distance.end());

            if (maxElement > breakValue)
            {
                Eigen::MatrixXf temp_1(2, 0);
                Eigen::MatrixXf temp_2(2, 0);
                for (int k = 0; k < maxElementIndex; k++)
                {
                    temp_1.conservativeResize(2, temp_1.cols() + 1);
                    temp_1(0, k) = objectClusters[index](0, k);
                    temp_1(1, k) = objectClusters[index](1, k);
                }
                int j = 0;
                for (int k = maxElementIndex; k < objectClusters[i].cols(); k++)
                {
                    temp_2.conservativeResize(2, temp_2.cols() + 1);
                    temp_2(0, j) = objectClusters[index](0, k);
                    temp_2(1, j) = objectClusters[index](1, k);
                    j++;
                }

                if (temp_1.cols() > 3)
                {
                    temp_objectClusters.push_back(temp_1);
                    // view_custer_set_and_start_end_points(temp_1, 0, 0, 0, 0);
                }

                if (temp_2.cols() > 3)
                {
                    temp_objectClusters.push_back(temp_2);
                    // view_custer_set_and_start_end_points(temp_2, 0, 0, 0, 0);
                }

                split.push_back(true);
                distance.clear();
            }
            else
            {
                split.push_back(false);

                temp_objectClusters.push_back(objectClusters[index]);
            }
            distance.clear();
        }

        executeSaM = false;
        for (int t = 0; t < split.size(); t++)
        {
            if (split[t] != false)
            {
                executeSaM = true;
                objectClusters = temp_objectClusters;
            }
        }
        split.clear();
    }

    // //filter clusters based on size
    // objectClusters = filter_clusters_from_size(objectClusters, 5);

    //create line link only to be filled if s&m does occur i.e. line is split
    // lineLink linelink;
    // for (int i = 0; i < objectClusters.size(); i++)
    // {
    //     dataClusters.push_back(objectClusters[i]);
    //     roughClusterLink.push_back(roughClusterID);

    //     if (objectClusters.size() > 1)
    //     {
    //         linelink.lineNumber.push_back(dataClusters.size() - 1);
    //     }
    // }

    // //only add linelink if lineLink is filled!
    // if (objectClusters.size() > 1)
    // {
    //     linkedLines.push_back(linelink);
    // }

    return objectClusters;
}

//Splits the given cluster once!
std::vector<Eigen::MatrixXf> split_and_merge_single_set(Eigen::MatrixXf objectCluster)
{
    Eigen::MatrixXf lineEquation(2, 0);
    std::vector<Eigen::MatrixXf> totalObjectCLusters;
    std::vector<Eigen::MatrixXf> temp_objectClusters;

    /*
    * Split-and-Merge (SaM)
    */

    //calculate line equation
    float slope, intercept;
    float x1, x2, y1, y2, dx, dy;
    float startPoint, endPoint;
    bool executeSaM = true;

    temp_objectClusters.resize(0);
    lineEquation.resize(2, 0);

    startPoint = 0;
    endPoint = objectCluster.cols() - 1;

    x1 = objectCluster(0, startPoint);
    x2 = objectCluster(0, endPoint);
    y1 = objectCluster(1, startPoint);
    y2 = objectCluster(1, endPoint);
    LineValues lVal = get_line(x1, x2, y1, y2);

    // calculate distance to line
    float breakValue = 0.21;
    float xPoint;
    float yPoint;
    std::vector<float> distance;
    int numAppended = 0;

    float a = lVal.slope;
    float b = -1;
    float c = lVal.intercept;

    for (int k = 0; k < objectCluster.cols(); k++)
    {
        xPoint = objectCluster(0, k);
        yPoint = objectCluster(1, k);
        distance.push_back(abs(a * xPoint + b * yPoint + c) / sqrt(pow(a, 2) + pow(b, 2)));
    }
    //get index of max element
    int maxElementIndex = std::max_element(distance.begin(), distance.end()) - distance.begin();
    float maxElement = *std::max_element(distance.begin(), distance.end());

    if (maxElement > breakValue)
    {
        Eigen::MatrixXf temp_1(2, 0);
        Eigen::MatrixXf temp_2(2, 0);
        for (int k = 0; k < maxElementIndex; k++)
        {
            temp_1.conservativeResize(2, temp_1.cols() + 1);
            temp_1(0, k) = objectCluster(0, k);
            temp_1(1, k) = objectCluster(1, k);
        }
        int j = 0;
        for (int k = maxElementIndex; k < objectCluster.cols(); k++)
        {
            temp_2.conservativeResize(2, temp_2.cols() + 1);
            temp_2(0, j) = objectCluster(0, k);
            temp_2(1, j) = objectCluster(1, k);
            j++;
        }

        if (temp_1.cols() > 3)
        {
            temp_objectClusters.push_back(temp_1);
        }

        if (temp_2.cols() > 3)
        {
            temp_objectClusters.push_back(temp_2);
        }

        distance.clear();
    }
    else
    {
        temp_objectClusters.push_back(objectCluster);
    }
    distance.clear();

    return temp_objectClusters;
}

//Further split clusters
std::vector<Eigen::MatrixXf> split_and_merge_total(std::vector<Eigen::MatrixXf> objectClusters)
{
    bool executeSaM = true;
    int vectorSize = 0;
    std::vector<Eigen::MatrixXf> samClusters = {};
    std::vector<Eigen::MatrixXf> returnCluster = {};

    for (int i = 0; i < objectClusters.size(); i++)
    {
        Eigen::MatrixXf objectCluster = objectClusters[i];

        samClusters = split_and_merge_single_set(objectCluster);

        while (executeSaM)
        {
            std::vector<Eigen::MatrixXf> sam = {};

            for (int i = 0; i < samClusters.size(); i++)
            {
                std::vector<Eigen::MatrixXf> temp_sam = {};

                temp_sam = split_and_merge_single_set(samClusters[i]);

                for (int k = 0; k < temp_sam.size(); k++)
                {
                    sam.push_back(temp_sam[k]);
                }
            }
            samClusters.clear();
            samClusters = sam;

            if (vectorSize != samClusters.size())
            {
                executeSaM = true;
                vectorSize = samClusters.size();
            }
            else
            {

                executeSaM = false;
            }
        }
        for (int j = 0; j < samClusters.size(); j++)
        {
            returnCluster.push_back(samClusters[j]);
        }
    }
    return returnCluster;
}

//calculate line based on linear regression + calculate error
LineValues get_line(Eigen::MatrixXf data)
{
    LineValues lineValues;
    int nPoints = data.cols();
    if (nPoints < 2)
    {
        return lineValues;
    }
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    for (int i = 0; i < nPoints; i++)
    {
        sumX += data(0, i);
        sumY += data(1, i);
        sumXY += data(0, i) * data(1, i);
        sumX2 += data(0, i) * data(0, i);
    }
    double xMean = sumX / nPoints;
    double yMean = sumY / nPoints;
    double denominator = sumX2 - sumX * xMean;
    // You can tune the eps (1e-7) below for your specific task
    if (std::fabs(denominator) < 1e-7)
    {
        // Fail: it seems a vertical line
        return lineValues;
    }
    float _slope = (sumXY - sumX * yMean) / denominator;
    float _yInt = yMean - _slope * xMean;

    //calculate regression
    float distance = 0;
    float xPoint, yPoint;
    for (int i = 0; i < nPoints; i++)
    {
        xPoint = data(0, i);
        yPoint = data(1, i);
        distance += abs(_slope * xPoint + -1 * yPoint + _yInt) / sqrt(pow(_slope, 2) + pow(-1, 2));
    }
    lineValues.slope = _slope;
    lineValues.intercept = _yInt;
    lineValues.regressionError = distance / nPoints;
    return lineValues;
}

cornerPoint get_intersection(LineValues line_1, LineValues line_2)
{
    cornerPoint cp;
    cp.x = (line_2.intercept - line_1.intercept) / (line_1.slope - line_2.slope);
    cp.y = line_1.slope * cp.x + line_1.intercept;

    return cp;
}

bool check_virtual_corners(std::vector<Eigen::MatrixXf> clusters, cornerPoint corner)
{
    bool virtualPoint = false;
    float jumpSize = 0.35;
    Eigen::MatrixXf cluster_1 = clusters[corner.lineID1];
    Eigen::MatrixXf cluster_2 = clusters[corner.lineID2];

    float d_points = sqrt(pow(cluster_1(0, (cluster_1.cols() - 1)) - cluster_2(0, 0), 2) + pow(cluster_1(1, (cluster_1.cols() - 1)) - cluster_2(1, 0), 2));
    float d_point1_corner = sqrt(pow(cluster_1(0, cluster_1.cols() - 1) - corner.x, 2) + pow(cluster_1(1, cluster_1.cols() - 1) - corner.y, 2));
    float d_point2_corner = sqrt(pow(corner.x - cluster_2(0, 0), 2) + pow(corner.y - cluster_2(1, 0), 2));

    if (d_points >= 1.5*d_point1_corner && d_points >= 1.5*d_point2_corner)
    {
        virtualPoint = true;
    }
    else if (d_points >= jumpSize || (d_point1_corner >= jumpSize && d_point2_corner >= jumpSize))
    {
        virtualPoint = true;
    }
    else
    {
        virtualPoint = false;
    }

    return virtualPoint;
}

LineValues calulate_ransac_full_regression(Eigen::MatrixXf data, LineValues line)
{
    float _slope = line.slope;
    float _yInt = line.intercept;
    float distance = 0;
    float xPoint = 0;
    float yPoint = 0;

    for (int i = 0; i < data.cols(); i++)
    {
        xPoint = data(0, i);
        yPoint = data(1, i);
        distance += abs(_slope * xPoint + -1 * yPoint + _yInt) / sqrt(pow(_slope, 2) + pow(-1, 2));
    }

    line.regressionError = distance / data.cols();
    return line;
}

//Ransac find line!
LineValues get_line_ransac(Eigen::MatrixXf data, float setSize, int iterations)
{
    int sampleSize = data.cols() * setSize;
    LineValues est_line;
    LineValues best_line;
    bool first_run = true;

    if (sampleSize < 2)
    {
        sampleSize = data.cols();
    }

    for (int k = 0; k < iterations; k++)
    {
        //create matrix with random chosen column indices
        Eigen::MatrixXf random = data;
        Eigen::VectorXi indices = Eigen::VectorXi::LinSpaced(data.cols(), 0, data.cols());
        Eigen::MatrixXf rand(2, sampleSize);

        // std::cout << indices << std::endl;
        std::random_shuffle(indices.data(), indices.data() + random.cols());
        // std::cout << indices << std::endl;
        random = random * indices.asPermutation(); //permute columns to permute rows -> indices.asPermutation * random

        for (int i = 0; i < sampleSize; i++)
        {
            rand(0, i) = random(0, i);
            rand(1, i) = random(1, i);
        }
        if (first_run == true)
        {
            best_line = get_line(rand);
            best_line = calulate_ransac_full_regression(data, best_line);
            first_run = false;
        }
        else
        {
            est_line = get_line(rand);
            est_line = calulate_ransac_full_regression(data, best_line);
            if (est_line.regressionError < best_line.regressionError)
            {
                best_line = est_line;
            }
        }
    }
    return best_line;
}

ClusteredData calculate_clusters_lines_corners(std::vector<double> rays, double view_angle)
{
    Eigen::MatrixXf carthesian = polar_to_carthesian(rays, view_angle);

    carthesian = remove_zeros(carthesian);

    std::vector<Eigen::MatrixXf> clusters = determine_rough_clusters(carthesian);

    clusters = remove_outlier_clusters(clusters);

    clusters = split_and_merge_total(clusters);

    //get all regression lines of each SaM cluster
    std::vector<LineValues> lines = {};
    for (int i = 0; i < clusters.size(); i++)
    {
        // lines.push_back(get_line(clusters[i]));
        lines.push_back(get_line_ransac(clusters[i], 0.25, 50));
    }

    //get all intersections of each line with the next line
    std::vector<cornerPoint> corners = {};
    for (int i = 0; i < lines.size() - 1; i++)
    {

        cornerPoint corner = get_intersection(lines.at(i), lines.at(i + 1));
        corner.lineID1 = i;
        corner.lineID2 = i + 1;
        corners.push_back(corner);
    }
    cornerPoint corner = get_intersection(lines.at(lines.size() - 1), lines.at(0));
    corner.lineID1 = lines.size() - 1;
    corner.lineID2 = 0;
    corners.push_back(corner);

    //Check if corner is an actual corner and not a virtual corner!
    std::vector<cornerPoint> realCorners;
    bool virtualCorner;
    for (int i = 0; i < corners.size(); i++)
    {
        virtualCorner = check_virtual_corners(clusters, corners[i]);
        if (virtualCorner == false)
        {
            realCorners.push_back(corners[i]);
        }
    }

    ClusteredData totalData;
    totalData.clusters = clusters;
    totalData.corners = realCorners;
    totalData.lines = lines;
    return totalData;
}
