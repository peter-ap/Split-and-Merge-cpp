#include <SaMlib/dataReader.h>

#include <iostream>
#include <fstream>
#include <iterator>

#include <Eigen/Dense>

#include <SaMlib/structs.h>

FileData get_file_size(std::string filename)
{
    int count_line = 0;
    int count_linewords = 0;

    std::string line;
    std::string linestring;
    /* Creating input filestream */
    std::ifstream file(filename);

    if (!file.is_open()){
        std::cout << "Failed to open file: " << filename << "\n";
        return FileData{};
    }

    while (getline(file, line))
    {
        count_line++;
        linestring = line;
    }

    for (int i = 0; linestring[i] != '\0'; i++)
    {
        if (linestring[i] == ' ')
        {
            count_linewords++;
        }
    }
    FileData data = {count_line, count_linewords};
    return data;
}

Laser get_laser_data(std::string laser_file, int i)
{
    FileData txt_laser = get_file_size(laser_file);

    std::string laser_data[txt_laser.linewords];
    std::string laserscan[(txt_laser.linewords - 4)];
    Laser scan;
    int count = 0;
    std::string line;
    /* Creating input filestream */
    std::ifstream file(laser_file);
    while (std::getline(file, line))
    {
        if (count == i)
        {
            
            std::istringstream buf(line);
            std::istream_iterator<std::string> beg(buf), end;

            std::vector<std::string> tokens(beg, end); // done!
            int j = 0;
            for (auto &s : tokens)
            {
                laser_data[j] = s;
                j++;
            }
        }
        count++;
    }
    int o = 0;
    for (int i = 4; i < txt_laser.linewords; i++)
    {
        scan.rays.push_back(stod(laser_data[i]));
        o++;
    }

    scan.angle = stod(laser_data[3]);
    return scan;
}
