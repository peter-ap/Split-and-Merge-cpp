#ifndef DATAREADER_H_
#define DATAREADER_H_

#include <string>

#include <SaMlib/structs.h>

FileData get_file_size(std::string filename);

Laser get_laser_data(std::string laser_file, int i);

#endif
