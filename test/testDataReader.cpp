#include <SaMlib/dataReader.h>

#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

bool try_to_open_file(const std::string& filename)
{
    bool result{false};
        std::ifstream file(filename);
    if (file.is_open())
    {
        result = true;
        file.close();
    }
    return result;
}

int main()
{
    std::string test_file_name{"../test/test_data.txt"};

    assert(try_to_open_file(test_file_name));

    FileData file_data = get_file_size(test_file_name);

    // the file should contain three lines with 7 spaces on every line
    assert(file_data.lines == 3);
    assert(file_data.linewords == 7);

    std::cout << "The test 'testDataReader' succeeded! Hooray!\n";

    return 0;
}
