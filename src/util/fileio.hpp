#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <unistd.h>
#include <limits.h>
#include <string>

class DataWriter
{
private:
    std::ofstream file;
    std::vector<std::string> headers;
    bool hasWrittenHeaders;
    std::string filename;

public:
    DataWriter(const std::string &filename, const std::vector<std::string> &headers);
    ~DataWriter();

    void writeHeaders();

    // Template function to write a row of data
    template <typename... Args>
    void addRow(const Args &...args)
    {
        file.open(filename, std::ios::out | std::ios::app);

        if (!file.is_open())
        {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            throw std::runtime_error("Unable to open file");
        }

        if (sizeof...(args) != headers.size())
        {
            throw std::runtime_error("Row data does not match header count");
        }

        writeHeaders();

        // Use a fold expression to write data without trailing comma
        writeData(args...);

        file << "\n";

        if (file.is_open())
        {
            file.close();
            std::cout << "File closed successfully." << std::endl;
        }
    }

private:
    // Helper function to write data without a trailing comma
    template <typename T, typename... Args>
    void writeData(const T &first, const Args &...rest)
    {
        file << first;
        ((file << ',' << rest), ...); // Fold expression to add commas between elements
    }

    // Overload for the case when there's only one argument (no comma needed)
    template <typename T>
    void writeData(const T &last)
    {
        file << last;
    }
};