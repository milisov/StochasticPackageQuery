#include "fileio.hpp"

// Constructor

DataWriter::DataWriter(const std::string &filename, const std::vector<std::string> &headers)
    : headers(headers), hasWrittenHeaders(false), filename(filename)
{
}

DataWriter::~DataWriter()
{
    if (file.is_open())
    {
        file.close();
        std::cout << "File closed successfully." << std::endl;
    }
}

void DataWriter::writeHeaders()
{
    std::cout << "I am writing headers" << std::endl;

    file.open(filename, std::ios::out | std::ios::app);

    if (!file.is_open())
    {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        throw std::runtime_error("Unable to open file");
    }

    for (size_t i = 0; i < headers.size(); ++i)
    {
        file << headers[i];
        if (i < headers.size() - 1)
        {
            file << ",";
        }
    }

    file << "\n";

    if (file.is_open())
    {
        file.close();
        std::cout << "File closed successfully." << std::endl;
    }
}
