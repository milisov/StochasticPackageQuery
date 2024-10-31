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
    if (!hasWrittenHeaders)
    {
        for (size_t i = 0; i < headers.size(); ++i)
        {
            file << headers[i];
            if (i < headers.size() - 1)
            {
                file << ",";
            }
        }
        file << "\n";
        hasWrittenHeaders = true;
    }
}
