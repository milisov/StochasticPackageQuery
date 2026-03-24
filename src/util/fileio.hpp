#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <unistd.h>
#include <limits.h>
#include <string>

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple> // Include this
#include <string>
#include <stdexcept>
#include <type_traits> // For std::enable_if

// ---------------------------------------------------------
// 1. FORWARD DECLARATIONS
// ---------------------------------------------------------

// Declare tuple operator<< so vector can see it (if you have vector<tuple>)
template <typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& tup);

// Declare vector operator<< so tuple can see it (for tuple<vector>)
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec);

// ---------------------------------------------------------
// 2. IMPLEMENTATIONS
// ---------------------------------------------------------

// --- Vector Implementation ---
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << '[';
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i + 1 < vec.size()) os << ';';
    }
    os << ']';
    return os;
}

// --- Tuple Implementation Helpers ---

template <std::size_t I = 0, typename... Ts>
typename std::enable_if<I == sizeof...(Ts), void>::type
printTupleElements(std::ostream&, const std::tuple<Ts...>&) {}

template <std::size_t I = 0, typename... Ts>
typename std::enable_if<I < sizeof...(Ts), void>::type
printTupleElements(std::ostream& os, const std::tuple<Ts...>& tup) {
    if constexpr (I > 0)
        os << ';';
    // This line caused the error previously because 
    // it couldn't see the vector overload:
    os << std::get<I>(tup); 
    printTupleElements<I + 1>(os, tup);
}

// --- Tuple Implementation ---
template <typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& tup) {
    os << '(';
    printTupleElements(os, tup);
    os << ')';
    return os;
}

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

        //writeHeaders();

        // Use a fold expression to write data without trailing comma
        writeData(args...);

        file << "\n";
        file.flush();

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