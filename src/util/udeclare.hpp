#ifndef UDECLARE_HPP
#define UDECLARE_HPP

#include <map>
#include <string>
#include <sstream>
#include <ostream>
#include <utility>
#include <vector>
#include <fmt/core.h>
#include <boost/variant.hpp>
#include <boost/preprocessor.hpp>

using std::string;
using std::ostringstream;
using std::pair;
using std::vector;
using std::ostream;
using std::map;

/**
 * @brief ENUM MACRO
 * 
 */ 

#define ENUM(name, ...)  \
    enum name { __VA_ARGS__ };\
    extern string str(name value); \
    extern name to_##name(size_t index); \
    extern int to_index(name value); \
    extern ostream& operator<<(ostream& os, const name& value);\
    extern map<string, name> to##name

ENUM(Ineq, lteq, gteq);
ENUM(ObjSense, maximize, minimize);
ENUM(Column, numeric_type, string_type, array_type, unsupported);

/**
 * @brief Boost Variant
 * 
 */

const double POS_INF = 1.0e30;
const double NEG_INF = -1.0e30;

using Bound = boost::variant<string, double>;

struct VariantToString : public boost::static_visitor<string> {
    template <typename T>
    string operator()(const T& val) const {
        ostringstream ss;
        ss << val;
        return ss.str();
    }
};

string str(const Bound& bound);
string str(const pair<Bound, Bound>& bounds);

#endif