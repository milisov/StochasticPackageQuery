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
#include <Highs.h>

using std::string;
using std::ostringstream;
using std::pair;
using std::vector;
using std::ostream;
using std::map;

const double POS_INF = 1.0e30;
const double NEG_INF = -1.0e30;

using SolType = map<size_t, double>;

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

ENUM(Inequality, lteq, gteq);
ENUM(ObjectiveSense, maximize, minimize);
ENUM(Column, numeric_type, string_type, array_type, unsupported);

bool sameSense(ObjectiveSense objectiveSense, ObjSense objSense);
bool sameSense(Inequality inequality, ObjSense objSense);

/**
 * @brief Boost Variant
 * 
 */

using Bound = boost::variant<string, double>;
using Option = boost::variant<bool, double, int>;

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