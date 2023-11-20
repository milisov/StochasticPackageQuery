#ifndef UDECLARE_HPP
#define UDECLARE_HPP

#include <map>
#include <string>
#include <sstream>
#include <ostream>
#include <utility>
#include <vector>
#include <boost/variant.hpp>
#include <fmt/core.h>
#include <boost/preprocessor.hpp>

using std::string;
using std::ostringstream;
using std::pair;
using std::vector;
using std::ostream;
using std::map;

using Bound = boost::variant<string, double>;

/**
 * @brief ENUM MACRO
 * 
 */ 

#define ENUM(name, ...)  \
    enum name { __VA_ARGS__ };\
    extern string to_string(name value); \
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

struct VariantToString : public boost::static_visitor<string> {
    template <typename T>
    string operator()(const T& val) const {
        ostringstream ss;
        ss << val;
        return ss.str();
    }
};

inline string to_string(const Bound& bound) {
    return boost::apply_visitor(VariantToString(), bound);
}

inline string to_string(const pair<Bound, Bound>& bounds){
    if (bounds.first.which() == 1 && boost::get<double>(bounds.first) == NEG_INF) return fmt::format("<= {}", to_string(bounds.second));
    if (bounds.second.which() == 1 && boost::get<double>(bounds.second) == POS_INF) return fmt::format(">= {}", to_string(bounds.first));
    return fmt::format("BETWEEN {} AND {}", to_string(bounds.first), to_string(bounds.second));
}

#endif