#include "udeclare.hpp"

#define PROCESS_ONE_ELEMENT(r, unused, idx, elem) \
  BOOST_PP_COMMA_IF(idx) BOOST_PP_STRINGIZE(elem)

#define ENUM_INIT(name, ...)\
    vector<string> from##name = { BOOST_PP_SEQ_FOR_EACH_I(PROCESS_ONE_ELEMENT, %%, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) };\
	string to_string(name value) { return from##name[static_cast<int>(value)];} \
	ostream& operator<<(ostream& os, const name& value){ os << to_string(value); return os;}\
    map<string, name> init##name(vector<string> vec) { map<string, name> res; int index = 0; for (const auto& str : vec) res[str] = static_cast<name>(index++); return res;} \
    map<string, name> to##name = init##name(from##name)

ENUM_INIT(Ineq, <=, >=);
ENUM_INIT(ObjSense, MAXIMIZE, MINIMIZE);
ENUM_INIT(Column, numeric_type, string_type, array_type, unsupported);