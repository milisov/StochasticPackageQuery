#include <cassert>

#include "udeclare.hpp"

#define PROCESS_ONE_ELEMENT(r, unused, idx, elem) \
	BOOST_PP_COMMA_IF(idx) BOOST_PP_STRINGIZE(elem)

#define ENUM_INIT(name, ...)\
	vector<string> from##name = { BOOST_PP_SEQ_FOR_EACH_I(PROCESS_ONE_ELEMENT, %%, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) };\
	string str(name value) { return from##name[static_cast<int>(value)];} \
	name to_##name(size_t index) { assert(index>=0 && index<from##name.size()); return static_cast<name>(index);}\
	int to_index(name value) { return static_cast<int>(value); }\
	ostream& operator<<(ostream& os, const name& value){ os << str(value); return os;}\
	map<string, name> init##name(vector<string> vec) { map<string, name> res; int index = 0; for (const auto& str : vec) res[str] = static_cast<name>(index++); return res;} \
	map<string, name> to##name = init##name(from##name)

ENUM_INIT(Ineq, <=, >=);
ENUM_INIT(ObjSense, MAXIMIZE, MINIMIZE);
ENUM_INIT(Column, numeric_type, string_type, array_type, unsupported);

string str(const Bound& bound) {
    return boost::apply_visitor(VariantToString(), bound);
}

string str(const pair<Bound, Bound>& bounds){
    if (bounds.first.which() == 1 && boost::get<double>(bounds.first) == NEG_INF) return fmt::format("<= {}", str(bounds.second));
    if (bounds.second.which() == 1 && boost::get<double>(bounds.second) == POS_INF) return fmt::format(">= {}", str(bounds.first));
    return fmt::format("BETWEEN {} AND {}", str(bounds.first), str(bounds.second));
}