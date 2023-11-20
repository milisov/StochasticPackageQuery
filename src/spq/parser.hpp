#ifndef PARSER_HPP
#define PARSER_HPP

#include <iostream>
#include <memory>
#include <boost/spirit/include/qi.hpp>
#include <boost/phoenix.hpp>

#include "spq.hpp"
#include "cons.hpp"
#include "obj.hpp"
#include "util/uio.hpp"

using std::shared_ptr;
using std::make_shared;
using std::cerr;

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

using qi::lit;
using qi::lexeme;
using qi::char_;
using qi::alpha;
using qi::alnum;
using qi::_1;
using qi::_2;
using qi::_3;
using qi::_4;
using qi::_5;
using qi::_val;
using qi::int_;
using qi::double_;
using qi::no_case;

template<typename Iterator, typename Skipper>
class SpaqlGrammar : public qi::grammar<Iterator, Skipper>{
public:
	shared_ptr<StochasticPackageQuery> spq;
	qi::rule<Iterator, Skipper> package_query;
	qi::rule<Iterator, string(), Skipper> name, ineq, obj_sense;
	qi::rule<Iterator, vector<string>(), Skipper> column_list;
	qi::rule<Iterator, shared_ptr<Constraint>(), Skipper> constraint, count_constraint, sum_constraint, expected_sum_constraint, var_constraint;
	qi::rule<Iterator, shared_ptr<Objective>(), Skipper> objective, count_objective, sum_objective, expected_sum_objective;
	qi::rule<Iterator, Bound(), Skipper> bound;
	SpaqlGrammar() : spq(make_shared<StochasticPackageQuery>()), SpaqlGrammar::base_type(package_query){
		name %= lexeme[alpha >> (*alnum)];
		ineq %= qi::string("<=") | qi::string(">=");
		obj_sense %= no_case["MAXIMIZE"] | no_case["MINIMIZE"];
		column_list = lit('*') [_val = vector<string>()] | (name % ',') [_val = _1];
		bound %= (name | double_);
		count_constraint = 
			no_case[lit("COUNT")] 
			>> lit('(') 
			>> lit('*') 
			>> lit(')')
			>> (				
				(no_case[lit("BETWEEN")]
				>> bound
				>> no_case[lit("AND")]
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<CountConstraint>(_1, _2))]
				| (lit("<=")
				>> bound
				) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<CountConstraint>(phx::construct<Bound>(NEG_INF), _1))]
				| (lit(">=")
				>> bound
				) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<CountConstraint>(_1, phx::construct<Bound>(POS_INF)))]
			);
		sum_constraint = 
			no_case[lit("SUM")]
			>> lit('(') 
			>> (
				(name
				>> lit(')')
				>>	no_case[lit("BETWEEN")]
				>> bound
				>> no_case[lit("AND")]
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<SumConstraint>(_1, _2, _3))]
				|(name
				>> lit(')')
				>> lit("<=")
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<SumConstraint>(_1, phx::construct<Bound>(NEG_INF), _2))]
				|(name
				>> lit(')')
				>> lit(">=")
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<SumConstraint>(_1, _2, phx::construct<Bound>(POS_INF)))]
			);
		expected_sum_constraint = 
			no_case[lit("EXPECTED")]
			>> no_case[lit("SUM")]
			>> lit('(') 
			>> (
				(name
				>> lit(')')
				>>	no_case[lit("BETWEEN")]
				>> bound
				>> no_case[lit("AND")]
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<ExpectedSumConstraint>(_1, _2, _3))]
				|(name
				>> lit(')')
				>> lit("<=")
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<ExpectedSumConstraint>(_1, phx::construct<Bound>(NEG_INF), _2))]
				|(name
				>> lit(')')
				>> lit(">=")
				>> bound) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<ExpectedSumConstraint>(_1, _2, phx::construct<Bound>(POS_INF)))]
			);
		var_constraint =
			no_case[lit("SUM")]
			>> lit('(') 
			>> (
				name
				>> lit(')')
				>> ineq
				>> bound
				>> no_case[lit("WITH")]
				>> no_case[lit("PROBABILITY")]
				>> ineq
				>> bound
			) [_val = phx::construct<shared_ptr<Constraint>>(phx::new_<VarConstraint>(_1, _3, _5, _2, _4))]
			;
		constraint %= count_constraint | var_constraint | expected_sum_constraint | sum_constraint;
		count_objective = 
			(obj_sense
			>> no_case[lit("COUNT")]
			>> lit('(')
			>> lit('*')
			>> lit(')')) [_val = phx::construct<shared_ptr<Objective>>(phx::new_<CountObjective>(_1))]
			;
		sum_objective =
			(obj_sense
			>> no_case[lit("SUM")]
			>> lit('(')
			>> name
			>> lit(')')) [_val = phx::construct<shared_ptr<Objective>>(phx::new_<SumObjective>(_1, _2))]
			;
		expected_sum_objective =
			(obj_sense
			>> no_case[lit("EXPECTED")]
			>> no_case[lit("SUM")]
			>> lit('(')
			>> name
			>> lit(')')) [_val = phx::construct<shared_ptr<Objective>>(phx::new_<ExpectedSumObjective>(_1, _2))]
			;
		objective %= count_objective | expected_sum_objective | sum_objective;
		package_query = no_case[lit("SELECT")]
			>> no_case[lit("PACKAGE")]
			>> lit('(') 
			>> column_list [phx::bind(&StochasticPackageQuery::setAttrList, spq, _1)]
			>> lit(')')
			>> no_case[lit("FROM")] 
			>> name [phx::bind(&StochasticPackageQuery::setTableName, spq, _1)]
			>> -(no_case[lit("REPEAT")] >> int_ [phx::bind(&StochasticPackageQuery::setRepeat, spq, _1)])
			>> -(
				no_case[lit("SUCH")]
				>> no_case[lit("THAT")]
				>> constraint [phx::bind(&StochasticPackageQuery::addConstraint, spq, _1)]
				>> *(no_case[lit("AND")] >> constraint [phx::bind(&StochasticPackageQuery::addConstraint, spq, _1)])
				)
			>> -(objective [phx::bind(&StochasticPackageQuery::setObjective, spq, _1)])
			;
	}
};

template<typename Iterator, typename Skipper>
shared_ptr<StochasticPackageQuery> parseSpaql(Iterator first, Iterator last){
	SpaqlGrammar<Iterator, Skipper> grammar;
	bool success = qi::phrase_parse(first, last, grammar, qi::space);
	if (success && first==last){
		return grammar.spq;
	}
	cerr << "Parsing query does not match from: '";
	for (auto it = first; it != last; ++it) cerr << *it;
	cerr << "'\n";
	grammar.spq.reset();
	return grammar.spq;
}

inline shared_ptr<StochasticPackageQuery> parseSpaql(string spaql){
	auto begin = spaql.begin();
	return parseSpaql<decltype(begin), qi::space_type>(begin, spaql.end());
}

#endif