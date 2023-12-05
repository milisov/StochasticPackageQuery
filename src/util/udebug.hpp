#ifndef UDEBUG_HPP
#define UDEBUG_HPP

#include <map>
#include <iostream>
#include <iterator>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <vector>
#include <gurobi_c.h>

using std::map;
using std::pair;
using std::string;
using std::make_pair;
using std::cout;
using std::endl;
using std::declval;
using std::begin;
using std::cerr;
using std::vector;
using std::istream_iterator;

extern const char* RED;
extern const char* GREEN;
extern const char* RESET;

template <typename T,typename U>                                                   
pair<T,U> operator+(const pair<T,U> & l,const pair<T,U> & r) {   
	return {l.first+r.first,l.second+r.second};                                    
}

#define SFINAE(x, ...)             \
	template <class, class = void> \
	struct x : std::false_type {}; \
	template <class T>             \
	struct x<T, std::void_t<__VA_ARGS__>> : std::true_type {}

SFINAE(DefaultIO, decltype(cout << declval<T &>()));
SFINAE(IsTuple, typename std::tuple_size<T>::type);
SFINAE(Iterable, decltype(begin(declval<T>())));

template <class T>
constexpr char Space(const T &) {
	return (Iterable<T>::value or IsTuple<T>::value) ? ' ' : ' ';
}

template <auto &os>
struct Writer {
	template <class T>
	void Impl(T const &t) const {
		if constexpr (DefaultIO<T>::value) os << t;
		else if constexpr (Iterable<T>::value) {
			int i = 0;
			os << "[";
			for (auto &&x : t) ((i++) ? (os << Space(x), Impl(x)) : Impl(x));
			os << "]";
		} else if constexpr (IsTuple<T>::value)
			std::apply([this](auto const &... args) {
				int i = 0;
				os << "{";
				(((i++) ? (os << ' ', Impl(args)) : Impl(args)), ...);
				os << "}";
			}, t);
		else static_assert(IsTuple<T>::value, "No matching type for print");
	}
	template <class F, class... Ts>
	auto &operator()(F const &f, Ts const &... ts) const {
		return Impl(f), ((os << ' ', Impl(ts)), ...), os <<'\n', *this;
	}
};

#ifdef DEBUG
	#define deb(args...)                                    \
		{                                                     \
			std::string _s = #args;                           \
			std::replace(_s.begin(), _s.end(), ',', ' ');     \
			std::stringstream _ss(_s);                        \
			istream_iterator<string> _it(_ss);      \
			cerr << RED << "File " << __FILE__  \
							<< ", Line " << __LINE__ << RESET << "\n"; \
			err(_it, args);                                   \
		}

inline void err(istream_iterator<string> it) {
	std::ignore = it;
}

template <typename T, typename... Args>
void err(istream_iterator<string> it, T a, Args... args) {
	cerr << *it << " = ";
	Writer<cerr>{}(a);
	err(++it, args...);
}

#define ASSERT(...) \
	if (not(__VA_ARGS__)) throw runtime_error(#__VA_ARGS__)
#else
	#define deb(...) 0
	#define ASSERT(...) 0
#endif

using Clock = std::chrono::high_resolution_clock;

class Profiler{
private:
	map<string, pair<double, int>> clocks;
	map<string, std::chrono::time_point<Clock>> timePoints;
public:
	Profiler();
	void clock(const string& label="");
	void stop(const string& label="");
	void add(const Profiler& pro);
	void print() const;
};

#ifdef DEBUG
	#define INIT(pro) Profiler pro
	#define CLK(pro, label) pro.clock(label)
	#define STP(pro, label) pro.stop(label)
	#define CLOCK(pro) pro.clock()
	#define STOP(pro) pro.stop()
	#define PRINT(pro) pro.print()
	#define ADD(pro, local)       \
	{                        \
		_Pragma("omp critical")\
		{											 \
			pro.add(local);			 \
		}											 \
	}
#else
	#define INIT(pro)
	#define CLK(pro, label)
	#define CLOCK(pro)
	#define STP(pro, label)
	#define STOP(pro)
	#define PRINT(pro)
	#define ADD(pro, local)
#endif

void checkGurobi(const bool& error, GRBenv* env, GRBmodel* model, const char* file, const int& line);
string getGurobiStatus(const int& status);

void getBasicVariables(GRBmodel* model, const int& numvars, vector<size_t>& basics);

#define ckg(error, env) checkGurobi(error, env, NULL, __FILE__, __LINE__)
#define ckgb(error, env, model) checkGurobi(error, env, model, __FILE__, __LINE__)

#endif