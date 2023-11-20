#ifndef UDEBUG_HPP
#define UDEBUG_HPP

#include <map>
#include <iostream>
#include <iterator>
#include <chrono>

using std::map;
using std::pair;
using std::string;
using std::make_pair;
using std::cout;
using std::endl;
using std::declval;
using std::begin;
using std::cerr;
using std::istream_iterator;

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
			cerr << "\033[0;31m" << "File " << __FILE__  \
							<< ", Line " << __LINE__ << "\033[0m" << "\n"; \
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
	void clock(string label="");
	void stop(string label="");
	void add(Profiler pro);
	void print();
};

#ifdef DEBUG
	#define INIT(pro) Profiler pro
	#define CLOCK(pro, label) pro.clock(label)
	#define STOP(pro, label) pro.stop(label)
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
	#define CLOCK(pro, label)
	#define STOP(pro, label)
	#define PRINT(pro)
	#define ADD(pro, local)
#endif

#endif