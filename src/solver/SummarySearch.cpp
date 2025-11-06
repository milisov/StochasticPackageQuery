#include "SummarySearch.hpp"
#include <boost/algorithm/string/join.hpp>
#include <fmt/ranges.h>


void SummarySearch::guessOptimalConservativeness(std::vector<std::vector<pair<double, double>>> &history, std::vector<double> &alpha)
{
    fitter.fit_and_predict(history, alpha);
}