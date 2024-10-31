#include "solver/CurveFit.hpp"


int main()
{
    CurveFitter fitter;
    std::vector<std::vector<std::pair<double, double>>> history(1); // initialize history with one row
    history[0] = {{0 -0.0206}, {0.4 -0.0206}};

    std::vector<double> alpha = {0.0}; // initialize alpha with one value

    fitter.fit_and_predict(history,alpha); 
}