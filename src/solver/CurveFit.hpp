#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <dlib/optimization.h>

typedef dlib::matrix<double, 1, 1> input_vector;
typedef dlib::matrix<double, 4, 1> parameter_vector;

class CurveFitter
{
public:
    double model(const input_vector &input, const parameter_vector &params);
    double residual(const std::pair<input_vector, double> &data, const parameter_vector &params);
    parameter_vector residual_derivative(const std::pair<input_vector, double> &data, const parameter_vector &params);
    void fit_and_predict(std::vector<std::vector<std::pair<double, double>>> &history, std::vector<double> &alpha);
    double predictAlphaK(parameter_vector learnedParams);
    CurveFitter();
};