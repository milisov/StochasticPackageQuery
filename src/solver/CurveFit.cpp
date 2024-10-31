#include "CurveFit.hpp"

using std::cout;
using std::endl;

CurveFitter::CurveFitter() {}

double sigmoid(double x) {
    double result = 1.0 / (1.0 + std::exp(-x));
    return std::round(result * 100.0) / 100.0;
}


double CurveFitter::model(const input_vector &input, const parameter_vector &params)
{
    const double A = params(0);
    const double B = params(1);
    const double C = params(2);
    const double D = params(3);

    const double alpha = input(0);

    const double temp = A * atan(B * alpha + C) + D;

    return temp;
}

double CurveFitter::residual(const std::pair<input_vector, double> &data, const parameter_vector &params)
{
    return model(data.first, params) - data.second;
}

parameter_vector CurveFitter::residual_derivative(const std::pair<input_vector, double> &data, const parameter_vector &params)
{
    parameter_vector der;

    const double A = params(0);
    const double B = params(1);
    const double C = params(2);
    const double D = params(3);

    const double alpha = data.first(0);
    const double temp = B * alpha + C;
    const double common = 1.0 / (1 + temp * temp);

    der(0) = atan(temp);         // Derivative w.r.t A
    der(1) = A * alpha * common; // Derivative w.r.t B
    der(2) = A * common;         // Derivative w.r.t C
    der(3) = 1;                  // Derivative w.r.t D
    return der;
}

double CurveFitter::predictAlphaK(parameter_vector learnedParams)
{
    double DbyA = -learnedParams(3) / learnedParams(0);
    std::cout << "DbyA" << " " << tan(DbyA) << std::endl;
    double alpha_k = (tan(DbyA) - learnedParams(2)) / learnedParams(1);
    std::cout << "Calculated New Alpha Value:" << " " << alpha_k << std::endl;
    return alpha_k;
}

void CurveFitter::fit_and_predict(std::vector<std::vector<std::pair<double, double>>> &history, std::vector<double> &alpha)
{
    for (int k = 0; k < history.size(); k++)
    {
        const parameter_vector params = 10 * dlib::randm(4, 1); 

        std::vector<std::pair<input_vector, double>> converted_history;
        for (const auto &dataPt : history[k])
        {
            input_vector input;
            input(0) = dataPt.first;
            converted_history.push_back(std::make_pair(input, dataPt.second));
        }

        auto residual_fn = [this](const std::pair<input_vector, double> &data, const parameter_vector &params)
        {
            return residual(data, params);
        };

        auto residual_derivative_fn = [this](const std::pair<input_vector, double> &data, const parameter_vector &params)
        {
            return residual_derivative(data, params);
        };

        std::cout << "derivative error: " << dlib::length(residual_derivative(converted_history[0], params) - dlib::derivative(residual_fn)(converted_history[0], params)) << std::endl;

        parameter_vector learnParams;
        learnParams = 1;

        std::cout << "Use Levenberg-Marquardt" <<std::endl;
        solve_least_squares_lm(dlib::objective_delta_stop_strategy(1e-5),
                               residual_fn,
                               residual_derivative_fn,
                               converted_history,
                               learnParams);

        std::cout << "inferred parameters: " << trans(learnParams) << std::endl;
        std::cout << "solution error:      " << length(learnParams - params) << std::endl;
        std::cout << std::endl;

        double alphaK = predictAlphaK(learnParams);
        std::cout<<"AlphaK no sigmoid "<<alphaK<<std::endl;
        if(alphaK > 1)
        {
            alpha[k] = std::max(0.40, sigmoid(alphaK));
        }else{
            alpha[k] = alpha[k] = std::max(0.40, alphaK);
        }
        std::cout<<"Alpha K "<<alpha[k]<<endl;
    }
}