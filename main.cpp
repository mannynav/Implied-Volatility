//Implied volatility we will be computed using the bisection, newton-raphson and secant methods


#include <iostream>
#include <ostream>
#include <numeric>
#include <cmath>
#include <boost/math/tools/roots.hpp>

#include "Solvers.h"

class TestFunction1
{
public:

	double eval(double x) const
	{
		return x * x - 1;
	}

	double derivative(double x) const
	{
		return 2 * x;
	}
};
class TestFunction2
{
public:

	explicit TestFunction2(double a) : a_{ a } {}

	double eval(double x) const
	{
		return x * x - a_;
	}

	double derivative(double x) const
	{
		return 2 * x;
	}

private:

	double a_{};
};

class EuropeanCall
{
public:

	explicit EuropeanCall(double T, double K,double S, double r) : K_{ K }, T_{ T }, S_{ S }, r_{ r }
	{}


	double Black_Scholes_Price(double S0, double sigma, double r)
	{
		return S0 * normalCDF(d_plus(S0,sigma,r)) - K_ * exp(-r * T_) * normalCDF(d_minus(S0,sigma,r));
	}

	double Vega(double S0, double sigma, double r)
	{
		double pi = 4 * std::atan(1);
		return S0 * sqrt(T_) * exp(-0.5 * d_plus(S0, sigma,r) * d_plus(S0,sigma,r)) / sqrt(2 * pi);
	}


	double normalCDF(double x) {

		const double t1 = 0.319381530;
		const double t2 = -0.356563782;
		const double t3 = 1.781477937;
		const double t4 = -1.821255978;
		const double t5 = 1.330274429;

		double pi = 4 * std::atan(1);

		double k = 1.0 / (1.0 + 0.2316419 * x);

		double poly = k * (t1 + k * (t2 + k * (t3 + k * (t4 + k * t5)))) * exp(-0.5 * x * x)/std::sqrt(2.0*pi);

		double approxCDF = 1.0 - poly;

		// Handle Symmetry of the standard normal distribution:
		if (x < 0) {
			approxCDF = 1.0 - normalCDF(-x);
		}

		return approxCDF;
	}

	double d_plus(double S0, double sigma, double r)
	{
		return (std::log(S0 / K_) + (r + 0.5 * std::pow(sigma, 2.0)) * T_) / (sigma * std::sqrt(T_));
	}

	double d_minus(double S0, double sigma, double r)
	{
		return d_plus(S0, sigma, r) - sigma * std::sqrt(T_);
	}

	double eval(double sigma)
	{
		return Black_Scholes_Price(S_, sigma, r_);
	}


	double derivative(double sigma)
	{
		return Vega(S_, sigma, r_);
	}

	double S_{};
	double T_{};
	double K_{};
	double r_{};

};

int main()
{

	double Acc = 0.000000001;
	double left = 0.01;
	double right = 1.0;
	
	double guess = 1.0;

	double target = 0.0; //of course we want to find the root. For the options, this is the market price!

	TestFunction1 f1;
	TestFunction2 f2(3);
	EuropeanCall c1(1, 20, 25, 0.05);

	std::cout << "Root of f1 using Bisection: " << BiSection(f1, left, right, Acc, target) << std::endl;
	std::cout << "Root of f2 using Bisection " << BiSection(f2, left, right, Acc, target) << std::endl;

	std::cout << "Root of f1 using Newton " << NewtonRaphson(f1,guess, Acc, target) << std::endl;
	std::cout << "Root of f2 using Newton " << NewtonRaphson(f2, guess, Acc, target) << std::endl;


	std::cout << "Root of BS using Bisection " << BiSection(c1, left, right, Acc, 7) << std::endl;
	std::cout << "Root of BS using Newton " << NewtonRaphson(c1, 0.25, Acc, 7) << std::endl;


}




