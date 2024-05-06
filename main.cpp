//Implied volatility we will be computed using the bisection, newton-raphson and secant methods


#include <iostream>
#include <ostream>

#include "Solvers.h"

class TestFunction1
{
public:

	double eval(double x) const
	{
		return x * x - 2;
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

int main()
{
	double Acc = 0.000001;
	double left = 0.0;
	double right = 2.0;
	
	double guess = 1.0;

	double target = 0.0; //of course we want to find the root

	TestFunction1 f1;
	TestFunction2 f2(3);

	std::cout << "Root of f1 using Bisection: " << BiSection(f1, left, right, Acc, target) << std::endl;
	std::cout << "Root of f2 using Bisection " << BiSection(f2, left, right, Acc, target) << std::endl;

	std::cout << "Root of f1 using Newton " << NewtonRaphson(f1,guess, Acc, target) << std::endl;
	std::cout << "Root of f2 using Newton " << NewtonRaphson(f2, guess, Acc, target) << std::endl;


}




