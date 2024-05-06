#pragma once

#ifndef Solvers_h
#define Solvers_h


template<typename T> double BiSection(T& function, double left, double right, double tol,double target)
{
	double left_ = left;
	double right_ = right;
	double middle_ = (left_ + right_) / 2;

	double leftVal = function.eval(left_)-target;
	double midVal = function.eval(middle_)-target;

	while(middle_-left_ > tol)
	{
		if((leftVal > 0 && midVal >0) || (leftVal<0 && midVal <0))
		{
			left_ = middle_;
			leftVal = midVal;
		}

		else
		{
			right_ = middle_;
		}

		middle_ = (left_ + right_) / 2;

		midVal = function.eval(middle_) - target;
	}
	return middle_;
}

template<typename T> double NewtonRaphson(T& function, double guess, double acc, double target)
{
	//we start with the guess
	double previous_value = guess;
	double next_value = previous_value - function.eval(previous_value) / function.derivative(previous_value);


	while (next_value - previous_value > acc || previous_value - next_value > acc)
	{
		previous_value = next_value;
		next_value = previous_value - function.eval(previous_value) / function.derivative(previous_value);
	}

	return next_value;
}

#endif