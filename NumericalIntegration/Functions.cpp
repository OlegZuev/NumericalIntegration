#include "Functions.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int count_of_calls_to_the_function = 0;

/*
 * func for task
 */
double func(double x)
{
	count_of_calls_to_the_function++;

	switch (VARIANT)
	{
	case 0:
		return pow(3, x - 1) + 4 - x;
	case 5:
		return pow(3, x - 1) + 2 - x;
	case 7:
		return exp(-2 * x) - 2 * x + 1;
	}

	throw std::runtime_error("Incorrect variant");
}

/*
 * max f^(number_of_derivative) on [a,b]
 */
double compute_m(int number_of_derivative)
{
	switch (VARIANT)
	{
	case 0:
	case 5:
		return fabs(3 * pow(log(3), number_of_derivative));
	case 7:
		return fabs(pow(-2, number_of_derivative) * exp(-2));
	}

	throw std::runtime_error("I don't know how to compute this derivative");
}

/*
 * compute first derivative func
 */
double compute_derivative_func(double x)
{
	switch (VARIANT)
	{
	case 0:
	case 5:
		return pow(3, x - 1) * log(3) - 1;
	case 7:
		return -2 * exp(-2 * x) - 2;
	}

	throw std::runtime_error("Incorrect variant");
}


/*
 * allocate matrix n * m
 * n - count of rows
 * m - count of columns
 */
void allocate_matrix(double**& matrix, int n, int m)
{
	matrix = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		matrix[i] = new double[m] { 0 };
	}
}

/*
 * delete matrix n * m
 * n - count of rows
 * m - count of columns
 */
void delete_matrix(double** matrix, int n, int m)
{
	for (int i = 0; i < n; ++i)
	{
		delete[] matrix[i];
	}

	delete[] matrix;
}

/*
 * print matrix n*m into ostr
 */
void print_matrix(double** matrix, int n, int m, std::ostream& ostr)
{
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			ostr << std::fixed << std::setprecision(6) << matrix[i][j] << " ";
		}

		ostr << std::endl;
	}

	ostr << std::endl;
}

void trapezium_method(double a, double b, std::ostream& ostr)
{
	ostr << "Trapezoid formula" << std::endl;
	ostr << "|N    |h     |Integral |Error estimation|  k   |" << std::endl;

	double error;
	int N = 1;
	double k;
	double inegral_h1;
	double inegral_h2;
	double inegral_h3;
	double h1;
	double h2;
	double h3;
	int begin_count_of_calls_to_the_function;
	int count_call1;
	int count_call2;
	int count_call3;

	begin_count_of_calls_to_the_function = count_of_calls_to_the_function;

	count_call2 = count_of_calls_to_the_function;
	comput_integral_trapezium_method(a, b, N, h2, inegral_h2);
	count_call2 = count_of_calls_to_the_function - count_call2;
	N *= 2;

	count_call3 = count_of_calls_to_the_function;
	comput_integral_trapezium_method(a, b, N, h3, inegral_h3);
	count_call3 = count_of_calls_to_the_function - count_call3;
	N *= 2;

	do
	{
		inegral_h1 = inegral_h2;
		inegral_h2 = inegral_h3;
		h1 = h2;
		h2 = h3;
		count_call1 = count_call2;
		count_call2 = count_call3;

		count_call3 = count_of_calls_to_the_function;
		comput_integral_trapezium_method(a, b, N, h3, inegral_h3);
		count_call3 = count_of_calls_to_the_function - count_call3;

		k = log((inegral_h3 - inegral_h1) / (inegral_h2 - inegral_h1) - 1) / log(0.5);
		error = ((inegral_h2 - inegral_h1) / (pow(h1, k) * (1 - pow(0.5, k)))) * pow(h1, k);

		ostr << std::fixed << std::setprecision(9) << "|" << N / 4 << "|" << h1 << "|" << inegral_h1 << "|" << error << "|" << k << "|" << std::endl;

		N *= 2;
	} while (fabs(error) > EPS);

	ostr << "Result " << inegral_h1 << std::endl;
	ostr << "Number of requests " << count_call1 << std::endl;
	ostr << "All Number of requests " << count_of_calls_to_the_function - begin_count_of_calls_to_the_function << std::endl;
}

void comput_integral_trapezium_method(double a, double b, int N, double& h, double& integral_value)
{
	h = fabs(b - a) / N;
	double sum = 0;
	for (int i = 1; i < N; i++)
	{
		sum += func(a + i * h);
	}
	integral_value = 0;
	integral_value = h * ((func(a) + func(b)) / 2.0 + sum);
}