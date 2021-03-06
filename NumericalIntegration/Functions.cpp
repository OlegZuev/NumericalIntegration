#include "Functions.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int count_of_calls_to_the_function = 0;

/*
 * func for task
 */
double func(double x) {
	count_of_calls_to_the_function++;

	switch (VARIANT) {
	case 0:
		return pow(2, x) - 3 * x + 2;
	case 5:
		return pow(3, x - 1) + 2 - x;
	case 7:
		return exp(-2 * x) - 2 * x + 1;
	}

	throw std::runtime_error("Incorrect variant");
}

/*
 * compute first derivative func
 */
double compute_derivative_func(double x) {
	switch (VARIANT) {
	case 0:
		return pow(2, x) * log(2) - 3;
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
void allocate_matrix(double**& matrix, int n, int m) {
	matrix = new double*[n];
	for (int i = 0; i < n; ++i) {
		matrix[i] = new double[m]{0};
	}
}

/*
 * delete matrix n * m
 * n - count of rows
 * m - count of columns
 */
void delete_matrix(double** matrix, int n, int m) {
	for (int i = 0; i < n; ++i) {
		delete[] matrix[i];
	}

	delete[] matrix;
}

/*
 * print matrix n*m into ostr
 */
void print_matrix(double** matrix, int n, int m, std::ostream& ostr) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j) {
			ostr << std::fixed << std::setprecision(6) << matrix[i][j] << " ";
		}

		ostr << std::endl;
	}

	ostr << std::endl;
}

/*
 * trapezoid method for finding the value of the quadrature form of the integral
 */
void trapezoid_method(double a, double b, bool modified, std::ostream& ostr) {
	void (*compute_integral)(double a, double b, int N, double& h, double& integral_value) = nullptr;
	double coef = 1 / 3.0;
	if (modified) {
		ostr << "Trapezoid modified formula" << std::endl;
		compute_integral = compute_integral_trapezium_modified_method;
	} else {
		ostr << "Trapezoid formula" << std::endl;
		compute_integral = compute_integral_trapezium_method;
	}

	ostr << "|N    |h     |Integral |Error estimation      |  k   |" << std::endl;

	double error;
	int N = 1;
	double k;
	double integral_h0;
	double integral_h1 = 0;
	double integral_h2;
	double integral_h3;
	double h1;
	double h2;
	double h3;
	int begin_count_of_calls_to_the_function;
	int count_call1;
	int count_call2;
	int count_call3;

	begin_count_of_calls_to_the_function = count_of_calls_to_the_function;

	count_call2 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h2, integral_h2);
	count_call2 = count_of_calls_to_the_function - count_call2;
	N *= 2;

	count_call3 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h3, integral_h3);
	count_call3 = count_of_calls_to_the_function - count_call3;
	N *= 2;

	do {
		integral_h0 = integral_h1;
		integral_h1 = integral_h2;
		integral_h2 = integral_h3;
		h1 = h2;
		h2 = h3;
		count_call1 = count_call2;
		count_call2 = count_call3;

		count_call3 = count_of_calls_to_the_function;
		compute_integral(a, b, N, h3, integral_h3);
		count_call3 = count_of_calls_to_the_function - count_call3;

		k = log((integral_h3 - integral_h1) / (integral_h2 - integral_h1) - 1) / log(0.5);
		error = (integral_h1 - integral_h0) * coef;

		if (N > 4) {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h1 << "|" <<
				std::scientific << std::setprecision(15) << error << "|" << std::fixed << std::setprecision(4) << k <<
				"|" << std::endl;
		} else {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h1 << "|" << std::endl;
		}


		N *= 2;
	} while (fabs(error / integral_h1) > EPS);

	ostr << "Result " << std::setprecision(15) << integral_h1 << std::endl;
	ostr << "Number of requests " << count_call1 << std::endl << std::endl;
}

/*
 * simpson method for finding the value of the quadrature form of the integral
 */
void simpson_method(double a, double b, std::ostream& ostr) {
	void (*compute_integral)(double a, double b, int N, double& h, double& integral_value) = compute_integral_simpson;
	double coef = 1 / 15.0;
	ostr << "Simpson formula" << std::endl;
	ostr << "|N    |h     |Integral |Error estimation      |  k   |" << std::endl;

	double error;
	int N = 1;
	double k;
	double integral_h1 = 0;
	double integral_h2;
	double integral_h3;
	double h1;
	double h2;
	double h3;
	int begin_count_of_calls_to_the_function;
	int count_call1;
	int count_call2;
	int count_call3;

	begin_count_of_calls_to_the_function = count_of_calls_to_the_function;

	count_call2 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h2, integral_h2);
	count_call2 = count_of_calls_to_the_function - count_call2;
	N *= 2;

	count_call3 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h3, integral_h3);
	count_call3 = count_of_calls_to_the_function - count_call3;
	N *= 2;

	do {
		integral_h1 = integral_h2;
		integral_h2 = integral_h3;
		h1 = h2;
		h2 = h3;
		count_call1 = count_call2;
		count_call2 = count_call3;

		count_call3 = count_of_calls_to_the_function;
		compute_integral(a, b, N, h3, integral_h3);
		count_call3 = count_of_calls_to_the_function - count_call3;

		k = log((integral_h3 - integral_h1) / (integral_h2 - integral_h1) - 1) / log(0.5);
		error = (integral_h2 - integral_h1) * coef;

		if (N > 4) {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h2 << "|" <<
				std::scientific << std::setprecision(15) << error << "|" << std::fixed << std::setprecision(4) << k <<
				"|" << std::endl;
		} else {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h2 << "|" << std::endl;
		}

		N *= 2;
	} while (fabs(error / integral_h2) > EPS);

	ostr << "Result " << std::setprecision(15) << integral_h2 << std::endl;
	ostr << "Number of requests " << count_call2 << std::endl << std::endl;
}

/*
 * gauss method for finding the value of the quadrature form of the integral
 */
void gauss_method(double a, double b, std::ostream& ostr) {
	void (*compute_integral)(double a, double b, int N, double& h, double& integral_value) = compute_integral_gauss;
	double coef = 1 / 63.0;
	ostr << "Gauss formula" << std::endl;
	ostr << "|N    |h     |Integral |Error estimation     |  k   |" << std::endl;

	double error;
	int N = 1;
	double k;
	double integral_h0;
	double integral_h1 = 0;
	double integral_h2;
	double integral_h3;
	double h1;
	double h2;
	double h3;
	int begin_count_of_calls_to_the_function;
	int count_call1;
	int count_call2;
	int count_call3;

	begin_count_of_calls_to_the_function = count_of_calls_to_the_function;

	count_call2 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h2, integral_h2);
	count_call2 = count_of_calls_to_the_function - count_call2;
	N *= 2;

	count_call3 = count_of_calls_to_the_function;
	compute_integral(a, b, N, h3, integral_h3);
	count_call3 = count_of_calls_to_the_function - count_call3;
	N *= 2;

	do {
		integral_h0 = integral_h1;
		integral_h1 = integral_h2;
		integral_h2 = integral_h3;
		h1 = h2;
		h2 = h3;
		count_call1 = count_call2;
		count_call2 = count_call3;

		count_call3 = count_of_calls_to_the_function;
		compute_integral(a, b, N, h3, integral_h3);
		count_call3 = count_of_calls_to_the_function - count_call3;

		k = log((integral_h3 - integral_h1) / (integral_h2 - integral_h1) - 1) / log(0.5);
		error = (integral_h1 - integral_h0) * coef;

		if (N > 4) {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h1 << "|" <<
				std::scientific << std::setprecision(15) << error << "|" << std::fixed << std::setprecision(4) << k <<
				"|" << std::endl;
		}
		else {
			ostr << std::fixed << "|" << std::setw(5) << N / 4 << "|" << std::setprecision(4) << h1 << "|" <<
				std::setw(9) << std::setprecision(6) << integral_h1 << "|" << std::endl;
		}

		N *= 2;
	} while (fabs(error / integral_h1) > EPS);

	ostr << "Result " << std::setprecision(15) << integral_h1 << std::endl;
	ostr << "Number of requests " << count_of_calls_to_the_function - count_call2 - count_call3 -
		begin_count_of_calls_to_the_function << std::endl << std::endl;
}

/*
 * finding the value of the quadrature form of the integral via trapezium method
 */
void compute_integral_trapezium_method(double a, double b, int N, double& h, double& integral_value) {
	h = (b - a) / N;
	double sum = 0;
	for (int i = 1; i < N; i++) {
		sum += func(a + i * h);
	}
	integral_value = h * ((func(a) + func(b)) / 2.0 + sum);
}

/*
 * finding the value of the quadrature form of the integral via trapezium modified method
 */
void compute_integral_trapezium_modified_method(double a, double b, int N, double& h, double& integral_value) {
	compute_integral_trapezium_method(a, b, N, h, integral_value);
	integral_value += (h * h) / 12.0 * (compute_derivative_func(a) - compute_derivative_func(b));
}

/*
 * finding the value of the quadrature form of the integral via simpson method
 */
void compute_integral_simpson(double a, double b, int N, double& h, double& integral_value) {
	h = (b - a) / N;
	double sum1 = 0;
	double sum2 = 0;
	for (int i = 1; i < N; i++) {
		if (i & 1) {
			sum1 += func(a + i * h);
		} else {
			sum2 += func(a + i * h);
		}
	}

	integral_value = h * ((func(a) + func(b)) + 4 * sum1 + 2 * sum2) / 3.0;
}

/*
 * finding the value of the quadrature form of the integral via gauss method
 */
void compute_integral_gauss(double a, double b, int N, double& h, double& integral_value) {
	h = (b - a) / N;
	double sum = 0;
	for (int i = 0; i < N; i++) {
		double x = a + (1.0 + 2.0 * i) * h / 2.0;
		sum += h / 2.0 * (
			(5.0 / 9.0) * func(x - sqrt(3.0 / 5.0) * h / 2.0)
			+ (8.0 / 9.0) * func(x)
			+ (5.0 / 9.0) * func(x + sqrt(3.0 / 5.0) * h / 2.0)
		);
	}

	integral_value = sum;
}
