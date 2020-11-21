#pragma once
#include <ostream>

const int VARIANT = 5;
const double EPS = 1E-8;//0.00000001;

enum TypeMethod {
	trapezoid, trapezoid_modified, simpson, gauss
};

double func(double x);

double compute_m(int number_of_derivative);

double compute_derivative_func(double x);

void allocate_matrix(double**& matrix, int n, int m);

void delete_matrix(double** matrix, int n, int m);

void print_matrix(double** matrix, int n, int m, std::ostream& ostr);

void trapezoid_method(double a, double b, bool modified, std::ostream& ostr);

void simpson_method(double a, double b, std::ostream& ostr);

void gauss_method(double a, double b, std::ostream& ostr);

void compute_integral_trapezium_method(double a, double b, int N, double& h, double& integral_value);

void compute_integral_trapezium_modified_method(double a, double b, int N, double& h, double& integral_value);

void compute_integral_simpson(double a, double b, int N, double& h, double& integral_value);

void compute_integral_gauss(double a, double b, int N, double& h, double& integral_value);