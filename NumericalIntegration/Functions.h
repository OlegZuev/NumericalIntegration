#pragma once
#include <ostream>

const int VARIANT = 0;
const double EPS = 10E-8;//0.00000001;

double func(double x);

double compute_m(int number_of_derivative);

double compute_derivative_func(double x);

void allocate_matrix(double**& matrix, int n, int m);

void delete_matrix(double** matrix, int n, int m);

void print_matrix(double** matrix, int n, int m, std::ostream& ostr);

void trapezium_method(double a, double b, std::ostream& ostr);

void comput_integral_trapezium_method(double a, double b, int N, double& h, double& integral_value);
