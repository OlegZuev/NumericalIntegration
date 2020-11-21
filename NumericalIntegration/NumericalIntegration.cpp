#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include "Functions.h"
using namespace std;

const bool PRINT_TO_FILE = false;


int main()
{
	ofstream fout("../" + to_string(VARIANT) + "_output" + ".txt");
	ostream& out = PRINT_TO_FILE ? fout : cout;

	out << "Variant: " << VARIANT << endl;

	// run_method(1, 2, trapezoid, out);
	// run_method(1, 2, trapezoid_modified, out);
	// run_method(1, 2, simpson, out);
	run_method(1, 2, gauss, out);
}