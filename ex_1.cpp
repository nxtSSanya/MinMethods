#include <iostream>
#include <algorithm>

using namespace std;

double f_test_1(double x) {
	return x * x - 1;
}

double f_test_2(double x) {
	return sin(2 * x) + 3 * cos(x - 1);
}

double f_test_3(double x) {
	if (x <= 0.0) {
		cout << "You entered the wrong input data." << endl;
	}
	return log(x) + 4 * x - 3 * sin(x);
}

double f_test_4(double x) {
	return x * x - 2 * cosh(x);
}
