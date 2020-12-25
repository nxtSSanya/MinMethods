#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

vector< vector <int> > Optimize_test_1(int n) {
	ifstream in("matrix1.txt");
	vector< vector<int> > test1;
	test1.resize(n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			int qwe;
			in >> qwe;
			test1[i].push_back(qwe);
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << test1[i][j] << "\t";
		}
		cout << endl;
	}
	return test1;
}

vector< vector <double> > Optimize_test_2(int n) {
	ifstream in("matrix2.txt");
	vector< vector<double> > test1;
	test1.resize(n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			double qwe;
			in >> qwe;
			test1[i].push_back(qwe);
		}
	}
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout << test1[i][j] << "\t";
		}
		cout << endl;
	}
	return test1;
}
