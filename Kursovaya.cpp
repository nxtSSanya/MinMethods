#pragma comment(linker, "/STACK:16777216")
#include <iomanip> 
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <Windows.h>
#include <ctime>
#include <climits>
#include <map>
#define MAX_HUNGARIAN
#define RELEASE_HUNGARIAN
#define FULL_COUNTING
#define Hungarian
#define ROSENBROCK

using namespace std;
#ifndef INPUT_DATA
int nq;
typedef pair<int, int> PInt;
typedef vector<int> VInt;
typedef vector<VInt> VVInt;
typedef vector<PInt> VPInt;
typedef vector<double> VDouble;
typedef vector<VDouble> VVDouble;

double f_test_1(double x);
double f_test_2(double x);
double f_test_3(double x);
double f_test_4(double x);

vector< vector <int> > Optimize_test_1(int n);
vector< vector <double> > Optimize_test_2(int n);
#endif

#ifndef MAX_HUNGARIAN
vector < vector<int> > a;      // Матрица эффективности a[разраб][задача]
vector<int> xy, yx;             // Паросочетания: xy[разраб], yx[задача]
vector<char> vx, vy;            // Альтернирующее дерево vx[разраб], vy[задача]
vector<int> maxrow, mincol, minrow;     // Способности

bool dotry(int i) {
	if (vx[i]) return false;
	vx[i] = true;
	for (int j = 0; j < n; ++j)
		if (a[i][j] - maxrow[i] - mincol[j] == 0)
			vy[j] = true;
	for (int j = 0; j < n; ++j)
		if (a[i][j] - maxrow[i] - mincol[j] == 0 && yx[j] == -1) {
			xy[i] = j;
			yx[j] = i;
			return true;
		}
	for (int j = 0; j < n; ++j)
		if (a[i][j] - maxrow[i] - mincol[j] == 0 && dotry(yx[j])) {
			xy[i] = j;
			yx[j] = i;
			return true;
		}
	return false;
}
#endif

#ifdef HUNGARIAN_ALGORITHM
vector<int> hugarian(vector<vector <int> > a) {
	vector<int> u(nq+ 1), v(nq + 1), p(nq + 1), way(nq + 1);
	for (int i = 1; i <= nq; ++i) {
		p[0] = i;
		int j0 = 0;
		vector<int> minv(nq + 1, 1 << 30);
		vector<char> used(nq + 1, false);
		do {
			used[j0] = true;
			int i0 = p[j0], delta = 1<<30, j1;
			for (int j = 1; j <= nq; ++j)
				if (!used[j]) {
					int cur = a[i0][j] - u[i0] - v[j];
					if (cur < minv[j])
						minv[j] = cur, way[j] = j0;
					if (minv[j] < delta)
						delta = minv[j], j1 = j;
				}
			for (int j = 0; j <= nq; ++j)
				if (used[j])
					u[p[j]] += delta, v[j] -= delta;
				else
					minv[j] -= delta;
			j0 = j1;
		} while (p[j0] != 0);
		do {
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;
		} while (j0);
	}
	vector<int> ans(nq + 1);
	for (int j = 1; j <= nq; ++j)
		ans[p[j]] = j;
	return ans;
}
#endif
int GetNumbersAfterDot(double eps) {
	return -1 * log10(eps);
}

void FullCounting(double eps, double a, double b, double f(double), int i) {
	vector<double> v;
	int n;
	double x_res = 0, f_min = 10000.0;
	cout << "Enter n(n>0): ";
	la:
	cin >> n;
	if (n <= 1) {
		cout << "This number of partition cant be solved."; exit(1);
	}
	double h = eps;//(b - a) / (double)n;
	for (i; i < n; ++i) {
		double x = a + h * i;
		v.push_back(f(x));
		if (f_min > f(x)) { x_res = x; f_min = f(x); }
	}
	cout.precision(GetNumbersAfterDot(eps));
	cout << fixed << "X* is: " << x_res << endl;
	cout << fixed << "F(x*) is: " << *min_element(v.begin(), v.end());
}
#ifndef RELEASE_HUNGARIAN
void VengerskiyAlgorithm() {
	cin >> n;
	a.resize(n);
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			int qq;
			cin >> qq;
			a[i].push_back(qq);
		}
	}
	mincol.assign(n, 0);
	minrow.assign(n, 0);
	maxrow.assign(n, 0);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			maxrow[i] = max(maxrow[i], a[i][j]);

	xy.assign(n, -1);
	yx.assign(n, -1);
	for (int c = 0; c < n; ) {
		vx.assign(n, 0);
		vy.assign(n, 0);
		int k = 0;
		for (int i = 0; i < n; ++i)
			if (xy[i] == -1 && dotry(i))
				++k;
		c += k;
		if (k == 0) {
			int z = 1 << 30; // INF
			for (int i = 0; i < n; ++i)
				if (vx[i])
					for (int j = 0; j < n; ++j)
						if (!vy[j])
							z = min(z, maxrow[i] + mincol[j] - a[i][j]);
			for (int i = 0; i < n; ++i) {
				if (vx[i]) maxrow[i] -= z;
				if (vy[i]) mincol[i] += z;
			}
		}
	}

	int ans = 0;
	for (int i = 0; i < n; ++i)
		ans += a[i][xy[i]];
	printf("%d\n", ans);
	for (int i = 0; i < n; ++i)
		printf("%d ", xy[i] + 1);
}    
#endif
VPInt hungarian(const VVInt& matrix) {

	int height = matrix.size(), width = matrix[0].size();


	VInt u(height, 0), v(width, 0);


	VInt markIndices(width, -1);

	for (int i = 0; i < height; i++) {
		VInt links(width, -1);
		VInt mins(width, 1 << 30);
		VInt visited(width, 0);
		int markedI = i, markedJ = -1, j;	
		while (markedI != -1) {
			j = -1;
			for (int j1 = 0; j1 < width; j1++)
				if (!visited[j1]) {
					if (matrix[markedI][j1] - u[markedI] - v[j1] < mins[j1]) {
						mins[j1] = matrix[markedI][j1] - u[markedI] - v[j1];
						links[j1] = markedJ;
					}
					if (j == -1 || mins[j1] < mins[j])
						j = j1;
				}
			int delta = mins[j];
			for (int j1 = 0; j1 < width; j1++)
				if (visited[j1]) {
					u[markIndices[j1]] += delta;
					v[j1] -= delta;
				}
				else {
					mins[j1] -= delta;
				}
			u[i] += delta;
			visited[j] = 1;
			markedJ = j;
			markedI = markIndices[j];
		}
		for (; links[j] != -1; j = links[j])
			markIndices[j] = markIndices[links[j]];
		markIndices[j] = i;
	}

	VPInt result;
	for (int j = 0; j < width; j++)
		if (markIndices[j] != -1)// a b
			result.push_back(PInt(markIndices[j], j));
	return result;
}

VPInt hungarian_double(const VVDouble& matrix) {

	int height = matrix.size(), width = matrix[0].size();


	VDouble u(height, 0), v(width, 0);


	VDouble markIndices(width, -1);

	for (int i = 0; i < height; i++) {
		VDouble links(width, -1);
		VDouble mins(width, static_cast<double>(1 << 30));
		VDouble visited(width, 0);
		int markedI = i, markedJ = -1, j;
		while (markedI != -1) {
			j = -1;
			for (int j1 = 0; j1 < width; j1++)
				if (!visited[j1]) {
					if (matrix[markedI][j1] - u[markedI] - v[j1] < mins[j1]) {
						mins[j1] = matrix[markedI][j1] - u[markedI] - v[j1];
						links[j1] = markedJ;
					}
					if (j == -1 || mins[j1] < mins[j])
						j = j1;
				}
			int delta = mins[j];
			for (int j1 = 0; j1 < width; j1++)
				if (visited[j1]) {
					u[markIndices[j1]] += delta;
					v[j1] -= delta;
				}
				else {
					mins[j1] -= delta;
				}
			u[i] += delta;
			visited[j] = 1;
			markedJ = j;
			markedI = markIndices[j];
		}
		for (; links[j] != -1; j = links[j])
			markIndices[j] = markIndices[links[j]];
		markIndices[j] = i;
	}

	VPInt result;
	for (int j = 0; j < width; j++)
		if (markIndices[j] != -1)
			result.push_back(PInt(markIndices[j], j));
	return result;
}

class Start_app {
public:
	Start_app() {	
		system("start C:\\Users\\SanyaBooster\\source\\repos\\Project78\\Release\\Project78.exe");
	}
};


int main() {
	cout << "What method do u wanna use?" << endl;
	cout << "1 - Full counting" << endl;
	cout << "2 - Hungarian Algorithm" << endl;
	cout << "3 - Rosenbrock`s Method" << endl;
	string inputStr;
	cin >> inputStr;
	if (inputStr == "1") {
#ifdef FULL_COUNTING
		double a, b, eps;
		cout << "Enter a, b ";
		cin >> a >> b;
		cout << endl << "Eps: ";
		cin >> eps;
		int k, i;
		cout << "Enter the number of your test function:" << endl;
		cout << "1: x^2 - 1" << endl;
		cout << "2: sin(2x) + 3cos(x-1)" << endl;
		cout << "3: log(x) + 4x - 3sin(x) [x > 0]" << endl;
		cout << "4: x^2 - 2cosh(x)" << endl;
		cout << endl;
		cin >> k;
		switch (k)
		{
		case 1: {
			FullCounting(eps, a, b, f_test_1, i = 0);
			exit(1);
		}
		case 2: {
			FullCounting(eps, a, b, f_test_2, i = 0);
			exit(1);
		}
		case 3: {
			FullCounting(eps, a, b, f_test_3, i = 1);
			exit(1);
		}
		case 4: {
			FullCounting(eps, a, b, f_test_4, i = 0);
			exit(1);
}
		default: {
			cout << "You entered the wrong number" << endl;
			exit(1);
		}
		}
#endif
	}
	if (inputStr == "2") {
#ifdef Hungarian
		vector< vector < int> > test;
		vector< vector <double> > test_1;
		cout << "Enter the matrix size(If you wanna use expample with int - enter 4, If example with double - enter 3, else enter any number and press 3 in next menu): " << endl;
		cin >> nq;
		test.resize(nq);
		test_1.resize(nq);
		int qq;
		cout << "What kind of input data do you want to use?" << endl;
		cout << "1 - int" << endl;
		cout << "2 - double" << endl;
		cout << "3 - I will enter matrix of integer" << endl;
		cin >> qq;
		bool f = 0, f1 = 0;
		switch (qq)
		{
		case 1: {
			test = Optimize_test_1(nq);
			f = 1;
			break;
		}
		case 2: {
			test_1 = Optimize_test_2(nq);
			f1 = 1;
			break;
		}
		case 4: {
			for (int i = 0; i < nq; ++i) {
				for (int j = 0; j < nq; ++j) {
					int qwe;
					cin >> qwe;
					test[i].push_back(qwe);
				}
			}
			f = 1;
			break;
		}
		case 5: {
			for (int i = 0; i < nq; ++i) {
				for (int j = 0; j < nq; ++j) {
					double qwe;
					cin >> qwe;
					test_1[i].push_back(qwe);
				}
			}
		}
		default: {
			cout << "You entered the wrong number" << endl;
			return 0;
		}
		}
		cout << endl;
		if (f == 1) {
			vector<pair<int, int> > ans = hungarian(test);
			int res = 0;
			cout << "Coordinates of optimized elements in matrix: " << endl;
			for (auto& i : ans) {
				cout << i.first << "	" << i.second << endl;
				res += test[i.first][i.second];
			}
			cout << endl << "ANSWER IS: " << res << endl;
		}
		if (f1 == 1) {
			vector<pair<int, int> > ans = hungarian_double(test_1);
			int res = 0;
			cout << "Coordinates of optimized elements in matrix: " << endl;
			for (auto& i : ans) {
				cout << i.first << "	" << i.second << endl;
				res += test_1[i.first][i.second];
			}
			cout << endl << "ANSWER IS: " << res << endl;
		}
#endif
	}
	if (inputStr == "3") {
#ifdef ROSENBROCK
		cout << "Enter the start point(be careful, this method uses a lot of memory and processor cache)." << endl << "PRESS ANY KEY TO CONTINUE." << endl;
		cin.get();
		Start_app();
#endif
	}
	else {
		cout << "You entered the wrong number of method" << endl;
		exit(1);
	}
	return 0;
}
