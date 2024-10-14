#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double u(double x, double t) {
	return -pow(x, 4) + x + t * x + t * t - t * exp(x);
}

double d2u_dx2(double x, double t) {
	return -12 * pow(x, 2) - t * exp(x);
}

double mu(double x) {
	return -pow(x, 4) + x;
}

double mu1(double t) {
	return pow(t, 2) - t;
}

double mu2(double t) {
	return t + pow(t, 2) - t * exp(1);
}

double f(double x, double t) {
	return 0.027 * 12 * x * x + x + 2 * t - exp(x) + 0.027 * t * exp(x);
}

void sweep_method(int N, double A, double B, double* vector_B, double* vector_X) {

	//straight stroke
	double* vector_V = new double[N];
	double* vector_U = new double[N];

	vector_V[1] = -(B) / (A);
	vector_U[1] = vector_B[1] / A;

	for (int i = 2; i < N - 1; i++) {
		vector_V[i] = (B) / (-A - (B * vector_V[i - 1]));
		vector_U[i] = (B * vector_U[i - 1] - vector_B[i]) / (-A - B * vector_V[i - 1]);
	}

	vector_V[N - 1] = 0;
	vector_U[N - 1] = (B * vector_U[N - 2] - vector_B[N - 1]) / (-A - B * vector_V[N - 2]);

	//reverse stroke

	vector_X[N - 1] = vector_U[N - 1];
	for (int i = N - 2; i >= 1; i--) {
		vector_X[i] = vector_V[i] * vector_X[i + 1] + vector_U[i];
	}
}

int main() {
	int N, T;
	ofstream file("value.txt");
	ofstream fout("value_function.txt");


	double tao = 0.1;
	double h = 0.1;
	double a = 0.027;

	N = 1 / h;
	T = 1 / tao;


	double A = (T + 2 * a * pow(N, 2));
	double B = (-a * pow(N, 2));
	double* vector_b = new double[N + 1];
	double delta = 0;
	double* array = new double[N + 1];

	for (int j = 0; j <= N; j++) {
		if (j == 0) {
			for (int i = 0; i <= N; i++) {
				array[i] = mu(i * h);
				file << array[i] << " ";
				fout << u(h * i, j * tao) << " ";
			}
			file << endl;
			fout << endl;
			continue;
		}
		array[0] = mu1(j * tao);
		array[N] = mu2(j * tao);

		//вычисляем вектор b: Ax=b
		vector_b[1] = f(h, tao * j) + T * array[1] + a * pow(N, 2) * array[0];
		for (int i = 2; i < N - 1; i++) {
			vector_b[i] = f(h * i, tao * j) + T * array[i];
		}
		vector_b[N - 1] = f(h * N - 1, tao * j) + T * array[N - 1] + a * pow(N, 2) * array[N];
		sweep_method(N, A, B, vector_b, array);
		double curr;
		double aprx;
		double var;
		for (int i = 0; i <= N; i++) {
			curr = u(h * i, j * tao);
			aprx = array[i];
			fout << curr << " ";
			file << aprx << " ";
			var = abs(curr - aprx);
			if (var > delta) {
				delta = var;
			}
		}
		file << endl;
		fout << endl;

	}


	file.close();
	fout.close();
	cout << delta << endl;
	cout << fixed << setprecision(10) << tao + pow(h, 2) << endl;
	return 0;
}