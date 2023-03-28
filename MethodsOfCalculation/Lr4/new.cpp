#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>

using namespace std;
const double pi = 3.14159265359;

double f1(double x) {
	if ((x <= 1) and (x >= -1))
		return x * x;
}

double f2(double x) {
	if ((x <= 1) and (x >= -1))
		return 1 / (1 + x * x);
}
double f3(double x) {
	if ((x <= 2) and (x >= 0))
		return exp(x);
}

//double f4(double x) {
//	if ((x <= 3) and (x >= -3))
//		return atan(1 + 10*x * x);
//}

//double f4(double x) {
//	return 1;
//}

//double f4(double x) {
//	return x * x;
//}

double f4(double x) {
	return x;
}


//double f4(double x) {
//	 return 1/(1+x*x);
//}



void UniformGrid(double a, double b, int n) {
	vector<double> x(n + 1, 0);
	double h = (b - a) / n;
	ofstream out;
	out.open("result.txt");

	for (int i = 0; i < n + 1; i++) {
		x[i] = a + h * i;
		out << x[i] << "   " << f4(x[i]) << endl;
	}

	out.close();
}


void ChebGrid(double a, double b, int n) {
	vector<double> x(n + 1, 0);

	ofstream out;
	out.open("result.txt");

	for (int i = 0; i < n + 1; i++) {
		x[i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * pi / (2 * (n + 1)));
		out << x[i] << "   " << f4(x[i]) << endl;
	}

	out.close();
}


double PolynomLagr(double X, int n) {
	vector<double> x(n + 1, 0);
	vector<double> y(n + 1, 0);
	vector<double> c(n + 1, 1);
	ifstream fin;
	//double h = 0.2;
	fin.open("result.txt");

	for (int i = 0; i < n + 1; i++) {
		fin >> x[i];
		fin >> y[i];
	}


	for (int k = 0; k < n + 1; k++) {
		for (int j = 0; j < n + 1; j++) {
			if (j != k)
				c[k] *= (X - x[j]) / (x[k] - x[j]);//определяем каждую функцию с
		}
	}

	double sum = 0;
	for (int k = 0; k < n + 1; k++)
		sum += c[k] * y[k];


	double max = 0;
	//for (int i = 0; i < n + 1; i++) {

		//if (y[i] > max)//исправить  тестовый файл
			//max = y[i];
	//}
	//max = max * pow(h, n);

	//cout <<"estimation of exp: "<< max << endl;
	fin.close();
	return sum;
}


double Splain3(double X, int n) {
	vector<double> x(n + 1, 1);
	vector<double> y(n + 1, 1);
	vector<double> c(n + 1, 1);
	vector<double> a(n + 1, 0);
	vector<double> b(n + 1, 0);
	vector<double> d(n + 1, 0);
	vector<double> g(n + 1, 0);
	vector<double> h(n + 1, 1);
	vector<double> alp(n + 1, 1);
	vector<double> bet(n + 1, 1);

	c[0] = 0;
	c[n] = 0;


	ifstream fin;
	fin.open("result.txt");
	double sum = 0;
	for (int i = 0; i < n + 1; i++) {
		fin >> x[i];
		fin >> y[i];
	}

	h[0] = 0;
	a[0] = y[0];
	g[0] = 0;

	for (int i = 1; i < n + 1; i++) {
		h[i] = x[i] - x[i - 1];
		//a[i] = y[i-1];
		g[i] = (y[i] - y[i - 1]) / h[i];



		//b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
		//d[i] = (c[i + 1] - c[i]) / 3*h[i];

	}


	alp[2] = h[2] / -(2 * (h[1] + h[2]));
	bet[1] = -3 * (g[1] - g[0]) / -(2 * (h[0] + h[1]));
	for (int i = 1; i < n; i++) {
		alp[i + 1] = h[i] / (-(2 * (h[i - 1] + h[i])) - h[i - 1] * alp[i]);
		bet[i + 1] = (-3 * (g[i] - g[i - 1]) + h[i - 1] * bet[i]) / (-(2 * (h[i - 1] + h[i])) - h[i - 1] * alp[i]);
	}

	//cout << alp[1] << ',' << bet[1];
	c[n] = 2;
	//c[n+1] = 0;
	c[0] = 0;

	for (int i = n-1 ; i > 0; i--) {
		c[i] = alp[i + 1] * c[i + 1] + bet[i + 1];
	}

	double A = 0, B = 0, C = 0, D = 0, z = 0;
	//возможно поменять счетчик цикла
	for (int i = 1; i < n; i++) {
		b[i] = g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
		a[i] = y[i - 1];
		if (X <= x[i]) {
			A = a[i];
			B = b[i];
			C = c[i];
			D = d[i];
			z = x[i - 1];
			break;
		}

	}

	

	sum = A + B * (X - z) + C * (X - z) * (X - z) + D * (X - z) * (X - z) * (X - z);

	//cout << a[7]<<' '<<b[7]<<' '<<c[7]<<' '<<d[7];

	fin.close();
	return sum;
}





int main()
{
	ifstream fin;
	fin.open("test.txt");
	double a, b, x;
	int n;
	fin >> a;
	fin >> b;
	fin >> n;

	UniformGrid(a, b, n);
	//ChebGrid(a, b, n);

	vector<double> SupportGrid;
	double xcur = a;
	for (int i = 0; i < 512; i++) {
		SupportGrid.push_back(xcur);
		xcur += (b - a) / 512;
	}

	//    double MaxDiff = -1e9;
	//    for (double h : SupportGrid){
	//        double interpolVal = PolynomLagr(h,n);
	//        double funcVal = f4(h);
	//        cout << "h " << h << " interpolVal " << interpolVal << endl;
	//        double diff = abs(interpolVal - funcVal);
	//        MaxDiff = (MaxDiff > diff)?MaxDiff:diff;
	//    }
	double MaxDiff = -1e9;
	for (double h : SupportGrid) {
		double interpolVal = Splain3(h, 3);
		double funcVal = f4(h);
		if (abs(interpolVal) < 0.00000001) break;
		cout << "h " << h << " interpolVal " << interpolVal << endl;
		double diff = abs(interpolVal - funcVal);
		MaxDiff = (MaxDiff > diff) ? MaxDiff : diff;
	}
	//double f = Splain3(0.33, 3);
	//cout << f;

	cout << "MaxDiff " << MaxDiff << endl;
	//ChebGrid(a, b, n);
	//double f=PolynomLagr(0.38,n);
	//cout << f;

	//double k = Splain3(-0.35, n);
	//cout << k;

	fin.close();
	return 0;
}

