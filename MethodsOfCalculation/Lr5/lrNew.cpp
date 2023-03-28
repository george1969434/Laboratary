#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>


using namespace std;
const double pi = 3.14159265359, eps1 = 0.01, eps2 = 0.000001;

double func1(double x) {
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

double func2(double x) {
	return sqrt(1 + x) - 1;
}

double func3(double x) {
	return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}

double funckrat(double x) {
	return (x-0.4)*(x-0.4) * (x - 0.4);
}






vector<double> Gausse(vector<vector<double>> a, vector<double> b, int N) {
	double c, sum;
	double max, temp;
	int main_el;
	bool flag = true;


	vector<double> x(N, 0);
	vector<vector<double>> aa(N,vector<double>(N));
	vector<double> bb(N, 1);
	//copy(b[0],end(b),)
	aa = a;
	bb = b;
	//cout << bb[0];

	for (int k = 0; k < N; k++) { //ищем строку с максимальльным по модулю элементу
		max = abs(aa[k][k]);
		main_el = k;
		for (int i = k + 1; i < N; i++) {
			if (abs(aa[i][k]) > max) {
				main_el = i;
				max = abs(aa[i][k]);
			}
		}
		if (max < 0.0000000001) {//если значение ведущего элемента близко к нулю, такое уравнение мы не решаем (решение тривиальное)
			flag = false;
		}
		for (int j = 0; j < N; j++) { //переставляем строки в матрице А и В
			temp = aa[k][j];
			aa[k][j] = aa[main_el][j];
			aa[main_el][j] = temp;
		}
		temp = bb[k];// и в 
		bb[k] = bb[main_el];
		bb[main_el] = temp;

	}


	for (int k = 0; k < N; k++) {//приводим к верхнетреугольному виду
		for (int i = k; i < N; i++) {
			double c = aa[i][k];
			if (abs(c) < 0.0000000001) continue;
			for (int j = k; j < N; j++)
				aa[i][j] = aa[i][j] / c;
			bb[i] = bb[i] / c;
			if (i == k) continue;
			for (int j = k; j < N; j++)
				aa[i][j] = aa[i][j] - aa[k][j];
			bb[i] = bb[i] - bb[k];
		}
	}

	for (int i = 0; i < N; i++) {
		if (abs(aa[i][i]) < 0.000000001)
			flag = false;
	}
	
	//for (unsigned int i = 0; i < N; i++) {
	 //   for (unsigned int j = 0; j < N; j++) {
	 //       cout << a[i][j]<<' ';
	 //   }
	 //   cout << endl;
	//}
	//cout << a[2][2] << endl;
	if (flag) {//находим корни
		x[N - 1] = bb[N - 1];
		for (int i = N - 2; i >= 0; i--) {
			sum = 0;
			for (int j = i + 1; j < N; j++) {
				sum += aa[i][j] * x[j];
			}
			x[i] = (bb[i] - sum) / aa[i][i];
		}
	}
	// for (int j = 0; j < N; j++) {
	  //   cout << x[j] << ' ';
	 //}
	// cout << endl;
 //}
	else {
		cout << "no solution except trivial" << endl;
		for (int i = 0; i < N; i++)
			x[i] = 0;
	}
	return x;
}

vector<double> SumVec(vector<double> x, vector<double> y, int n) {
	vector<double> ans(n, 0);
	for (int i = 0; i < n; i++) {
		ans[i] = x[i] + y[i];
	}
	return ans;
}

vector<double> DifVec(vector<double> x, vector<double> y, int n) {
	vector<double> ans(n, 0);
	for (int i = 0; i < n; i++) {
		ans[i] = x[i] - y[i];
	}
	return ans;
}

vector<double> MultiplyMatrixes(vector<double>  m1, vector<double> m2) {
	int n = m1.size();
	vector<double> ans(n, 0);
	for (int i = 0; i < m1.size(); i++) {	
			ans[i] += m1[i] * m2[i];
	}
	return ans;
}


vector<double> func4(double x1, double x2) {
	int n = 2;
	vector<double> ans(n, 0);
	ans[0] = x1 * x1 - x2 * x2 - 15;
	ans[1] = x1 * x2 + 4;
	return ans;
}

vector<double> func5(double x1, double x2) {
	int n = 2;
	vector<double> ans(n, 0);
	ans[0] = x1 * x1 + x2 * x2+x1+x2-8;
	ans[1] = x1 * x1+ x2 * x2 +x1*x2 -7;
	return ans;
}


double Derive_Analitic(double (*f)(double), double x, double eps) {
	return (f(x + eps) - f(x)) / eps;
}

vector<vector<double>> Jacoby(vector<double> x,vector<double>(*f)(double x1,double x2)) {
	int n = x.size();
	vector<vector<double> > v(n,vector<double>(n));

	v[0][0] = (f(x[0]+eps2,x[1])[0]-f(x[0],x[1])[0])/eps2;



	v[0][1] = (f(x[0], x[1]+eps2)[0] - f(x[0], x[1])[0]) / eps2;
	v[1][0] = (f(x[0] + eps2, x[1])[1] - f(x[0], x[1])[1]) / eps2;
	v[1][1] = (f(x[0], x[1]+eps2)[1] - f(x[0], x[1])[1]) / eps2;
	
	return v;
}

vector<vector<double>> JacobyAnal(vector<double> x, vector<double>(*f)(double x1, double x2)) {
	int n = x.size();
	vector<vector<double> > v(n, vector<double>(n));

	v[0][0] = 2*x[0];



	v[0][1] = -2*x[1];
	v[1][0] = x[1];
	v[1][1] = x[0];

	return v;
}

void Localization(double A, double B) {


	ofstream fout;
	fout.open("intervals.txt");
	double h = (B - A) / 31;
	double a1 = A;
	double b1 = B;
	fout << a1 << endl;
	while (a1 < b1) {
		a1 += h;
		fout << a1 << endl;
	}

	fout.close();
}


void Besek(double (*f)(double)) {
	ifstream fin;
	fin.open("intervals.txt");
	double a, b, x;
	double eps = eps2;
	int k = 0;
	fin >> a;
	cout << "Besek:" << endl;
	while (!fin.eof()) {
		fin >> b;
		//cout << b<<endl;
		if ((f(a) * f(b) < 0) || (abs(f(b)) < eps)) {
			x = (a + b) / 2;
			while ((abs(f(x)) >= eps) && (abs((b - a) / 2) >= eps)) {
				k++;
				x = (a + b) / 2;
				if (f(a) * f(x) <= 0) {
					b = x;
				}
				else {
					a = x;
				}

			}
			if (abs(f(a)) < eps)
				x = a;
			cout << "x=" << x << ' ' << "f=" << f(x) << endl;
		}

		a = b;
	}
	cout << "kol:" << k;
	cout << endl;
	fin.close();
}


void Newton(double (*f)(double)) {
	ifstream fin;
	fin.open("intervals.txt");
	double a, b, x;
	double eps = eps2;
	vector<double> xk;
	vector<double> p;
	fin >> a;
	//xk.push_back(a);
	cout << "Newton:" << endl;
	int k = 0;
	double xprev = 0;
	while (!fin.eof()) {
		fin >> b;
		
		x = 0.9;
		if ((f(a) * f(b) < 0) || (abs(f(b)) < eps)) {
			while ((abs(f(x)) >= eps)&&(abs(x-xprev)>eps)) {
				 xprev = x;
				// cout << xprev << ' ';
				 if ((f(x) * Derive_Analitic(f, Derive_Analitic(f, x, eps), eps)) / (Derive_Analitic(f, x, eps) * Derive_Analitic(f, x, eps)) >= 1) {
					 x -= f(x) * (b - x) / (f(b) - f(x));
				 }
				k++;
				
				xk.push_back(x);
				x -= (f(x) / Derive_Analitic(f, x, eps));
				//cout << x << ' ';
			}
			if (abs(f(a)) < eps)
				x = a;
			cout << "x=" << x << ' ' << "f=" << f(x) << endl;
		}
		a = b;
	}
	//for (int i = 0; i < xk.size(); i++) {
		//cout << xk[i] << ' ';
	//}
	
	//cout << abs(xk[0] - x);
	for (int i = 1; i < xk.size()-1; i++) {
		if ((log(abs(xk[i] - x) / abs(xk[i - 1] - x))) >=eps) {
			double f = (log(abs(xk[i + 1] - x) / abs(xk[i] - x))) / (log(abs(xk[i] - x) / abs(xk[i - 1] - x)));
			p.push_back(f);
			
		}
	}
	cout<<"p:";
	for (int i = 0; i <p.size(); i++) {
		cout << p[i] << ' ' ;
	}
	cout << endl;
	//cout << "p=" << f << endl;

	cout << "kol:" << k;
	cout << endl;
	fin.close();

}

void NewtonMatrNum(vector<double>(*f)(double x1, double x2)) {
	vector<double> x(2, 0);
	vector<double> y(2, 1);
	vector<vector<double>> J(2,vector<double>(2));
	vector<double> F(2);
	x[0] = 5;
	x[1] = 5;
	int k = 0;
	while(abs(y[0]+y[1])>eps2) {
		if ((abs(x[1]) > 10) || (abs(x[0]) > 10)) {

			vector<double> n(2, 1);
			x[0] = (x[0] + n[0]) / 2;
			x[1] = (x[1] + n[1]) / 2;

		}
		k++;
		J = Jacoby(x,f);
		F = f(x[0], x[1]);
		
		F[0] = -F[0];
		F[1] = -F[1];
		//cout << J[1][1];
		y=Gausse(J, F, 2);
		//cout << y[0]<<endl;
		x = SumVec(x, y, 2);
	}

	
	cout << "Newton for matrix(numerical):" << endl;
	cout <<"(X): "<< x[0] << " " << x[1]<<endl;
	cout << "kol iter:" << k << endl;

}


void NewtonMatrAnalit(vector<double>(*f)(double x1, double x2)) {
	vector<double> x(2, 0);
	vector<double> y(2, 1);
	vector<vector<double>> J(2, vector<double>(2));
	vector<double> F(2);
	x[0] = 1;
	x[1] = 1;
	int k = 0;
	while (abs(y[0] + y[1]) > eps2) {
		if ((abs(x[1]) > 10) || (abs(x[0]) > 10)) {
			
			vector<double> n(2, 1);
			
			x [0]= (x[0]+n[0])/2;
			x[1] = (x[1] + n[1]) / 2;

		}
		k++;
		J = JacobyAnal(x, f);
		F = f(x[0], x[1]);

		F[0] = -F[0];
		F[1] = -F[1];
		//cout << J[1][1];
		y = Gausse(J, F, 2);
		//cout << y[0]<<endl;
		x = SumVec(x, y, 2);
	}


	cout << "Newton for matrix(analitic):" << endl;
	cout << "(X): " << x[0] << " " << x[1] << endl;
	cout << "kol iter:" << k << endl;

}


int main()
{
	ifstream fin;
	fin.open("test3.txt");
	double a, b, n;

	fin >> a;
	fin >> b;
	Localization(a, b);
	Besek(func3);
	Newton(func3);
	NewtonMatrNum(func5);
	NewtonMatrAnalit(func4);

	fin.close();
	//cout << 0;
}