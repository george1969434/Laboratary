#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


using namespace std;
typedef double T1;
const double eps = 10e-4;
//const ofstream fout;



double Neviaska(T1** a, T1* b, T1* x, int N) {

    T1 sum, otv;
    T1* nev = new T1[N];
    otv = 0;
    for (int i = 0; i < N; i++) {
        sum = 0;
        for (int j = 0; j < N; j++) {
            sum += a[i][j] * x[j];
        }
        nev[i] = abs(b[i] - sum);
    }

    for (int i = 0; i < N; i++) {
        otv += nev[i] * nev[i];
    }

    delete[] nev;
    if (abs(sum) > 0.000000001)
        return sqrt(otv);
    else
        return -1;
}

T1 Norm(T1** c, int N) {
    T1 max = 0, sum;
    bool flag = true;
    for (int i = 0; i < N; i++) {
        sum = 0;
        for (int j = 0; j < N; j++) {// вычисляем кубическую норму
            sum += abs(c[i][j]);
        }
        if (max < sum)
            max = sum;
    }
    return max;
}

T1 Norm_x(T1* x1, T1* x2, int N) {
    vector<T1> v;
    for (int k = 0; k < N; k++) {
        v.push_back(abs(x1[k] - x2[k]));
    }
    return *max_element(v.begin(), v.end());
}

void SimpleIteration(T1** a, T1* b, T1* x, int N) {
    //ofstream fout;
    //fout.open("result.txt");
    bool flag = true, flag2 = true;
    T1 sum=0,  norm_x = 100, kol = 0, max = 0;
    T1 xk;
    T1 tau = 0.0001;
    T1** c = new T1 * [N];
    T1* y = new T1[N];
    T1* x1 = new T1[N];

    vector<T1> v;

    for (int i = 0; i < N; i++) {
        c[i] = new T1[N];
    }
    

    for (int j = 0; j < N; j++) {
        sum = 0;
        for (int i = 0; i < N; i++) {
            if (i!=j) {
                sum+=a[i][j];
            }
            if (max < sum)
                max = sum;
        }
    }

    //tau = 2 / 200;
    //cout << max;
    for (int i = 0; i < N; i++) {
        tau = 1 / a[i][i];
        for (int j = 0; j < N; j++) {//заполняем матрицу С
            if (i == j)
                c[i][j] = -tau * a[i][j] + 1;
            else
                c[i][j] = -tau * a[i][j];
        }
    }
    cout << "Matrix C" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << c[i][j] << " ";
        }
        cout << endl;
    }

    for (int i = 0; i < N; i++) {
        tau = 1 / a[i][i];
        y[i] = tau * b[i];
    }
    cout << "Vector y" << endl;
    for (int i = 0; i < N; i++) {

        cout << y[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < N; i++) {
        x[i] = 5;
        x1[i] = 10000;

    }





    if (Norm(c, N) >= 1)
        flag = false;
        T1 norm = Norm(c, N);
    cout << norm << endl;
    if (flag) {
        while (kol<13) {
            if (Norm_x(x, x1, N) < (1-norm)/norm *eps)
                flag2 = false;
            for (int i = 0; i < N; i++)
                x[i] = x1[i];
            for (int i = 0; i < N; i++) {
                x1[i] = 0;
                for (int j = 0; j < N; j++) {
                    x1[i] += c[i][j] * x[j];
                }
                x1[i] += y[i];
            }
            if (kol == 0)
                max = Norm_x(x, x1, N);
            kol++;
        }

        cout << "Norm C: " << endl;
        // cout << x[N-1]- (c[N-1][N-1] * x[N-1] + y[N-1]) << endl;
        cout << norm << endl;
        cout << "solve of matrix equation of simple iteration:" << endl;
        for (int i = 0; i < N; i++)
            cout << x1[i] << " ";
        cout << "kol=" << kol << endl;
        cout << "kest: " << endl;
        cout << log((1-norm)/ max *eps)/log(norm)<< endl;
        cout << " Neviaska: "<<Neviaska(a,b,x1,N) << endl;
        for (int i = 0; i < N; i++)
            x[i] = 10000;
       // cout << log((1 - norm) * eps / Norm_x(x, x1, N))/log(norm) << endl;
    }
    else
        cout << "no solution" << endl;

    //cout << Neviaska(a, b, x1, N) << endl;
    for (int i = 0; i < N; i++) {
        delete[] c[i];
    }

    delete[] c;

    // fout.close();
    delete[] x1;
    delete[] y;
}



void Jacoby(T1** a, T1* b, int N) {
    T1* x1 = new T1[N];
    T1* x2 = new T1[N];
    T1 max=0;
    T1** c = new T1 * [N];
    T1* y = new T1[N];
    bool flag = true, flag2 = true;
    T1 kol = 0;

    for (int i = 0; i < N; i++) {
        c[i] = new T1[N];
    }
    for (int i = 0; i < N; i++) {
        x1[i] = 0;
        x2[i] = b[i];
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c[i][j] = -a[i][j] / a[i][i];
            y[i] = b[i] / a[i][i];
        }
    }

    for (int i = 0; i < N; i++)
        c[i][i] = 0;

    cout << "Matrix C" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << c[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Vector Y" << endl;
    for (int i = 0; i < N; i++) {
        cout << y[i] << ' ';
    }
    cout << endl;

    T1 norm = Norm(c, N);
    if (norm > 1)
        flag = false;
    cout << norm << endl;

    int iterCount = 0;
    while (kol<8) {
        //cout << Norm_x(x2, x1, N) << endl;
        iterCount++;
        if (Norm_x(x2, x1, N) < eps)
        flag2 = false;
        for (int i = 0; i < N; i++)
            x2[i] = x1[i];
        for (int i = 0; i < N; i++) {
            x1[i] = 0;
            for (int j = 0; j < N; j++) {
                x1[i] += c[i][j] * x2[j];
            }
            x1[i] += y[i];
            //cout << x2[i] << ' ';
//if (Neviaska(a, b, x1, N) < eps)
            //((1 - norm) / norm)*

        }

        kol++;
        if (kol == 1)
            max = Norm_x(x2, x1, N);
            cout << "kest: " << endl;
        cout << log((1 - norm) / max * eps) / log(norm) << endl;
    }

    cout << "Jacoby:" << endl;
    for (int i = 0; i < N; i++)
        cout << x1[i] << " ";
    cout << " kol=" << kol;
    cout << endl;
    cout << " Neviaska: " << Neviaska(a, b, x1, N) << endl;
    cout << "kest: " << endl;
    cout << log((1 - norm) / max * eps) / log(norm) << endl;
    for (int i = 0; i < N; i++) {
        delete[] c[i];
    }
    delete[] c;

    delete[] x2;
    delete[] x1;
}


int main()
{

    ifstream fin;
    fin.open("test2.txt");

    int N;
    fin >> N;
    T1* x = new T1[N];
    T1** a = new T1 * [N];
    T1* b = new T1[N];
    for (int i = 0; i < N; i++) {
        a[i] = new T1[N];
    }




    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fin >> a[i][j];
        }
    }



    for (int j = 0; j < N; j++) {
        fin >> b[j];
    }


    SimpleIteration(a, b, x, N);
    cout << endl;
    Jacoby(a, b, N);


    for (int i = 0; i < N; i++) {
        delete[] a[i];
    }
    delete[] a;




    delete[] b;
    delete[] x;


}

