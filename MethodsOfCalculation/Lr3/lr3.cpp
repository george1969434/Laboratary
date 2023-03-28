#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


using namespace std;
typedef double T1;
const double eps = 10e-4;


void Gausse(T1** a, T1* b, T1* x, int N) {
    T1 c, sum;
    T1 max, temp;
    int main_el;
    bool flag = true;
    T1** aa = new T1 * [N];
    T1* bb = new T1[N];

    for (int i = 0; i < N; i++) {
        aa[i] = new T1[N];
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aa[i][j] = a[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        bb[i] = b[i];
    }

    //copy(b[0],end(b),)

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
            T1 c = aa[i][k];
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
        //cout << "no solution except trivial" << endl;
        for (int i = 0; i < N; i++)
            x[i] = 0;
    }

    for (int i = 0; i < N; i++) {
        delete[] aa[i];
    }
    delete[] aa;


    delete[] bb;
}

T1 Norm_x(T1* x1, int N) {
    T1 sum = 0;
    for (int i = 0; i < N; i++) {
        sum += x1[i] * x1[i];
    }
    return sqrt(sum);
}

T1 Mult_xx(T1* x1, T1* x2, int N) {
    T1 mult = 0;
    for (int i = 0; i < N; i++)
        mult += x1[i] * x2[i];
    return mult;
}

void Mult_Ax(T1** a, T1* x2, T1* y, int N) {


    for (int j = 0; j < N; j++)
        y[j] = 0;


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            y[i] += a[i][j] * x2[j];
        }
    }


}

void RevIter(T1** a, T1* x, int N) {
    T1 lambda1 = 0, norm;
    T1** aa = new T1 * [N];
    T1* y = new T1[N];
    T1* y2 = new T1[N];
    
    for (int i = 0; i < N; i++) {
        aa[i] = new T1[N];
    }

    

    for (int i = 0; i < N; i++) {
        y2[i] = 0;
    }

    x[0] = 0.01;
    x[1] = 0.7;
    x[2] = -0.6;
    x[3] = 0.34;

    //for (int i = 0; i < N; i++) {
    //    x[i] = 1/sqrt(N);
    //}
    for (int k = 0; k < 10; k++) {
        Mult_Ax(a, x, y2, N);
        
            lambda1 = Mult_xx(y2, x, N);// рассмотреть. неправильно считает
            //cout<<Mult_xx(Mult_Ax(a, x, N), x, N)<<endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j)
                    aa[i][j] = a[i][j] - lambda1;
                else
                    aa[i][j] = a[i][j];
            }
        }

        Gausse(aa, x, y, N);
        norm = Norm_x(y, N);
        if (norm < eps)
            break;
        for (int i = 0; i < N; i++) {
            x[i] = y[i] / norm;
            //cout << x[i] << ' ';
        }
        //cout << lambda1 << endl;
    }
    cout << "self vector: " << endl;

    for (int i = 0; i < N; i++)
        cout << x[i] << ' ';
    cout << endl;
    cout << "self value: " << lambda1 << endl;
    for (int i = 0; i < N; i++) {
        delete[] aa[i];
    }
    delete[] aa;


    delete[] y2;
    delete[] y;
    delete[] x;
}




int main()
{

    ifstream fin;
    fin.open("test.txt");

    int N;
    fin >> N;
    T1* x = new T1[N];
    T1** a = new T1 * [N];

    for (int i = 0; i < N; i++) {
        a[i] = new T1[N];
    }

    for (int i = 0; i < N; i++) {
        x[i] = 0;
    }




    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fin >> a[i][j];
        }
    }

    RevIter(a, x, N);
    




    for (int i = 0; i < N; i++) {
        delete[] a[i];
    }
    delete[] a;

   
cout << 0 << endl;

    delete[] x;
    
    fin.close();
   

}