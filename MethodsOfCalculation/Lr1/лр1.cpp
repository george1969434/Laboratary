#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;
typedef float T1;




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
     cout << "no solution except trivial" << endl;
     for (int i = 0; i < N; i++)
         x[i] = 0;
 }

 for (int i = 0; i < N; i++) {
     delete[] aa[i];
 }
 delete[] aa;


 delete[] bb;
}


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


T1 Cond_Number(T1** a,  int N) {
    T1 norm, rev_norm,sum,maximum,maximum_,maximum_2,maximum2,max,temp,temp2;
    bool flag;
    int main_el;
    norm = 0;
    sum = 0;
    maximum = 0;
    maximum_ = 0;
    maximum2 = 0;
    maximum_2 = 0;
    rev_norm = 0;
    maximum = a[0][0];
    T1* e = new T1[N];
    T1* x = new T1[N];
    T1** aa = new T1 * [N];
    T1** rev_a = new T1 * [N];
    T1** ed = new T1 * [N];


    for (int i = 0; i < N; i++) {
        aa[i] = new T1[N];
    }

    for (int i = 0; i < N; i++) {
        ed[i] = new T1[N];
    }

    for (unsigned int i = 0; i < N; i++) {
        rev_a[i] = new T1[N];
    }

    for (int i = 0; i < N; i++) {
        sum = 0;
        for (int j = 0; j < N; j++) {// вычисляем кубическую норму
            sum += abs(a[i][j]);
        }
        if (maximum < sum)
            maximum = sum;
    }

    for (int j = 0; j < N; j++) {
        sum = 0;
        for (int i = 0; i < N; i++) {// вычисляем октаэдерическую норму
            sum += abs(a[i][j]);
        }
        if (maximum_ < sum)
            maximum_ = sum;
    }



    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            aa[i][j] = a[i][j];
        }
    }



    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j == i)
                e[j] = 1;
            else
                e[j] = 0;
        }

        Gausse(aa, e, x, N);
        for (int j = 0; j < N; j++) {
            rev_a[j][i] = x[j];
        }




    }

    maximum2 = 0;

    for (int i = 0; i < N; i++) {
        sum = 0;
        for (int j = 0; j < N; j++) {// вычисляем кубическую норму обратной матрицы
            sum += abs(rev_a[i][j]);
        }
        if (maximum2 < sum)
            maximum2 = sum;
    }

    for (int j = 0; j < N; j++) {
        sum = 0;
        for (int i = 0; i < N; i++) {// вычисляем октаэдерическую  норму обратной матрицы
            sum += abs(rev_a[i][j]);
        }
        if (maximum_2 < sum)
            maximum_2 = sum;
    }



        for (int k = 0; k < N; k++) {
            for (int i = 0; i < N; i++) {
                ed[k][i] = 0;
                for (int j = 0; j < N; j++) {
                    ed[k][i] += (rev_a[k][j] * a[j][i]);
                }

            }

        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (ed[i][j] < 0.00000000001)
                    ed[i][j] = 0;
                cout << ed[i][j] << ' ';
            }
            cout << endl;
        }


        for (int i = 0; i < N; i++) {
            delete[] aa[i];
        }
        delete[] aa;

        for (int i = 0; i < N; i++) {
            delete[] rev_a[i];
        }
        delete[] rev_a;


        for (int i = 0; i < N; i++) {
            delete[] ed[i];
        }
        delete[] ed;

        delete[] e;
        delete[] x;

        cout << "Conditional number for octaider norm: " << endl;
        cout << maximum_ * maximum_2 << endl;

        cout << "Conditional number for cubic norm :" << endl;
    return maximum * maximum2;
}



int main() {

    ifstream fin;
    fin.open("test1.txt");
    int N;
    fin >> N;
    T1* x = new T1[N];
    T1** a = new T1 * [N];

    for (int i = 0; i < N; i++) {
        a[i] = new T1[N];
    }
    T1* b = new T1[N];



    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fin >> a[i][j];
        }
    }



    for (int j = 0; j < N; j++) {
        fin >> b[j];
    }

    cout << "Our matrix:" << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }


    Gausse(a, b, x, N);
    cout << "solve of matrix equation:" << endl;
    for (int i = 0; i < N; i++)
        cout << x[i] << " ";
    cout << endl;
    cout << "vector of neviaska: " << endl;
    cout << Neviaska(a, b, x, N) << endl;

    cout << "Multipline of matrix: " << endl;
    cout << Cond_Number(a, N) << endl;

    //cout << a[0][0] << endl;

    fin.close();
    
    for (int i = 0; i < N; i++) {
        delete[] a[i];
    }
    delete[] a;




    delete[] b;


    return 0;

}



//  
//for (int j = k; j < N; j++) {
 //   for (int i = k + 1; j < N; j++) {
  //      c = a[i][k] / a[k][k];
   //     a[i][j] = a[i][j] - c * a[k][j];
   //     b[i] = b[i] - c * b[k];
 //   }
//}


//for (int i = k; i < N; i++) {
//double c = a[i][k];
//if (abs(c) > 10 ^ (-10)) continue;
//for (int j = 0; j < N; j++)
  //  a[i][j] = a[i][j] / c;
//b[i] = b[i] / c;
//if (i == k) continue;
//for (int j = 0; j < N; j++)
  //  a[i][j] = a[i][j] - a[k][j];
//b[i] = b[i] - b[k];     }
