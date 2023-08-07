#include <iostream>
#include <locale.h>
#include <complex>

using namespace std;

typedef complex<double> cpx; // определение типа комплексного числа

cpx** createCPXMatrix(int& rows, int& cols) { //функция, создающая динамическую матрицу
    std::cout << "Введите количество строк: ";
    std::cin >> rows;
    std::cout << "Введите количество столбцов: ";
    std::cin >> cols;
    cpx** matrix = new cpx * [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new cpx[cols];
        cout << "Введите действительную и мнимую часть элементов строки " << i << " : ";
        for (int j = 0; j < cols; j++) {
             double real, imag;
            cin >> real >> imag;
            matrix[i][j] = cpx(real, imag);
        }
    }
    return matrix;
}
cpx** transpose_matrix(cpx** matrix, int n, int m) { //функция для транспонирования матрицы

    // Создание новой матрицы для результата
    cpx** result = new cpx* [m];
    for (int i = 0; i < m; i++) {
        result[i] = new cpx[n];
    }

    // Транспонирование матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}
cpx** multiply_matrix(cpx** matrix1, cpx** matrix2, int rows1, int cols1, int cols2) { // функкция для перемножения матриц
    // Создание новой матрицы для результата
    cpx** result = new cpx* [rows1];
    for (int i = 0; i < rows1; i++) {
        result[i] = new cpx[cols2];
    }

    // Вычисление произведения матриц
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols2; j++) {
            result[i][j] = 0;
            for (int k = 0; k < cols1; k++) {
                result[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return result;
}

cpx** concatenate_matrix(cpx** matrix1, int rows1, int cols1, cpx** matrix2, int rows2, int cols2) { // функция для склеивания матриц

    int total_cols = cols1 + cols2;
    cpx** result_matrix = new cpx* [rows1];
    // Выделяем память для результирующей матрицы
    for (int i = 0; i < rows1; i++) {
        result_matrix[i] = new cpx[total_cols];
    }

    // Копируем элементы из первой матрицы
    for (int i = 0; i < rows1; i++) {
        for (int j = 0; j < cols1; j++) {
            result_matrix[i][j] = matrix1[i][j];
        }
    }

    // Копируем элементы из второй матрицы
    for (int i = 0; i < rows2; i++) {
        for (int j = 0; j < cols2; j++) {
            result_matrix[i][j + cols1] = matrix2[i][j];
        }
    }
    return result_matrix;
}

cpx* Gauss(cpx** G, int n_A) { // метод Гаусса для решение системы
    int n = n_A;

    cpx* x = new cpx[n]; // выделение памяти под вектор x

    // приведение расширенной матрицы к ступенчатому виду, прямой ход Гаусса
    int rank = 0;  // переменная для хранения ранга матрицы
    for (int j = 0; j < n; j++) {
        int k = rank;
        while (k < n && G[k][j].real() == 0 && G[k][j].imag() == 0) {
            k++;
        }
        if (k == n) {
            continue; // все коэффициенты в столбце j равны нулю, пропускаем этот столбец
        }
        for (int i = j; i <= n; i++) {
            swap(G[rank][i], G[k][i]); // меняем местами k-ю и rank-ю строки
        }
        for (int i = rank + 1; i < n; i++) {
            std::complex<double> coef = G[i][j] / G[rank][j];
            for (int k = j; k <= n; k++) {
                G[i][k] -= coef * G[rank][k];
            }
        }
        rank++; // увеличиваем ранг матрицы, так как нашли очередную ненулевую строку
    }
    // проверка на отсутствие решений
    for (int i = rank; i < n; i++) {
        if (G[i][n].real() != 0 && G[i][n].imag() != 0) {
            cout << "Система не имеет решений\n";
            return 0;
        }
    }
    // обратный ход метода Гаусса
    for (int i = n - 1; i >= rank; i--) {
        if (G[i][n].real() != 0 && G[i][n].imag() != 0) {
            cout << "Система не имеет решений\n";
            return 0;
        }
    }
    // вывод решения системы 
    if (rank < n) {
        cout << "Система имеет бесконечное число решений\n";
        return 0;
    }

    else {
        // обратный ход метода Гаусса
        for (int i = n - 1; i >= 0; i--) {
            std::complex<double> sum = 0;
            for (int j = i + 1; j < n; j++) {
                sum += G[i][j] * x[j];
            }
            x[i] = (G[i][n] - sum) / G[i][i];
        }
    }
    // вывод решения системы
    cout << "Решение системы:\n";
    for (int i = 0; i < n; i++) {
        cout << "I[" << i << "] = " << x[i] << "\n";
    }
}

int main() {
    setlocale(LC_ALL, "ru");
    int n = 0, m = 0;
    std::cout << "Ввод матрицы контуров C" << endl;
    cpx** C = createCPXMatrix(n,m);
    int n_C = n, m_C = m, n_CT = m, m_CT = n; 
    cpx** CT = transpose_matrix(C, n, m);
    std::cout << "Транспонированная матрица контуров C:" << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << CT[i][j].real() << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Ввод матрицы сопротивлений Z" << endl;
    cpx** Z = createCPXMatrix(n,m);
    int n_Z = n, m_Z = m;
    std::cout << "Ввод матрицы ЭДС E" << endl;
    cpx** E = createCPXMatrix(n,m);
    int n_E = n, m_E = m;
    cpx** B = multiply_matrix(C, E, n_C, m_C, m_E);
    int n_B = n_C, m_B = m_E;
    cpx** A1 = multiply_matrix(C, Z, n_C, m_C, m_Z);
    int n_A1 = n_C, m_A1 = m_Z;
    cpx** A = multiply_matrix(A1, CT, n_A1, m_A1, m_CT);
    int n_A = n_A1, m_A = m_CT;
    cpx** G = concatenate_matrix(A, n_A, m_A, B, n_B, m_B);
       
    Gauss(G, n_A);
  
    return 0;
}