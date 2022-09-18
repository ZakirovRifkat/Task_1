#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

//void getMatrixWithoutRowAndCol(vector<vector<double>> &matrix, int size, int row, int col, vector<vector<double>> &newMatrix) 
//{
//	int offsetRow = 0;
//	int offsetCol = 0;
//	for (int i = 0; i < size - 1; i++) 
//	{
//		if (i == row) 
//		{
//			offsetRow = 1;
//		}
//		offsetCol = 0;
//		for (int j = 0; j < size - 1; j++) 
//		{
//			if (j == col) 
//			{
//				offsetCol = 1;
//			}
//			newMatrix[i][j] = matrix[i + offsetRow][j + offsetCol];
//		}
//	}
//}
//int Determination( vector<vector<double>>&matrix, int size) 
//{
//	double det = 0;
//	int degree = 1;
//
//	if (size == 1) {
//		return matrix[0][0];
//	}
//
//	if (size == 2) {
//		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
//	}
//
//	vector<vector<double>> newMatrix(size - 1, vector<double>(size-1));
//
//	for (int j = 0; j < size; j++) 
//	{
//		getMatrixWithoutRowAndCol(matrix, size, 0, j, newMatrix);
//		det = det + (degree * matrix[0][j] * Determination(newMatrix, size - 1));
//		degree = -degree;
//	}
//	return det;
//}
void show(vector<vector <double>> A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << A[i][j] << "\t";
		}
		cout << endl;
	}
}
vector<double> multiplyMatrixVector(vector<vector<double>>& A, vector<double>& b, int n)
{
	vector<double>Ab(n);
	for (int i = 0; i < n; i++)
	{
		Ab[i] = 0;
		for (int j = 0; j < n; j++)
			Ab[i] += A[i][j] * b[j];
	}
	return Ab;
}
void multiplyMatrixMatrix(vector <vector <double>> A, vector <vector <double>> B,vector <vector <double>>& AB, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				AB[i][j] += A[i][k] * B[k][j];
}
vector<double> singleGauss(vector<vector<double>>matrix_A, vector<double>vector_b, int n)
{
	double d, s;
	vector<double>x(n),b(n);
	vector<vector<double>>a(n, vector<double>(n));
	a = matrix_A;
	b = vector_b;
	for (int k = 0; k < n; k++) 
	{
		for (int j = k + 1; j < n; j++)
		{
			d = a[j][k] / a[k][k]; 
			for (int i = k; i < n; i++)
			{
				a[j][i] = a[j][i] - d * a[k][i]; 
			}
			b[j] = b[j] - d * b[k]; 
		}
	}
	for (int k = n - 1; k >= 0; k--) 
	{
		d = 0;
		for (int j = k; j < n; j++)
		{
			s = a[k][j] * x[j]; 
			d = d + s; 
		}
		x[k] = (b[k] - d) / a[k][k]; 
	}
	return x;
}
vector<double> modGauss (vector<vector<double>> matrix_A,  vector<double> vector_b, int n)
{
	int p;
	double r,c,s;
	vector<double>x(n), b(n);
	vector<vector<double>>a(n, vector<double>(n));
	a = matrix_A;
	b = vector_b;
	for (int k = 0; k < n; k++)
	{
		p = k;
		for (int m = k + 1; m < n; m++)
		{
			if (abs(a[p][k]) < abs(a[m][k])) //поиск максимального ведущего элемента
			{
				p = m;
			}
		}
		for (int j = k; j < n; j++)
		{
			r = a[k][j];
			a[k][j] = a[p][j];   //перестановка строк
			a[p][j] = r;
		}
		r = b[k];
		b[k] = b[p];   //перестановка свободных членов
		b[p] = r;
		for (int m = k+1; m < n; m++)
		{
			c = a[m][k] / a[k][k];
			b[m] = b[m] - c * b[k]; //приведение матрицы к верхнетреугольному виду
			for (int i = k; i < n; i++)
			{
				a[m][i] = a[m][i] - c * a[k][i];
			}
		}
	}
	x[n-1] = b[n-1] / a[n-1][n-1];
	for (int k = n - 1; k >= 0; k--)
	{
		s = 0;
		for (int i = k + 1; i < n; i++)				//обратный ход метода Гаусса
		{												
			s = s + a[k][i] * x[i];
		}
		x[k] = (b[k] - s) / a[k][k];
	}
	return x;
}
void LU(vector <vector <double>> A, vector <vector <double>>& L,vector <vector <double>>& U, int n)
{
	U = A;

	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			L[j][i] = U[j][i] / U[i][i];

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++)
				L[j][i] = U[j][i] / U[i][i];

		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i][j] = U[i][j] - L[i][k - 1] * U[k - 1][j];
	}

}
double normMatrix(vector<vector<double>>& A, int n)
{
	vector<double>sums(n);
	for (int j = 0; j < n; j++)
	{
		double sum = 0;
		for (int i = 0; i < n; i++)
		{
			sum += abs(A[i][j]);
		}
		sums[j] = sum;
	}
	double max = *max_element(sums.begin(), sums.end());
	return max;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	int n;
	cout << "Задание 1. Прямые методы решения линейных систем\n"
		<< "Введите размерность матрицы (А):\n"
		<< "n=";
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double>(n)),
		L(n, vector<double>(n)),
		U(n, vector<double>(n)),
		matrix_LU(n, vector<double>(n)),
		reverse_A(n, vector<double>(n));
	vector<double> vector_b(n),
		solution_Gauss(n),
		solution_modGauss(n),
		solution_LU(n),
		R(n),
		Ax(n),
		y(n);
	//R(n) - вектор невязки

	cout << "\nЗаполняем матрицу (А):\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			{
				cout << "[" << i + 1 << "][" << j + 1 << "] = ";
				cin >> matrix_A[i][j];
			}
		}
	}
	cout << "\nЗаполняем вектор (b):\n";
	for (int i = 0; i < n; i++)
	{
		cout << "b[" << i + 1 << "] = ";
		cin >> vector_b[i];
	}
	cout << "\nРасширенная матрица (А|b):\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << "\t" << matrix_A[i][j];
		cout << "\t" << vector_b[i] << endl;
	}

	//singleGauss
	{
		solution_Gauss = singleGauss(matrix_A, vector_b, n);
		cout << "\nРешение методом Гаусса единственного деления:\n";
		for (int i = 0; i < n; i++)
			cout << "x[" << i + 1 << "] = " << solution_Gauss[i] << "\n";
		Ax = multiplyMatrixVector(matrix_A, solution_Gauss, n);
		cout << "Вектор невязки метода Гаусса единственного деления R = ( ";
		for (int i = 0; i < n; i++)
		{
			if (i < n - 1)
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i] << " ";
			}
			else
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i];
			}
		}
		cout << " )\n";
	}
	//modGauss
	{
		solution_modGauss = modGauss(matrix_A, vector_b, n);
		cout << "\nРешение методом Гаусса с выбором главного элемента:\n";
		for (int i = 0; i < n; i++)
			cout << "x[" << i + 1 << "] = " << solution_modGauss[i] << "\n";
		Ax = multiplyMatrixVector(matrix_A, solution_modGauss, n);
		cout << "Вектор невязки метода Гаусса c выбором главного элемента R = ( ";
		for (int i = 0; i < n; i++)
		{
			if (i < n - 1)
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i] << " ";
			}
			else
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i];
			}
		}
		cout << " )\n";
	}

	//LU
	{
		LU(matrix_A, L, U, n);
		cout << "\n\tМатрица U\n";
		show(U, n);
		cout << "\n\tМатрица L\n";
		show(L, n);
		multiplyMatrixMatrix(L, U, matrix_LU, n);
		cout << "\n\tL*U matrix" << endl;
		show(matrix_LU, n);
		y = singleGauss(L, vector_b, n);
		solution_LU = singleGauss(U, y, n);
		cout << "\nРешение методом LU-разложения:\n";
		for (int i = 0; i < n; i++)
			cout << "x[" << i + 1 << "] = " << solution_LU[i] << "\n";
		Ax = multiplyMatrixVector(matrix_A, solution_LU, n);
		cout << "Вектор невязки метода LU-разложения R = ( ";
		for (int i = 0; i < n; i++)
		{
			if (i < n - 1)
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i] << " ";
			}
			else
			{
				R[i] = vector_b[i] - Ax[i];
				cout << R[i];
			}
		}
		cout << " )\n";
	}
		vector<double>temp_b(n);
		vector<double>solve(n);
		for (int i = 0; i < n; i++)
		{
			temp_b[i] = 1;
			solve = singleGauss(matrix_A, temp_b, n);
			for (int j = 0; j < n; j++)
			{
				reverse_A[j][i] = solve[j];
			}
			temp_b[i] = 0;
		}
		cout << "\n\tОбратная матрица А^-1:\n";
		show(reverse_A, n);
		
		cout << "Число обусловленности cond(A) = " << normMatrix(reverse_A, n) * normMatrix(matrix_A, n)<<"\n";
}


