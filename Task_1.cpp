#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void fill_matrix(vector<vector<double>>&matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "Элемент [" << i + 1 << "][" << j + 1 << "] = ";
			cin >> matrix[i][j];
		}
	}
}
void print_matrix(vector<vector<double>>&matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << matrix[i][j] << "\t";
		}
		cout << endl;
	}
}

int main()
{
	setlocale(LC_ALL, "Russian");	
	int n;
	double element;
	cout<<"Задание 1. Прямые методы решения линейных систем\n";
	cout << "Введите размерность матрицы А\n"; 
	cin >> n;
	vector<vector<double>> matrix_A(n, vector<double> (n));
	vector<double> vector_b(n);
	cout << "Заполняем матрицу (А):\n";
	fill_matrix(matrix_A, n); 
	cout << "Ваша матрица А:\n";
	print_matrix(matrix_A, n); 
	cout << "Заполняем вектор (b):\n";
	for (int i = 0; i < n; i++)
	{
		cin >> vector_b[i];
	}
	cout << "Расширенная матрица (A|b):\n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "\t" << matrix_A[i][j];
		}
		cout <<"\t" << vector_b[i] << endl;
    }

	
	
}

