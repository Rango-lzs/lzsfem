// func_test.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include<vector>

#include <math.h>
#include <iomanip>

using namespace std;

template<class T>
using Vector = vector<T>;

void cholesky(Vector<double>& kv, const Vector<int>& kdiag)
{
	int n = kdiag.size();
	kv[0] = sqrt(kv[0]);
	double x = 0;
	for (int i = 1; i < n; ++i)
	{		
		int ki = kdiag[i] - i;
		int col = kdiag[i - 1] - ki +1;  // ith row first non - zero element's col 
		for (int j = col; j <= i; ++j)  //注意这里是 <=
		{
			x = kv[ki + j];  // 第i行 第j个元素非0元素
			int kj = kdiag[j] - j;
			if (j != 0)
			{
				int col_1 = kdiag[j - 1] - kj +1; // jth row first non - zero element's col
				col_1 = std::max(col, col_1);
				if (col_1 != j)
				{
					int m = j - 1;
					for (int k = col_1; k <= m; ++k) //注意这里是 <=
					{
						x = x - kv[ki + k] * kv[kj + k];  //K[i,k]*K[j,k]
					}
				}
			}
			kv[ki + j] = x / kv[kj + j];
		}
		kv[ki + i] = sqrt(x);
	}
}

//L*LT*x = b
//LT*x =y
//L*y = b
void spabac(Vector<double>& kv, const Vector<int>& kdiag, Vector<double>& loads)
{
	//This subroutine performs Cholesky forward and back-substitution
	// on a symmetric skyline global matrix.

	double x = 0.0;
	int n = kdiag.size();
	loads[0] = loads[0] / kv[0];
	for (int i = 1; i < n; ++i)
	{
		int ki = kdiag[i] - i;
		int col = kdiag[i - 1] - ki + 1; //first non-zero of row i
		x = loads[i];
		if (col != i)
		{
			int m = i - 1;
			for (int j = col; j<= m; ++j)
			{
				x = x - kv[ki + j] * loads[j];
			}
		}
		loads[i] = x / kv[ki + i];
	}

	for (int it = 1; it < n; ++it)
	{
		int i = n - it;
		int ki = kdiag[i] - i;
		x = loads[i] / kv[ki + i];
		loads[i] = x;
		int col = kdiag[i - 1] - ki + 1;
		if (col != i)
		{
			int m = i - 1;
			for (int k = col; k <= m; ++k)
			{
				loads[k] = loads[k] - x * kv[ki + k];
			}
		}
	}
	loads[0] = loads[0] / kv[0];
}

//Cholesky 分解
void Cholesky(double a[4][4], double L[4][4], int m, int n) {
	for (int j = 0; j < n; j++)
		for (int i = j; i < m; i++) {
			if (i == j) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum = sum + L[j][k] * L[j][k];
				}
				L[i][j] = sqrt(a[i][j] - sum);
			}
			else {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum = sum + L[i][k] * L[j][k];
				}
				L[i][j] = (a[i][j] - sum) / L[j][j];
			}
		}
}
//初始化矩阵 L
void InitL(double L[4][4], int m, int n) {
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			if (i < j) {
				L[i][j] = 0;
			}
		}
}
//显示 L 和 LT 矩阵
void Display(double L[4][4], double LT[4][4], int m, int n) {
	cout << "Cholesky 分解的 L 矩阵为：" << endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << setprecision(8) << setw(12) << L[i][j] << "";
		}
		cout << endl;
	}
	cout << "Cholesky 分解的 LT 矩阵为：" << endl;
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << setprecision(8) << setw(12) << L[j][i] << "";
			LT[i][j] = L[j][i];
		}
		cout << endl;
	}
}
//Ax=b 的解答 中间步骤求 y
double* SolveOne(double L[4][4], double b[4], int m, int n) {
	static double y[4];
	y[0] = b[0] / L[0][0];
	for (int i = 1; i < m; i++) {
		double sum = 0;
		for (int k = 0; k < i; k++) {
			sum = sum + L[i][k] * y[k];
		}
		y[i] = (b[i] - sum) / L[i][i];
	}
	cout << "解答的中间结果 y 为：" << endl;
	for (int i = 0; i < m; i++) {
		cout << setprecision(8) << setw(12) << y[i] << endl;
	}
	cout << endl;
	return y;
}
//Ax=b 的解答 最终结果
void SolveTwo(double L[4][4], double y[4], int m, int n) {
	double x[4];
	x[m - 1] = y[m - 1] / L[n - 1][n - 1];
	for (int i = m - 2; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k < m; k++) {
			sum = sum + L[k][i] * x[k];
		}
		x[i] = (y[i] - sum) / L[i][i];
	}
	cout << "解答的最终结果 x 为：" << endl;
	for (int i = 0; i < m; i++) {
		cout << setprecision(8) << setw(12) << x[i] << endl;
	}
	cout << endl;
}
//主函数
int main() {
	double a[4][4] = { {2,1,-1,1},{1,5,2,7},{-1,2,10,1},{1,7,1,11} };
	double b[4] = { 13,-9,6,0 };
	double L[4][4];//L 矩阵
	double LT[4][4];//LT 矩阵
	double* y;//中间矩阵 y
	InitL(L, 4, 4);//初始化矩阵 L
	Cholesky(a, L, 4, 4);//Cholesky 分解矩阵 A
	Display(L, LT, 4, 4);//显示矩阵 L 和 LT
	//求解原方程 Ax=b 中间步骤：L*y = b; LT*x = y;
	y = SolveOne(L, b, 4, 4);
	SolveTwo(L, y, 4, 4);


	Vector<double> kv{
		2,1,5,-1,2,10,1,7,1,11
	};

	Vector<int> kdiag{
		0,2,5,9
	};

	Vector<double> vb{ 13,-9,6,0 };

	cholesky(kv, kdiag);
	for (auto v : kv)
		cout << v << " ";
	cout << endl;

	spabac(kv, kdiag, vb);
	for (auto v : vb)
		cout << v << " ";
	cout << endl;
	
	return 0;
}


// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
