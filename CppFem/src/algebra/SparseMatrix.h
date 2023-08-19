/*****************************************************************//**
 * \file   SparseMatrix.h
 * \brief
 *
 * \author Leizs
 * \date   August 2023
 *********************************************************************/

#ifndef _SPARSE_MATRIX_HH
#define _SPARSE_MATRIX_HH

#include "Vector.h"
#include "Matrix.h"

 /**
  * @ingroup xxx
  * @~English
  * @brief This class define the sparse matrix.
  *
  * @~Chinese
  * @brief 简要描述。
  */
class SparseMatrix
{

public:
	/**
	 * @~Chinese
	 * @brief 此函数用于计算Sky Line存储格式的稀疏矩阵的结构。
	 * @param[in,out]	colHeight	At inital colHeight[i] = 1， 计算完后存储每一列的非零元素的高度。
	 * @param[in]		index		单元局部自由度对应的全局自由度
	 * @return void。
	 *
	 */
	void skyLineHeight(Vector<int>& colHeight, Vector<int>& index)
	{
		int ndof = index.size();
		for (int i = 0; i < ndof; ++i)
		{
			int maxH = 1;
			for (int j = 0; j < ndof; ++j)
			{
				maxH = std::max(index[i] - index[j] + 1, maxH);
			}

			int k = index[i];
			colHeight[k] = std::max(colHeight[k], maxH);
		}
	}

	/**
	 * @~Chinese
	 * @brief brief-description-about-skyLineProfile .
	 * @param[??] colHeight 列高
	 * @param[??] diagIndex 对角线元素的索引，前缀和数组。
	 * @return void brief-description-about-void .
	 * @eg
	 *	 *  *  *  *  0  0
	 *	 *  *  *  0  0  0
	 *	 *  *  *  0  *  0
	 *	 *  0  0  *  *  *
	 *   0  0  *  *  *  *
	 *   0  0  0  *  *  *
	 * @ the colHeight is [1 2 3 4 3 3]
	 * @ the diagIndex is [0 2 5 9 12 15] 一共有16个非0元素，前缀和数组第
	 * @ i列的非0元素 diagIndex[i] - diagIndex[i-1]
	 *
	 */
	void skyLineProfile(Vector<int>& colHeight, Vector<int>& diagIndex)
	{
		diagIndex.resize(colHeight.size(), 0);
		diagIndex[0] = colHeight[0];

		// diagIndex[i] restore the dial element index int the valus array.
		for (int i = 1; i < colHeight.size(); ++i)
		{
			diagIndex[i] = diagIndex[i - 1] + colHeight[i];
		}
	}

	/**
	 * @~Chinese
	 * This subroutine performs Cholesky factorisation on a symmetric
	   skyline global matrix.
	 * @brief 对稀疏矩阵进行cholesky分解。
	 * @param[??] kv brief-description-about-kv .
	 * @param[??] kdiag brief-description-about-kdiag .
	 * @return void brief-description-about-void .
	 * k[i][j] = kv[ki+j] ,ki = kdiag[i]-i
	 */
	void cholesky(Vector<double>& kv, const Vector<int>& kdiag)
	{
		int n = kdiag.size();
		kv[0] = sqrt(kv[0]);
		double x = 0;
		for (int i = 1; i < n; ++i)
		{
			int ki = kdiag[i] - i;
			int col = kdiag[i - 1] - ki + 1;  // ith row first non - zero element's col 
			for (int j = col; j <= i; ++j)  //注意这里是 <=
			{
				x = kv[ki + j];  // 第i行 第j个元素非0元素
				int kj = kdiag[j] - j;
				if (j != 0)
				{
					int col_1 = kdiag[j - 1] - kj + 1; // jth row first non - zero element's col
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

	/**
	 * @~Chinese
	 * @brief 将单元刚度矩阵组装到稀疏矩阵中去
	 * @param[??] kv brief-description-about-kv .
	 * @param[??] kdiag brief-description-about-kdiag .
	 * @param[??] g brief-description-about-g .
	 * @param[??] km brief-description-about-km .
	 * @return void brief-description-about-void .
	 *
	 */
	void assembleToSparse(Vector<double>& kv, const Vector<int>& kdiag, const Vector<int>& g, Matrix<double> km)
	{
		int idof = g.size();
		for (int i = 0; i < idof; ++i)
		{
			int gi = g[i];
			if (gi != -1)  // -1表示 自由度被约束
			{
				for (int j = 0; j < idof; ++j)
				{
					int gj = g[j];  // gi gj 为 i, j自由度对应的全局自由度 都是从0开始
					if (gj != 0)
					{
						int iw = gi - gj;
						if (iw >= 0) //只需要 下三角， 所有 i >= j
						{
							int ival = kdiag[gi] - iw;  // k[i,j] -> k[gi,gj]= kv[kdiag[gi] - (gi-gj)] ;
							kv[ival] += km(i, j);
						}
					}
				}
			}
		}
	}

	void spabac(Vector<double>& kv, const Vector<int>& kdiag, Vector<double>& loads)
	{
		//This subroutine performs Cholesky forward and back-substitution
	    // on a symmetric skyline global matrix.
		
		double x = 0;
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
				for(int j = col;i<= m;++j)
				{
					x = x - kv[ki + j] * loads[j];
				}
			}
			loads[i] = x / kv[ki + i];
		}

		for (int it = 1;it< n;++it)
		{
			int i = n - it;
			int ki = kdiag[i] - i;
			x = loads[i] / kv[ki + i];
			loads[i] = x;
			int col = kdiag[i - 1] - ki + 1;
			if (col != i)
			{
				int m = i - 1;
				for(int k = col;k<= m;++k)
				{
					loads[k] = loads[k] - x * kv[ki + k];
				}

			}
		}
		loads[1] = loads[1] / kv[1];
	}		

};



//ToDo:: 检查Programing the fem 书上的稀疏矩阵和自己实现的是不是一致的, 书上其实计算的是行的高度，不过和列的高度是一样的

#endif