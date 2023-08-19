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
  * @brief ��Ҫ������
  */
class SparseMatrix
{

public:
	/**
	 * @~Chinese
	 * @brief �˺������ڼ���Sky Line�洢��ʽ��ϡ�����Ľṹ��
	 * @param[in,out]	colHeight	At inital colHeight[i] = 1�� �������洢ÿһ�еķ���Ԫ�صĸ߶ȡ�
	 * @param[in]		index		��Ԫ�ֲ����ɶȶ�Ӧ��ȫ�����ɶ�
	 * @return void��
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
	 * @param[??] colHeight �и�
	 * @param[??] diagIndex �Խ���Ԫ�ص�������ǰ׺�����顣
	 * @return void brief-description-about-void .
	 * @eg
	 *	 *  *  *  *  0  0
	 *	 *  *  *  0  0  0
	 *	 *  *  *  0  *  0
	 *	 *  0  0  *  *  *
	 *   0  0  *  *  *  *
	 *   0  0  0  *  *  *
	 * @ the colHeight is [1 2 3 4 3 3]
	 * @ the diagIndex is [0 2 5 9 12 15] һ����16����0Ԫ�أ�ǰ׺�������
	 * @ i�еķ�0Ԫ�� diagIndex[i] - diagIndex[i-1]
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
	 * @brief ��ϡ��������cholesky�ֽ⡣
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
			for (int j = col; j <= i; ++j)  //ע�������� <=
			{
				x = kv[ki + j];  // ��i�� ��j��Ԫ�ط�0Ԫ��
				int kj = kdiag[j] - j;
				if (j != 0)
				{
					int col_1 = kdiag[j - 1] - kj + 1; // jth row first non - zero element's col
					col_1 = std::max(col, col_1);
					if (col_1 != j)
					{
						int m = j - 1;
						for (int k = col_1; k <= m; ++k) //ע�������� <=
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
	 * @brief ����Ԫ�նȾ�����װ��ϡ�������ȥ
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
			if (gi != -1)  // -1��ʾ ���ɶȱ�Լ��
			{
				for (int j = 0; j < idof; ++j)
				{
					int gj = g[j];  // gi gj Ϊ i, j���ɶȶ�Ӧ��ȫ�����ɶ� ���Ǵ�0��ʼ
					if (gj != 0)
					{
						int iw = gi - gj;
						if (iw >= 0) //ֻ��Ҫ �����ǣ� ���� i >= j
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



//ToDo:: ���Programing the fem ���ϵ�ϡ�������Լ�ʵ�ֵ��ǲ���һ�µ�, ������ʵ��������еĸ߶ȣ��������еĸ߶���һ����

#endif