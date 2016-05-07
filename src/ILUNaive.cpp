#include "ILUNaive.h"
#include <vector>
#include <time.h>

using namespace std;

void ILUNaive::Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p)
{
	SparseMatrixCRS L, U(A);
	L.mN = A.mN;
	L.mValues.reserve(A.mNZ);
	L.mCol.reserve(A.mNZ);
	L.mRowIndex.resize(A.mN + 1);
	size_t LNZ = 0;
	L.mRowIndex[0] = LNZ;
	for (size_t i = 0; i < A.mN; i++)
	{
		size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i+1];
		size_t LNZRow = 0, UNZRow = 0;
		for (size_t j = j1; j < j2; j++)
		{
			if (i > A.mCol[j])
			{
				L.mValues.push_back(A.mValues[j]);
				L.mCol.push_back(A.mCol[j]);
				LNZRow++;
			}
			else if (i == A.mCol[j])
			{
				L.mValues.push_back(1.0);
				L.mCol.push_back(A.mCol[j]);
				LNZRow++;
			}
		}
		L.mRowIndex[i+1] = L.mRowIndex[i] + LNZRow;
		LNZ += LNZRow;
	}

	L.mNZ = LNZ;
	vector<size_t> posInRowL(L.mN);
	for (size_t i = 0; i < L.mN; i++)
	{
		posInRowL[i] = L.mRowIndex[i];
	}
	for (size_t i = 0; i < U.mN; i++)
	{
		for (size_t j = i + 1; j < U.mN; j++)
		{
			if (U.IsNonZero(j, i))
			{
				size_t k1 = U.mRowIndex[i], k2 = U.mRowIndex[j]; 
				while (i > U.mCol[k1] && k1 < U.mRowIndex[i+1]) k1++;
				while (i > U.mCol[k2] && k2 < U.mRowIndex[j+1]) k2++;
				if (k1 <= k2 && k1 < U.mRowIndex[i+1] && k2 < U.mRowIndex[j+1])
				{
					double mn = U.mValues[k2] / U.mValues[k1];
					//if (fabs(mn) > ZERO_IN_CRS)
					{
						L.mCol[posInRowL[j]] = i;
						L.mValues[posInRowL[j]] = mn;
						posInRowL[j]++;

						while (k1 < U.mRowIndex[i+1] && k2 < U.mRowIndex[j+1])
						{
							if (U.mCol[k1] == U.mCol[k2])
							{
								U.mValues[k2] = U.mValues[k2] - mn * U.mValues[k1];
								
								k1++; k2++;
							}
							else if (U.mCol[k1] < U.mCol[k2])
							{
								k1++;
							}
							else if (U.mCol[k1] > U.mCol[k2])
							{
								k2++;
							}
						}
					}
				}
			}
		}
	}

	SparseMatrixCRS UT = U.Transpose();

	L.MultiplyNaive(L, UT, M);
}

bool ILUNaive::isCorrectMatrix(SparseMatrixCRS &A)
{
	bool result = true;

	for (size_t i = 0; i < A.mN && result; i++)
	{
		size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i + 1];
		for (size_t j = j1; j < j2 && result; j++)
		{
			if (A.mCol[j] == i && fabs(A.mValues[j]) < ZERO_IN_CRS)
			{
				result = false;
			}
		}
	}

	return result;
}

bool ILUNaive::CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M)
{
	bool result = true;

	for (size_t i = 0; i < A.mN && result; i++)
	{
		for (size_t j = 0; j < A.mN && result; j++)
		{
			double a = A.Get(i, j);
			double m = M.Get(i, j);
			if (!(fabs(a - m) < ZERO_IN_CRS || (fabs(a) < ZERO_IN_CRS && fabs(m) > ZERO_IN_CRS)))
			{
				printf("%lf %lf - (%d, %d)\n", a, m, i, j);
				result = false;
				break;
			}
		}
	}

	return result;
}

bool ILUNaive::CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M)
{
	
	bool result = false;
	vector<vector<double>> full_A;
	vector<vector<double>> inverse_full_A;
	vector<vector<double>> full_M;
	vector<vector<double>> inverse_full_M;
	vector<vector<double>> full_M_A;
	vector<vector<double>> inverse_full_M_A;
	full_A.resize(A.mN);
	inverse_full_A.resize(A.mN);
	full_M.resize(A.mN);
	inverse_full_M.resize(A.mN);
	full_M_A.resize(A.mN);
	inverse_full_M_A.resize(A.mN);
	for (int i = 0; i < A.mN; i++)
	{
		full_A[i].resize(A.mN);
		inverse_full_A[i].resize(A.mN);
		full_M[i].resize(A.mN);
		inverse_full_M[i].resize(A.mN);
		full_M_A[i].resize(A.mN);
		inverse_full_M_A[i].resize(A.mN);
	}
	for (int i = 0; i < A.mN; i++)
	{
		for (int j = 0; j < A.mN; j++)
		{
			full_A[i][j] = 0;
			full_M[i][j] = 0;
			inverse_full_A[i][j] = 0;
			inverse_full_M[i][j] = 0;
			full_M_A[i][j] = 0;
			inverse_full_M_A[i][j] = 0;
			if (i == j)
			{
				inverse_full_A[i][j] = 1;
				inverse_full_M[i][j] = 1;
				inverse_full_M_A[i][j] = 1;
			}
		}
	}


	A.recovery_matrix(A.mValues, A.mCol, A.mRowIndex, A.mN + 1, full_A);
	double norm_A = A.norm_of_matrix(full_A,A.mN);

	M.recovery_matrix(M.mValues, M.mCol, M.mRowIndex, M.mN + 1, full_M);
	double norm_M = M.norm_of_matrix(full_M, M.mN);

	M.inverse_Matrix(full_M, M.mN, inverse_full_M);
	double norm_inv_M = M.norm_of_matrix(inverse_full_M, M.mN);

	A.MultiplyFullMatrix(inverse_full_M, full_A, full_M_A, A.mN);
	double norm_M_A = A.norm_of_matrix(full_M_A, A.mN);

	A.inverse_Matrix(full_A, A.mN, inverse_full_A);	
	double norm_inv_A = A.norm_of_matrix(inverse_full_A, A.mN);

	double ma = norm_A * norm_inv_A;

	A.inverse_Matrix(full_M_A, A.mN, inverse_full_M_A);
	double norm_inverse_M_A = A.norm_of_matrix(inverse_full_M_A, A.mN);

	double mmm = norm_inverse_M_A * norm_M_A;
	
	result = ma > mmm;
	printf("%lf %lf \n", norm_inverse_M_A, norm_M_A);
	printf("ma = %lf, mmm = %lf\n", ma, mmm);

	return result;
}