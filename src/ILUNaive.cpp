#include "ILUNaive.h"
#include <vector>

using namespace std;

void ILUNaive::Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p)
{
	SparseMatrixCRS L, U;
	L.mN = A.mN;
	L.mValues.reserve(A.mNZ);
	L.mCol.reserve(A.mNZ);
	L.mRowIndex.resize(A.mN + 1);

	U.mN = A.mN;
	U.mValues.reserve(A.mNZ);
	U.mCol.reserve(A.mNZ);
	U.mRowIndex.resize(A.mN + 1);

	size_t LNZ = 0, UNZ = 0;
	L.mRowIndex[0] = LNZ;
	U.mRowIndex[0] = UNZ;
	for (size_t i = 0; i < A.mN; i++)
	{
		size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i+1];
		size_t LNZRow = 0, UNZRow = 0;
		for (size_t j = j1; j < j2; j++)
		{
			if (i < A.mCol[j])
			{
				U.mValues.push_back(A.mValues[j]);
				U.mCol.push_back(A.mCol[j]);
				UNZRow++;
			}
			else if (i > A.mCol[j])
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

				U.mValues.push_back(A.mValues[j]);
				U.mCol.push_back(A.mCol[j]);
				UNZRow++;
			}
		}
		L.mRowIndex[i+1] = L.mRowIndex[i] + LNZRow;
		U.mRowIndex[i+1] = U.mRowIndex[i] + UNZRow;
		LNZ += LNZRow;
		UNZ += UNZRow;
	}

	L.mNZ = LNZ;
	U.mNZ = UNZ;
	vector<size_t> posInRowL(L.mN);
	for (size_t i = 0; i < L.mN; i++)
	{
		posInRowL[i] = L.mRowIndex[i];
	}
	for (size_t i = 0; i < A.mN; i++)
	{
		for (size_t j = i + 1; j < A.mN; j++)
		{
			size_t u = U.mRowIndex[j];
			if (A.IsNonZero(j, i))
			{
				size_t k1 = A.mRowIndex[i], k2 = A.mRowIndex[j]; 
				while (A.mValues[k1] == 0) k1++;
				while (A.mValues[k2] == 0) k2++;
				double mn = A.mValues[k2] / A.mValues[k1];
				if (mn != 0)
				{
					L.mCol[posInRowL[j]] = i;
					L.mValues[posInRowL[j]] = mn;
					posInRowL[j]++;
				}

				for (size_t k = i; k < A.mN; k++)
				{
					if (A.IsNonZero(j, k) && A.IsNonZero(i, k))
					{
						A.mValues[k2] = A.mValues[k2] - mn * A.mValues[k1];
						if (k >= j)
						{
							U.mCol[u] = k;
							U.mValues[u] = A.mValues[k2];
							u++;
						}
						k1++; k2++;
					}
					else if (A.IsNonZero(j, k))
					{
						k2++;
					}
					else if (A.IsNonZero(i, k))
					{
						k1++;
					}
				}
			}
		}
	}

	SparseMatrixCRS UT = U.Transpose();
	L.Multiply(L, UT, M);
}

bool ILUNaive::CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M)
{
	bool result = true;

	for (size_t i = 0; i < A.mN && result; i++)
	{
		for (size_t j = 0; j < A.mN && result; j++)
		{
			if (!(A.Get(i, j) == M.Get(i, j) || (A.Get(i, j) == 0 && M.Get(i, j) != 0)))
			{
				result = false;
			}
		}
	}

	return result;
}

bool ILUNaive::CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M)
{
	
	bool result = false;
	vector<vector<double>>full_A;
	vector<vector<double>>inverse_full_A;
	vector<vector<double>>full_M;
	vector<vector<double>>inverse_full_M;
	full_A.resize(A.mN);
	inverse_full_A.resize(A.mN);
	full_M.resize(A.mN);
	inverse_full_M.resize(A.mN);
	for (int i = 0; i < A.mN; i++)
	{
		full_A[i].resize(A.mN);
		inverse_full_A[i].resize(A.mN);
		full_M[i].resize(A.mN);
		inverse_full_M[i].resize(A.mN);

	}
	for (int i = 0; i < A.mN; i++)
	{
		for (int j = 0; j < A.mN; j++)
		{
			full_A[i][j] = 0;
			full_M[i][j] = 0;
			inverse_full_A[i][j] = 0;
			inverse_full_M[i][j] = 0;
			if (i == j)
			{
				inverse_full_A[i][j] = 1;
				inverse_full_M[i][j] = 1;
			}
		}
	}
	A.recovery_matrix(A.mValues, A.mCol, A.mRowIndex, A.mN + 1, full_A);
	M.recovery_matrix(M.mValues, M.mCol, M.mRowIndex, M.mN + 1, full_M);
	double norm_A=A.norm_of_matrix(full_A,A.mN);
	double norm_M=M.norm_of_matrix(full_M, M.mN);
	A.inverse_Matrix(full_A, A.mN, inverse_full_A);
	M.inverse_Matrix(full_M, M.mN, inverse_full_M);
	double norm_inv_A = A.norm_of_matrix(inverse_full_A, A.mN);
	double norm_inv_M = M.norm_of_matrix(inverse_full_M, M.mN);

	double ma = norm_A * norm_inv_A;
	double mm = norm_M * norm_inv_M;
	
	result = ma > mm;

	return result;
}