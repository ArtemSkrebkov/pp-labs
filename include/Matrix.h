#ifndef _MATRIX_
#define _MATRIX_


#include <vector>
#include <fstream>
#include <cmath>
#include "SparseMatrix.h"

class Matrix
{
public:
    Matrix(size_t N = 0);
	Matrix(const Matrix &c);
	Matrix Transpose();
    double Get(int i, int j) const;
    void Set(int i, int j, double val) ;
    void Print( size_t outSize = 0, size_t start_i = 0, size_t start_j = 0);
    ~Matrix();

    Matrix GetInverseMatrix();
    double NormMatrix();
    void fromCRS(SparseMatrixCRS &matr);

	inline bool operator==(const Matrix& rhs)
	{
		bool result = true;
		if (this->mN != rhs.mN)
		{
			result = false;
		}
		else
		{
			for (size_t i = 0; i < rhs.mN && result; i++)
			{
				for (size_t j = 0; j < rhs.mN && result; j++)
				{
					if (fabs(this->Get(i, j) - rhs.Get(i, j)) > ZERO_IN_CRS)
					{
						result = false;
						break;
					}
				}
			}
		}

		return result;
	}
	inline bool operator!=(const Matrix& rhs){ return !(*this == rhs); }

	friend inline bool operator==(const Matrix& lhs, const Matrix& rhs);
public:
    size_t mN;

    std::vector<std::vector<double>> mValues;
};

inline bool operator==(const Matrix& lhs, const Matrix& rhs)
{
	bool result = true;
	if (lhs.mN != rhs.mN)
	{
		result = false;
	}
	else
	{
		for (size_t i = 0; i < rhs.mN; i++)
		{
			for (size_t j = 0; j < rhs.mN; j++)
			{
				if (fabs(lhs.mValues[i][j] - rhs.mValues[i][j]) > ZERO_IN_CRS)
				{
					result = false;
					break;
				}
			}
		}
	}
	return result;
}
	

#endif