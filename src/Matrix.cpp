#include "Matrix.h"
#include <iostream>
#include <climits>

using namespace std;

//CRS implementation
Matrix::Matrix(size_t N)
{
    mN = N;
    mValues.resize(mN);
    for (size_t i = 0; i < mN; i++)
    {
        mValues[i].resize(mN);
        for (size_t j = 0; j < mN; j++)
        {
            mValues[i][j] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix &c)
{
    mN = c.mN;
	mValues.resize(mN);

    for (size_t i = 0; i < mN; i++)
    {
    	mValues[i].resize(mN);
        for (size_t j = 0; j < mN; j++)
        {
            mValues[i][j] = c.mValues[i][j];
        }
    }
}

double Matrix::Get(int i, int j) const
{
    double result = 0.0;

    result = mValues[i][j];

    return result;
}

void Matrix::Print(size_t outSize, size_t start_i, size_t start_j)
{
	size_t size = outSize && outSize < mN ? outSize : mN;
    for (size_t i = start_i; i < start_i + size; i++)
    {
        for (size_t j = start_j; j < start_j + size; j++)
        { 
			cout.width(13);
			cout.fill(' ');
            cout << Get(i, j) << " ";
        }
        cout << endl;
    }
}

void Matrix::Set(int i, int j, double val)
{
    mValues[i][j] = val;
}

double Matrix::NormMatrix()
{
	double summ = 0, max_summ = INT_MIN;
	for (size_t i = 0; i < mN; i++)
	{
		summ = 0;
		for (size_t j = 0; j < mN; j++)
		{
			summ = summ + mValues[i][j];
		}
		if (summ > max_summ)
		{
			max_summ = summ;
		}
	}
	return max_summ;
}

Matrix Matrix::GetInverseMatrix()
{
    Matrix matr(*this);
    Matrix invMatr(mN);
    for (size_t i = 0; i < mN; i++)
    {
        invMatr.mValues[i][i] = 1.0;
    }
	for (size_t i = 0; i < mN; i++)
	{
		double mult = matr.mValues[i][i];
		for (size_t j = 0; j < mN; j++)
		{
			matr.mValues[i][j] = matr.mValues[i][j] / mult;
			invMatr.mValues[i][j] = invMatr.mValues[i][j] / mult;

		}
		for (size_t j = 0; j < mN; j++)
		{
			if (j != i)
			{
				double multiplier = matr.mValues[j][i];
				for (size_t k = 0; k < mN; k++)
				{
					matr.mValues[j][k] = matr.mValues[j][k] - matr.mValues[i][k] * multiplier;
					invMatr.mValues[j][k] = invMatr.mValues[j][k] - invMatr.mValues[i][k] * multiplier;
				}
			}
		}
	}

    return invMatr;
}

void Matrix::fromCRS(SparseMatrixCRS &matr)
{
	for (size_t i = 0; i < matr.mN; i++)
	{
        size_t j1 = matr.mRowIndex[i], j2 = matr.mRowIndex[i + 1];
		for (size_t j = j1; j < j2; j++)
		{
			this->mValues[i][matr.mCol[j]] = matr.mValues[j];
		}
	}
}

Matrix::~Matrix()
{

}