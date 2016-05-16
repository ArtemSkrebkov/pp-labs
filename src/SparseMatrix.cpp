#include "SparseMatrix.h"
#include <iostream>
#include <string>
#include <climits>

using namespace std;
//CRS implementation
SparseMatrixCRS::SparseMatrixCRS()
{

}

SparseMatrixCRS::SparseMatrixCRS(const string &filename)
{
    ifstream in(filename);

	ReadFromMtx(filename);

    in.close();
}

SparseMatrixCRS::SparseMatrixCRS(const SparseMatrixCRS &c)
{
	mValues = c.mValues;
	mCol = c.mCol;
	mRowIndex = c.mRowIndex;
	mNZ = c.mNZ;
	mN = c.mN;
}

SparseMatrixCRS SparseMatrixCRS::Transpose()
{
	SparseMatrixCRS result;
	result.mValues = mValues;
	result.mCol = mCol;
	result.mRowIndex = mRowIndex;
	result.mN = mN;
	result.mNZ = mNZ;

	for (size_t i = 0; i < mN + 1; i++)
	{
		result.mRowIndex[i] = 0;
	}
	for (size_t i = 0; i < mNZ; i++)
	{
		result.mCol[i] = 0;
		result.mValues[i] = 0.0;
	}
	for (size_t i = 0; i < mNZ; i++)
	{
		result.mRowIndex[mCol[i] + 1]++;
	}

	size_t S = 0;
	for (size_t i = 1; i <= mN; i++)
	{
		size_t tmp = result.mRowIndex[i];
		result.mRowIndex[i] = S;
		S = S + tmp;
	}

	for (size_t i = 0; i < mN; i++)
	{
		size_t j1 = mRowIndex[i], j2 = mRowIndex[i + 1];
		size_t Col = i;
		for (size_t j = j1; j < j2; j++)
		{
			double V = mValues[j];
			size_t RIndex = mCol[j];
			size_t IIndex = result.mRowIndex[RIndex + 1];
			result.mValues[IIndex] = V;
			result.mCol[IIndex] = Col;
			result.mRowIndex[RIndex + 1]++;
		}
	}

	return result;
}

void SparseMatrixCRS::InitializeMatrix(size_t N, size_t NZ) 
{ 
	mN = N; 
	mNZ = NZ; 
	mValues.resize(NZ); 
	mCol.resize(NZ); 
    for (int i = 0; i < NZ; i++)
    {
        mValues[i] = 0.0;
        mCol[i] = 0.0;
    }
	
    mRowIndex.resize(N + 1);
    for (int i = 0; i < (N + 1); i++)
    {
    	mRowIndex[i] = 0;
    }
} 


bool SparseMatrixCRS::IsNonZero(size_t i, size_t j)
{
	bool result = false;

	for (size_t k = mRowIndex[i]; k < mRowIndex[i + 1]; k++)
	{
		if (mCol[k] == j && mValues[k] != 0.0)
		{
			result = true;
			break;
		}
	}

	return result;
}

void SparseMatrixCRS::ReadFromMtx(const std::string filename)
{
    mValues.clear();
    mCol.clear();
    mRowIndex.clear();
	
	ifstream in(filename);
	size_t M, N, NZ;
	in >> M >> N >> NZ;
	vector<size_t> is(NZ);
	vector<size_t> js(NZ);
	vector<double> vals(NZ);
	for (size_t k = 0; k < NZ; k++)
	{
		in >> is[k] >> js[k] >> vals[k];
	}

	mN = N;
	mNZ = NZ;
    mValues.resize(NZ);
    mCol.resize(NZ);
    mRowIndex.resize(N + 1);

    size_t k = 0;
    size_t countNonZeroInRow = 0;
    mRowIndex[0] = countNonZeroInRow;
    for (size_t i = 0; i < mN; i++) 
    {
        for (size_t j = 0; j < mN; j++)
        {
			double val = INT_MAX;
			for (size_t ii = 0; ii < is.size(); ii++)
			{
				if (is[ii] - 1 == i && js[ii] - 1 == j)
				{
					val = vals[ii];
					break;
				}
			} 

            if (val != INT_MAX)
            {
                mValues[k] = val;
                mCol[k] = j;
				k++;
                countNonZeroInRow++;
            }
        }
        mRowIndex[i + 1] = mRowIndex[i] + countNonZeroInRow;
		countNonZeroInRow = 0;
    }
	
	in.close();
}

void SparseMatrixCRS::Print(size_t outSize, size_t start_i, size_t start_j)
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

double SparseMatrixCRS::Get(int i, int j) const
{
    double result = 0.0;

    size_t startI = mRowIndex[i];
    for (size_t ii = startI; ii < mRowIndex[i + 1]; ii++)
    {
        if (mCol[ii] == j)
        {
            result = mValues[ii];
            break;
        }
    }

    return result;
}

void SparseMatrixCRS::Set(int i, int j, double val)
{
    //TO DO
}

SparseMatrixCRS::~SparseMatrixCRS()
{

}