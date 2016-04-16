#include "SparseMatrix.h"
#include <iostream>

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

double SparseMatrixCRS::GershgorinConditionNumber()
{
	double minEigVal = INT_MAX;
	double maxEigVal = INT_MIN;
	for (size_t i = 0; i < mN; i++)
	{
		size_t j1 = mRowIndex[i], j2 = mRowIndex[i + 1];
		double aii = 0.0;
		double ri = 0.0;
		for (size_t j = j1; j < j2; j++)
		{
			if (mCol[j] == i)
			{
				aii = mValues[j];
			}
			else
			{
				ri += fabs(mValues[j]);
			}
		}
		double leftEigVal = aii - ri;
		double rightEigVal = aii + ri;
		if (leftEigVal < minEigVal)
		{
			minEigVal = leftEigVal;
		}
		if (rightEigVal > maxEigVal)
		{
			maxEigVal = rightEigVal;
		}
	}
	return fabs(maxEigVal) / fabs(minEigVal);
}

SparseMatrixCRS SparseMatrixCRS::Transpose()
{
	SparseMatrixCRS result;
	result.mValues = mValues;
	result.mCol = mCol;
	result.mRowIndex = mRowIndex;
	result.mN = mN;
	result.mNZ = mNZ;

	for (size_t i = 0; i < mN; i++)
	{
		result.mRowIndex[i] = 0;
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

void SparseMatrixCRS::InitializeMatrix(int N, int NZ, SparseMatrixCRS &mtx) 
{ 
	mtx.mN = N; 
	mtx.mNZ = NZ; 
	mtx.mValues.resize(NZ); 
	mtx.mCol.resize(NZ); 
	mtx.mRowIndex.resize(N + 1);
} 

void SparseMatrixCRS::Multiply(SparseMatrixCRS &A, SparseMatrixCRS &BT, SparseMatrixCRS &C)
{
	size_t N = this->mN; 
	const double ZERO_IN_CRS = 0.00000001;
	
	vector<int> columns; 
	vector<double> values; 
	vector<int> row_index; 
 
	size_t rowNZ = 0; 
	row_index.push_back(0); 
	for (size_t i = 0; i < N; i++) 
	{ 
		rowNZ = 0; 
		for (size_t j = 0; j < N; j++) 
		{
			double sum = 0; 
			//скалярное произведение
			for (size_t k = A.mRowIndex[i]; k < A.mRowIndex[i + 1]; k++) 
			{
				for (size_t l = BT.mRowIndex[j]; l < BT.mRowIndex[j + 1]; l++) 
				{
					if (A.mCol[k] == BT.mCol[l]) 
					{ 
						sum += A.mValues[k] * BT.mValues[l]; 
						break; 
					} 				}			}
 
			if (fabs(sum) > ZERO_IN_CRS) 
			{ 
				columns.push_back(j); 
				values.push_back(sum); 
				rowNZ++; 
			} 
		} 
		row_index.push_back(rowNZ + row_index[i]); 
	} 
 
	InitializeMatrix(N, columns.size(), C); 
 
	for (size_t j = 0; j < columns.size(); j++) 
	{ 
		C.mCol[j] = columns[j]; 
		C.mValues[j] = values[j]; 
	} 

	for (size_t i = 0; i <= N; i++) 
	{
		C.mRowIndex[i] = row_index[i]; 
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
			double val = 0.0;
			for (size_t ii = 0; ii < is.size(); ii++)
			{
				if (is[ii] - 1 == i && js[ii] - 1 == j)
				{
					val = vals[ii];
					break;
				}
			} 

            if (val != 0.0)
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

void SparseMatrixCRS::Print()
{
    for (size_t i = 0; i < mN; i++)
    {
        for (size_t j = 0; j < mN; j++)
        { 
            cout << Get(i, j) << " ";
        }
        cout << endl;
    }
}

void SparseMatrixCRS::Set(int i, int j, double val)
{
    //TO DO
}

SparseMatrixCRS::~SparseMatrixCRS()
{

}
double SparseMatrixCRS::norm_of_matrix(vector<vector<double>> &a, size_t nn)
{
	double summ = 0, max_summ = 0;
	int n = nn;
	for (int i = 0; i < n; i++)
	{
		summ = 0;
		for (int j = 0; j < n; j++)
		{
			summ = summ + a[i][j];
		}
		if (summ > max_summ)
		{
			max_summ = summ;
		}
	}
	return max_summ;
}
void SparseMatrixCRS::inverse_Matrix(std::vector<std::vector<double>>&a, size_t nn, std::vector<std::vector<double>>&e)
{
	int n = nn;
	for (int i = 0; i < n; i++)
	{
		double mult = a[i][i];
		for (int j = 0; j < n; j++)
		{
			a[i][j] = a[i][j] / mult;
			e[i][j] = e[i][j] / mult;

		}
		for (int j = 0; j < n; j++)
		{
			if (j != i)
			{
				double multiplier = a[j][i];
				for (int k = 0; k < n; k++)
				{
					a[j][k] = a[j][k] - a[i][k] * multiplier;
					e[j][k] = e[j][k] - e[i][k] * multiplier;
				}
			}
		}
	}
}
void SparseMatrixCRS::recovery_matrix(vector<double>val, vector<size_t>colum, vector<size_t>rind, size_t srow, vector<vector<double>>&A)
{
	for (int i = 0; i < srow - 1; i++)
	{
		int temp = rind[i + 1] - rind[i];
		if (temp != 0)
		{
			for (int j = rind[i]; j < rind[i] + temp; j++)
			{
				A[i][colum[j]] = val[j];
			}
		}
	}
}