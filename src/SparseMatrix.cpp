#include "SparseMatrix.h"
#include <iostream>
#include <omp.h>
#include <tbb\tbb.h>

using namespace std;
using namespace tbb;

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

void SparseMatrixCRS::InitializeMatrix(int N, int NZ, SparseMatrixCRS &mtx) 
{ 
	mtx.mN = N; 
	mtx.mNZ = NZ; 
	mtx.mValues.resize(NZ); 
	mtx.mCol.resize(NZ); 
	mtx.mRowIndex.resize(N + 1);
} 

void SparseMatrixCRS::MultiplyNaive(SparseMatrixCRS &A, SparseMatrixCRS &BT, SparseMatrixCRS &C)
{
	size_t N = this->mN; 
	
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
					} 
				}
			}

 
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

void SparseMatrixCRS::MultiplyOpenMP(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
	int N = A.mN; 
	int i, j, k; 
 
  vector<int>* columns = new vector<int>[N]; 
  vector<double> *values = new vector<double>[N]; 
 
  int* row_index = new int[N + 1]; 
  memset(row_index, 0, sizeof(int) * N); 
 
#pragma omp parallel 
  { 
    int *temp = new int[N]; 
 #pragma omp for private(j, k) schedule(static, chunk)
    for (i = 0; i < N; i++) 
    { 
      memset(temp, -1, N * sizeof(int)); 
      int ind1 = A.mRowIndex[i], ind2 = A.mRowIndex[i + 1]; 
      for (j = ind1; j < ind2; j++) 
      { 
        int col = A.mCol[j]; 
        temp[col] = j; // Значит, что a[i, НОМЕР] лежит  
        // в ячейке массива Value с номером temp[НОМЕР] 
      } 
      // Построен индекс строки i матрицы A 
      // Теперь необходимо умножить ее на каждую из строк  
      // матрицы BT 
      for (j = 0; j < N; j++) 
      { 
        // j-я строка матрицы B 
        double sum = 0; 
        int ind3 = B.mRowIndex[j], ind4 = B.mRowIndex[j + 1]; 
        // Все ненулевые элементы строки j матрицы B 
        for (k = ind3; k < ind4; k++) 
        { 
          int bcol = B.mCol[k]; 
          int aind = temp[bcol]; 
          if (aind != -1) 
            sum += A.mValues[aind] * B.mValues[k]; 
        } 
        if (fabs(sum) > ZERO_IN_CRS) 
        { 
          columns[i].push_back(j); 
		            values[i].push_back(sum); 
          row_index[i]++; 
        } 
      } 
    } 
    delete [] temp; 
  } 
 
  int NZ = 0; 
  for(i = 0; i < N; i++) 
  { 
    int tmp = row_index[i]; 
    row_index[i] = NZ; 
    NZ += tmp; 
  } 
  row_index[N] = NZ; 
 
  InitializeMatrix(N, NZ, C); 
 
  int count = 0; 
  for (i = 0; i < N; i++) 
  { 
    int size = columns[i].size(); 
    memcpy(&C.mCol[count], &columns[i][0],  
           size * sizeof(int)); 
    memcpy(&C.mValues[count], &values[i][0],  
           size * sizeof(double)); 
    count += size; 
  } 

	for (size_t i = 0; i <= N; i++) 
	{
		C.mRowIndex[i] = row_index[i]; 
	}
 
  delete [] row_index; 
  delete [] columns; 
  delete [] values; 
}

class Multiplicator 
{ 
  SparseMatrixCRS A, B; 
  vector<int>* columns; 
  vector<double>* values; 
  int *row_index; 
public: 
  Multiplicator(SparseMatrixCRS& _A, SparseMatrixCRS& _B, 
    vector<int>* &_columns, vector<double>* &_values, 
    int *_row_index) : A(_A), B(_B), columns(_columns), 
    values(_values), row_index(_row_index) 
  {} 
 
  void operator()(const blocked_range<int>& r) const 
  { 
    int begin = r.begin(); 
    int end = r.end(); 
    int N = A.mN; 
 
    int i, j, k; 
    int *temp = new int[N]; 
 
    for (i = begin; i < end; i++) 
    { 
      memset(temp, -1, N * sizeof(int)); 
      int ind1 = A.mRowIndex[i], ind2 = A.mRowIndex[i + 1]; 
      for (j = ind1; j < ind2; j++) 
      { 
        int col = A.mCol[j]; 
        temp[col] = j;
      } 
      for (j = 0; j < N; j++) 
      { 
        double sum = 0; 
        int ind3 = B.mRowIndex[j], ind4 = B.mRowIndex[j + 1]; 
        for (k = ind3; k < ind4; k++) 
        { 
          int bcol = B.mCol[k]; 
          int aind = temp[bcol]; 
          if (aind != -1) 
            sum += A.mValues[aind] * B.mValues[k]; 
        } 
        if (fabs(sum) > ZERO_IN_CRS) 
        { 
          columns[i].push_back(j); 
          values[i].push_back(sum); 
          row_index[i]++; 
        } 
      } 
    } 
    delete [] temp; 

  } 
}; 


void SparseMatrixCRS::MultiplyTBB(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
 
  int N = A.mN; 
  int i; 
 
  task_scheduler_init init(); 
 
  vector<int>* columns = new vector<int>[N]; 
  vector<double> *values = new vector<double>[N]; 
  int* row_index = new int[N + 1]; 
  memset(row_index, 0, sizeof(int) * N); 
 
  int grainsize = 10; 
 
  parallel_for(blocked_range<int>(0, A.mN, grainsize), 
    Multiplicator(A, B, columns, values, row_index)); 
 
  int NZ = 0; 
  for(i = 0; i < N; i++) 
  { 
    int tmp = row_index[i]; 
    row_index[i] = NZ; 
    NZ += tmp; 
  } 
  row_index[N] = NZ; 
 
  InitializeMatrix(N, NZ, C); 
 
  int count = 0; 
  for (i = 0; i < N; i++) 
  { 
    int size = columns[i].size(); 
    memcpy(&C.mCol[count], &columns[i][0],  
           size * sizeof(int)); 
    memcpy(&C.mValues[count], &values[i][0],  
           size * sizeof(double)); 
    count += size; 
  } 
	for (size_t i = 0; i <= N; i++) 
	{
		C.mRowIndex[i] = row_index[i]; 
	}
}

void SparseMatrixCRS::MultiplyCilk(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
	int N = A.mN; 
	int i, j, k; 
 
  vector<int>* columns = new vector<int>[N]; 
  vector<double> *values = new vector<double>[N]; 
 
  int* row_index = new int[N + 1]; 
  memset(row_index, 0, sizeof(int) * N); 
 
#pragma omp parallel 
  { 
    int *temp = new int[N]; 
 #pragma omp for private(j, k) schedule(static, chunk)
    for (i = 0; i < N; i++) 
    { 
      memset(temp, -1, N * sizeof(int)); 
      int ind1 = A.mRowIndex[i], ind2 = A.mRowIndex[i + 1]; 
      for (j = ind1; j < ind2; j++) 
      { 
        int col = A.mCol[j]; 
        temp[col] = j; // Значит, что a[i, НОМЕР] лежит  
        // в ячейке массива Value с номером temp[НОМЕР] 
      } 
      // Построен индекс строки i матрицы A 
      // Теперь необходимо умножить ее на каждую из строк  
      // матрицы BT 
      for (j = 0; j < N; j++) 
      { 
        // j-я строка матрицы B 
        double sum = 0; 
        int ind3 = B.mRowIndex[j], ind4 = B.mRowIndex[j + 1]; 
        // Все ненулевые элементы строки j матрицы B 
        for (k = ind3; k < ind4; k++) 
        { 
          int bcol = B.mCol[k]; 
          int aind = temp[bcol]; 
          if (aind != -1) 
            sum += A.mValues[aind] * B.mValues[k]; 
        } 
        if (fabs(sum) > ZERO_IN_CRS) 
        { 
          columns[i].push_back(j); 
		            values[i].push_back(sum); 
          row_index[i]++; 
        } 
      } 
    } 
    delete [] temp; 
  } 
 
  int NZ = 0; 
  for(i = 0; i < N; i++) 
  { 
    int tmp = row_index[i]; 
    row_index[i] = NZ; 
    NZ += tmp; 
  } 
  row_index[N] = NZ; 
 
  InitializeMatrix(N, NZ, C); 
 
  int count = 0; 
  for (i = 0; i < N; i++) 
  { 
    int size = columns[i].size(); 
    memcpy(&C.mCol[count], &columns[i][0],  
           size * sizeof(int)); 
    memcpy(&C.mValues[count], &values[i][0],  
           size * sizeof(double)); 
    count += size; 
  } 

	for (size_t i = 0; i <= N; i++) 
	{
		C.mRowIndex[i] = row_index[i]; 
	}
 
  delete [] row_index; 
  delete [] columns; 
  delete [] values; 
}

bool SparseMatrixCRS::IsNonZero(size_t i, size_t j)
{
	bool result = false;

	for (size_t k = mRowIndex[i]; k < mRowIndex[i + 1]; k++)
	{
		if (mCol[k] == j && fabs(mValues[k]) != 0.0/*ZERO_IN_CRS*/)
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

void SparseMatrixCRS::Set(int i, int j, double val)
{
    //TO DO
}

SparseMatrixCRS::~SparseMatrixCRS()
{

}

void SparseMatrixCRS::MultiplyFullMatrix(vector<vector<double>> &a, vector<vector<double>> &b, vector<vector<double>> &c, size_t n)
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			double sum = 0.0;
			for (size_t k = 0; k < n; k++)
			{
				sum += a[i][k] * b[k][j];
			}
			c[i][j] = sum;
		}
	}
}

double SparseMatrixCRS::norm_of_matrix(vector<vector<double>> &a, size_t nn)
{
	double summ = 0, max_summ = INT_MIN;
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
		for (int j = rind[i]; j < rind[i + 1]; j++)
		{
			A[i][colum[j]] = val[j];
		}
	}
}