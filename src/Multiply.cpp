#include "Multiply.h"

#include <omp.h>
#include <tbb/tbb.h>

using namespace std;
using namespace tbb;

void MultiplyFullMatrix(Matrix &a, Matrix &b, Matrix &c)
{
    for (size_t i = 0; i < a.mN; i++)
    {
        for (size_t j = 0; j < a.mN; j++)
        {
            double sum = 0.0;
            for (size_t k = 0; k < a.mN; k++)
            {
                sum += a.mValues[i][k] * b.mValues[k][j];
            }
            c.mValues[i][j] = sum;
        }
    }
}

void MultiplyNaive(SparseMatrixCRS &A, SparseMatrixCRS &BT, SparseMatrixCRS &C)
{
	size_t N = A.mN; 
	
	vector<int> columns; 
	vector<double> values; 
	vector<int> rowIndex; 
 
	size_t rowNZ = 0; 
	rowIndex.push_back(0); 
	for (size_t i = 0; i < N; i++) 
	{ 
		rowNZ = 0; 
		for (size_t j = 0; j < N; j++) 
		{
			double sum = 0; 
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
		rowIndex.push_back(rowNZ + rowIndex[i]); 
	} 
 
	C.InitializeMatrix(N, columns.size()); 
 
	for (size_t j = 0; j < columns.size(); j++) 
	{ 
		C.mCol[j] = columns[j]; 
		C.mValues[j] = values[j]; 
	} 

    copy(rowIndex.begin(), rowIndex.end(), C.mRowIndex.begin());
} 

void MultiplyOpenMP(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
	size_t N = A.mN; 
 
    vector<size_t> *columns = new vector<size_t>[N]; 
    vector<double> *values = new vector<double>[N]; 
 
    vector<size_t> rowIndex(N + 1, 0);
    #pragma omp parallel 
    { 
        int *temp = new int[N]; 
        #pragma omp for private(j, k) schedule(static, chunk)
        for (size_t i = 0; i < N; i++) 
        { 
            memset(temp, -1, N * sizeof(int)); 
            size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i + 1]; 
            for (size_t j = j1; j < j2; j++) 
            { 
                size_t col = A.mCol[j]; 
                temp[col] = j; 
            } 
            for (size_t j = 0; j < N; j++) 
            { 
                double sum = 0; 
                size_t k1 = B.mRowIndex[j], k2 = B.mRowIndex[j + 1]; 
                for (size_t k = k1; k < k2; k++) 
                { 
                    size_t bcol = B.mCol[k]; 
                    int aind = temp[bcol]; 
                    if (aind != -1) 
                    {
                        sum += A.mValues[aind] * B.mValues[k]; 
                    }
                }
                if (fabs(sum) > ZERO_IN_CRS) 
                { 
                    columns[i].push_back(j); 
                    values[i].push_back(sum); 
                    rowIndex[i]++; 
                } 
            } 
        } 
        delete[] temp; 
    } 
 
    size_t NZ = 0; 
    for(size_t i = 0; i < N; i++) 
    { 
        int tmp = rowIndex[i]; 
        rowIndex[i] = NZ; 
        NZ += tmp; 
    } 
    rowIndex[N] = NZ; 
 
    C.InitializeMatrix(N, NZ);
 
    size_t count = 0; 
    for (size_t i = 0; i < N; i++) 
    { 
        size_t size = columns[i].size(); 
        for (size_t j = 0; j < size; j++)
        {
            C.mCol[count + j] = columns[i][j];
            C.mValues[count + j] = values[i][j];
        }
        count += size; 
    } 

    copy(rowIndex.begin(), rowIndex.end(), C.mRowIndex.begin());
 
    delete[] columns; 
    delete[] values; 
}

class Multiplicator 
{
private:
    SparseMatrixCRS A, B; 
    vector<size_t> *columns; 
    vector<double> *values; 
    vector<size_t> *rowIndex; 
public: 
    Multiplicator(SparseMatrixCRS &_A, SparseMatrixCRS &_B, 
      vector<size_t> *&_columns, vector<double> *&_values, 
      vector<size_t> *_rowIndex) : A(_A), B(_B), columns(_columns), 
      values(_values), rowIndex(_rowIndex) 
    {} 
 
    void operator()(const blocked_range<size_t> &r) const 
    { 
        size_t begin = r.begin(); 
        size_t end = r.end(); 
        size_t N = A.mN; 
 
        int *temp = new int[N]; 
        for (size_t i = begin; i < end; i++) 
        { 
            memset(temp, -1, N * sizeof(int)); 
            size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i + 1]; 
            for (size_t j = j1; j < j2; j++) 
            { 
                size_t col = A.mCol[j]; 
                temp[col] = j;
            } 
            for (size_t j = 0; j < N; j++) 
            { 
                double sum = 0; 
                size_t k1 = B.mRowIndex[j], k2 = B.mRowIndex[j + 1]; 
                for (size_t k = k1; k < k2; k++) 
                { 
                    size_t bcol = B.mCol[k]; 
                    int aind = temp[bcol]; 
                    if (aind != -1) 
                    {
                        sum += A.mValues[aind] * B.mValues[k]; 
                    }
                } 
                if (fabs(sum) > ZERO_IN_CRS) 
                { 
                    columns[i].push_back(j); 
                    values[i].push_back(sum); 
                    (*rowIndex)[i]++; 
                } 
            } 
        } 
        delete[] temp; 
    } 
}; 

void MultiplyTBB(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
    size_t N = A.mN;
 
    task_scheduler_init init(); 
 
    vector<size_t> *columns = new vector<size_t>[N]; 
    vector<double> *values = new vector<double>[N]; 
    vector<size_t> rowIndex(N + 1, 0); 
 
    size_t grainSize = 10; 
 
    parallel_for(blocked_range<size_t>(0, A.mN, grainSize), 
                 Multiplicator(A, B, columns, values, &rowIndex)); 
 
    size_t NZ = 0; 
    for(size_t i = 0; i < N; i++) 
    { 
        size_t tmp = rowIndex[i]; 
        rowIndex[i] = NZ; 
        NZ += tmp; 
    } 
    rowIndex[N] = NZ; 
 
    C.InitializeMatrix(N, NZ);
 
    size_t count = 0; 
    for (size_t i = 0; i < N; i++) 
    { 
        size_t size = columns[i].size(); 
        for (size_t j = 0; j < size; j++)
        {
            C.mCol[count + j] = columns[i][j];
            C.mValues[count + j] = values[i][j];
        }
        count += size; 
    }

    copy(rowIndex.begin(), rowIndex.end(), C.mRowIndex.begin());
}

void MultiplyCilk(SparseMatrixCRS &A, SparseMatrixCRS &B, SparseMatrixCRS &C)
{
	size_t N = A.mN; 
 
    vector<size_t>* columns = new vector<size_t>[N]; 
    vector<double> *values = new vector<double>[N]; 
    vector<size_t> rowIndex(N + 1, 0);   
    cilk_for (size_t i = 0; i < N; i++) 
    { 
        int temp[N];
        memset(temp, -1, N * sizeof(int)); 
        size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i + 1]; 
        for (size_t j = j1; j < j2; j++) 
        { 
            size_t col = A.mCol[j]; 
            temp[col] = j;
        } 

        for (size_t j = 0; j < N; j++) 
        { 
            double sum = 0; 
            size_t k1 = B.mRowIndex[j], k2 = B.mRowIndex[j + 1]; 
            for (size_t k = k1; k < k2; k++) 
            { 
                size_t bcol = B.mCol[k]; 
                int aind = temp[bcol]; 
                if (aind != -1)
                {
                    sum += A.mValues[aind] * B.mValues[k]; 
                }
            } 
            if (fabs(sum) > ZERO_IN_CRS) 
            { 
                columns[i].push_back(j); 
                values[i].push_back(sum); 
                rowIndex[i]++; 
            } 
        } 
    } 
 
    size_t NZ = 0; 
    for(size_t i = 0; i < N; i++) 
    { 
        size_t tmp = rowIndex[i]; 
        rowIndex[i] = NZ; 
        NZ += tmp; 
    } 
    rowIndex[N] = NZ; 
 
    C.InitializeMatrix(N, NZ); 
 
    size_t count = 0; 
    for (size_t i = 0; i < N; i++) 
    { 
        size_t size = columns[i].size(); 
        for (size_t j = 0; j < size; j++)
        {
            C.mCol[count + j] = columns[i][j];
            C.mValues[count + j] = values[i][j];
        }
        count += size; 
    } 

    copy(rowIndex.begin(), rowIndex.end(), C.mRowIndex.begin());
  
    delete[] columns; 
    delete[] values; 
}
