#include "ILU.h"
#include <vector>
#include "Multiply.h"
#include <iostream>
#include <iomanip>

using namespace std;

void ILU::Compute(SparseMatrixCRS &A, SparseMatrixCRS &M, int p, int flags)
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

    switch (flags)
    {
        case 0:
            MultiplyNaive(L, UT, M);
            break;
        case 1:
            MultiplyOpenMP(L, UT, M);
            break;
        case 2:
            MultiplyTBB(L, UT, M);
            break;
        case 3:
            MultiplyCilk(L, UT, M);
            break;
        default:
            MultiplyNaive(L, UT, M);
    };
}

bool ILU::isCorrectMatrix(SparseMatrixCRS &A)
{
    bool result = true;

    for (size_t i = 0; i < A.mN && result; i++)
    {
        size_t j1 = A.mRowIndex[i], j2 = A.mRowIndex[i + 1];
        for (size_t j = j1; j < j2 && result; j++)
        {
            if (A.mCol[j] == i && fabs(A.mValues[j]) == 0.0)
            {
                result = false;
            }
        }
    }

    return result;
}

bool ILU::CheckAinM(SparseMatrixCRS &A, SparseMatrixCRS &M)
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
                
                cout << fixed << setprecision(13) << a << " " << m << " - (" << i <<", " << j << ")" << endl;;
                result = false;
                break;
            }
        }
    }

    return result;
}

bool ILU::CheckInverse(SparseMatrixCRS &A, SparseMatrixCRS &M)
{
    bool result = false;
    Matrix full_A(A.mN);
    Matrix inverse_full_A(A.mN);
    Matrix full_M(A.mN);
    Matrix inverse_full_M(A.mN);
    Matrix full_M_A(A.mN);
    Matrix inverse_full_M_A(A.mN);

    full_A.fromCRS(A);
    double norm_A = full_A.NormMatrix();

    full_M.fromCRS(M);
    double norm_M = full_M.NormMatrix();

    inverse_full_M = full_M.GetInverseMatrix();
    double norm_inv_M = inverse_full_M.NormMatrix();

    MultiplyFullMatrix(inverse_full_M, full_A, full_M_A);
    double norm_M_A = full_M_A.NormMatrix();

    inverse_full_A = full_A.GetInverseMatrix(); 
    double norm_inv_A = inverse_full_A.NormMatrix();

    double ma = norm_A * norm_inv_A;

    inverse_full_M_A = full_M_A.GetInverseMatrix();
    double norm_inverse_M_A = inverse_full_M_A.NormMatrix();

    double mmm = norm_inverse_M_A * norm_M_A;
    
    result = ma > mmm;
    printf("%lf %lf \n", norm_inverse_M_A, norm_M_A);
    printf("ma = %lf, mmm = %lf\n", ma, mmm);

    return result;
}

