#include "SparseMatrix.h"


SparseMatrixMtx::SparseMatrixMtx()
{

}

SparseMatrixMtx::SparseMatrixMtx(const string &filename)
{
    ifstream in(filename);

    int m, n;
    in >> m >> n >> mCountNonZero;

    mIs.resize(mCountNonZero);
    mJs.resize(mCountNonZero);
    mValues.resize(mCountNonZero);

    for (int i = 0; i < mCountNonZero; i++)
    {
        in >> mIs[i] >> mJs[i] >> mValues[i];
    }

    in.close();
}

double SparseMatrixMtx::Get(int i, int j)
{
    double result = 0.0;

    for (int k = 0; k < mCountNonZero; k++)
    {
        if (mIs[k] == i && mJs[k] == j)
        {
            result = mValues[k];
        }
    }

    return result;
}


void SparseMatrixMtx::Set(int i, int j, double val)
{
    //TO DO
}

SparseMatrixMtx::~SparseMatrixMtx()
{

}


//

SparseMatrixCRS::SparseMatrixCRS()
{

}

SparseMatrixCRS::SparseMatrixCRS(const string &filename)
{
    ifstream in(filename);

    //TO DO

    in.close();
}

double SparseMatrixCRS::Get(int i, int j)
{
    double result = 0.0;

    for (int k = 0; k < mCountNonZero; k++)
    {
        if (mIs[k] == i && mJs[k] == j)
        {
            result = mValues[k];
        }
    }

    return result;
}

void SparseMatrixCRS::fromMtx(SparseMatrixMtx &matr)
{

}

void SparseMatrixCRS::Set(int i, int j, double val)
{
    //TO DO
}

SparseMatrixCRS::~SparseMatrixCRS()
{

}