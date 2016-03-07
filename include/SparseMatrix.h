#include <vector>
#include <fstream>

class SparseMatrixBase
{
public:
    SparseMatrixBase();
    SparseMatrixBase(int N);
    SparseMatrixBase(const string &filename);
    virtual double Get(int i, int j) = 0;
    virtual void Set(int i, int j, double val) = 0;
    virtual ~SparseMatrixBase();

    int GetSize() { return mN; }
protected:
    int mN;
};

class SparseMatrixMtx : public SparseMatrixBase
{
public:
    SparseMatrixMtx();
    SparseMatrixMtx(const string &filename);

    virtual double Get(int i, int j);
    virtual void Set(int i, int j, double val) ;
    virtual ~SparseMatrixMtx();
    
private:
    std::vector<int> mIs;
    std::vector<int> mJs;
    std::vector<double> mValues;
};

class SparseMatrixCRS :  : public SparseMatrixBase
{
public:
    SparseMatrixCRS();
    SparseMatrixCRS(const string &filename);

    virtual double Get(int i, int j);
    virtual void Set(int i, int j, double val) ;
    virtual ~SparseMatrixCRS();
    
    void fromMtx(SparseMatrixMtx &matr);
private:
    int mCountNonZero;
    std::vector<double> mValues;
    std::vector<int> mCol;
    std::vector<int> mRowIndex;
};

