#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <cassert>
#include <cmath>

using namespace std;

template<typename T>
class Mat{
private:
    unique_ptr<T[]> data;
    int _size;
public:
    Mat(int size);
    Mat(const Mat<T> &m);
    Mat<T> & operator=(const Mat<T> &m);
    Mat<T> & operator*=(const Mat<T> &m);
    T & operator()(unsigned int i, unsigned int j);
    T operator()(unsigned int i, unsigned int j) const;
    Mat<T> transpose();
    bool isDiagonal();
    void print();
};

template<typename T>
Mat<T>::Mat(const int size) : _size(size)
{
    data = make_unique<T[]>(size * size);
}
template<typename T>
Mat<T>::Mat(const Mat<T> &m) : _size(m._size)
{
    data = make_unique<T[]>(_size * _size);
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            data[i * _size + j] = m(i, j);
        }
    }
}
template<typename T>
inline T & Mat<T>::operator()(unsigned int i, unsigned int j)
{
    return data[i * _size + j];
}
template<typename T>
inline T Mat<T>::operator()(unsigned int i, unsigned int j) const
{
    return data[i * _size + j];
}
template<typename T>
Mat<T> & Mat<T>::operator=(const Mat& m)
{
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            data[i * _size + j] = m(i, j);
        }
    }
    return *this;
}
template<typename T>
Mat<T> & Mat<T>::operator*=(const Mat& m)
{
    Mat tmp(_size);
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            for(int k = 0; k < _size; ++k) {
                tmp(i, j) += data[i * _size + k] * m(k, j);
            }
        }
    }
    return (*this = tmp);
}
template<typename T>
Mat<T> Mat<T>::transpose()
{
    Mat<T> tmp(_size);
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            tmp(j, i) = data[i * _size + j];
        }
    }
    return tmp;
}
template<typename T>
bool Mat<T>::isDiagonal()
{
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            if(i != j && data[i * _size + j] != 0) return false;
        }
    }
    return true;
}
template<typename T>
void Mat<T>::print()
{
    for(int i = 0; i < _size; ++i) {
        for(int j = 0; j < _size; ++j) {
            cout << /*setprecision(3) << */data[i * _size + j] << " ";
        }
        cout << endl;
    }
}

template<typename T>
Mat<T> operator*(const Mat<T>& m1, const Mat<T>& m2)
{
    Mat<T> tmp(m1);
    return (tmp *= m2);
}

template<int N>
class Jacob{
private:
    const int FTRDIM = N;
    ifstream fpcov;
    ifstream fpout;
    Mat<double> cov;
    Mat<double> A;
    Mat<double> Eigen;

public:
    Jacob(string input_file, string output_file);
    void readCov(void);
    tuple<int,int,double> searchMaxIdx(Mat<double> &m);
    void run(void);
};

template<int N>
Jacob<N>::Jacob(string input_file, string output_file) : cov(FTRDIM), A(FTRDIM), Eigen(FTRDIM)
{
    fpcov.open(input_file, ios::in);
    fpout.open(output_file, ios::out); 
}

template<int N>
void
Jacob<N>::readCov(void)
{
    for(int i = 0; i < FTRDIM; ++i)
        for(int j = 0; j < FTRDIM; ++j) fpcov >> cov(i, j);
    A = cov;

    // Make Identity 
    for(int i = 0; i < FTRDIM; ++i)
        for(int j = 0; j < FTRDIM; ++j) Eigen(i, j) = (i == j) ? 1 : 0;
}

template<int N>
tuple<int,int,double>
Jacob<N>::searchMaxIdx(Mat<double> &m)
{
    tuple<int,int,double> res;
    double maximum = 0;

    for(int i = 0; i < FTRDIM; ++i) {
        for(int j = 0; j < FTRDIM; ++j) {
            if(i == j) continue;
            if(maximum < abs(m(i, j))) {
                maximum = abs(m(i, j));
                res = make_tuple(i, j, maximum);
            }
        }
    }
    return res;
}

template<int N>
void
Jacob<N>::run(void)
{
    Mat<double> U(FTRDIM);
    tuple<int,int,double> loc;
    int n = 0;

    while(1) {
        loc = searchMaxIdx(A);
        cout << "n = " << n << ", max = " << get<2>(loc) << endl;
        if(get<0>(loc) > get<1>(loc)) swap(get<0>(loc), get<1>(loc));
        double theta = 0.5 * atan(2 * A(get<0>(loc), get<1>(loc)) / (A(get<1>(loc), get<1>(loc)) - A(get<0>(loc), get<0>(loc))));
        for(int i = 0; i < FTRDIM; ++i) {
            for(int j = 0; j < FTRDIM; ++j) {
                if(i == get<0>(loc) && j == get<0>(loc)) U(i, j) = cos(theta);
                else if(i == get<0>(loc) && j == get<1>(loc)) U(i, j) = sin(theta);
                else if(i == get<1>(loc) && j == get<0>(loc)) U(i, j) = -1.0 * sin(theta);
                else if(i == get<1>(loc) && j == get<1>(loc)) U(i, j) = cos(theta);
                else if(i == j) U(i, j) = 1;
                else U(i, j) = 0;
            }
        }
        A = U.transpose() * A * U;
        Eigen = Eigen * U;
        if(A.isDiagonal()) return;
        n++;
    }
}

int
main(void)
{
    constexpr int S = 196;
    unique_ptr<Jacob<S>> j = make_unique<Jacob<S>>("cov.dic", "eigenvv.dat");

    j->readCov();
    j->run();
}
