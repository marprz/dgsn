#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <stdlib.h>

#include "Def.h"

class GaussianMatrix
{
  public:
    GaussianMatrix(){}
    GaussianMatrix( std::vector< std::vector< double > >  aData );
    void setData( std::vector< std::vector< double > >  aData );
    void makeGaussian();
    void printData();
    double operator() ( int r, int c ){ return (mData.at(r-1)).at(c-1); } 
    int getRowsNb();
    int getColsNb();

  private:
    double get( int row, int col ){ return (mData.at(row)).at(col); }
    void swapRows( int first, int second );
    int findMaxRow( int col ); 
    void multiply( int col );
    void subtractRow( int row );
    void makeDiagonalOnes();
    Matrix mData;
};

#endif
