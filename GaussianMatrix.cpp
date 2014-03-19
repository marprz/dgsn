#include <iostream>
#include <vector>
#include <stdlib.h>

#include "GaussianMatrix.h"

GaussianMatrix::GaussianMatrix( std::vector< std::vector< double > > aData )
: mData( aData )
{
    makeGaussian();
}

void GaussianMatrix::setData( std::vector< std::vector< double > >  aData )
{
    mData = aData;
}
 
void GaussianMatrix::swapRows( int first, int second )
{
    Row aRow = mData.at(first);
    mData.at(first) = mData.at(second);
    mData.at(second) = aRow;
}

int GaussianMatrix::findMaxRow( int col )
{
    int ret = col; 
    for( int i = col; i<getRowsNb(); ++i )
    {
	if( abs(get(i,col)) > abs(get(ret,col)) )
	{
	  ret = i;
	}
    }
    return ret;
}

void GaussianMatrix::multiply( int col )
{
    double coefficient;
    for( int i = 0; i<getRowsNb(); ++i )
    {
      coefficient = get(i,col);
      if( coefficient != 1 && coefficient != 0 )
      {
        for( int j=0; j<getColsNb(); ++j )
        {
          (mData.at(i)).at(j) = get(i,j)/coefficient;
        }
      }
    }
}

void GaussianMatrix::subtractRow( int row )
{
    for( int i=0; i<getRowsNb(); ++i )
    {
	if( i!=row && get(i,row)!=0 )
	{
	  for( int j=0; j<getColsNb(); ++j )
	  {
	    (mData.at(i)).at(j) = get(i,j) - get(row,j);
	  }
	}
    }
}

void GaussianMatrix::makeDiagonalOnes()
{
    double coefficient;
    for( int i = 0; i < mData.size(); ++i )
    {
	if( get(i,i) != 0 && get(i,i) !=1 )
	{
	  coefficient = get(i,i);
	  for( int j = 0; j < mData.at(0).size(); ++j )
	  {
	    (mData.at(i)).at(j) = get(i,j)/coefficient;
	  }
	}
    }
}

int GaussianMatrix::getRowsNb()
{
    return mData.size();
}

int GaussianMatrix::getColsNb()
{
    return (mData.at(0)).size();
}

void GaussianMatrix::makeGaussian()
{
    int max_index;
    for( int i=0; i<getRowsNb(); ++i )
    {
	max_index = findMaxRow( i );
        if( max_index != i )
	{
	  swapRows( max_index, i );
	}
	multiply( i );
	subtractRow(i);
    }
    makeDiagonalOnes();
}

void GaussianMatrix::printData()
{
    std::cout << std::endl;
    for( int i=0; i<getRowsNb(); ++i )
    {
	for( int j=0; j<getColsNb(); ++j )
	{
	  std::cout << (mData.at(i)).at(j) << " " ;
	}
	std::cout << std::endl;
    }
}

