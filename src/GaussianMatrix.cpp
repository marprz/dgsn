#include <iostream>
#include <vector>
#include <stdlib.h>

#include "GaussianMatrix.h"

GaussianMatrix::GaussianMatrix()
{
    if( mData.size() > 4 )
    {
        std::cout << "overdetermination!" << std::endl;
        overdetermined();
    }
}
    
GaussianMatrix::GaussianMatrix( std::vector< std::vector< double > > aData )
: mData( aData )
{
    if( mData.size() > 4 )
    {
        std::cout << "overdetermination!" << std::endl;
        overdetermined();
    }
    makeGaussian();
}

void GaussianMatrix::setData( std::vector< std::vector< double > >  aData )
{
    mData = aData;
}
 
void GaussianMatrix::swapRows( int first, int second )
{
    if( first == second )
        return;
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

// very dirty code:
void GaussianMatrix::overdetermined()
{
    Matrix A;
    int N = getRowsNb();
    for( int i=0; i<N; ++i )
    {
        Row row( 4, 0 );
        A.push_back( row );
    }
    Row b( N, 0 );

    for( int i=0; i<N; ++i )
    {
        for( int j=0; j<N-1; ++j )
        {
            A[i][j] = mData[i][j];
        }
        b[i] = mData[i][4];
    }
    
    Matrix ATA;
    for( int i=0; i<4; ++i )
    {
        Row row( 4,0 );
        ATA.push_back( row );
    }
    Row ATb( 4, 0 );
    
    for( int i=0; i<4; ++i )
    {
        for( int j=0; j<4; ++j )
        {
            for( int k=0; k<N; ++k )
            {
                ATA[i][j] = ATA[i][j] + A[k][i]*A[k][j];
            }
        }
    }

    for( int i=0; i<4; ++i )
    {
        for( int k=0; k<N; ++k )
        {
            ATb[i] = ATb[i] + A[k][i]*b[k];
        }
    }

    mData.clear();
    mData.resize( 4 );
    for( int i=0; i<4; ++i )
    {
        mData[i].resize( 5 );
    }

    for( int i=0; i<4; ++i )
    {
        for( int j=0; j<4; ++j )
        {
            mData[i][j] = ATA[i][j];
        }
        mData[i][4] = ATb[i];
    }

}

void GaussianMatrix::makeGaussian()
{
   for( int k=0; k<getRowsNb(); ++k )
   {
       // Find pivot for column k:
       int iMaxPos = k;
       double iMax = (mData.at(k)).at(k);
       for( int j=k; j<getRowsNb(); ++j )
       {
           if( iMax < abs((mData.at(j)).at(k)) )
           {
               iMax = abs((mData.at(j)).at(k));
               iMaxPos = j;
           }
       }
        if( iMax == 0 )
            std::cout << "ERROR! MATRIX IS SINGULAR!" << std::endl;
        swapRows( iMaxPos, k );
        
        for( int i=k+1; i<getRowsNb(); ++i )
        {
            for( int j=k+1; j<getColsNb(); ++j )
            {
                (mData.at(i)).at(j) -= (mData.at(k)).at(j)*( (mData.at(i)).at(k) / (mData.at(k)).at(k));
            }
            (mData.at(i)).at(k) = 0;
        }
   }

   // '1' na diagonali
   for( int k=0; k<getRowsNb(); ++k )
   {
        if( (mData.at(k)).at(k) != 0 )
        {
            double temp = (mData.at(k)).at(k);
            for( int j=0; j<getColsNb(); ++j )
            {
                (mData.at(k)).at(j) = (mData.at(k)).at(j)/temp;
            }
        }
   }

   for( int k=getRowsNb()-2; k>=0; --k ) // start: przedostatni wiersz
   {
        for( int j=k+1; j<getRowsNb(); ++j )
        {
            if( (mData.at(k)).at(j) != 0 )
            {
                // odjecie j-tego wiersza pomnozonego przez (mData.at(k)).at(j)
                if( (mData.at(j)).at(j) != 0 )
                {
                    double temp = (mData.at(k)).at(j);
                    for( int i=k; i<getColsNb(); ++i )
                    {
                        (mData.at(k)).at(i) = (mData.at(k)).at(i) - (mData.at(j)).at(i)*temp;
                    }
                }
            }
        }
   }

   printData();
}

void GaussianMatrix::makeGaussian2()
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
      std::cout.precision(15);
	  std::cout << (mData.at(i)).at(j) << " " ;
	}
	std::cout << std::endl;
    }
}

