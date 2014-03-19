#include <iostream>
#include <stdint.h>
#include <vector>

uint64_t initialBitCombination( int n, int k )
{
  uint64_t result = 1;
  uint64_t temp = 1;
  for( int i = 1; i<k; ++i )
  {
    temp = temp<<1;
    result = result + temp;
  }
  return result;
}

// algorithm: http://www-graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
// t = v | (v - 1); // t gets v's least significant 0 bits set to 1
// Next set to 1 the most significant bit to change,
// set to 0 the least significant ones, and add the necessary 1 bits. 
uint64_t nextBitCombination( uint64_t currentBitCombination )
{
  uint64_t t = ( currentBitCombination | ( currentBitCombination -1 )) + 1;
  return t | ((((t & -t) / ( currentBitCombination & -currentBitCombination )) >> 1) - 1);

//  linux version (faster):
//  uint64_t t = currentBitCombination | ( currentBitCombination -1 );
//  return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz( currentBitCombination ) + 1));
}

std::vector< std::vector< int > > getStationsCombinations( int n, int k )
{
  uint64_t v = initialBitCombination( n, k ); // current permutation of bits
  uint64_t lastCombination = ( v<<(n-k) );  
  uint64_t w = v; 

  std::vector< uint64_t > bitCombinations;
  std::vector< std::vector< int > > stationsCombinations;

  int counter = 1;
  while( v <= lastCombination )
  { 
    bitCombinations.push_back( v );
    w = nextBitCombination( v );
    v = w;
    ++counter;
  }

  std::vector< uint64_t >::iterator iter; 
  for( iter = bitCombinations.begin(); iter != bitCombinations.end(); ++iter )
 {
    uint64_t curr = *iter;
    uint64_t temp = 1;
    std::vector< int > currentStations;
    for( int i=0; i<n; ++i )
    {
        if( (temp & curr) == temp )
	{
	  currentStations.push_back(i); // i in 0,1,...,n-1
	}
	temp = temp<<1;
    } 
    stationsCombinations.push_back( currentStations );
  }
  return stationsCombinations;
}

