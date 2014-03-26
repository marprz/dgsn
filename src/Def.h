#ifndef DEF_H
#define DEF_H

#include <iostream>
#include <vector>
#include <array> 
#include <string>

typedef std::vector< double > Row;
typedef std::vector< Row > Matrix;

class Station
{
  public:
    Station( double ax, double ay, double az, double at )
    : x( ax )
    , y( ay )
    , z( az )
    , t( at )
    {
    }

    void setR( double aT0 ) 
    {
        double clight = 299792458; // [m/s]
	dt = t-aT0;
	r = dt*clight;
    }
    double getX() { return x; }
    double getY() { return y; }
    double getZ() { return z; }
    double getR() { return r; }

  private:
    double x;
    double y;
    double z;
    double t;
    double r;
    double dt; 
};

class Stations
{
  public:
    Stations(){}
    void addStation( Station aStation )
    {
	aStation.setR( mT );
        mStations.push_back( aStation );
        std::vector< Station >::iterator iter;
    }

    void addStations( double ax, double ay, double az, double at )
    {
	Station aStation( ax, ay, az, at );
	mStations.push_back( aStation );
    }

    void setTime( double aT )
    {
	mT = aT;
    }

    int size() { return mStations.size(); }

    Station getStation( int iStationNb )
    {
      return mStations.at(iStationNb ); 
    }

    void printStations()
    {
        std::vector< Station >::iterator iter;
	std::cout << std::endl;
        for( iter = mStations.begin(); iter != mStations.end(); ++iter )
        {
		std::cout << "***" << (*iter).getX() << ", " << (*iter).getY() << ", " << (*iter).getZ() << ", " << (*iter).getR() << std::endl;
        }
	std::cout << std::endl;
    }
//  private:
    std::vector< Station > mStations;  
    double mT;
};

typedef std::array< double, 4 > Position; // xs,ys,zs,rs from Apollonius 

class PositionsList
{
  public:
    PositionsList(){}
    PositionsList( Position aPosition ) { mPositions.push_back( aPosition ); }

    void addPosition( Position aPosition ) { mPositions.push_back( aPosition ); }
    void addPosition( std::vector< double > aPosVec ) 
    {
	Position aPosition = { aPosVec.at(0), aPosVec.at(1), aPosVec.at(2), aPosVec.at(3) };
        mPositions.push_back( aPosition ); 
    }
    void addPositions( std::vector< Position > aPositions )
    {
        std::vector< Position >::iterator iter;
        for( iter = aPositions.begin(); iter != aPositions.end(); ++iter )
        {
            mPositions.push_back( *iter );
        }
    }

    Position getPosition( int aIndex ) { return mPositions.at(aIndex); }
    void addPositions( PositionsList aList )
    {
        for( int i = 0; i<aList.size(); ++ i )
        {
	    mPositions.push_back( aList.getPosition(i) );
	}
    }
    int size() { return mPositions.size(); }
    double getX( int aIndex ) { return mPositions.at(aIndex).at(0); }
    double getY( int aIndex ) { return mPositions.at(aIndex).at(1); }
    double getZ( int aIndex ) { return mPositions.at(aIndex).at(2); }
    double getR( int aIndex ) { return mPositions.at(aIndex).at(3); }

    void printPositions()
    {
	std::cout << std::endl << "CALCULATED POSITIONS: " << std::endl;
        std::vector< Position >::iterator iter;
	for( iter = mPositions.begin(); iter != mPositions.end(); ++iter )
	{
	    std::cout << (*iter).at(0) << ", " << (*iter).at(1) << ", " << (*iter).at(2) << std::endl;
	}
    }
    
    void printAveragePosition()
    {
        double x = 0, y = 0, z = 0;
	int positionsNb = mPositions.size();
	std::vector< Position >::iterator iter;
        for( iter = mPositions.begin(); iter != mPositions.end(); ++iter )
        {
	    x += (*iter).at(0)/positionsNb;
	    y += (*iter).at(1)/positionsNb;
	    z += (*iter).at(2)/positionsNb;
        }
	std::cout << "Average position: " << x << ", " << y << ", " << z << std::endl;
    }

  private:
    std::vector< Position > mPositions;
};
#endif
