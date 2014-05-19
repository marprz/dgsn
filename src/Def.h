#ifndef DEF_H
#define DEF_H

#include <iostream>
#include <vector>
#include <array> 
#include <string>

/// Row type
typedef std::vector< double > Row;
/// Matrix type
typedef std::vector< Row > Matrix;

/**
 * @brief Station class
 */
class Station
{
  public:
    /**
     * @brief Constructor.
     * @param ax coordinate x
     * @param ay coordinate y
     * @param az coordinate z
     * @param at time
     */
    Station( double ax, double ay, double az, double at )
    : x( ax )
    , y( ay )
    , z( az )
    , t( at )
    {
    }

    /**
     * @brief Sets radius.
     * @param aT0 time of sending signal
     */
    void setR( double aT0 ) 
    {
        double clight = 299792458; // [m/s]
	dt = t-aT0;
	r = dt*clight;
    }

    /**
     * @brief Returns coordinate x.
     * @return coordinate x
     */
    double getX() const { return x; }

    /**
     * @brief Returns coordinate y.
     * @return coordinate y
     */
    double getY() const { return y; }

    /**
     * @brief Returns coordinate z.
     * @return coordinate z
     */
    double getZ() const { return z; }

    /**
     * @brief Returns radius.
     * @return radius r
     */
    double getR() const { return r; }

  private:
    double x;  ///< Coordinate x.
    double y;  ///< Coordinate y.
    double z;  ///< Coordinate z.
    double t;  ///< Time of receiving signal.
    double r;  ///< Distance of vehicle from the station at time t-dt.
    double dt; ///< Difference between time of receiving signal by ground station and sending by the vehicle.
};

/**
 * @brief Stations class which stores stations.
 */
class Stations
{
  public:
    /**
     * @brief Constructor.
     */
    Stations(){}

    /**
     * @brief Adding a single station.
     * @param aStation a single station added to the storage
     */
    void addStation( Station aStation )
    {
    	aStation.setR( mT );
        mStations.push_back( aStation );
        std::vector< Station >::iterator iter;
    }

    /**
     * @brief Adding of station using its position and time of receiving signal.
     * @param ax coordinate x
     * @param ay coordinate y
     * @param az coordinate z
     * @param at time of receiving signal
     */
    void addStations( double ax, double ay, double az, double at )
    {
	    Station aStation( ax, ay, az, at );
	    mStations.push_back( aStation );
    }

    /**
     * @brief Sets time of receiving signal.
     * @param aT time of receiving signal
     */
    void setTime( double aT )
    {
	    mT = aT;
    }

    /**
     * @brief Getting number of ground stations.
     * @return Number of stations.
     */
    int size() const { return mStations.size(); }

    /**
     * @brief Getting stations of given index.
     * @param iStationNb index of station
     * @return Station with given index.
     */
    Station getStation( int iStationNb )
    {
        return mStations.at(iStationNb ); 
    }

    /**
     * @brief Printing information about all stations.
     */
    void const printStations()
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

/**
 * @brief Class which stores list of positions.
 */
class PositionsList
{
  public:
    /**
     * @brief Default constructor.
     */
    PositionsList(){}

    /**
     * @brief Constructor.
     * @param aPosition first position which will be on the list
     */
    PositionsList( Position aPosition ) { mPositions.push_back( aPosition ); }

    /**
     * @brief Adding position to the list.
     * @param aPosition position which is added
     */
    void addPosition( Position aPosition ) { mPositions.push_back( aPosition ); }

    /**
     * @brief Adding position using its coordinates.
     * @param aPosVec vector of coordinates.
     */
    void addPosition( std::vector< double > aPosVec ) 
    {
    	Position aPosition = { aPosVec.at(0), aPosVec.at(1), aPosVec.at(2), aPosVec.at(3) };
        mPositions.push_back( aPosition ); 
    }

    /**
     * @brief Adding vector of positions.
     * @param aPosVec vector of positions
     */
    void addPositions( std::vector< Position > aPositions )
    {
        std::vector< Position >::iterator iter;
        for( iter = aPositions.begin(); iter != aPositions.end(); ++iter )
        {
            mPositions.push_back( *iter );
        }
    }

    /**
     * @brief Getting position with given index.
     * @param aIndex index of position which is returned
     */
    Position getPosition( int aIndex ) const { return mPositions.at(aIndex); }

    /**
     * @brief Adding on the end of the PositionsList another PositionsList.
     * @param aList PositionList which is added 
     */
    void addPositions( PositionsList aList )
    {
        for( int i = 0; i<aList.size(); ++ i )
        {
    	    mPositions.push_back( aList.getPosition(i) );
	    }
    }

    /**
     * @brief Getting number of stored positions.
     */
    int size() const { return mPositions.size(); }

    /**
     * @brief Getting coordinate x of given position.
     * @param aIndex index of position which coordinate is returned
     */
    double getX( int aIndex ) { return mPositions.at(aIndex).at(0); }

    /**
     * @brief Getting coordinate y of given position.
     * @param aIndex index of position which coordinate is returned
     */
    double getY( int aIndex ) { return mPositions.at(aIndex).at(1); }

    /**
     * @brief Getting coordinate z of given position.
     * @param aIndex index of position which coordinate is returned
     */
    double getZ( int aIndex ) { return mPositions.at(aIndex).at(2); }

    /**
     * @brief Getting distance between vehicle and ground station.
     * @param aIndex index of position
     */
    double getR( int aIndex ) { return mPositions.at(aIndex).at(3); }

    /**
     * @brief Printing data of all positions on the list.
     */
    void printPositions()
    {
    	std::cout << std::endl << "CALCULATED POSITIONS: " << std::endl;
        std::vector< Position >::iterator iter;
    	for( iter = mPositions.begin(); iter != mPositions.end(); ++iter )
	    {
	        std::cout << (*iter).at(0) << ", " << (*iter).at(1) << ", " << (*iter).at(2) << std::endl;
    	}
    }
    
    /**
     * @brief Printing average of all stored positions.
     */
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
    std::vector< Position > mPositions; ///< Container of positions.
};
#endif
