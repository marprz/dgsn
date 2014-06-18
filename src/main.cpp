#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>
#include <unistd.h>
#include <tuple>

#include <cmath> // square root
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "Def.h"
#include "Combinations.h"
#include "GaussianMatrix.h"
#include "Apollonius.h"

PositionsList solveApol( std::vector< Station > aStations );

// int - timestamp, Stations - positions of gs and distances from satellite
//typedef std::tuple< int, Stations > Signal;
//typedef std::vector< Signal > Signals;

void processSignalData()
{
    std::cout << "********* PROCESSING SIGNALS *****************************" << std::endl;
    std::vector< Signal >::iterator iter;
    for( iter = mSignals.begin(); iter != mSignals.end(); ++iter )
    {
        int satId = (*iter).getSatId();
        long double timestamp = (*iter).getTimestamp();

        if( (*iter).getSize() > 3 )
        {
            std::cout << "For satellite " << (*iter).getSatId() << " and timestamp " << (*iter).getTimestamp() << " there are " << (*iter).getSize() << " GS" << std::endl;
            (*iter).printSignal();
            Stations mStations;
            
            (*iter).convertSignalToStation( mStations );

            std::vector< std::vector< int > > stationsComb;

            std::cout << "Number of ground stations: " << mStations.size() << std::endl;

            int aNbStationsToOneMatrix;
//            aNbStationsToOneMatrix = 4; // TODO: temporary
            if( mStations.size() > 4 )
                aNbStationsToOneMatrix = 5;
            else
                aNbStationsToOneMatrix = 4;
            stationsComb = getStationsCombinations( mStations.size(), aNbStationsToOneMatrix );
	        std::cout << "Number of combinations: " << stationsComb.size() << std::endl;
            PositionsList xyzr;
            //PositionsList xyzr2;
        

            int it = 0;
//            std::vector< Signal >::iterator iter2;
            for( auto iter2 = stationsComb.begin(); it < aNbStationsToOneMatrix; ++iter2 )
            {
                std::cout << mStations.getStation(it).getX() << ", " << mStations.getStation(it).getY() << ", " << mStations.getStation(it).getZ() << ". R:" << mStations.getStation(it).getR() << std::endl;
                ++it;
            }
        	for( auto iter2 = stationsComb.begin(); iter2 != stationsComb.end(); ++iter2 )
	        {
                std::vector< Station > takenStations;
	            for( int j=0; j<aNbStationsToOneMatrix; ++j )
         	    {
                    takenStations.push_back(  mStations.getStation( (*iter2).at(j) ));
                }
                xyzr.addPositions( solveApol( satId, timestamp, takenStations ) );
            }
            xyzr.printPositions();
            mStations.printStations();
            xyzr.printAveragePosition();
/*
            std::vector< Station > takenStations2;
            for( int i=0; i<mStations.size(); ++i )
            {
                std::cout << "adding station: " << i << std::endl;
                takenStations2.push_back( mStations.getStation( i ) );
            }
            xyzr2.addPositions( solveApol( satId, timestamp, takenStations2 ) );
            xyzr2.printPositions();
*/
        }
        else
        {
            std::cout << "For satellite " << (*iter).getSatId() << " with sending time " << (*iter).getTimestamp() << " only " << (*iter).getSize() << " GS found" << std::endl;
        }
    }

}
void loadGSData( const char* lFileName )//, Signals& mSignals )
{
    Stations mStations;
    std::fstream lFile;
    lFile.open( lFileName, std::ios::in );
    if( !lFile.is_open() )
    	std::cout << "ERROR: problem with file" << std::endl;

    double ax, ay, az;
    long double at0, adt;
    int satId;
    std::string op1;
   // lFile >> at;
  //  mStations.setTime( at );
    while( lFile >> ax >> ay >> az >> adt >> at0 >> satId >> op1  )
    {
    //    at0 = 0; // tymczasowo
        bool satKnown = false;
        std::vector< Signal >::iterator iter;
//        std::cout << "mSignals.size() == " << mSignals.size() << std::endl;
/*        if( !mSignals.empty() )
        {
            std::cout << "first signal for satellite " << satId << " and timestamp " << at0 << std::endl;
        }*/

        for( iter = mSignals.begin(); !satKnown && iter != mSignals.end(); ++iter )
        {
//            std::cout << "   Comparing timestamps: " << at0 << " and " << (*iter).getTimestamp() << std::endl;
            if( (*iter).getSatId() == satId && (*iter).getTimestamp() == at0 )
            {
                if (!(*iter).positionKnown( ax, ay, az ))
                {
//                std::cout << "Known satellite with id: " << satId << std::endl;
                    (*iter).addGroundStation( ax, ay, az, at0, adt );
                }
                satKnown = true; // bez tego tez sie da
            }
        }

        if( !satKnown )
        {
            Signal mSignal;
            mSignal.setSatId( satId );
            mSignal.setTimestamp( at0 );
//        std::cout << "Adding new satellite: " << ax << " " << ay << " " <<  az << " " <<  at0 << " " <<  adt << " " <<  satId << " " << op1 << std::endl;
            mSignal.addGroundStation( ax, ay, az, at0, adt ); 
            mSignals.push_back( mSignal );
          /*  double clight = 299792458;
            double ar = clight*(atR-at0);
    	    mStations.addStation( ax, ay, az, at0, atR ); // zostanie tylko signal
          */
        }
    }
   /* std::cout << "first signal: " << mSignals.begin()->getSize() << std::endl;

    // printing:
    std::cout << "Known signals (" << mSignals.size() << "): " << std::endl;
    std::vector< Signal >::iterator it;
    for( it = mSignals.begin(); it != mSignals.end(); ++it )
    {
        (*it).printSignal();
    }
    std::cout << "****************" << std::endl << std::endl;
*/
    lFile.close();
}

void loadFromDirectory( char* lDirName )//, Signals& mSignals )
{
    std::string file;
    DIR *dir;
    struct dirent *dirEnt;
    struct stat filestat;
    dir = opendir( lDirName );
    if( dir == NULL )
    {
        std::cout << "ERROR: Problem with directory" << std::endl;
    }
    else
    {
        while( dirEnt = readdir( dir ) )
        {
            file = std::string(lDirName) + "/" + dirEnt->d_name;

            if( stat( file.c_str(), &filestat )) continue;
            if( S_ISDIR( filestat.st_mode ) ) continue;

            std::cout << "Processing file: " << file << std::endl;
            loadGSData( file.c_str() ); //, mSignals );            

        }
        closedir( dir );
    }
}

/** loading data from file "stations.txt"
 *  @lFileName a name of file which contains data
 *  @param mStations a contener of stations 
 */
void loadStations( char* lFileName, Stations& mStations )
{
//    std::string lFileName("stations.txt");
    std::fstream lFile;
    lFile.open( lFileName, std::ios::in );
    //lFile.open( lFileName.c_str(), std::ios::in );
    if( !lFile.is_open() )
	std::cout << "ERROR: problem with file" << std::endl;

    double ax, ay, az, at;
    lFile >> at;
    mStations.setTime( at );
    while( lFile >> ax >> ay >> az >> at )
    {
    	mStations.addStation( Station( ax, ay, az, at ) );
    }
    lFile.close();
}

/** solving Apollonius problem for set of stations
 *  @aStations a vector which contains stations
 */
/*
PositionsList solveApol( std::vector< Station > aStations ) 
{
    PositionsList calculatedPositions;
    std::vector< std::vector< double > >  matrix;
    double M,N,P,Q,R,S,a,b,c,rs,xs,ys,zs;

    double x1 = aStations.at(0).getX();
    double y1 = aStations.at(0).getY();
    double z1 = aStations.at(0).getZ();
    double r1 = aStations.at(0).getR();

    double el1, el2, el3, el4;
    el2 = (aStations.at(0).getX())*(aStations.at(0).getX()) + (aStations.at(0).getY())*(aStations.at(0).getY()) + (aStations.at(0).getZ())*(aStations.at(0).getZ());
	el3 = (aStations.at(0).getR())*(aStations.at(0).getR());

    for( int i=1; i<=3; ++i )
    {
	el1 = (aStations.at(i).getX())*(aStations.at(i).getX()) + (aStations.at(i).getY())*(aStations.at(i).getY()) + (aStations.at(i).getZ())*(aStations.at(i).getZ());
	el4 = (aStations.at(i).getR())*(aStations.at(i).getR());

	Row row = { ((aStations.at(i).getX())-x1)*2,
		     ((aStations.at(i).getY())-y1)*2,
		     ((aStations.at(i).getZ())-z1)*2,
		     (aStations.at(i).getR()-r1)*2,
		     el1-el2+el3-el4 };
        matrix.push_back( row );
    }
    GaussianMatrix gaussMatrix( matrix );
    M = gaussMatrix(1,5);
    N = -(gaussMatrix(1,4));
    P = gaussMatrix(2,5);
    Q = -(gaussMatrix(2,4));
    R = gaussMatrix(3,5);
    S = -(gaussMatrix(3,4));

    a = N*N+Q*Q+S*S-1;
    b = 2*(M-x1)*N+2*(P-y1)*Q+2*(R-z1)*S-2*r1;
    c = (M-x1)*(M-x1)+(P-y1)*(P-y1)+(R-z1)*(R-z1)-r1*r1;

    double p1 = (-b-sqrt(b*b-4*a*c))/(2*a);
    double p2 = (-b+sqrt(b*b-4*a*c))/(2*a);
    rs = p1>p2 ? p1 : p2 ;
    xs = M+N*rs;
    ys = P+Q*rs;
    zs = R+S*rs;
    std::vector< double > solution = { xs, ys, zs, rs };
    calculatedPositions.addPosition( solution );

     
/*  commented for possible future changes:
    std::array< double, 2 > r = { p1, p2 };

    for( int i=0; i<2; ++i )
    {
        rs = r.at(i);
	std::cout << "r = " << rs << std::endl;
        xs = M+N*rs;
        ys = P+Q*rs;
        zs = R+S*rs;
        std::vector< double > solution = { xs, ys, zs, rs };
        calculatedPositions.addPosition( solution );
    }
*/
//    return calculatedPositions;
//}

int main( int argc, char* argv[] )
{
    Stations mStations;
//    Signals mSignals;
    bool test = false;
    int c = getopt( argc, argv, "fd:" );
    switch( c ){
    case 'd':
        loadFromDirectory( optarg );
        processSignalData();
        //loadFromDirectory( optarg, mSignals );
        break;
	case 'f':
	    loadStations( optarg, mStations );
        test = true;
	    break;
    default:
	    std::cout << "Missing file name! Use parameter -f filename." << std::endl;
	    return -1;
    }

    if( test )
    {
        std::vector< std::vector< int > > stationsComb;
        std::vector< Station > takenStations;

        if( mStations.size() > 3 )
        {
            std::cout << "Number of ground stations: " << mStations.size() << std::endl;
            stationsComb = getStationsCombinations( mStations.size(), 4 );
	        std::cout << "Number of combinations: " << stationsComb.size() << std::endl;
            PositionsList xyzr;
        
            std::vector< std::vector< int > >::iterator iter; 
        	for( iter = stationsComb.begin(); iter != stationsComb.end(); ++iter )
	        {
	            for( int j=0; j<4; ++j )
         	    {
                    takenStations.push_back(  mStations.getStation( (*iter).at(j) ));
                }
                xyzr.addPositions( solveApol( 0, 0, takenStations ) );
            }
            xyzr.printAveragePosition();
        }
        else
        {
            std::cout << "File stations.txt have to contain more ground stations!" << std::endl;
        }
    }
    return 0;
}
