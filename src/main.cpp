#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>
#include <unistd.h>

#include <cmath> // square root


#include "Def.h"
#include "Combinations.h"
#include "GaussianMatrix.h"

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
    return calculatedPositions;
}

int main( int argc, char* argv[] )
{
    Stations mStations;
    int c = getopt( argc, argv, "f:" );
    switch( c ){
	case 'f':
	    loadStations( optarg, mStations );
	    break;
        default:
	    std::cout << "Missing file name! Use parameter -f filename." << std::endl;
	    return -1;
    }
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
            xyzr.addPositions( solveApol( takenStations ) );
        }

        xyzr.printAveragePosition();
    }
    else
    {
        std::cout << "File stations.txt have to contain more ground stations!" << std::endl;
    }
    
    return 0;

}
