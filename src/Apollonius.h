#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>

#include <cmath> // square root

#include "Def.h"
#include "GaussianMatrix.h"

std::vector< Signal > mSignals;



/** solving Apollonius problem for set of stations
 *  @aStations a vector which contains stations
 */
PositionsList solveApol( int aSatId, long double aTimestamp, Stations aStations ) 
{
    std::cout << std::endl << "Solving Apollonius problem for positions: " << std::endl;
/*    std::vector< Station >::iterator iter;
    for( iter = aStations.begin(); iter != aStations.end(); ++iter )
    {
        std::cout << "      " << (*iter).getX() << ", " << (*iter).getY() << ", " << (*iter).getZ() << std::endl;
    }
*/
    std::vector< std::tuple< int, int, int, int > > Si;
    Si.push_back( std::make_tuple( 1, 1, 1, 1 ) );
    Si.push_back( std::make_tuple(-1,-1,-1,-1 ) );


PositionsList tempcalculatedPositions;
std::vector< std::tuple< int, int, int, int > >::iterator it;
it = Si.begin();
for( it = Si.begin(); it != Si.end(); ++it )
{
    int s[4];
    for( int i=0; i<4; ++i )
    { 
        s[0]=std::get<0>(*it);
        s[1]=std::get<1>(*it);
        s[2]=std::get<2>(*it);
        s[3]=std::get<3>(*it);
    }

    PositionsList calculatedPositions;
    std::vector< std::vector< double > >  matrix;
    double M,N,P,Q,R,S,a,b,c,rs,xs,ys,zs;

    double x1 = aStations.getStation(0).getX();
    double y1 = aStations.getStation(0).getY();
    double z1 = aStations.getStation(0).getZ();
    double r1 = aStations.getStation(0).getR(); // XXX to all times

    double el1, el2, el3, el4;
    el2 = (aStations.getStation(0).getX())*(aStations.getStation(0).getX()) + (aStations.getStation(0).getY())*(aStations.getStation(0).getY()) + (aStations.getStation(0).getZ())*(aStations.getStation(0).getZ());
	el3 = (aStations.getStation(0).getR())*(aStations.getStation(0).getR());

    std::cout << "Macierz: " << std::endl;
    for( int i=1; i<aStations.size() ; ++i )
    {
	    el1 = (aStations.getStation(i).getX())*(aStations.getStation(i).getX()) + (aStations.getStation(i).getY())*(aStations.getStation(i).getY()) + (aStations.getStation(i).getZ())*(aStations.getStation(i).getZ());
	    el4 = (aStations.getStation(i).getR())*(aStations.getStation(i).getR());

    	Row row = { ((aStations.getStation(i).getX())-x1)*2,
		     ((aStations.getStation(i).getY())-y1)*2,
		     ((aStations.getStation(i).getZ())-z1)*2,
		     (s[0]*r1-s[i]*aStations.getStation(i).getR())*2,
		     el1-el2+el3-el4 };
        matrix.push_back( row );
        std::cout << row.at(0) << " " << row.at(1) << " " << row.at(2) << " " << row.at(3) << std::endl; 
    }

    if( aStations.size() > 4 )
    {
            GaussianMatrix tempMatrix( matrix );
            /*GaussianMatrix tempMatrix;
            for( int i=0; i<aStations.size(); ++i )
            {
                tempMatrix.addRow( (aStations.getStation(i)).stationToVector() );
            }*/
            tempMatrix.printData();
            std::cout << "overdetermination: " << std::endl;
            tempMatrix.overdetermined();
/*            aStations.clear();
            for( int i=0; i<tempMatrix.getRowsNb(); ++i )
            {
                aStations.addStation( Station(tempMatrix.getRow(i)));
            }*/
            std::cout << "after overdetermination: " << std::endl;
            tempMatrix.printData();
            tempMatrix.makeGaussian();
            tempMatrix.printData();
            double xs = tempMatrix.get(0,4);
            double ys = tempMatrix.get(1,4);
            double zs = tempMatrix.get(2,4);
            double rs = tempMatrix.get(3,4);
            std::vector< double > solution = { xs, ys, zs, rs };
            calculatedPositions.addPosition( solution );
            std::cout << /*"satId=" << aSatId <<*/ "timestamp=" << aTimestamp << " " << xs << " " << ys << " " << std::setprecision(20) <<  zs/* << " " << rs*/ << std::endl << std::endl;
    }
    else
    {
        GaussianMatrix gaussMatrix( matrix );
        std::cout << "gaussian: " << std::endl;
        gaussMatrix.makeGaussian();
        M = gaussMatrix(1,5);
        N = -(gaussMatrix(1,4));
        P = gaussMatrix(2,5);
        Q = -(gaussMatrix(2,4));
        R = gaussMatrix(3,5);
        S = -(gaussMatrix(3,4));

        std::cout << "ignored calculations " << aTimestamp << ": " << gaussMatrix(1,1) << " " << gaussMatrix(2,2) << " " << gaussMatrix(3,3) << std::endl;
        a = N*N+Q*Q+S*S-1;
        b = 2*(M-x1)*N+2*(P-y1)*Q+2*(R-z1)*S-2*r1;
        c = (M-x1)*(M-x1)+(P-y1)*(P-y1)+(R-z1)*(R-z1)-r1*r1;

        long double p1 = (-b-sqrt(b*b-4*a*c))/(2*a);
        long double p2 = (-b+sqrt(b*b-4*a*c))/(2*a);

        if( b*b < 4*a*c )
        {
            p1 = (-b)/(2*a);
            p2 = (-b)/(2*a);

        }

        int rEarth = 6371000;
        rs = p1 ;
        xs = M+N*rs;
        ys = P+Q*rs;
        zs = R+S*rs;
    
        if( sqrt( xs*xs+ys*ys+zs*zs ) > rEarth )
        {
            std::vector< double > solution = { xs, ys, zs, rs };
            calculatedPositions.addPosition( solution );

            std::cout << /*"satId=" << aSatId <<*/ "timestamp=" << aTimestamp << " " << xs << " " << ys << " " << zs/* << " " << rs*/ << std::endl << std::endl;
        }

        rs = p2 ;
        xs = M+N*rs;
        ys = P+Q*rs;
        zs = R+S*rs;

        if( sqrt( xs*xs+ys*ys+zs*zs ) > rEarth )
        {
            std::vector< double > solution2 = { xs, ys, zs, rs };
            calculatedPositions.addPosition( solution2 );
    
            std::cout << /*"satId=" << aSatId <<*/ "timestamp=" << aTimestamp << " " << xs << " " << ys << " " << zs/* << " " << rs*/ << std::endl << std::endl;
        }
        tempcalculatedPositions=calculatedPositions;
    } 
}

    return tempcalculatedPositions;
}

