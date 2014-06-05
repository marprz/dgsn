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
PositionsList solveApol( int aSatId, long double aTimestamp, std::vector< Station > aStations ) 
{
    std::cout << std::endl << "Solving Apollonius problem for positions: " << std::endl;
    std::vector< Station >::iterator iter;
    for( iter = aStations.begin(); iter != aStations.end(); ++iter )
    {
        std::cout << "      " << (*iter).getX() << ", " << (*iter).getY() << ", " << (*iter).getZ() << std::endl;
    }

    std::vector< std::tuple< int, int, int, int > > Si;
    Si.push_back( std::make_tuple( 1, 1, 1, 1 ) );
/*    Si.push_back( std::make_tuple(-1, 1, 1, 1 ) );
    Si.push_back( std::make_tuple( 1,-1, 1, 1 ) );
    Si.push_back( std::make_tuple( 1, 1,-1, 1 ) );
    Si.push_back( std::make_tuple( 1, 1, 1,-1 ) );
    Si.push_back( std::make_tuple(-1,-1, 1, 1 ) );
    Si.push_back( std::make_tuple(-1, 1,-1, 1 ) );
    Si.push_back( std::make_tuple(-1, 1, 1,-1 ) );
    Si.push_back( std::make_tuple( 1,-1,-1, 1 ) );
    Si.push_back( std::make_tuple( 1,-1, 1,-1 ) );
    Si.push_back( std::make_tuple( 1, 1,-1,-1 ) );
    Si.push_back( std::make_tuple( 1,-1,-1,-1 ) );
    Si.push_back( std::make_tuple(-1, 1,-1,-1 ) );
    Si.push_back( std::make_tuple(-1,-1, 1,-1 ) );
    Si.push_back( std::make_tuple(-1,-1,-1, 1 ) );*/
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

    double x1 = aStations.at(0).getX();
    double y1 = aStations.at(0).getY();
    double z1 = aStations.at(0).getZ();
    double r1 = aStations.at(0).getR();

    aStations.at(1).addToZ( 1 );
    aStations.at(2).addToZ( 2 );
    aStations.at(3).addToZ( 3 );

//    std::cout << "solveApol(): " << x1 << ", " << y1 << ", " << z1 << ", " << r1 << std::endl;
    double el1, el2, el3, el4;
    el2 = (aStations.at(0).getX())*(aStations.at(0).getX()) + (aStations.at(0).getY())*(aStations.at(0).getY()) + (aStations.at(0).getZ())*(aStations.at(0).getZ());
	el3 = (aStations.at(0).getR())*(aStations.at(0).getR());

//    std::cout << "Macierz: " << std::endl;
    for( int i=1; i<=3; ++i )
    {
	    el1 = (aStations.at(i).getX())*(aStations.at(i).getX()) + (aStations.at(i).getY())*(aStations.at(i).getY()) + (aStations.at(i).getZ())*(aStations.at(i).getZ());
	    el4 = (aStations.at(i).getR())*(aStations.at(i).getR());

    	Row row = { ((aStations.at(i).getX())-x1)*2,
		     ((aStations.at(i).getY())-y1)*2,
		     ((aStations.at(i).getZ())-z1)*2,
		     (s[0]*r1-s[i]*aStations.at(i).getR())*2,
		     el1-el2+el3-el4 };
        matrix.push_back( row );
//        std::cout << row.at(0) << " " << row.at(1) << " " << row.at(2) << " " << row.at(3) << std::endl; 
    }
    GaussianMatrix gaussMatrix( matrix );
    M = gaussMatrix(1,5);
    N = -(gaussMatrix(1,4));
    P = gaussMatrix(2,5);
    Q = -(gaussMatrix(2,4));
    R = gaussMatrix(3,5);
    S = -(gaussMatrix(3,4));

//    std::cout << "M, N, P, Q, R, S: " << M << ", " << N << ", " << P << ", " << Q << ", " << R << ", " << S << std::endl;

    a = N*N+Q*Q+S*S-1;
    b = 2*(M-x1)*N+2*(P-y1)*Q+2*(R-z1)*S-2*r1;
    c = (M-x1)*(M-x1)+(P-y1)*(P-y1)+(R-z1)*(R-z1)-r1*r1;

//    std::cout << "a, b, c: " << a << ", " << b << ", " << c << std::endl;
    long double p1 = (-b-sqrt(b*b-4*a*c))/(2*a);
    long double p2 = (-b+sqrt(b*b-4*a*c))/(2*a);
/*
    std::cout << "b^2=" << b*b << ", 4ac=" << 4*a*c << ", sqrt() =" << sqrt(b*b-4*a*c) << std::endl;
    if( sqrt(b*b-4*a*c) > 0 ) 
        std::cout << std::endl << "DOBRZE!" << std::endl << std::endl;
    if( b*b == 4*a*c ) 
        std::cout << "b^2 = 4ac" << std::endl;
    else 
        std::cout << "b^2-4ac = " << b*b-4*a*c << std::endl;
*/
    if( b*b < 4*a*c )
    {
        p1 = (-b)/(2*a);
        p2 = (-b)/(2*a);

  //      std::cout << "b^2 < 4ac. p1=" << p1 << ", p2=" << p2 << std::endl;
    }

 //   std::cout << "p1, p2: " << p1 << ", " << p2 << std::endl;
    //rs = p1>p2 ? p1 : p2 ;
    rs = p1 ;
    xs = M+N*rs;
    ys = P+Q*rs;
    zs = R+S*rs;
    std::vector< double > solution = { xs, ys, zs, rs };
    calculatedPositions.addPosition( solution );

    std::cout << "satId=" << aSatId << ", timestamp=" << aTimestamp << ": " << xs << ", " << ys << ", " << zs << ", " << rs << std::endl << std::endl;
    
    
    rs = p2 ;
    xs = M+N*rs;
    ys = P+Q*rs;
    zs = R+S*rs;
    std::vector< double > solution2 = { xs, ys, zs, rs };
    calculatedPositions.addPosition( solution2 );

    std::cout << "satId=" << aSatId << ", timestamp=" << aTimestamp << ": " << xs << ", " << ys << ", " << zs << ", " << rs << std::endl << std::endl;

    tempcalculatedPositions=calculatedPositions;
} 
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
    return tempcalculatedPositions;
}

