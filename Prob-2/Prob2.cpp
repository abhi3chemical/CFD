/*
 * Prob2.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: abhishek
 */
/*
 * Assign2-Prob2.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: Abhishek Bagusetty
 */

////////////////////////////////////////
//                                    //
//       FTCS SCHEME       //
//                                    //
////////////////////////////////////////

// Boost includes
//#include <boost/program_options.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include <cmath>

//#define ENABLE_DEBUG

using std::cout;
using std::endl;

float d;                 // thermal diffusivity
float mu = 0.000217;     // kinematic viscosity
float delY = 5e-4;       // delta Y
float h = 0.04;          // gap between the plates
float delt = 0.0001;     // time step
float t = 0;             // current time
float tf = 1.0;          // final time
int ny, nt=10000;              // Spatial & temporal points

boost::posix_time::ptime start, stop;
boost::posix_time::time_duration elapsed;

int main( int argc, char* argv[] ) {

  start = boost::posix_time::microsec_clock::universal_time();
  delt = (tf - t)/nt;
  ny = h/delY;
  nt = tf/delt;

//  {
//    boost::program_options::options_description desc( "Supported Options");
//    desc.add_options()
//        ("help", "print help message )")
//        ("ntime", boost::program_options::value<int>(&nt)->default_value(10000), " Number of time steps " );
//
//    boost::program_options::variables_map args;
//    boost::program_options::store( boost::program_options::parse_command_line( argc, argv, desc), args );
//    boost::program_options::notify( args );
//
//    if( args.count("help")) {
//      std::cout << desc << std::endl;
//      return 1;
//    }
//  }

  // velocity Array
  float u[ny], uold[ny];

  d = (delt*mu)/(delY*delY);
  // Stability criteria
  if( d < 0.5 ) std::cout << "Solution is Expected to be STABLE... \n";
  else          {
    std::ostringstream msg;
    msg << "Error solution is UNSTABLE, at " << __FILE__ << " : " << __LINE__
        << std::endl;
    throw( std::runtime_error(msg.str()));
  }

  // initial & boundary conditions
  for( int j = 0; j < ny; j++ ){
    if( j == 0 ){
      uold[0] = 40.0; // Lower Boundary Condition
      u[0] = 40.0;
    }else{
      uold[j] = 0.0;  // Upper  Boundary Condition
      u[j] = 0.0;
    }
  }

  for( int k = 0; k < nt; k++ ){ // time loop

    u[0] = 40.0;
    u[ny-1] = 0.0;

#ifdef ENABLE_DEBUG
    std::cout << "Time Step : " << k << std::endl;
#endif

    for( int j=1; j < ny-1; j++ ){ // Interior Points
        u[j] = d*uold[j-1] + (1-2*d)*uold[j] + d*uold[j+1];
    }

    if( k != 0 ){
      // Steady State Condition
      float TV = 0, tmp = 0;
      for( int j = 1; j < ny-1; j++){
        tmp += std::abs( u[j] - uold[j] );
      }
      TV = tmp/float(ny);
      if( TV < 0.0001 ) {
        std::cout << "Steady State Attained....\n";
        std::cout << "Time Step (" << k << ") : " << t << std::endl;
        stop = boost::posix_time::microsec_clock::universal_time();
        elapsed = stop - start;
        std::cout << "Time Taken : " << (elapsed.total_microseconds()*1e-6) <<
            " sec " << std::endl;
        return 0;
      }

      // Preparing for next time step
      t = t + delt;                        // current time
      for( int l=1; l < ny-1; l++){
        uold[l] = u[l];
      }
    }


  } // time loop

  stop = boost::posix_time::microsec_clock::universal_time();
  elapsed = stop - start;
  std::cout << "Time Taken : " << (elapsed.total_microseconds()*1e-6) <<
      " sec " << std::endl;

#ifdef ENABLE_DEBUG
  std::cout << "Final Temperatures \n";
  for( int l = 0; l < ny; l++ ){
    std::cout << "l : " << l << " , " << u[l] << std::endl;
  }
#endif

}




