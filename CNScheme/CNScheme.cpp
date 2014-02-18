/*
 * CNScheme.cpp
 *
 *  Created on: Feb 17, 2014
 *      Author: abhishek
 */

////////////////////////////////////////
//                                    //
//       CRANK-NICHOLSON SCHEME       //
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
using std::cout;
using std::endl;

// Debug flag
//#define ENABLE_DEBUG

float d;                 // thermal diffusivity
float mu = 0.000217;     // kinematic viscosity
float delY = 5e-4;       // delta Y
float h = 0.04;          // gap between the plates
float delt = 0.0001;     // time step
float t = 0;             // current time
float tf = 1.0;          // final time
int ny, nt=10000;        // Spatial & temporal points
float alpha;

boost::posix_time::ptime start, stop;
boost::posix_time::time_duration elapsed;

int main( int argc, char* argv[] ) {

  start = boost::posix_time::microsec_clock::universal_time();
  delt = (tf - t)/nt;
  ny = 10;//h/delY;
  nt = 10;//tf/delt;

  // Diagonal vectors
  float a[ny], b[ny], c[ny], rhs[ny];

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
    float unew[ny], uold[ny];

    // Initializing data structures
    for( int k = 0; k < ny; k++ ){
      unew[k] = 0.0;
      uold[k] = 0.0;
      a[k]    = 0.0;
      b[k]    = 0.0;
      c[k]    = 0.0;
      rhs[k]  = 0.0;
    }

    alpha = (delt*mu)/(2*delY*delY);

    for( int k = 0; k < nt; k++ ){ // time loop

#ifdef ENABLE_DEBUG
      std::cout << "time step and time : " << k << "," << t << std::endl;
#endif

      t += delt;                      // current time

      for( int j = 0; j < ny; j++ ){ // spatial looping

        // Framing the RHS vector
        if( j == 0 ){
          rhs[j] = (1+(2*alpha))*uold[j] + (alpha)*uold[j+1];
        }
        else if( j == ny-1 ){
          rhs[j] = (alpha)*uold[j-1] + (1+(2*alpha))*uold[j];
        }
        else{
          rhs[j] = (alpha)*uold[j-1] + (1+(2*alpha))*uold[j] + (alpha)*uold[j+1];
        }

        // Framing the LHS Matrix
        b[j] = (1+(2*alpha));    // main  diagonal

        if( j == 0 ){
          a[j] = 0;              // sub   diagonal
          c[j] = -alpha;         // super diagonal
        }
        else if( j == ny-1 ){
          a[j] = -alpha;         // sub    diagonal
          c[j] = 0;              // super  diagonal
        }
        else{
          a[j] = -alpha;         // sub   diagonal
          c[j] = -alpha;         // super diagonal
        }

      } // spatial loop

    // TDMA / Thomas Algorithm at each time step
    float* cprime = (float*)malloc( sizeof(float) * ny );

    if( !cprime ){
      std::ostringstream msg;
      msg << "Error occurred in THOMAS ALGORITHM, at " << __FILE__ << " : " << __LINE__
          << std::endl;
      throw( std::runtime_error(msg.str()));
    }

    // Note : Make sure we are not dividing by 0
    if( b[0] != 0 ){
      cprime[0] = c[0]/b[0];
      rhs[0]    = rhs[0]/b[0];
    }else{
      std::ostringstream msg;
      msg << "Error occurred in THOMAS ALGORITHM, at " << __FILE__ << " : " << __LINE__
          << std::endl;
      throw( std::runtime_error(msg.str()));
    }

    for( int k=1; k < ny; k++ ){
      float m  = 1.0/(b[k] - a[k] * cprime[k-1]);
      cprime[k] = c[k] * m;
      rhs[k] = (rhs[k] - a[k] * rhs[k-1]) * m;
    }

    // output vector at new time step
    for( int k = ny - 1; k-- > 1; ){
      unew[k] = rhs[k] - cprime[k] * rhs[k+1];
    }

    free(cprime);

    // TODO : do the steady state check here
    for( int k = 0; k < ny; k++ ){
      uold[k] = unew[k];
    }

    // Apply Boundary Conditions
    unew[0]    = 40.0; uold[0] = 40.0;   // Lower  B.C
    unew[ny-1] =  0.0; uold[ny-1] = 0.0; // Upper  B.C

#ifdef ENABLE_DEBUG
    // Debug : Print the matrices
    for( int k = 0; k < ny; k++ ){
      std::cout << "unew[" << k << "]" << unew[k] << std::endl;
    }
#endif

    } // time loop



  } // main

