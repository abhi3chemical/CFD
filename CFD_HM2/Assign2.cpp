/*
 * Assign2.cpp
 *
 *  Created on: Feb 6, 2014
 *      Author: Abhishek Bagusetty
 */

// Boost includes
#include <boost/date_time/posix_time/posix_time.hpp>

#include <iostream>
#include <sstream>
#include <stdlib>
#include <stdexcept>
#include <cmath>
using std::cout;
using std::endl;

  int i,j;
  float alpha = 1.1234e-4;   // const thermal diffusivity
  float Lx = 0.3, Ly = 0.4;  // Length of Grid in x & y direction
  int Nrows = 81, Ncols = 61;      // No of Grid Points in x & y direction
  float delx = Lx/Ncols;
  float dely = Ly/Nrows;
  float delt = 0.1; // seconds

  int main( int argc, char* argv[] ){

  float T[Nrows][Ncols];

  // Initial Condition
  for( j = 0; j<Nrows; j++){
    for( i = 0; i < Ncols; i++ ){
      T[j][i] = 0.0;
    }
  }

  // For boundary conditions
  for( i = 0; i < Ncols; i++ ){
    T[0][i]       = 10.0; // Top
    T[Nrows-1][i] = 40.0; // bottom
  }
  for( j = 1; j < Nrows-1; j++ ){
    T[j][0]       = 0.0;  // Left
    T[j][Ncols-1] = 0.0;  // Right
  }

#ifdef NDEBUG
  for( j = 0; j<Nrows; j++){
    for( i = 0; i < Ncols; i++ ){
      std::cout << T[j][i] << " ";
    }
    cout << endl;
  }
#endif

  float rhs[Nrows-1], rhs1[Ncols-1];
  float a[Nrows-1], a1[Ncols-1];
  float b[Nrows-1], b1[Ncols-1];
  float c[Nrows-1], c1[Ncols-1];
  float output[Nrows-1], output1[Nrows-1];
  float t;
  float TV = 999999.99; // Tolerance Check variable
  float err;

  for( j = 0; j < Nrows; j++){
    rhs[j] = 0.0;
    a[j]   = 0.0;
    b[j]   = 0.0;
    c[j]   = 0.0;
  }
  for( i = 0; i < Ncols; i++){
    rhs1[i] = 0.0;
    a1[i]   = 0.0;
    b1[i]   = 0.0;
    c1[i]   = 0.0;
  }

  float ky = (alpha*delt)/(2*dely*dely);
  float kx = (alpha*delt)/(2*delx*delx);

for( t = 0; t < 11; t+= delt ){

//  //Tolerance Check
//  if( TV > 0.0001 ){
//    err = 0;
//    for( i = 0; i < Nrows; i++ ){
//      for(j = 0; j < Ncols; j++ ){
//        err = err + std::abs( Tnew[i][j] - Told[i][j] );
//      }
//    }
//    TV = (1 /(Nrows*Ncols))*err;
//  }
//  else{
//    std::cout << "Steady State ... \n";
//    std::cout << "Current TimeStep : " << t << std::endl;
//    return 0;
//  }

  // ---------------------- Y sweep -------------------------------//

  for( i = 0; i < Nrows; i++ ){   // (rows)
    for( j = 0; j < Ncols; j++ ){ // (columns)

      // Framing the RHS vector
      if( i == 0 ){
        rhs[j] = (1-(2*kx))*T[i][j] + (kx)*T[i+1][j];
      }
      else if( i == Nrows-1 ){
        rhs[j] = (kx)*T[i-1][j] + (1-(2*kx))*T[i][j];
      }
      else{
        rhs[j] = (kx)*T[i-1][j] + (1-(2*kx))*T[i][j] + (kx)*T[i+1][j];
      }

      // Framing the LHS Matrix
      b[j] = (1+(2*ky));  // main  diagonal

      if( j == 0 ){
        a[j] = 0;           // sub   diagonal
        c[j] = -ky;         // super diagonal
      }
      else if( j == Ncols -1 ){
        a[j] = -ky;         // sub    diagonal
        c[j] = 0;           // super  diagonal
      }
      else{
        a[j] = -ky;         // sub   diagonal
        c[j] = -ky;         // super diagonal
      }
    }

    float* cprime = (float*)malloc( sizeof(float) * Ncols );

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

    for( int k=1; k < Ncols; k++ ){
      float m  = 1.0/(b[k] - a[k] * cprime[k-1]);
      cprime[k] = c[k] * m;
      rhs[k] = (rhs[k] - a[k] * rhs[k-1]) * m;
    }

    for( int k = Ncols - 1; k-- > 0; ){
      output[k] = rhs[k] - cprime[k] * rhs[k+1];
    }

    free(cprime);

    for( int k=0; k < Ncols; k++ ){
      T[i][k] = output[k];
    }
  } // i loop

  // ---------------------- X sweep -------------------------------//

  for( j = 0; j < Ncols; j++ ){   // x sweep (columns)
    for( i = 0; i < Nrows; i++ ){ // y sweep (rows)

      // Framing the RHS vector
      if( j == 0 ){
        rhs1[i] = (1-(2*ky))*T[i][j] + (ky)*T[i][j+1];
      }
      else if( j == Ncols - 1 ){
        rhs1[i] = (ky)*T[i][j-1] + (1-(2*ky))*T[i][j];
      }
      else{
        rhs1[i] = (ky)*T[i][j-1] + (1-(2*ky))*T[i][j] + (ky)*T[i][j+1];
      }

      // Framing the LHS Matrix
      b1[i] = (1+(2*kx));  // main  diagonal
      if( i == 0 ){
        a1[i] = 0;         // sub diagonal
        c1[i] = -kx;       // super diagonal
      }
      else if( i == Nrows - 1 ){
        a1[i] = -kx;       // sub diagonal
        c1[i] = 0;         // super diagonal
      }
      else{
        a1[i] = -kx;         // sub   diagonal
        c1[i] = -kx;         // super diagonal
      }
    }

    float* cprime = (float*)malloc( sizeof(float) * Nrows );

    if( !cprime ){
      std::ostringstream msg;
      msg << "Error occurred in THOMAS ALGORITHM, at " << __FILE__ << " : " << __LINE__
          << std::endl;
      throw( std::runtime_error(msg.str()));
    }

    // Note : Make sure we are not dividing by 0
    if( b1[0] != 0 ){
      cprime[0] = c1[0]/b1[0];
      rhs1[0]    = rhs1[0]/b1[0];
    }else{
      std::ostringstream msg;
      msg << "Error occurred in THOMAS ALGORITHM, at " << __FILE__ << " : " << __LINE__
          << std::endl;
      throw( std::runtime_error(msg.str()));
    }

    for( int k = 1; k < Nrows; k++ ){
      float m  = 1.0/(b1[k] - a1[k] * cprime[k-1]);
      cprime[k] = c1[k] * m;
      rhs1[k] = (rhs1[k] - a1[k] * rhs1[k-1]) * m;
    }

    for( int k = Nrows - 1; k-- > 0; ){
      output1[k] = rhs1[k] - cprime[k] * rhs1[k+1];
    }

    free(cprime);

    for( int k=0; k < Nrows; k++ ){
      T[k][j] = output1[k];
    }
  }

  // For boundary conditions
  for( i = 0; i < Ncols; i++ ){
    T[0][i]       = 10.0; // Top
    T[Nrows-1][i] = 40.0; // bottom
  }
  for( j = 1; j < Nrows-1; j++ ){
    T[j][0]       = 0.0;  // Left
    T[j][Ncols-1] = 0.0;  // Right
  }

} // time loop

  // Print the Temperature Matrix
  for( j = 0; j<Nrows; j++){
    for( i = 0; i < Ncols; i++ ){
      std::cout << T[j][i] << " ";
    }
    cout << endl;
  }

} // main
//}
