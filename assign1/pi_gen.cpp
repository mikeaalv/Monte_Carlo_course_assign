// this program approximate pi through random drawing points in square
// pi is calculated by 4*Q, Q is the ratio of sample in the quarter circile within the square
// Compile:
// module load icc/2018.1.163-GCC-6.4.0-2.28
// icc -O2 pi_gen.cpp -o pi_gen
// Run
// ./pi_gen N_sample randseed
// N_sample is the number of drawn points
// randseed is the start random seed
// there will be 10 (N_INDEP_CHAIN) independent run each with length N_sample and with random seed randseed..randseed+9
// Return: the program will output 10 rows, and each row is an estimated pi
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
using namespace std;
#define N_INDEP_CHAIN 10
int main(int argc, char* argv[])
{
  long Nsample=stol(argv[1]);
  long randseed=stol(argv[2]);
  double pi_est[10];
  // sample for each sampling chain
  for(long i=0;i<N_INDEP_CHAIN;i++){
    //different random seed for each sampling chain
    long randseednew=randseed+i;
    srand48(randseednew);
    double rand_axis[2];
    long Ncircle=0;
    // sample within one sampling chain
    for(long j=0;j<Nsample;j++){
      rand_axis[0]=drand48();
      rand_axis[1]=drand48();
      // within circle
      if(pow(rand_axis[0],2)+pow(rand_axis[1],2) < 1){
        Ncircle++;
      }
    }
    pi_est[i]=4.0*(double)Ncircle/(double)Nsample;
    printf("%.7f\n",pi_est[i]);
  }
}
