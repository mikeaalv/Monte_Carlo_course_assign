//this program will use simple random sampling to study LXL Ising model
// in the calculation, both k_b and J will be assumed with value 1 for convinence.
// temperature is assumed to be infinite
// Compile:
// module load icc/2018.1.163-GCC-6.4.0-2.28
// icc -O2 randomsampling.cpp -o randomsampling
// run:
// ./randomsampling L N randseeds inip gap
// Argument:
//        L: lattice size. The number of spine will be L^2
//        N: the length of Markove chain. L^2 spin flip is one Monte Caro step
//        randseed: random seed
//        inip: initial proportion of 1
//        gap: the number of skipping samples(if gap=1, no scape, collect samples after every Monte Carlo step)
//it will return:
//        a table with microstate H(first column) and M(second column) for each Monte Caro step
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <array>
#include <iterator>
using namespace std;
int delH_cal(int **lattice,int i,int j,int Lsize){
  int lasti=i-1;
  int lastj=j-1;
  int nexti=i+1;
  int nextj=j+1;
  lasti=lasti<0?(Lsize-1):lasti;
  lastj=lastj<0?(Lsize-1):lastj;
  nexti=nexti>=Lsize?0:nexti;
  nextj=nextj>=Lsize?0:nextj;
  int center=lattice[i][j];
  int delH= -(center*lattice[lasti][j]+center*lattice[i][lastj]+center*lattice[nexti][j]+center*lattice[i][nextj]);
  return delH;
}
int main(int argc, char* argv[])
{
  int Lsize=stoi(argv[1]);
  // double temperat=stod(argv[2]);
  long Nsamp=stol(argv[2]);
  long rseed=stol(argv[3]);
  double inip=stod(argv[4]);
  long gap=stod(argv[5]);
  int sizeL=Lsize*Lsize;
  int **lattice;
  lattice= new int *[Lsize];
  for(int i=0;i<Lsize;i++){
    lattice[i]= new int[Lsize];
    for(int j=0;j<Lsize;j++){
      lattice[i][j]= -1;
    }
  }
  //initial state
  /// random_shuffle implementaion
  /// reshuffle index vector randomly and use the first nexit index(uniform)
  long nexit=long(inip*(double)sizeL);//0.5 as the initial state of lattice with random 50% spin flipped
  std::vector<int> locindvec(sizeL);
  for(int i=0;i<sizeL;i++){
    locindvec[i]=i;
  }
  std::mt19937 g(rseed);
  std::shuffle(locindvec.begin(),locindvec.end(),g);
  for(long i=0;i<nexit;i++){
    long number=locindvec[i];
    int loci=number/Lsize;
    int locj=number%Lsize;
    if(lattice[loci][locj]== -1){
      lattice[loci][locj]=1;
    }else{
      printf("problem on random number generator (not enough period)\n");
      exit(1);
    }
  }
  //calculating initial H
  long H_micro=0;
  for(int i=0;i<Lsize;i++){
    for(int j=0;j<Lsize;j++){
      H_micro+=delH_cal(lattice,i,j,Lsize);
    }
  }
  H_micro=H_micro/2;// as all bounds are calculated twice
  long M_micro=0;
  for(int mi=0;mi<Lsize;mi++){
    for(int mj=0;mj<Lsize;mj++){
      M_micro+=lattice[mi][mj];
    }
  }
  srand(rseed);
  srand48(rseed);
  //initial value
  printf("%ld\t%ld\n",H_micro,M_micro);
  for(long i=0;i<Nsamp;i++){
    // int accepn=0;
    for(int gapi=0;gapi<gap;gapi++){
      for(int j=0;j<sizeL;j++){
        int loc=rand()%sizeL;
        int loci=loc/Lsize;
        int locj=loc%Lsize;
        long delH=0;
        delH= -2*delH_cal(lattice,loci,locj,Lsize);
        // double r=drand48();
        // // printf("prob:%f\t",exp(-delH/temperat));
        // if(r<exp(-delH/temperat)){
        H_micro=H_micro+delH;
        lattice[loci][locj]= -lattice[loci][locj];
        // accepn++;
        // }
      }
    }
    
    M_micro=0;
    for(int mi=0;mi<Lsize;mi++){
      for(int mj=0;mj<Lsize;mj++){
        M_micro+=lattice[mi][mj];
        // printf("%d,",lattice[mi][mj]);
      }
      // printf(";");
    }
    // printf("accep: %f ",(double)accepn/(double)sizeL);
    printf("%ld\t%ld\n",H_micro,M_micro);
  }
}
