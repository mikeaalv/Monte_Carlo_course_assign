//this program will use Wang-Landau sampling to study LXL Ising model
// in the calculation, both k_b and J will be assumed with value 1 for convinence.
// Compile:
// module load icc/2018.1.163-GCC-6.4.0-2.28
// icc -O2 wanglandau.cpp -o wanglandau
// run:
// ./wanglandau L hflat hf randseeds inip
// Argument:
//        L: lattice size. The number of spin will be L^2
//        hflat: the threshold for flat in histogram. The histogram is flat when min(H)>=hflat*mean(H)
//        fini: the initial value for f.
//        hf: the threshold for iterations. The iteration is finished when f<=hf
//        randseed: random seed
//        inip: initial proportion of 1
//it will return:
//        a density function table: energy(1st column), log density g(2nd column), last histogram (3rd column)
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
long locind(long *locvec,long locele,long nlen){
  for(long i=0;i<nlen;i++){
    if(locele==locvec[i]){
      return i;
    }
  }
  return -1;
}
double mean(long *vec,long len){
  double sum=0.0;
  for(long i=0;i<len;i++){
    sum+=(double)vec[i];
  }
  return sum/((double)len);
}
long min(long *vec,long len){
  long min=vec[0];
  for(long i=0;i<len;i++){
    if(min>vec[i]){
      min=vec[i];
    }
  }
  return min;
}
int main(int argc, char* argv[])
{
  //input parameters
  int Lsize=stoi(argv[1]);
  double hflat=stod(argv[2]);
  double fini=stod(argv[3]);
  double fmin=stod(argv[4]);
  long rseed=stol(argv[5]);
  double inip=stod(argv[6]);
  int sizeL=Lsize*Lsize;
  //initilize the lattice
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
  //initialize vectors according to energy space
  long nElevels=sizeL-1;
  long Evec[nElevels];
  double gveclog[nElevels];
  long startE= -sizeL*2;
  for(long i=0;i<nElevels;i++){
    Evec[i]=startE;
    // printf("%d ",startE);
    gveclog[i]=0.0;
    if(i==0||i==(nElevels-2)){
      startE=startE+8;
    }else{
      startE=startE+4;
    }
  }
  // printf("\n");
  //calculating initial E
  long E_micro=0;
  for(int i=0;i<Lsize;i++){
    for(int j=0;j<Lsize;j++){
      E_micro+=delH_cal(lattice,i,j,Lsize);
    }
  }
  E_micro=E_micro/2;
  //initilize f
  double flog=log(fini);//exp(1);
  //the density estimation loop
  srand(rseed);
  srand48(rseed);
  long Hvec[nElevels];
  do{
    for(long i=0;i<nElevels;i++){
      Hvec[i]=0;
    }
    int flatflag=0;
    long counti=0;
    do{
      counti++;
      int loc=rand()%sizeL;
      int loci=loc/Lsize;
      int locj=loc%Lsize;
      long delE=0;
      // printf("{%d}\n",loc);
      delE= -2*delH_cal(lattice,loci,locj,Lsize);
      long E_micro2=E_micro+delE;
      // printf("loc %d e1 %ld e2 %ld dele %ld ",loc,E_micro,E_micro2,delE);
      // printf("e1 %ld e2 %ld dele %ld\n",E_micro,E_micro2,delE);
      long gE1ind=locind(Evec,E_micro,nElevels);
      long gE2ind=locind(Evec,E_micro2,nElevels);
      double gE1=gveclog[gE1ind];
      double gE2=gveclog[gE2ind];
      double r=drand48();
      // if(1){
      //   printf("\n");
      //   for(int locherei=0;locherei<Lsize;locherei++){
      //     for(int locherej=0;locherej<Lsize;locherej++){
      //       printf("%d",lattice[locherei][locherej]);
      //     }
      //     printf("\n");
      //   }
      // }
      // printf("ind1 %d, ind2 %d g1: %lf g2: %lf ",gE1ind,gE2ind,gE1,gE2);
      //acceptance
      if(r<exp((double)gE1-(double)gE2)){
        // printf("**");
        E_micro=E_micro2;
        gE1ind=gE2ind;
        lattice[loci][locj]= -lattice[loci][locj];
      }
      // if(1){
      //   printf("\n");
      //   for(int locherei=0;locherei<Lsize;locherei++){
      //     for(int locherej=0;locherej<Lsize;locherej++){
      //       if(loci==locherei&&locherej==locj&&r>=(double)gE1/(double)gE2){
      //         printf("%d",-lattice[locherei][locherej]);
      //       }else{
      //         printf("%d",lattice[locherei][locherej]);
      //       }
      //     }
      //     printf("\n");
      //   }
      // }
      gveclog[gE1ind]=gveclog[gE1ind]+flog;
      Hvec[gE1ind]=Hvec[gE1ind]+1;
      if(counti%100==0){
        double ratio=((double)min(Hvec,nElevels))/mean(Hvec,nElevels);
        // long minH=min(Hvec,nElevels);
        // long minHind=locind(Hvec,minH,nElevels);
        // printf("ratio:%lf H:%ld minH %ld minHind %ld\n",ratio,Hvec[gE1ind],minH,minHind);
        if(ratio>hflat){
          flatflag=1;
        }
      }
    }while(flatflag!=1);
    flog=flog/2;
    // printf("%f \n",flog);
  }while(exp(flog)>fmin);
  for(long i=0;i<nElevels;i++){
    printf("%ld\t%lf\t%ld\n",Evec[i],gveclog[i],Hvec[i]);
  }
}
