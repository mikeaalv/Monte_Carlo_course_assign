//this program will simulate percolation on LXL sqaure with proportion p.
// Compile:
// module load icc/2018.1.163-GCC-6.4.0-2.28
// icc -O2 percolation.cpp -o percolation
// run:
// ./percolation L p randseed matflag sizeflag
// Argument:
//          L length of square side
//          p proportion of points
//          randseed initial random seed
//          matflag: if true output the random generated square array
//          sizeflag: if true output the cluster size vector
//it will return:
//               whether have infinite cluster,
//               fraction of points belong to the infinite cluster
//               the random generated square (if matflag=1),
//               vector of cluster size (if sizeflag=1)
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
int main(int argc, char* argv[])
{
  long Lsize=stol(argv[1]);
  double prop=stod(argv[2]);
  long rseed=stol(argv[3]);
  int matflag=stoi(argv[4]);
  int sizeflag=stoi(argv[5]);
  int percolmat[Lsize][Lsize];
  long sizeL=Lsize*Lsize;
  long nexit=long(prop*(double)sizeL);
  // printf("%d %f %d\n",sizeL,prop,nexit);
  for(int i=0;i<Lsize;i++){
    for(int j=0;j<Lsize;j++){
      percolmat[i][j]=0;
    }
  }
  //generate randome samples
  ///Naive random number implementaion
  // // srand(rseed);
  // std::default_random_engine generator(rseed);
  // std::uniform_int_distribution<int> distribution(0,sizeL-1);
  // for(long i=0;i<nexit;i++){
  //   long number=distribution(generator);
  //   // long number=rand()%(sizeL-1);
  //   printf("%d\n",number);
  //   if( (&percolmat[0][0])[number]==0 ){
  //     (&percolmat[0][0])[number]=1;
  //   }else{
  //     // it turns out pretty easy to repeat random number(short period)
  //     printf("problem on random number generator (not enough period)\n");
  //     exit(1);
  //   }
  // }
  
  /// random_shuffle implementaion
  /// reshuffle index vector randomly and use the first nexit index(uniform)
  std::vector<int> locindvec(sizeL);
  for(int i=0;i<sizeL;i++){
    locindvec[i]=i;
  }
  std::mt19937 g(rseed);
  std::shuffle(locindvec.begin(),locindvec.end(),g);
  for(long i=0;i<nexit;i++){
    long number=locindvec[i];
    // printf("%d %d\n",i,number);
    if( (&percolmat[0][0])[number]==0 ){
      (&percolmat[0][0])[number]=1;
    }else{
      printf("problem on random number generator (not enough period)\n");
      exit(1);
    }
  }
  
  //count clusters by H-K algorithm
  int clustlable[Lsize][Lsize];//first layer labels
  long clustlabellabel[sizeL];//label of label{also do count}
  int infilabel[sizeL][4];//table for labeling inifinite clusters(4 bounds)
  int labellargest=1;
  for(int i=0;i<Lsize;i++){
    for(int j=0;j<Lsize;j++){
      clustlable[i][j]=0;
    }
  }
  for(int i=0;i<sizeL;i++){
    clustlabellabel[i]=0;
    for(int j=0;j<4;j++){
      infilabel[i][j]=0;
    }
  }
  for(int i=0;i<Lsize;i++){
    for(int j=0;j<Lsize;j++){
      if(percolmat[i][j]==1){
        int left=0;
        int up=0;
        int label;
        //bound index condtion
        if(i!=0){
          up=percolmat[i-1][j];
        }
        if(j!=0){
          left=percolmat[i][j-1];
        }
        if(up==0&&left==0){
          label=labellargest;
          clustlable[i][j]=label;
          clustlabellabel[label]=1;
          labellargest++;
        }else if(up==0&&left==1){
          label=clustlable[i][j-1];
          //fetch real cluster
          if(clustlabellabel[label]<0){
            label= -clustlabellabel[label];
          }
          clustlable[i][j]=label;
          clustlabellabel[label]++;
        }else if(up==1&&left==0){
          label=clustlable[i-1][j];
          //fetch real cluster
          if(clustlabellabel[label]<0){
            label= -clustlabellabel[label];
          }
          clustlable[i][j]=label;
          clustlabellabel[label]++;
        }else if(up==1&&left==1){//deal with disagreement between left and up sides
          int label1=clustlable[i][j-1];
          int label2=clustlable[i-1][j];
          //fetch real cluster
          if(clustlabellabel[label1]<0){
            label1= -clustlabellabel[label1];
          }
          if(clustlabellabel[label2]<0){
            label2= -clustlabellabel[label2];
          }
          if(label1!=label2){
            //small value as new label
            if(label1>label2){
              int temp=label1;
              label1=label2;
              label2=temp;
            }
            clustlabellabel[label1]=clustlabellabel[label1]+clustlabellabel[label2]+1;
            clustlabellabel[label2]= -label1;
            //cluster combine
            for(int k=1;k<labellargest;k++){
              if(clustlabellabel[k]== -label2){
                clustlabellabel[k]= -label1;
              }
            }
          }else{
            clustlabellabel[label1]++;
          }
          // if(label1!=label2&&clustlabellabel[label2]<0&&clustlabellabel[label2]!= -label1){
          //   printf("%d %d\n",i,j);
          // }
          // printf("%d ",clustlabellabel[label1]);
          label=label1;
          clustlable[i][j]=label;
          //combine infinite information
          for(int j=0;j<4;j++){
            infilabel[label][j]=(infilabel[label][j]+infilabel[label2][j])>0;
          }
        }
        //connect to the four bounds
        if(i==0){
          infilabel[label][0]=((infilabel[label][0]+1)>0);
        }
        if(i==Lsize-1){
          infilabel[label][1]=((infilabel[label][1]+1)>0);
        }
        if(j==0){
          infilabel[label][2]=((infilabel[label][2]+1)>0);
        }
        if(j==Lsize-1){
          infilabel[label][3]=((infilabel[label][3]+1)>0);
        }
      }
    }
  }
  // for(long i=0;i<Lsize;i++){
  //   for(long j=0;j<Lsize-1;j++){
  //     printf("%d ",clustlable[i][j]);
  //   }
  //   printf("%d\n",clustlable[i][Lsize-1]);
  // }
  // printf("%d",labellargest);
  int exitinfi=0;
  long infini_ele_count=0;
  for(int i=1;i<labellargest;i++){
    if(clustlabellabel[i]>0){
      if(infilabel[i][0]>0&&infilabel[i][1]>0&&infilabel[i][2]>0&&infilabel[i][3]>0){//touch four bounds and so infinite
        // printf("%d",clustlabellabel[i]);
        exitinfi=1;
        infini_ele_count=clustlabellabel[i];
      }
    }
  }
  printf("%d %f %d %f\n",Lsize,prop,exitinfi,((double)infini_ele_count)/((double)sizeL));
  if(matflag==1){
    for(int i=0;i<Lsize;i++){
      for(int j=0;j<Lsize-1;j++){
        printf("%d ",percolmat[i][j]);
      }
      printf("%d\n",percolmat[i][Lsize-1]);
    }
  }
  if(sizeflag==1){
    for(int i=0;i<sizeL;i++){
      if(clustlabellabel[i]>0){
        printf("%d ",clustlabellabel[i]);
      }
    }
    printf("\n");
  }
}
