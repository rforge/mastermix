//This function calculates the likelihood over all genotype combinations for each loci for given parameters
//R CMD SHLIB contBProbCpp.cpp -llapack -lblas
#include <cmath>
#include <RcppArmadillo.h> //require RcppArmadillopackage and Namespaced defined
using namespace std;
using namespace arma;
//#define DEBUG
//helpfunctions:
void cumsum0(int* x, int *sz, int *y);

/************************/
/****Helpfunctions*******/
/************************/

//function calculating cumulative values in x
//x: steps, sz: vectorlength, y: vector
void cumsum0(int* x, int *sz, int *y) {
 int i;
 y[0] = 0;
 for(i=0; i<*sz; i++) {
  y[i+1] = y[i] + x[i];
 }
}


extern "C" {

//function for calculating inner product of integrale of theta
void calcLthetaC(double *PE,int *nA, int *nC, int *nL, double *allY,  double *allX,  double *allO,  double *alliW,  int *cdfX, int *nQi, double *pG, double *theta) {

 int i,j,cc;
 int *CnA = new int[*nL+1]; //cumulative index for vector
 int *CnA2 = new int[*nL+1]; //cumulative index for matrix (vectorized matrix)
 int *nAsq = new int[*nL+1]; //squared numbers of nA
 int *CnQ = new int[*nL+1]; //cumulative index for matrix (vectorized matrix)

 Mat<double> *Yi; 
 Mat<double> *iWi; 
 Mat<double> *Xij; 
 Mat<double> *Oij; 
 Mat<double> *omega; //mixture proportion 
 omega = new Mat<double>( theta, *nC-1, 1 ,false); //insert mixture proportion

 //temporary variables:
  Mat<double> mui;
  Mat<double> ri;
  Mat<double> Di;
  double bigsum;
  double sigma;
  sigma = sqrt(theta[*nC-1]);

 for(i=0; i<*nL; i++) { 
  nAsq[i] = nA[i]*nA[i]; 
 } //square numbers
 cumsum0(nA, nL, CnA); //cumulative numbers of nA
 cumsum0(nAsq, nL, CnA2); //cumulative squared number of alleles nA
 cumsum0(nQi, nL, CnQ); //cumulative of numbers of genotypes. Used in pG

#ifdef DEBUG
 for(i=0; i<*nL; i++) {
  printf("nA[%i]=%d\n",i,nA[i]);
  printf("CnA[%i]=%d\n",i,CnA[i]);
  printf("CnA2[%i]=%d\n",i,CnA2[i]);
  printf("nQi[%i]=%d\n",i,nQi[i]);
  for(j=0; j<nA[i]; j++) { printf("Y[%d,%d]=%f\n",i,j,allY[CnA[i]+j]); }
  printf("\n"); 
 }
#endif

 //calculating probabilities on the move:
 cc = 0; //counter for X-matrix and O-matrix
 for(i=0; i<*nL; i++) {
  Yi =  new Mat<double>( &allY[CnA[i]], nA[i], 1,false); //init datavector
  iWi = new Mat<double>( &alliW[CnA2[i]],nA[i],nA[i],false); //init. covariance-matrix
  bigsum = 0.0; //sum over all genotypes

  for(j=0; j<nQi[i]; j++) { //for all genotypes: Init space and insert relevant matrice directly
   Xij = new Mat<double>( &allX[cdfX[cc]-1] ,nA[i],*nC-1,false); //init. Xij-matrix
   Oij = new Mat<double>( &allO[cdfX[cc]-1] ,nA[i],*nC-1,false); //init. Pi-matrices as blocked   
   mui = (*Xij)*(*omega) + *Oij; //mean
   ri = *Yi - mui; //residual
   Di = trans(ri)*(*iWi)*(ri); //mahalonobis distance
   bigsum = bigsum + exp(-0.5/theta[*nC-1]*(Di.at(0,0)))/pow(sigma,nA[i])*pG[CnQ[i]+j]; //#likelihood value   
   cc = cc + 1; //update counter
   delete Xij;
   delete Oij;
  }//end for each j: genotype combination
  *PE = (*PE)*bigsum;
  delete Yi;
  delete iWi;
 }//end for each loci i:

 delete[] CnA;
 delete[] CnA2;
 delete[] CnQ;
 delete[] nAsq;
 delete omega;
}

} 
