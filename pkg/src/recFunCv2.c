//This script includes helpfunctions used combined with R
//Compile code: R CMD SHLIB recProd.c
#include <stdio.h>
#include <stdlib.h>
//#define DEBUG

//function calculating cumulative values in x
void cumsum(int* x, int *sz, int *y) {
 int i;
 y[0] = 0;
 for(i=0; i<*sz; i++) {
  y[i+1] = y[i] + x[i];
 }
}

void sum(double* x, int sz, double *y) {
 int i;
 *y=0;
 for(i=0; i<sz; i++) {
  (*y) = (*y) + x[i];
 }
}

//Recursion function used
double recfun(int marker, double score, double pvalue, double *LRsuspect, double *LRlist, double *Plist, int *M, int *G, int *cG, double *inflationfactor) {
  int j;
  double tmp;
  double pvals[G[marker]]; 
  for(j=0; j<G[marker];j++) { pvals[j]=0; }
  if( marker == (*M-1) ) { 
    int j;
    for(j=0; j<G[marker]; j++) { 
      if( (float) (score * LRlist[cG[marker]+j]) >= (float) *LRsuspect ){
        pvals[j] = pvalue * Plist[cG[marker]+j];
     } else { break; }
    }
    sum(pvals,G[marker],&tmp);
//    fprintf(stderr,"tmp=%lf\n",*tmp);
    return(tmp);
  } else { 
    double scoreNew;
    double pvalueNew;
    for(j=0; j<G[marker]; j++) { 
      scoreNew = score * LRlist[cG[marker]+j];
      pvalueNew = pvalue * Plist[cG[marker]+j];
      if( (float) (scoreNew * inflationfactor[marker+1] ) >= (float) *LRsuspect ) { 
      pvals[j] = recfun(marker+1, scoreNew, pvalueNew, LRsuspect, LRlist, Plist, M, G, cG, inflationfactor);
//     fprintf(stderr,"marker=%d and geno=%d\n",marker,j);
      } else {  break;  }
    }
    sum(pvals,G[marker],&tmp);
    return(tmp);
  }
}




//X,P must be given at logscale
//wrapper function for recursion methods:
void prodrecfunction(double *pvalue, double *LRsuspect, double *LRlist, double *Plist, int *M, int *G,double *inflationfactor) {

 int *cG = malloc(sizeof(int)*((*M)+1));
 cumsum(G, M, cG);

#ifdef DEBUG 
 printf("pval=%lf\n",*pvalue);
 printf("LRsuspect=%lf\n",*LRsuspect);
 printf("#markers=%d\n",*M);
 int i,j;
 for(i=0; i<*M; i++) {
  printf("nA[%i]=%d\n",i,G[i]);
  printf("CnA[%i]=%d\n",i,cG[i]);
  printf("Inflator[%i]=%lf\n",i,inflationfactor[i]);
  for(j=0; j<G[i]; j++) { printf("X[%d,%d]=%lf\n",i,j,LRlist[cG[i]+j]); }
  for(j=0; j<G[i]; j++) { printf("P[%d,%d]=%lf\n",i,j,Plist[cG[i]+j]); }
  printf("\n"); 
 }
#endif

 int marker = 0;
 double score = 1.0;
 double pval = 1;
 *pvalue = recfun(marker, score, pval, LRsuspect, LRlist, Plist, M, G, cG, inflationfactor);
 free(cG);
}



