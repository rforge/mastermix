//This function permutates and fits a lin. gauss. model (Tvedebrink 2011).
//cd ~/Dropbox/Forensic/MixtureProj/myDev/mastermix/R
//R CMD SHLIB contProbCpp.cpp -llapack -lblas
#include <cmath> //math expressions
#include <vector> //list-structure for loci
#include <armadillo> //
using namespace std;
using namespace arma;
const double PI = std::acos(0.0)*2;
//#define DEBUG
//helpfunctions:
void cumsum(int* x, int *sz, int *y);

/*********************************/
/*Class defining for Loci-objects*/
/*********************************/
class Loci {
 public: 
  /*Array-structure*/
  int locinr;
  int nC; //number of contributors
  int nA; //number alleles at loci
  int nCombs; //number of combined genotypes
  int combsel; //index of selected genotypecomb
  int* cumPind; //cumulative row-index for genotypecomb in Pmatrix
  double *PG; //probability of genotype

  Mat<double> *Y; //pointer to data vector
  Mat<double> *iW; //pointer to inverse weight-matrix
  Mat<double> *Pmatrix; //pointer to matrix used to pull out Pi
  Mat<double> *Pi; //points to matrix used in model 
  Mat<double> *PiO; //points to matrix used in model 
  Mat<double> *XTWXi; //used to store matrix temporary in recursion 
  Mat<double> *XTWYi; //used to store vector temporary in recursion 
  Mat<double> *YTWYi; //used to store vector temporary in recursion 
  Mat<double> *Plist ; //points to list of Pi-matrices created by genPlist
  double S_Y; //sum of Y-values


//constructor:
  Loci(int i, int nA, double* data,int nC, double* weight, int nCombs, double* subPlist, double *PG) {
  //i - locisnr
  //n - number of alleles, data - values of data
  //C - number of contributors
   int j;
   locinr = i; //insert loci-number
   combsel = 0; //first combination selected as default.
   Y =  new Mat<double>(data,nA,1,false); //init datavector
   iW = new Mat<double>(weight,nA,nA,false); //init. covariance-matrix
   Pmatrix = new Mat<double>(subPlist,nA*nCombs,nC,false); //init. Pi-matrices as blocked
   Pi = new Mat<double>(nA,nC); //init. Pi-matrix
   PiO = new Mat<double>(nA,nC-1); //init PiC-matrix at each column
   XTWXi = new Mat<double>(nC-1,nC-1); //init
   XTWYi = new Mat<double>(nC-1,1); //init
   YTWYi = new Mat<double>(1,1); //init
   cumPind = new int[nCombs];
   for(j=0; j<nCombs; j++) { cumPind[j] = j*nA; } //cumulative index for Pmatrix 
   S_Y = accu(*Y);
   //if(nA>2*nC) printf("Too few contributors!!");
   this->nC = nC;
   this->nA = nA;
   this->nCombs = nCombs;
   this->PG = PG;
#ifdef DEBUG 
   Y->print();
   iW->print();
   Pmatrix->print();
   for(j=0; j<nCombs; j++) { 
    printf("%i",cumPind[j]);
    printf("%f",PG[j]);
   } //cumulative index for Pmatrix 
#endif
  }
  ~Loci() {
   delete Y,iW,Pmatrix,Pi,PiO,XTWXi,XTWYi,YTWYi;
   delete[] cumPind;
  }

   //function for updating combmatrix Pi for current loci
  void updateCombMatrix() {
   int k;
   Pi->submat(0,0,nA-1,nC-1)  = Pmatrix->submat(cumPind[combsel],0,cumPind[combsel]+nA-1,nC-1); 
   for(k=0;k<(nC-1); k++) {
    PiO->submat(0,0,nA-1,k) = Pi->submat(0,nC-1,nA-1,nC-1); //insert last column
   }
   combsel++; //update combsel for next extraction
  }

   //function for updating a specific combmatrix Pi for current loci
  void selCombMatrix(int csel) {
   int k;
   Pi->submat(0,0,nA-1,nC-1)  = Pmatrix->submat(cumPind[csel],0,cumPind[csel]+nA-1,nC-1); 
   for(k=0;k<(nC-1); k++) {
    PiO->submat(0,0,nA-1,k) = Pi->submat(0,nC-1,nA-1,nC-1); //insert last column
   }
  }

 };


/*************************************/
/*Class defining for Search strategy**/
/*************************************/
//MasterMix0 goes through all combined combinations and fit model
class MasterMix0 {
 public: 
  int nL; //number of loci
  int nC; //number of contr
  int* CnA; //cdf index of alleles
  int n; //total length of data
  double scale; //used as scale
  //double LAsum,pGprod; //sum of probabilities and prod of geno-probabilities
  Mat<double> LAsum;
  Mat<double> tmpLA;
  vector<Loci*> *lociList;
  Mat<double> mhat;
  Mat<double> MD;
  Mat<double> *Y; 
  Mat<double> *X; 
  Mat<double> *O; 
  Mat<double> tmpX; 
  Mat<double> tmpO; 
  Mat<double> tmpW; 
  Mat<double> tmpYtilde; 
  Mat<double> Ytilde; //Y-O
  Mat<double> sumXtWX; //used to sum matrix
  Mat<double> sumXtWY; //used to sum matrix
  Mat<double> sumYtWY; //used to sum matrix
  Mat<double> ginvXtWX; //used to sum matrix
  Mat<double> prodPG; //used to prod PG

  //variables used in importance sample:
  Mat<uword> *simX; //simulated input-matrix
  Mat<uword> rowrange; //ranges of rows
  int M; //number of samples

  //constructor in search method:
  MasterMix0(int nL, int nC,int *CnA, vector<Loci*> *lociList, uword *simXtmp, int M) {
   this->nL = nL;
   this->nC = nC;
   this->CnA = CnA;
   this->lociList = lociList;
   n = CnA[nL]; //number of alleles
   LAsum.zeros(1,1); //init
   //pGprod=1;//.ones(1,1); //init
   Y = new Mat<double>(n,1);
   O = new Mat<double>(n,1);
   X = new Mat<double>(n,nC-1);
   int i,S=1; //S is number of genotype-combinations
   for(i=0; i<nL; i++) {
    Y->submat(CnA[i],0,CnA[i+1]-1,0) =  ((*lociList)[i]->Y)->col(0);
    S = S*((*lociList)[i]->nCombs);
   }
   //importance part: simX is a non-zero matrix with genotype combination indices
   this->M = M;
   simX = new Mat<uword>(simXtmp,M,nL,false); //store sampling-matrix
   rowrange.zeros(2,nL); //init matrix

   #ifdef DEBUG
    simX->print();
    Y->print();
    if(M==0) {  printf("Number of combinations to fit: %i",S);  }
   #endif
  }
  ~MasterMix0() {
   delete Y,O,X,simX;
  }
  
  //recursive function for going through all combinations to fit
  void recfit(int i) {
   //xi is locus visited
   int j,k;
   for(j=0;j<(*lociList)[i]->nCombs;j++) {  //number 2 means after cleared combinations
    scale = ((*lociList)[i]->S_Y)/2;
    (*lociList)[i]->updateCombMatrix(); //first update in loci
    tmpO = scale*(  ((*lociList)[i]->Pi)->col(nC-1));   
    tmpX = scale*(  ((*lociList)[i]->Pi)->cols(0,nC-2) - ((*lociList)[i]->PiO)->cols(0,nC-2) );
    tmpW = *((*lociList)[i]->iW); 
    tmpYtilde = *((*lociList)[i]->Y)-tmpO; //Response-offset
    (*lociList)[i]->XTWXi->submat(0,0,nC-2,nC-2) = trans(tmpX)*tmpW*(tmpX); //store matrix
    (*lociList)[i]->XTWYi->submat(0,0,nC-2,0) = trans(tmpX)*tmpW*tmpYtilde;
    (*lociList)[i]->YTWYi->submat(0,0,0,0) = trans(tmpYtilde)*tmpW*tmpYtilde;
    if(i==0) { 
     sumXtWX = ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); //add to matrix
     sumXtWY = ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0);
     sumYtWY = ((*lociList)[i]->YTWYi)->submat(0,0,0,0);
     prodPG = ((*lociList)[i]->PG)[j];
    } else { 
     sumXtWX = sumXtWX + ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); 
     sumXtWY = sumXtWY + ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0);
     sumYtWY = sumYtWY + ((*lociList)[i]->YTWYi)->submat(0,0,0,0);
     prodPG = prodPG*((*lociList)[i]->PG)[j];
    }//add to matrix
    if((i+1)==nL) {   //modelfit if in last loci:
     ginvXtWX = pinv(sumXtWX); //generalized inverse
     mhat = ginvXtWX*sumXtWY;
     MD = sumYtWY - 2*trans(sumXtWY)*mhat + trans(mhat)*sumXtWX*mhat; //MAHALONOBIS DISTANCE
     tmpLA = 1/sqrt(pow(2*PI*MD/n,n))*exp(-0.5*n)*prodPG;
     LAsum = LAsum + tmpLA;
    } else {
     recfit(i+1);
    }
    sumXtWX = sumXtWX - ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); //subtract from matrixsum
    sumXtWY = sumXtWY - ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0); //subtract from matrixsum
    sumYtWY = sumYtWY - ((*lociList)[i]->YTWYi)->submat(0,0,0,0); //subtract from matrixsum
    prodPG = prodPG/((*lociList)[i]->PG)[j];
   } //end combination:j
   (*lociList)[i]->combsel = 0; //important to reset to zero again after looping once
  } //end recurse


  //recursive function to go through all combination in a sorted matrix X:
  void recfitX(int i) {
   //X is a matrix with sorted combination-elements.
   int j,k;
   Col<uword> unXi = unique(simX->col(i)); //get unique values
   Col<uword> equalind,filt1,filt2;
   int unlen = unXi.n_elem;
   for(j=0;j<unlen;j++) {  //number 2 means after cleared combinations
     equalind = find(simX->col(i)==unXi[j]); //get equal indices
     rowrange(0,i) = equalind(0);
     rowrange(1,i) = equalind(equalind.n_elem-1);
     if(i>0) {
      filt1 = find(equalind>=rowrange(0,i-1)); 
      filt2 = find(equalind<=rowrange(1,i-1)); 
      if(filt1.n_elem==0 || filt2.n_elem==0) { 
       continue; //continoue loop if some are not found
      }
      rowrange(0,i) = max(rowrange(0,i),equalind(filt1(0)));
      rowrange(1,i) = min(rowrange(1,i),equalind(filt2(filt2.n_elem-1)));
     } //finished finding rows and combination

     //assign    
     scale = ((*lociList)[i]->S_Y)/2;
     (*lociList)[i]->selCombMatrix(unXi[j]-1); //update Pi-matrix with selected combination:
     tmpO = scale*(  ((*lociList)[i]->Pi)->col(nC-1));   
     tmpX = scale*(  ((*lociList)[i]->Pi)->cols(0,nC-2) - ((*lociList)[i]->PiO)->cols(0,nC-2) );
     tmpW = *((*lociList)[i]->iW); 
     tmpYtilde = *((*lociList)[i]->Y)-tmpO; //Response-offset
     (*lociList)[i]->XTWXi->submat(0,0,nC-2,nC-2) = trans(tmpX)*tmpW*(tmpX); //store matrix
     (*lociList)[i]->XTWYi->submat(0,0,nC-2,0) = trans(tmpX)*tmpW*tmpYtilde;
     (*lociList)[i]->YTWYi->submat(0,0,0,0) = trans(tmpYtilde)*tmpW*tmpYtilde; 
     if(i==0) { 
      sumXtWX = ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); //add to matrix
      sumXtWY = ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0);
      sumYtWY = ((*lociList)[i]->YTWYi)->submat(0,0,0,0);
      prodPG = ((*lociList)[i]->PG)[unXi[j]-1];
     } else { 
      sumXtWX = sumXtWX + ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); 
      sumXtWY = sumXtWY + ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0);
      sumYtWY = sumYtWY + ((*lociList)[i]->YTWYi)->submat(0,0,0,0);
      prodPG = prodPG*((*lociList)[i]->PG)[unXi[j]-1];
     }//add to matrix

     if((i+1)==nL) {   //modelfit if in last loci:
      ginvXtWX = pinv(sumXtWX); //generalized inverse
      mhat = ginvXtWX*sumXtWY;
      MD = sumYtWY - 2*trans(sumXtWY)*mhat + trans(mhat)*sumXtWX*mhat; //MAHALONOBIS DISTANCE
      tmpLA = 1/sqrt(pow(2*PI*MD/n,n))*exp(-0.5*n)*prodPG;
      LAsum = LAsum + (rowrange(1,i)-rowrange(0,i)+1)*tmpLA/M;
     } else {
      recfitX(i+1);
     }
     sumXtWX = sumXtWX - ((*lociList)[i]->XTWXi)->submat(0,0,nC-2,nC-2); //subtract from matrixsum
     sumXtWY = sumXtWY - ((*lociList)[i]->XTWYi)->submat(0,0,nC-2,0); //subtract from matrixsum
     sumYtWY = sumYtWY - ((*lociList)[i]->YTWYi)->submat(0,0,0,0); //subtract from matrixsum
     prodPG = prodPG/((*lociList)[i]->PG)[unXi[j]-1];
   } //end combination:j
  }

};


/************************/
/****Helpfunctions*******/
/************************/

//function calculating cumulative values in x
//x: steps, sz: vectorlength, y: vector
void cumsum(int* x, int *sz, int *y) {
 int i;
 y[0] = 0;
 for(i=0; i<*sz; i++) {
  y[i+1] = y[i] + x[i];
 }
}


extern "C" {


void doAnalysisC(double *LA,int *nA, int *nC, int *nL, double *Y, double *iW, int *nCombs, double *Plist, double *pGvec,uword *simX,int *M) {
 //simX is 0 if not importance sample are applied. M is number of samples in simX
 int i,j,nA2;
 int *CnA = new int[*nL+1]; //cumulative index for vector
 int *CnA2 = new int[*nL+1]; //cumulative index for matrix (vectorized matrix)
 int *nAsq = new int[*nL+1]; //squared numbers of nA
 int *nCG = new int[*nL+1]; //#combined genotypes*#contributors*#allele
 int *CGind = new int[*nL+1]; //cumulative index for Punlisted (loci-start)
 int *CPGind = new int[*nL+1]; //cumulative index for pGvec (loci-start)
 for(i=0; i<*nL; i++) { 
  nAsq[i] = nA[i]*nA[i]; 
  nCG[i] = nA[i]*(*nC)*nCombs[i];
 } //square numbers
 cumsum(nA, nL, CnA);
 cumsum(nAsq, nL, CnA2);
 cumsum(nCG, nL, CGind);
 cumsum(nCombs, nL, CPGind);

#ifdef DEBUG
 for(i=0; i<*nL; i++) {
  printf("nA[%i]=%d\n",i,nA[i]);
  printf("CnA[%i]=%d\n",i,CnA[i]);
  printf("CnA2[%i]=%d\n",i,CnA2[i]);
  printf("nCombs[%i]=%d\n",i,nCombs[i]);
  printf("nCG[%i]=%d\n",i,nCG[i]);
  printf("CGind[%i]=%d\n",i,CGind[i]);
  for(j=0; j<nA[i]; j++) { printf("Y[%d,%d]=%f\n",i,j,Y[CnA[i]+j]); }
  printf("\n"); 
 }
#endif

 //Create Loci-object and insert them into array.
 vector<Loci*> lociList;
 for(i=0; i<*nL; i++) {
  lociList.push_back(new Loci(i+1,nA[i], &Y[CnA[i]],*nC, &iW[CnA2[i]], nCombs[i] ,&Plist[CGind[i]], &pGvec[CPGind[i]]));
 }
#ifdef DEBUG
 printf("Check elements:");
 for(i=0; i<*nL; i++) {
  printf("%i",lociList[i]->locinr);
//  (lociList[i]->Y)->print();
//  (lociList[i]->Pmatrix)->print();
 }
#endif

  MasterMix0 *modelfit = new MasterMix0(*nL,*nC,CnA,&lociList, simX,*M);
  if(*M==0) { modelfit->recfit(0); } //if exact method
  if(*M>0) { modelfit->recfitX(0); } //if sampling  method
 *LA = modelfit->LAsum(0,0);
 //deallocation:
 delete modelfit;
 for(i=0; i<*nL; i++) {
  delete lociList[i]; //delete objects
 }
 delete[] CnA,CnA2,CGind,CPGind,nCG,nAsq;
}

} 
