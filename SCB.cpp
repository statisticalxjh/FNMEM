// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>  
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double Kh(double mu, double h) {

double output = (1-pow(mu/h, 2))*0.75*(abs(mu/h)<1)/h;
return output;

}

// [[Rcpp::export]]
mat spG(int n, int M, int p, int q, int ni, vec t, mat X, mat Y, vec Seed, double h1, mat betahat, mat bhat, mat Phihat, vec sigmahat, int s, vec tau){

int i, m;

//mat Sxx(p, p), XX(ni, p), phihat(p, p), XXi(ni, p), ddeltainv(q, q), Sigmai(ni, ni);
//colvec Sxy(p), YY(ni), wi(ni); 
//uvec indexi(ni); 
//
//mat sp = zeros<mat>(p, 1);
//vec K = zeros<vec>(M);

mat::fixed<2, 2> Sxx;
mat::fixed<5, 2> XX;
mat::fixed<2, 2> phihat;
mat::fixed<5, 2> XXi;
mat::fixed<2, 2> ddeltainv;
mat::fixed<5, 5> Sigmai;
mat::fixed<2, 1> Sxy;
mat::fixed<5, 1> YY;
mat::fixed<5, 1> wi;
mat::fixed<2, 1> sp; sp.zeros();
vec::fixed<25> K; K.zeros();
uvec::fixed<5> indexi; 

double Kh1;
for (m=0;m<M;m++){
  if (abs(t[m] - t[s])<h1){
      Kh1 = Kh(t[m] - t[s], h1);
      K[m] = Kh1;
      Sxx.zeros();
      Sxy.zeros(); 
      for (i=1;i<=n;i++){    
        indexi = find(Seed==i);
        XX = X.rows(indexi);
        YY = Y.submat(indexi[0], m, indexi[ni-1], m);
         phihat.zeros();
    	  phihat(0, 0) = betahat(0, m) + bhat(q*(i-1), m);
    	  phihat(1, 1) = betahat(1, m) + bhat(q*i-1, m);
    	  XXi = XX%exp(XX*phihat);
      	wi = YY - exp(XX.col(0)*phihat(0, 0) + XX.col(1)*phihat(1, 1)) + XXi*betahat.col(m) + XXi*bhat.submat(q*(i-1), m,q*i-1, m);
        ddeltainv = Phihat.rows(q*m, q*m+1)/sigmahat[m];
        Sigmai = inv(eye<mat>(ni, ni) + XXi*ddeltainv*XXi.t());
        Sxx += XXi.t()*Sigmai*XXi;
        Sxy += XXi.t()*Sigmai*(wi-XXi*betahat.col(m))*tau[i-1];
      }    
  sp += Kh1*solve(Sxx, Sxy);
  }
}  
sp = sp/sum(K);
return sp;
}

// [[Rcpp::export]]
mat SCB(int n, int M, int p, int q, int ni, int G, double h1, vec t, mat X, mat Y, vec Seed, mat betahat, mat bhat, mat Phihat, vec sigmahat){

int i, g, s;
//vec tau(n);
//mat sp(p, M);
//mat allC = zeros<mat>(p, G);

vec::fixed<50> tau;
mat::fixed<2, 25> sp;
mat::fixed<2, 500> allC; allC.zeros();


for (g=0;g<G;g++){   
    tau.randn();
    sp.zeros();
    for (s=0;s<M;s++){
        sp.col(s) = spG(n, M, p, q, ni, t, X, Y, Seed, h1, betahat, bhat, Phihat, sigmahat, s, tau);
    }
    for (i=0;i<p;i++){
        allC(i, g) = max(abs(sp.row(i)));        
    }
}
return allC;
}  
