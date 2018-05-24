# include <RcppArmadillo.h>
// [[ Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector EYgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z,
                        NumericVector se, NumericVector sp, 
                        int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj, u;
  NumericVector WW(N);
  
  for(g=0;g<GI;g++){
  for(i=0;i<N;i++){
    pi1=p(i);
    pi2=1-p(i);
    np=Y(i,1);
    for(l=0;l<np;l++){
      j=Y(i,(l+2));
      Zj=Z(j-1,0);
      cj=Z(j-1,1);
      tid=Z(j-1,2);
      sej=se(tid-1);
      spj=sp(tid-1);
      ybar=0;
      Y(i,0)=0;
      for(t=0;t<cj;t++){
        id=Z(j-1,(t+3));
        ybar=ybar+Y(id-1,0);
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
      if(ybar > 0){
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
      }
    }
    pistar=(pi1/(pi1+pi2));
    u = R::runif(0,1);
//    u=rand() / double(RAND_MAX);
    if(u<pistar){
      Y(i,0)=1;
    }else{Y(i,0)=0;}
    WW(i)=WW(i)+Y(i,0);
  }}  

  return WW;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix EYiYjgibbs_slow(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z,
                        NumericVector se, NumericVector sp, 
                        int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t;
  float pi1, pi2, pistar, sej, spj, u;
  NumericMatrix EYiYj(N,N);
  
  for(g=0;g<GI;g++){
  for(i=0;i<N;i++){
    pi1=p(i);
    pi2=1-p(i);
    np=Y(i,1);
    for(l=0;l<np;l++){
      j=Y(i,(l+2));
      Zj=Z(j-1,0);
      cj=Z(j-1,1);
      tid=Z(j-1,2);
      sej=se(tid-1);
      spj=sp(tid-1);
      ybar=0;
      Y(i,0)=0;
      for(t=0;t<cj;t++){
        id=Z(j-1,(t+3));
        ybar=ybar+Y(id-1,0);
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj));
      if(ybar > 0){
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj));
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));
      }
    }
    pistar=(pi1/(pi1+pi2));
    u = R::runif(0,1);
    //u=rand() / double(RAND_MAX);
    if(u<pistar){
      Y(i,0)=1;
    }else{Y(i,0)=0;}
    // WW(i)=WW(i)+Y(i,0);
    for(j=std::max(0,i - Z.nrow()*Z.nrow());j<=i;j++){
	    // If individuals are ordered according to the groups they are in,  
	    // then we only need to compute EYiYj for |i-j|< cj, where cj is the
	    // size of the group. This speeds things MASSIVELY.
	EYiYj(i,j) = EYiYj(i,j) + Y(i,0)*Y(j,0);
    }}}  

  return EYiYj;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix CovYiYjgibbs(int N, NumericVector p, NumericMatrix Y, NumericMatrix Z, NumericMatrix W,
                        NumericVector se, NumericVector sp, NumericVector EY,
                        int na, int GI) {
  int g, i, np, l, j, Zj, cj, ybar, tid, id, t, k;
  float pi1, pi2, pistar, sej, spj, u;
  NumericMatrix EYiYj(N,N);
  NumericMatrix CovYiYj(N,N);
  
  for(g=0;g<GI;g++){
  for(i=0;i<N;i++){
    pi1=p(i); // the value P_{\alpha,\beta}(Y_i=1|X_i)
    pi2=1-p(i);// the value P_{\alpha,\beta}(Y_i=0|X_i)
    np=Y(i,1); //number of assays in which individual i participated
    k=0;
    for(l=0;l<np;l++){
      j=Y(i,(l+2)); // number of the lth assay in which individual i participated
      Zj=Z(j-1,0); // outcome of the lth assay in which individual i participated
      cj=Z(j-1,1); // number of individuals participating in the lth assay in which individual i participated
      tid=Z(j-1,2); // an index to the sens or spec of the lth assay in which individual i participated
      sej=se(tid-1); // sensitivity of the lth assay in which individual i participated
      spj=sp(tid-1); // specificity of the same.
      ybar=0;
      Y(i,0)=0;
      for(t=0;t<cj;t++){  
        id=Z(j-1,(t+3)); // number of the tth individual in the lth assay in which individual i participated
        ybar=ybar+Y(id-1,0); // ybar > 0 gives the true disease status of the group of individuals in the lth assay of individual i
      }
      pi1=pi1*(sej*Zj + (1-sej)*(1-Zj)); // product of independent assay probabilities, given Y_i = 1 (group is positive), times P_{\alpha,\beta}(Y_i=1|X_i)
      if(ybar > 0){// since Y(i,0)=0, ybar > 0 only if OTHER individuals in the group are positive.
        pi2=pi2*(sej*Zj + (1-sej)*(1-Zj)); // product of independent assay probabilities, given Y_i = 0 but group is positive, times P_{\alpha,\beta}(Y_i=0|X_i)
      }else{
        pi2=pi2*((1-spj)*Zj + spj*(1-Zj));// product of independent assay probabilities, given Y_i = 0 and group is negative, times P_{\alpha,\beta}(Y_i=0|X_i)
      }
    }
    pistar=(pi1/(pi1+pi2));// probability that Y_i = 1 given all other disease statuses and outcomes of all assays in which individual i participated
    u = R::runif(0,1);
    //u=rand() / double(RAND_MAX);
    if(u<pistar){
      Y(i,0)=1;
    }else{Y(i,0)=0;}
    k = W(i,0);
    for(j=0;j<k;j++){
    	t = W(i,1+j)-1;
		EYiYj(i,t) = EYiYj(i,t) + Y(i,0)*Y(t,0);
    }}}

    // Now compute the covariance matrix
    for(i=0;i<N;i++){
    	k = W(i,0);
    	for(j=0;j<k;j++){
    		t = W(i,1+j)-1;
		CovYiYj(i,t) = EYiYj(i,t)/GI - EY(i)*EY(t) ;
    	}
    }
    								
	return CovYiYj;
}


// [[Rcpp::export()]]

Rcpp::List logistic_enet(Rcpp::NumericVector Yr, 
			Rcpp::NumericMatrix Xr,
			float lambda,
			Rcpp::NumericVector gammar,
			float theta,
			Rcpp::NumericVector binitr,
			float delta){
				
int n = Xr.nrow(), p = Xr.ncol()-1, i, k;
float uj, vj, wj, sj;
	
arma::mat X(Xr.begin(), n, p+1, false); 
arma::colvec Y(Yr.begin(),Yr.size(), false);
arma::colvec gamma(gammar.begin(),gammar.size(),false);
arma::colvec binit(binitr.begin(),binitr.size(),false);
	
arma::colvec b1 = binit;
arma::colvec b0 = binit;
arma::colvec diff = arma::ones(p+1);
	
arma::colvec b11 = binit;
arma::colvec b00 = binit;
arma::colvec diff2 = arma::ones(p+1);
	
arma::colvec px(n);
arma::colvec w(n);
arma::colvec r(n);
	
arma::colvec rj(n);
arma::mat Xj = X;
arma::colvec b11j(n);

i = 0;

while( (i < 500) & (diff.max() > delta))
{
	b0 = b1;
	
	px = exp(X * b0) / (1 + exp(X * b0));
	w = px % (1 - px);
	r = X * b0 + ( Y - px) / w;

	k = 0;
	diff2 = arma::ones(p+1);
	
	while( (k < 500) & (diff2.max() > delta))
	{
		
		b00 = b11; 
		
		b11(0) = sum( w % ( r - X.cols(1,p) * b00.rows(1,p))) / sum(w);
		
		
		for(int j=1; j < (p+1) ; j++)
		{
			Xj = X;
			Xj.shed_col(j);
		
			b11j = b11;
			b11j.shed_row(j);
			
			rj = Xj * b11j;
															
			uj = sum( w % ( r - rj) % X.col(j) );

			vj = theta * lambda * arma::as_scalar(gamma(j-1));//why do I have j-1 here ?
			
			wj = sum( w % pow(X.col(j),2)) + lambda * (1 - theta); 
			
			if( (uj > 0) & (vj < std::abs(uj)) ) // soft-thresholding
			{
				sj = uj - vj;
				
			} else if( (uj < 0) & (vj < std::abs(uj)))
			{
				
				sj = uj + vj;
				
			} else {
				
				sj = 0;
				
			}
			
			
			b11(j) = sj / wj;
			
		}
		
		
		diff2 = abs(b11 - b00);
		k++;
		
	}
	
	b1 = b11;


	diff = abs(b1 - b0);
	i++;
}

if(b1.has_nan()) 
{				
	Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
	
}

return Rcpp::List::create(Named("b") = b1,
					Named("lambda") = lambda,
					Named("theta") = theta,
					Named("gamma") = gamma
					);
}
		
// [[Rcpp::export()]]

Rcpp::List llj_array(	Rcpp::IntegerVector Zjr, 
			Rcpp::IntegerVector Zjc,
			Rcpp::IntegerVector Yji,
			Rcpp::IntegerVector whichjretest,
			Rcpp::NumericVector pxji,
			Rcpp::NumericVector Se,
			Rcpp::NumericVector Sp,
			int B){

arma::ivec Zjrow(Zjr.begin(),Zjr.size(),false);
arma::ivec Zjcol(Zjc.begin(),Zjc.size(),false);
arma::ivec Yj_retested(Yji.begin(),Yji.size(),false);			
arma::ivec retestint(whichjretest.begin(), whichjretest.size(),false);
arma::uvec retest = arma::conv_to<arma::uvec>::from(retestint);	
retest = retest - 1;
arma::colvec px(pxji.begin(),pxji.size(),false);
arma::colvec Sens(Se.begin(),Se.size(),false);
arma::colvec Spec(Sp.begin(),Sp.size(),false);				


int cj = Zjr.size(),
	M = whichjretest.size();

arma::uword cjsq = cj*cj;

float pZjrZjcYji_tYj, pZjrow_tYj, pZjcol_tYj, pYji_tYj, p1, p0, llj, Lj;

arma::colvec U(cjsq);
arma::vec tYj = arma::zeros(cjsq);			
arma::vec tYj_retested(M);
				
arma::umat Array(cj,cj);
arma::uvec tZjrow(cj);
arma::uvec tZjcol(cj);

Lj = 0;

for(int b=0;b<B;b++)								
{
	
	U.randu();
	
	for(arma::uword i=0;i < cjsq ;i++)
	{ 	
		if(U(i) <= px(i)) 
		{
			tYj(i) = 1;
			
		} else {
			
			tYj(i) = 0 ; 
		}
	
	}
	
	for(int row=0;row<cj;row++)
		for(int col = 0 ; col < cj ; col++)
			{
				
				Array(row,col) = tYj( col*cj + row);
				
			}
								
	pZjrow_tYj = 1;
	pZjcol_tYj = 1;
	
	for(int k=0;k<cj;k++)
	{
		
		tZjrow(k) = max(Array.row(k));
		tZjcol(k) = max(Array.col(k));
		
		p0 = (Sens(0)*Zjrow(k) + (1-Sens(0))*(1 - Zjrow(k))) * tZjrow(k) + ((1-Spec(0))*Zjrow(k) + Spec(0)*(1 - Zjrow(k))) * (1-tZjrow(k));
		p1 = (Sens(0)*Zjcol(k) + (1-Sens(0))*(1 - Zjcol(k))) * tZjcol(k) + ((1-Spec(0))*Zjcol(k) + Spec(0)*(1 - Zjcol(k))) * (1-tZjcol(k));
		
		pZjrow_tYj = pZjrow_tYj * p0 ;
		pZjcol_tYj = pZjcol_tYj * p1	 ;					
		
	}
		
	tYj_retested = tYj.rows(retest);
	
	pYji_tYj = 1;

	for(arma::uword i=0;i<tYj_retested.n_elem ;i++)
	{
		
		p1 = Sens(1)*Yj_retested(i) + (1-Sens(1))*(1-Yj_retested(i));
		p0 = (1-Spec(1))*Yj_retested(i)+ Spec(1)*(1-Yj_retested(i));
		
		pYji_tYj = pYji_tYj *  ( p1 *  tYj_retested(i) + p0 * ( 1 -  tYj_retested(i)) ) ;
				
	}
	
	
		pZjrZjcYji_tYj = pZjrow_tYj * pZjcol_tYj * pYji_tYj	;						
	
		Lj = Lj + pZjrZjcYji_tYj ;
						
}	

llj = log(Lj/B);
				
return Rcpp::List::create(Named("retest") = retest,
					Named("Yj_retested") = Yj_retested,
					Named("tYj_retested") = tYj_retested,
					Named("tYj") = tYj,
					Named("Se") = Sens,
					Named("Sp") = Spec,
					Named("Array") = Array,
					Named("tZjrow") = tZjrow,
					Named("tZjcol") = tZjcol,
					Named("pYji_tYj") = pYji_tYj,
					Named("pZjrZjcYji_tYj") = pZjrZjcYji_tYj,
					Named("Lj") = Lj,
					Named("llj") = llj,
					Named("U") = U
					);
			}
			

