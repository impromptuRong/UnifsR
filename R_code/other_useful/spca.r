################################################################################
#####    Type		: R functions
#####   Title		: Sparse Principal Component Analysis
#####	Version		: 1.0
#####	Date		: 2011-03-17
#####	Author		: Seokho Lee; Michael Epstein; Richard Duncan; Xihong Lin
#####	Reference	: Seokho Lee, Micheal Epstein, Richard Duncan, and Xihong Lin (2011) Sparse Principal Component Analysis for Identifying Ancestry-Informative Markers in Genome Wide Association Studies. Submitted.
#####	Maintainer	: Seokho Lee (leesh12@gmail.com; lees@hufs.ac.kr) 
#####	Depends		: 
#####	Description	: Sparse PCA using Lasso and Adaptive Lasso
#####	R topics	: SOFT.THRESHOLD(), S.PCA()
################################################################################

SOFT.TRESHOLD<-function(vec,lam){
    # Soft threshold function
    #
    # ARGUMENTS
    #	vec	: vector that soft tresholding is applied
    #	lam	: non-negative scalar- or vector-valued soft tresholding parameter. If lam is a vector, then the length of lam should be the same as that of vec
    #
    # VALUES
    #	res	: resulting vector of the same size as vec
    #
    if ( length(lam)>1 & length(lam)!=length(vec) ) {
        cat('\n ERROR: THE SIZE OF THE SECOND ARGUMENT SHOULD BE 1 OR THE SAME AS THE SIZE OF THE FIRST ARGUMENT.\n')
        return ( 0 )
    }
    idx.1<-which(vec < -lam)
    idx.2<-which(vec > lam)
    res<-rep(0,length(vec))
    if (length(lam)==1){
        if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam
        if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam
    } else if (length(lam)>1) {
        if ( length(idx.1)>0 ) res[idx.1]<- vec[idx.1]+lam[idx.1]
        if ( length(idx.2)>0 ) res[idx.2]<- vec[idx.2]-lam[idx.2]
    }
    return( res )
}

S.PCA<-function(mat,lambda,method='lasso',prt=FALSE){
    # Sparse PCA function.
    # Dependence : SOFT.THRESHOLD()
    #
    # ARGUMENTS
    #	mat		: column-centered data matrix (n times d)
    #	lambda	: scalar- or vector-valued penalty parameter
    #	method	: sparse PCA method. Either 'lasso' or 'alasso'
    #	prt		: logical. Iteration status is printed if TRUE
    #
    # VALUES
    #	u		: vector of size n. It is principal component score of unit size
    #	v		: vector of size d. It is the estimated principal component.
    #	s		: scalar. It is a scale factor (or variance) of principal component score. s*u is the vector a in Lee et al. (2011)
    #	bic		: vector. Bic values at each lambda grid are returned
    #	lam		: same as lambda.
    #	nzr		: scalar. The number of nonzero loadings on v
    #
    # Note : Initial estimate of v is given by the standard pc
    #
    X<-mat;	sv<-svd(X)
    
    bic<-NULL; bic.min<-Inf; bic.u<-NULL; bic.v<-NULL; bic.s<-NULL; nzr<-NULL
    for(lam in lambda){
        if (prt==TRUE) cat('===== lam =', lam, 'is starting.... =====\n')
        u<-sv$u[,1]; t.v<-sv$v[,1]*sv$d[1]
        obj.old<-n*d; iter<-1
        repeat{
            Xu<-colSums(X*u)
            if (method=='lasso') {
                v<-SOFT.TRESHOLD(Xu,lam)
                if ( sum(v)!=0 ) v<-v/sqrt(sum(v^2))
                s<-as.vector(t(u)%*%X%*%v)
                obj.new<-mean((X-s*outer(u,v))^2)+2*lam*sum(abs(v))
            }
            if (method=='alasso') {
                v<-SOFT.TRESHOLD(Xu,lam/abs(t.v))
                if ( sum(v)!=0 ) v<-v/sqrt(sum(v^2))
                s<-as.vector(t(u)%*%X%*%v)
                t.v[t.v==0]<-1E-12
                obj.new<-mean((X-s*outer(u,v))^2)+2*lam*sum(abs(v)/abs(t.v))
            }
            dif<-obj.old-obj.new
            if (prt==TRUE) cat('( iter, obj, dif ) = (', iter, ',', obj.new, ',', dif, ')\n')
            if ( abs(dif)<1E-4 | iter>100 ) break
            u<-as.vector(X%*%v)
            if ( sum(u)!=0 ) u<-u/sqrt(sum(u^2))
            s<-as.vector(t(u)%*%X%*%v) 
            iter<-iter+1; obj.old<-obj.new
        }
        
        #### BIC computation
        bic.tmp <- ifelse( sum(u^2)==0, Inf, n*d*log(mean((X-s*outer(u,v))^2)) + n*log(mean(u^2)*s^2) + log(n*d)*(n+sum(v!=0)) )
        bic<-append(bic, bic.tmp)
        nzr<-append(nzr, sum(v!=0))
        if (bic.min>bic.tmp) {
            bic.min<-bic.tmp
            bic.u<-u; bic.v<-v; bic.s<-s
        }
        
    }
    
    return ( list(u=bic.u, v=bic.v, s=bic.s, bic=bic, lambda=lambda, nzr=nzr) )
    
}