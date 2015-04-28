

grid_solver <- function(data_A, group, branch, r, alpha, D){
    A <- data.matrix(data_A)
    npar <- nrow(A)
    nfea <- ncol(A)
    y <- (as.numeric(group)-1.5)*2
    dA <- A %*% diag(sqrt(branch*abs(sin(r*pi)/pi)))
    
    #######   Interior Point Method for alpha=0   #######
    #	Kappa <- kernelMatrix(kerfun, data_B, data_B)
    Hess0 <- (dA*y)%*%t(dA*y)
    Grad0 <- function(u){	return(t(u)%*%Hess0)	}
    fn0 <- function(u){	return(1/2*t(u)%*%Hess0%*%u)	}
    ####   Solve Quandratic Formula   ####
    if(any(diag(Hess0)==0)){
        method <- "LU"
    } else{
        method <- "CHOL"
    }
    ans0 <- LowRankQP(Vmat=Hess0, dvec=rep(0,npar), Amat=rbind(rep(1,npar),y), bvec=c(200,0), u=rep(D*100,npar), method=method)
    
    opt_u0 <- as.vector(ans0$alpha)
    opt_p0 <- t(dA)%*%(opt_u0*0)
    
    #######   Linear Programming for alpha=1   #######
    fn1 <- function(u){	return(sum(abs(t(dA)%*%(u*y))))	}
    ####   Solve Linear Formula  ####
    E_mat <- cbind(rbind(rep(1,npar),y),matrix(0,2,nfea))
    G_mat <- cbind(rbind(-t(dA)%*%diag(y),t(dA)%*%diag(y),diag(-1,npar,npar)),rbind(diag(1,nfea,nfea),diag(1,nfea,nfea),matrix(0,npar,nfea)))
    ans1 <- linp(E=E_mat, F=c(200,0), G=G_mat, H=c(rep(0,2*nfea),rep(-D*100,npar)), Cost=c(rep(0,npar),rep(1,nfea)), ispos=TRUE)
    opt_u1 <- ans1$X[1:npar]
    opt_p1 <- t(dA)%*%(opt_u1*y)
    
    opt_u <- opt_u1
    opt_p <- opt_p1
    
    #######   non-linear programming   #######
    cost=c(rep(0,npar),rep(1,nfea))
    fn <- function(x){	alpha*cost%*%x+sqrt(alpha^2*(cost%*%x)^2+4*(1-alpha)*cost%*%x^2)	}
    grad <- function(x){	alpha*cost+(alpha^2*(cost%*%x)*cost+4*(1-alpha)*cost*x)/sqrt(alpha^2*(cost%*%x)^2+4*(1-alpha)*cost%*%x^2)	}
        
    #######   Calculate maximum Distance   #######
    #####        0<opt_u[i]<D,	(a-b)/2+y[i]*(a+b)/2 = y[i]*t(w)%*%A[i,],		SV==0
    #####        opt_u[i]=D,	(a-b)/2+y[i]*(a+b)/2 > y[i]*t(w)%*%A[i,],		SV==1
    #####        opt_u[i]=0,	(a-b)/2+y[i]*(a+b)/2 < y[i]*t(w)%*%A[i,],		SV==-1
    tol <- 10^-4
    SV <- rep(0,npar)
    SV[opt_u<=tol] <- -1
    SV[opt_u>=D*100-tol] <- 1
    rescale <- abs(opt_p)<tol
    
    k_fac <- sqrt((alpha/2*sum(abs(opt_p)))^2+(1-alpha)*sum(opt_p^2))
    opt_w <- (1-alpha)*opt_p/k_fac + (alpha/2+alpha^2/4*sum(abs(opt_p))/k_fac)*ifelse(rescale,0,sign(opt_p))
    dA_ori <- dA%*%opt_w*y
    dA_res <- cbind(-(1+y)/2,(1-y)/2,dA[,rescale]*y*(alpha/2+alpha^2/4*sum(abs(opt_p))/k_fac))
    
    res <- varranges(E=dA_res[SV==0,], F=-dA_ori[SV==0,], G=rbind(-dA_res[SV!=0,]*SV[SV!=0],cbind(0,0,rbind(diag(1,sum(rescale)),diag(-1,sum(rescale))))), H=c(dA_ori[SV!=0,]*SV[SV!=0],rep(-1,2*sum(rescale))), EqA=c(1,-1,rep(0,sum(rescale))), EqB=NULL, ispos=FALSE, tol=1e-8)
    xranges(E=dA_res[SV==0,], F=-dA_ori[SV==0,], G=rbind(-dA_res[SV!=0,]*SV[SV!=0],cbind(0,0,rbind(diag(1,sum(rescale)),diag(-1,sum(rescale))))), H=c(dA_ori[SV!=0,]*SV[SV!=0],rep(-1,2*sum(rescale))), ispos=FALSE, tol=1e-8)
    d <- res[2]
    return(d)
}


##########   Sequential Quandratic Programming   ###############
SQP <- function(x0, f, grad, hess, Amat, bvec, meq=0, maxiter=100) {
    require(quadprog)
    for(i in 1:maxiter) {
        Dmat <- hess(x0)
        # If the problem is not convex, 
        # the hessian is not always positive semi-definite: 
        # force it to be so.
        e <- eigen(Dmat)
        stopifnot(any(e$values>0))
        if(any(e$values<-1e-6)){	warning("Non-convex problem: Hessian not positive definite")	}
        Dmat <- e$vectors %*% diag(pmax(1e-6,e$values)) %*% t(e$vectors)
        dvec <- - grad(x0) + t(x0) %*% hess(x0)
        r <- solve.QP(Dmat, dvec, Amat, bvec, meq=meq)
        #print(f(r$solution))
        if(sum(abs(r$solution-x0)) < 1e-8){	break	}
        x0 <- r$solution
    }
    r$sqp_iterations  <- i
    r$sqp_convergence <- i < maxiter
    r
}

####   Initial a p0 Value s.t.: sum(p)=200,sum(p*y)=0,0<=p<=D ####
init_p0 <- function(group){
    p0 <- (as.numeric(group)-1.5)*2
    n1 <- length(p0[p0==1])
    n2 <- length(p0[p0==-1])
    tmp1 <- rep(floor(100/n1),n1)
    if(100-sum(tmp1)){	tmp1[c(1:(100-sum(tmp1)))] <- floor(100/n1)+1	}
    p0[p0==1] <- tmp1
    tmp2 <- rep(floor(100/n2),n2)
    if(100-sum(tmp2)){	tmp2[c(1:(100-sum(tmp2)))] <- floor(100/n2)+1	}
    p0[p0==-1] <- tmp2
    return(p0)
}
###########################

####   Random Generate r vector based on original vector  ####
par_rand <- function(par){
    tmp <- runif(length(par),-80,80)
    while(abs(sum(tmp))>10){	tmp <- runif(length(par),-80,80)	}
    if(sum(tmp)>0){	tmp[which.max(tmp)] <- max(tmp)-sum(tmp)	}
    if(sum(tmp)<0){	tmp[which.min(tmp)] <- min(tmp)-sum(tmp)	}
    x <- par + tmp
    while(any(abs(x)>=100)){
        d <- max(c(-99-min(x),max(x)-99))
        x[which.max(x)] <- max(x)-d
        x[which.min(x)] <- min(x)+d
    }
    return(x)
}
###########################

#########   take a sample species matrix and a phylogenetic tree    #############
fastUnifrac <- function(phylo, tree=NULL, seq.depth=TRUE, parallel=FALSE, method=c("Edge","dPCoA","non_w","w.non_nor","w.nor")){
    if(class(phylo)=="phyloseq"){
        otu_table <- otu_table(phylo)@.Data
        if(otu_table(phylo)@taxa_are_rows){otu_table <- t(otu_table)}
        tree <- phy_tree(phylo)
    } else {
        otu_table <- as.matrix(phylo)
    }
    
    if(class(tree)!="phylo" || class(otu_table)!="matrix"){ stop("Input Data or Tree are Invalid")  }
    if(is.null(tree$edge.length)){  stop("Tree has no branch lengths, cannot compute UniFrac")  }
    Nsample <- nrow(otu_table)
    Ntaxa <- length(tree$tip.label)
    Nnode <- tree$Nnode
    Nedge <- nrow(tree$edge)
    edge_length <- tree$edge.length
    if((Ntaxa+Nnode-1)!=Nedge){ stop("Tree structure error")    }
    
    otu_table <- otu_table[, tree$tip.label]
    if(seq.depth==TRUE){
        seq.depth <- rowSums(otu_table)
    } else if(seq.depth==FALSE){
        seq.depth <- 1
    }
    otu_table <- otu_table/seq.depth
    
    extend_par <- function(edge, tree=tree, Ntaxa=Ntaxa){
        int.nodes <- tree$edge[edge, 2]
        if(int.nodes<=Ntaxa){	return(int.nodes)	}
        sons <- c()
        repeat{
            int.nodes <- tree$edge[which(tree$edge[, 1]%in%int.nodes), 2]
            if(length(int.nodes)==0){	break	}
            sons <- c(sons, int.nodes)
        }
        sons <- sons[sons<=Ntaxa]
        return(sons)
    }
    extend <- function(edge){extend_par(edge, tree, Ntaxa)}
    
    ###   Non-parallel and foreach   ###
    if(parallel){
        edge_list <- try(foreach(edge=1:Nedge) %dopar% extend(edge))
        if(isTRUE(all.equal(class(edge_list), "try-error"))){
            sfInit(parallel=TRUE, cpus=8)
            # sfClusterSetupRNG(type="SPRNG", seed=12345)
            edge_list <- sfLapply(1:Nedge, extend_par, tree=tree, Ntaxa=Ntaxa)
        }
    } else{	system.time(edge_list <- lapply(1:Nedge, extend))	}
    edge_name <- c(tree$tip.label, tree$node.label)[tree$edge[,2]]
    if(is.rooted(tree)){
        edge_list[[Nedge+1]] <- c(1:Nedge)
        edge_name <- c(edge_name, tree$node.label[which(table(c(tree$edge))==2)-Ntaxa])
        edge_length <- c(edge_length, 0)
    }
    
    edge_bool <- sapply(edge_list, function(x){1:Ntaxa %in% x})*1;
    rownames(edge_bool) <- tree$tip.label
    colnames(edge_bool) <- edge_name;
    edge_matrix <- otu_table %*% edge_bool;
    
    ###   Organize Return List   ###
    method = c("Edge","dPCoA","non_w","w.non_nor","w.nor") %in% method;
    if(!sum(method)){
        warning("Input Methods is Invalid!\nReturn edge_list, edge_matrix, edge_bool Only.\nProper Methods are: Edge; dPCoA; non_w; w.non_nor; w.nor");
        method = c(1,0,0,0,0)
    }
    
    result <- list();
    if(method[1]){
        result[["edge_list"]] <- edge_list
        result[["edge_bool"]] <- edge_bool
        result[["edge_matrix"]] <- edge_matrix
    }
    if(method[2]){
        D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
        for(i in 1:Nsample){
            for(j in 1:Nsample){
                D[i,j] <- (edge_matrix[i,]-edge_matrix[j,])^2 %*% edge_length
            }
        }
        result[["dist"]][["dPCoA"]] = as.dist(D);
    }
    if(method[3]){
        D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
        for(i in 1:Nsample){
            for(j in 1:Nsample){
                D[i,j] <- (abs((edge_matrix[i,]>0)-(edge_matrix[j,]>0)) %*% edge_length)/(abs((edge_matrix[i,]>0)|(edge_matrix[j,]>0)) %*% edge_length)
            }
        }
        result[["dist"]][["non_w"]] = as.dist(D);
    }
    if(method[4]){
        D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
        for(i in 1:Nsample){
            for(j in 1:Nsample){
                D[i,j] <- abs(edge_matrix[i,]-edge_matrix[j,]) %*% edge_length
            }
        }
        result[["dist"]][["w.non_nor"]] = as.dist(D);
    }
    if(method[5]){
        D <- matrix(0, Nsample, Nsample, dimnames=list(rownames(otu_table), rownames(otu_table)))
        leaf_depth <- setNames(node.depth.edgelength(tree), c(tree$tip.label, tree$node.label))[colnames(otu_table)]
        for(i in 1:Nsample){
            for(j in 1:Nsample){
                D[i,j] <- (abs(edge_matrix[i,]-edge_matrix[j,]) %*% edge_length)/((otu_table[i,]+otu_table[j,]) %*% leaf_depth)
                }
        }
        result[["dist"]][["w.nor"]] = as.dist(D);
    }
    return(result)
}


sim_group <- function(tree=NULL, Nfeature=NULL, Nnoise=NULL, srange=c(4,10)){
    if(class(tree)!="phylo"){
        nOTU <- try(round(tree)[1], silent=TRUE)
        if(class(nOTU)=="try-error"){nOTU <- 1000}
        tree <- rtree(nOTU, rooted = TRUE)
        digit <- floor(log10(nOTU))+1
        tree$tip.label <- sprintf(paste("Otu%0",digit,"d",sep=""), 1:nOTU)
    }
    if(is.null(Nfeature)){Nfeature <- round((tree$Nnode+1)/20)}
    if(is.null(Nnoise)){Nnoise <- round((tree$Nnode+1)/40)}
    
    extend_par <- function(edge, tree=tree, Ntaxa=Ntaxa){
        int.nodes <- tree$edge[edge, 2]
        if(int.nodes<=Ntaxa){    return(int.nodes)	}
        sons <- c()
        repeat{
            int.nodes <- tree$edge[which(tree$edge[, 1]%in%int.nodes), 2]
            if(length(int.nodes)==0){	break	}
            sons <- c(sons, int.nodes)
        }
        sons <- sons[sons<=Ntaxa]
        return(sons)
    }
    
    nOTU <- tree$Nnode+1
    tree$edge.length <- ifelse(tree$edge.length>0.03, tree$edge.length, tree$edge.length+0.03)
    edgeorder <- order(tree$edge.length)
    edge_list <- lapply(edgeorder, function(x){extend_par(x, tree, nOTU)})
    edge_node <- sapply(edge_list, length)
    edge_dist <- tree$edge.length[edgeorder]
    edge_cand <- edge_node>srange[1] & edge_node<srange[2]
    edge_rest <- edge_cand
    sr.nei <- c(1:nOTU)
    
    #####  Feature Group: OTUs with Large grouped branch length   #####
    feature_list <- list()
    for(i in 1:Nfeature){
        s <- length(edge_cand)-which.max(rev(edge_rest))+1
        s.nei <- edge_list[[s]]
        s.inner <- sapply(edge_list, function(x){length(x)<length(s.nei)&&sum(x%in%s.nei)==length(x)})
        edge_dist[s] <- max(edge_dist[s], max(edge_dist[s.inner])+0.1)
        reduce <- (edge_dist>=edge_dist[s]-0.5) & s.inner
        edge_dist[s] <- edge_dist[s]*1.5
        edge_dist[reduce] <- sqrt(edge_dist[reduce])/2
        dis <- sapply(edge_list, function(x){sum(x%in%s.nei)==0})
        edge_rest <- edge_rest & dis
        sr.nei <- setdiff(sr.nei, s.nei)
        feature_list[[i]] <- s.nei
    }
    
    #####   Noise Group: OTUs with Small grouped branch length   #####
    noise_list <- list()
    for(i in 1:Nnoise){
        s <- which.max(edge_rest)
        s.nei <- edge_list[[s]]
        s.inner <- sapply(edge_list, function(x){length(x)<length(s.nei)&&sum(x%in%s.nei)==length(x)})
        reduce <- (edge_dist>=edge_dist[s]+0.05) & s.inner
        edge_dist[reduce] <- runif(sum(reduce), edge_dist[s], edge_dist[s]+0.1)
        dis <- sapply(edge_list, function(x){sum(x%in%s.nei)==0})
        edge_rest <- edge_rest & dis
        sr.nei <- setdiff(sr.nei, s.nei)
        noise_list[[i]] <- s.nei
    }
    
    #####    Single OTU    #####
    tree$edge.length <- edge_dist[order(edgeorder)]
    tip.color <- rep("black", nOTU)
    names(tip.color) <- tree$tip.label
    tip.color[unlist(feature_list)] <- "green"
    tip.color[unlist(noise_list)] <- "red"
    plot(tree, tip.color=tip.color, font=1, cex=1, direction="downwards")
    return(list(tree=tree, Ntaxa=length(feature_list)+length(noise_list)+length(sr.nei), feature=feature_list, noise=noise_list, rest=sr.nei))
}


#######   Fit Abundance Model   #######
sim_abufit <- function(data, group, plot=TRUE){
    n <- nlevels(group)
    if(plot){par(mfrow=c(floor(sqrt(n)), ceiling(n/floor(sqrt(n)))))}
    
    parlist <- list()
    for(type in levels(group)){
        sub <- data[group==type, ]
        sub <- sub[, colSums(sub)!=0]
        ntaxa <- ncol(sub)
        prob <- sort(colSums(sub), decreasing=TRUE)/sum(sub)
        
        #######   Fit (Truncated) Log-Normal Distribution   #######
        sub.mod <- radfit(colSums(sub))
        par <- sub.mod$models$Lognormal$coefficients
        err <- ntaxa
        k <- 0
        repeat{
            p <- dlnorm(1:ntaxa+k, meanlog=par[1], sdlog=par[2])
            se <- sum(abs(p/sum(p)-prob))
            if(se>err){break}
            err <- se
            k <- k+1
        }
        par["k"] <- k-1
        parlist[[type]] <- par
        
        if(plot){
            abu_p <- dlnorm(1:ntaxa+par[3], meanlog=par[1], sdlog=par[2])
            abu_p <- abu_p/sum(abu_p)
            plot(prob, xlim=c(1, 100), col="red", ylim=c(0, max(prob, abu_p)+0.03), type='l', xlab=type, ylab="Abu")
            points(abu_p, col="blue")
        }
    }
    return(parlist)
    
    #         prestonfit(colSums(sub), tiesplit=TRUE)
    #         prestondistr(colSums(sub), truncate=-1)
    #         par <- prestonfit(colSums(sub), tiesplit=TRUE)$coefficients
    
    #         sub.sep <- radfit(sub)
    #         par_ind <- sapply(sub.sep, function(x){x$models$Lognormal$coefficients})
    #         abu_p <- dlnorm(1:nAbu+1, meanlog=par[1], sdlog=par[2])
    #         plot(abu_p/sum(abu_p), col="blue", ylim=c(0,0.5))
    #         abu_pave <- dlnorm(1:nAbu+2, meanlog=mean(parlist[1,]), sdlog=mean(parlist[2,]))
    #         points(abu_pave/sum(abu_pave), col="red")
    #         
    #         for(i in 1:ncol(parlist)){
    #             p <- dlnorm(1:nAbu+2, meanlog=parlist[1,i], sdlog=parlist[2,i])
    #             points(p/sum(p), type='l', col="yellow")
    #             p <- c(unlist(sort(sub[i,], decreasing=TRUE)/sum(sub[i,])))
    #             points(p/sum(p), type='l', col="green")
    #         }
    #         points(abu_p/sum(abu_p), col="blue")
    #         points(abu_pave/sum(abu_pave), col="red")
    
    #######   Neg-Geometric Distribution   #######
    # a <- sort(dnbinom(1:100, 30, prob=0.75), decreasing=TRUE)
    # points(1000/sum(a)*a, col="purple")
    # b <- sort(dnbinom(1:100, 10, prob=0.7), decreasing=TRUE)
    # points(1000/sum(b)*b, col="green")
    # k1 = 30
    # p1 = 0.75
    # mu1 = 10
    # 
    # k2 = 5
    # p2 = 0.6
    # mu2 = 10/3
    # 
    # abu_p1 <- sort(colSums(nor.g1)/nrow(nor.g1), decreasing=TRUE)
    # abu_p2 <- sort(colSums(nor.g2)/nrow(nor.g2), decreasing=TRUE)
    # abu_p1 <- as.vector(abu_p1[1:(nOTU-length(s1.nei)-length(s2.nei)-length(s3.nei)-length(s4.nei)+4)])
    # abu_p2 <- as.vector(abu_p2[1:(nOTU-length(s1.nei)-length(s2.nei)-length(s3.nei)-length(s4.nei)+4)])
}


#######   Simulate grouped OTUs   #######
sim_matrix <- function(sim.pack, abu_mod, N=c(50,50), cutoff=NULL){
    sim.tree <- sim.pack$tree
    otuName <- sim.tree$tip.label
    nOTU <- length(otuName)
    Ntaxa <- sim.pack$Ntaxa
    
    par_p1 <- abu_mod[[1]]
    abu_p1 <- dlnorm(1:Ntaxa+par_p1[3], meanlog=par_p1[1], sdlog=par_p1[2])
    abu_p1 <- abu_p1/sum(abu_p1)
    par_p2 <- abu_mod[[2]]
    abu_p2 <- dlnorm(1:Ntaxa+par_p2[3], meanlog=par_p2[1], sdlog=par_p2[2])
    abu_p2 <- abu_p2/sum(abu_p2)
    
    maxdiff <- abu_p2[1]-abu_p1[1]
    if(maxdiff<0){
        maxdiff <- -maxdiff
        gname <- rev(gname)
        N <- rev(N)
        tmp <- abu_p1
        abu_p1 <- abu_p2
        abu_p2 <- tmp
    }
    print(maxdiff)
    
    gname <- names(abu_mod)[c(1,2)]
    if(is.null(gname)){gname <- c("A","B")}
    nGroup1 <- N[1]
    nGroup2 <- N[2]
    p_mat1 <- matrix(0, nGroup1, nOTU)
    p_mat2 <- matrix(0, nGroup2, nOTU)
    colnames(p_mat1) <- otuName
    colnames(p_mat2) <- otuName
    
    #######   Generate dirichlet distribution for features   #######
    if(is.null(cutoff)){cutoff <- c(0.010, 0.010, min(0.03,maxdiff/2), max(0.06,maxdiff))}
    diff_mat1 <- matrix(abu_p1, length(abu_p1), length(abu_p2))
    diff_mat2 <- matrix(abu_p2, length(abu_p1), length(abu_p2), byrow=TRUE)
    diff_mat0 <- abs(diff_mat1-diff_mat2)
    
    repeat{
        diff_feat <- diff_mat0>=cutoff[3] & diff_mat0<=cutoff[4]
        index1 <- which(abu_p1>=cutoff[1])
        index2 <- which(abu_p2>=cutoff[2])
        index <- rep(0, length(index2))
        for(j in index2){
            for(i in index1){
                if(diff_feat[i,j]){
                    adj <- TRUE
                    index_tmp <- index
                    slimit <- index_tmp %in% i
                    while(sum(slimit) && adj){
                        adj <- adj && diff_feat[index_tmp[slimit]+1, which(slimit)]
                        index_tmp[slimit] <- index_tmp[slimit]+1
                        slimit <- index %in% index_tmp[slimit]
                    }
                    if(adj){
                        index <- index_tmp
                        index[j] <- i
                    }
                }
                if(index[j]){break}
            }
        }
        index
        
        if(sum(index!=0)>=length(sim.pack$feature)){
            index2 <- sample(which(index!=0), length(sim.pack$feature))
            index1 <- index[index2]
            index <- cbind(index1, index2)
            break;
        }
        cutoff <- cutoff - c(0.001,0.001,0.001,0)
    }
    
    for(i in 1:nrow(index)){
        s <- sim.pack$feature[[i]]
        p_mat1[, s] <- abu_p1[index[i,1]] * rdirichlet(nGroup1, rep(1,length(s)))
        p_mat2[, s] <- abu_p2[index[i,2]] * rdirichlet(nGroup2, rep(1,length(s)))
        print(c(abu_p1[index[i,1]], abu_p2[index[i,2]], abs(abu_p1[index[i,1]]-abu_p2[index[i,2]])))
    }
    abu_p1 <- abu_p1[-index1]
    abu_p2 <- abu_p2[-index2]
    
    #######   Generate exponential distribution for noises   #######
    for(i in 1:length(sim.pack$noise)){
        diff_mat1 <- matrix(abu_p1, length(abu_p1), length(abu_p1))
        diff_mat2 <- matrix(abu_p2, length(abu_p2), length(abu_p2), byrow=TRUE)
        s <- sim.pack$noise[[i]]
        print(i)
        index <- arrayInd(sample(which(abs(diff_mat1-diff_mat2)<0.01 & abs(diff_mat1-0.01)<0.005 & abs(diff_mat2-0.01)<0.005), 1), .dim=dim(diff_mat1))
        p_mat1[, s] <- diff_mat1[index] * rdirichlet(nGroup1, 2^(1:length(s))/2)
        p_mat2[, s] <- diff_mat2[index] * rdirichlet(nGroup2, 2^(length(s):1)/2)
        abu_p1 <- abu_p1[-index[1]]
        abu_p2 <- abu_p2[-index[2]]
    }
    
    #######   Other OTUs   #######
    nRest <- length(sim.pack$rest)
    rpick <- sample(1:nRest)
    abu_p1[rpick[1:round(nRest/10)]] <- sort(abu_p1[rpick[1:round(nRest/10)]])
    p_mat1[, sim.pack$rest] <- t(matrix(abu_p1[rpick], nRest, nGroup1))
    p_mat2[, sim.pack$rest] <- t(matrix(abu_p2[rpick], nRest, nGroup2))
    
    #######   Random Sample based on probability   #######
#    r_mat1 <- t(sapply(1:nGroup1, function(x){table(factor(sample(otuName, 1000, replace=TRUE, prob=p_mat1[x,]),levels=otuName))}))
#    r_mat2 <- t(sapply(1:nGroup2, function(x){table(factor(sample(otuName, 1000, replace=TRUE, prob=p_mat2[x,]),levels=otuName))}))
    rownames(p_mat1) <- sprintf(paste(gname[1],"%0",floor(log10(nGroup1))+1,"d",sep=""), 1:nGroup1)
    rownames(p_mat2) <- sprintf(paste(gname[2],"%0",floor(log10(nGroup2))+1,"d",sep=""), 1:nGroup2)
    sim.mat <- rbind(p_mat1, p_mat2)
    sim.pack$abu_mod <- abu_mod
    sim.pack$prob.mat <- sim.mat
    return(sim.pack)
}


enet_cir <- function(r, alpha){
    plotseq <- seq(-1,1,by=0.001)
    x1 <- c(plotseq,rev(plotseq))*r
    y1 <- c(1-abs(plotseq),abs(plotseq)-1)*r
    x2 <- c(sin(plotseq*pi/2), rev(sin(plotseq*pi/2)))*r
    y2 <- c(cos(plotseq*pi/2),-rev(cos(plotseq*pi/2)))*r
    x3 <- c(plotseq, rep(1,length(plotseq)), rev(plotseq), rep(-1,length(plotseq)))*r
    y3 <- c(rep(1,length(plotseq)), rev(plotseq), rep(-1,length(plotseq)), plotseq)*r
    
    #####  (1-alpha)*(x^2/r^2+y^2/r^2)+alpha*(abs(x/r)+abs(y/r))=1
    a <- (1-alpha)
    b <- alpha
    c <- (1-alpha)*plotseq^2+alpha*abs(plotseq)-1
    y4 <- (-b+sqrt(b^2-4*a*c))/2/a
    x4 <- c(plotseq,rev(plotseq),y4,-y4)
    y4 <- c(y4,-y4,plotseq,rev(plotseq))
    
#     #####  (1-alpha)*(x^2/r^2+y^2/r^2)+alpha*max(abs(x/r)+abs(y/r))=1
#     x_norm5 <- seq(-1,1,by=0.01)*r
#     ytmp1 <- sqrt((1-alpha*abs(x_norm5)/r-(1-alpha)*x_norm5^2/r^2)/(1-alpha)*r^2)
#     a=(1-alpha)
#     b=alpha*r
#     c=(1-alpha)*x_norm4^2-r^2
#     ytmp2 <- (-b+sqrt(b^2-4*a*c))/2/a
#     y_norm5 <- ifelse(ytmp1<ytmp2,ytmp1,ytmp2)
#     a=2*(1-alpha)
#     b=alpha*r
#     c=-r^2
#     x0 = (-b+sqrt(b^2-4*a*c))/2/a
#     x_norm5 <- c(x_norm5,-x0,x0)
#     y_norm5 <- c(y_norm5,x0,x0)
    
    #####  (1-alpha)*(x^2/r^2+y^2/r^2)+alpha*max(abs(x/r)+abs(y/r))=1
    x5 <- x2*(1-alpha) + alpha*r*sign(x2)
    y5 <- y2*(1-alpha) + alpha*r*sign(y2)
    x5 <- c(x5, x3[abs(x3)<alpha*r | abs(y3)<alpha*r])
    y5 <- c(y5, y3[abs(x3)<alpha*r | abs(y3)<alpha*r])

    plot(x1,y1,col="red",cex=0.7)
    points(x2,y2,col="green",cex=0.7)
    points(x3,y3,col="cyan",cex=0.7)
    points(x4,y4,col="purple",cex=0.7)
    points(x5,y5,col="blue",cex=0.7)
}

feature_pallete <- function(coef, alpha_list){
    m = ceiling(sqrt(length(alpha_list)))
    n = floor(sqrt(m^2-length(alpha_list)))
    png("feature_pallete.png",wid=100*(m+n),hei=100*(m-n),res=300)
    par(mfrow=c(m-n,m+n))
    for(k in 1:length(alpha_list)){
        alpha <- alpha_list[k]
        weight <- abs(coef[,k])
        weight <- sort(weight/max(weight))
        weight.1 <- log10(1+weight*9)/1
        weight.2 <- log10(1+weight*99)/2
        weight.3 <- log10(1+weight*999)/3
        
        col.edge <- rev(topo.colors(30))[as.numeric(cut(weight, breaks=30))]
        col.edge.1 <- rev(topo.colors(30))[as.numeric(cut(weight.1, breaks=30))]
        col.edge.2 <- rev(topo.colors(30))[as.numeric(cut(weight.2, breaks=30))]
        col.edge.3 <- rev(topo.colors(30))[as.numeric(cut(weight.3, breaks=30))]
        
        plot(weight[weight>0], col=col.edge[weight>0])
        points(1:sum(weight.1>0), weight.1[weight.1>0], col=col.edge.1[weight>0])
        points(1:sum(weight.2>0), weight.2[weight.2>0], col=col.edge.2[weight>0])
        points(1:sum(weight.3>0), weight.3[weight.3>0], col=col.edge.3[weight>0])
    }
    dev.off()
}

# 
# unifs_dis <- function(x, y, alpha, beta1, beta2, npar){
#     p = abs(x-y)
#     if(alpha == 0){
#         result = beta2 %*% p
#     }
#     p = abs(x-y)
#     Rseq = unique(c(sort(p*beta1, decreasing=TRUE),0))
#     R = Rseq(1)
#     M = 0
#     while()
#     
#     rho = [par(npar+1:npar+nfea);0];
#     bk = [bw;1];
#     [Rseq, index] = sort(rho.*bk, 'descend');
#     R = Rseq(1);
#     M = 0;
#     i = 1;
#     while(lambda^2*R^2 > M)
#         i = i + 1;
#     while(Rseq(i) > R-10^-4)
#         i = i + 1;
#     end
#     R = Rseq(i);
#     M = sum((Rseq(1:i)-R).^2./bk(index(1:i)));
#     end
#     list = index(1:i-1);
#     c_rho = rho(list);
#     c_bk = bk(list);
#     tel = (sum(1./c_bk)-lambda^2);
#     sq = sqrt(sum(c_rho)^2-tel*c_bk'*(c_rho.^2));
#         fval = (sum(c_rho)-sq)/tel;
#     
#     
#     
#     % sum(bw*(rho-(alpha*R)./bw).^2) = lambda^2*(alpha*R)^2
#     % Rseq = alpha*R = (sum(rho)-sqrt(rho'*rho-(bw*rho.^2)*(sum(1./bw)-lambda^2)))/(sum(1./bw)-lambda^2);
#     
#     
    
    

