



pspline.setting <- function(x,knots=select.knots(x,35),
                            p=3,m=2,weight=NULL,type="full",
                            knots.option="equally-spaced"){
  
  # x: the marginal data points
  # knots: the list of interior knots or the numbers of interior knots
  # p: degrees for B-splines, with defaults values 3
  # m: orders of difference penalty, with default values 2
  # knots.option: type of knots placement, with default values "equally-spaced"
  
  #require(splines)
  #require(Matrix)
  
  ### design matrix 
  K = length(knots)-2*p-1
  B = create.fourier.basis(rangeval = range(mod$argvals.new), nbasis = 5)
  
  bs = "ps"
  
  if(knots.option == "quantile"){
    bs = "bs"
  }
  
  s.object = s(x=x, bs=bs, k=K+p,m=c(p-1,2), sp=NULL)
  object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
  P =  object$S[[1]]
  if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
  
  if(is.null(weight)) weight <- rep(1,length(x))
  
  if(type=="full"){
    
    Sig = crossprod(matrix.multiply(B,weight,option=2),B)
    eSig = eigen(Sig)
    V = eSig$vectors
    E = eSig$values
    if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
      #cat("A small identity matrix is added!\n");
      E <- E + 0.000001;
      
    }
    Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
    
    tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
    Esig = eigen(tUPU,symmetric=TRUE)
    U = Esig$vectors
    s = Esig$values
    s[(K+p-m+1):(K+p)]=0
    A = B%*%(Sigi_sqrt%*%U)
  }
  
  if(type=="simple"){
    
    A = NULL
    s = NULL
    Sigi_sqrt = NULL
    U = NULL
  }
  List = list(
    "A" = A,
    "B" = B,
    "s" = s,
    "Sigi.sqrt" = Sigi_sqrt,
    "U" = U,
    "P" = P)
  
  return(List)
}





###################################################################
###################################################################




pspline <- function(data,argvals.new= NULL,knots=35,
                    # knots.option="equally-spaced",
                    p=3,m=2,lambda=NULL,
                    search.length = 100,
                    lower=-20,upper=20){
  
  
  #require(splines)  
  #source("pspline.setting.weight.R")
  
  ## read in data
  check.data(data)
  Y <- data$y
  t <- data$argvals
  subj <- data$subj
  tnew <- argvals.new
  if(is.null(tnew)){tnew <- t}
  ###
  J <- length(Y)
  if(is.null(t))  t <- (1:J)/J-1/2/J ## if NULL, assume equally spaced
  subj_unique <- unique(subj)
  n <- length(subj_unique)
  WY <- Y
  N <- rep(NA,n)
  W <- 0*Y
  for(i in 1:n){
    seq <- (1:J)[subj==subj_unique[i]]
    N[i] <- length(seq)
    WY[seq] <- Y[seq]/N[i]
    W[seq] <- 1/N[i]
  }
  
  p.p <- p
  m.p <- m
  
  ######### precalculation for smoothing ############
  knots.option="equally-spaced"
  knots <- construct.knots(t,knots,knots.option,p)
  
  
  List <- pspline.setting(t,knots=knots,p.p,m.p,weight=W,knots.option=knots.option)
  AS <- as.matrix(List$A)
  s <- List$s
  Sigi.sqrt <- List$Sigi.sqrt
  U <- List$U
  B <- List$B
  Bnew <- create.fourier.basis(rangeval = range(mod$argvals.new), nbasis = 5)
  
  AStY <- as.vector(t(AS)%*%Y)
  AStWY <- as.vector(t(AS)%*%WY)
  AStAS <- as.matrix(t(AS)%*%AS)
  AStAS_eigen <- eigen(AStAS)
  AStAS_eigen$values[AStAS_eigen$values<=0] <- 0
  AStAS_half <- AStAS_eigen$vectors%*%diag(sqrt(AStAS_eigen$values))%*%t(AStAS_eigen$vectors)
  
  ##calculate Term1
  AStAS_N <- matrix(NA,n,length(s)^2)
  ASitYi <- matrix(NA,n,length(s))
  Term1 <- matrix(0,length(s),length(s))
  Term0 <- rep(0,length(s))
  for(i in 1:n){
    seq <- (1:J)[subj==subj_unique[i]]
    ASi <- matrix(AS[seq,],N[i])
    ASitASi <- t(ASi)%*%ASi
    AStAS_N[i,] <- as.vector(ASitASi)
    temp <- as.vector(t(ASi)%*%Y[seq]/N[i])
    ASitYi[i,] <- temp
    Term1 <- Term1 + diag(temp)%*%as.matrix(ASitASi%*%diag(AStWY))
    Term0 <- Term0 + temp^2*N[i] 
  }
  
  
  pspline_gcv <- function(x){
    lambda <- exp(x)
    lambda_s <- 1/(1 + lambda*s)
    gcv <- sum((AStAS_half%*%(lambda_s*AStWY))^2)- 2*sum(lambda_s*(AStY*AStWY))
    gcv <- gcv + 2*sum(lambda_s*Term0)
    gcv <- gcv - 4*sum(lambda_s*(Term1%*%lambda_s))
    for(i in 1:n){
      gcv <- gcv + 2*sum((sqrt(lambda_s)*(matrix(AStAS_N[i,],length(s),length(s))%*%(lambda_s*AStWY)))^2)/N[i]
    }
    return(gcv) 
  }
  
  if(is.null(lambda)){
    Lambda <- seq(lower,upper,length=search.length)
    Length <- length(Lambda)
    Gcv <- rep(0,Length)
    for(i in 1:Length) 
      Gcv[i] <- pspline_gcv(Lambda[i])
    i0 <- which.min(Gcv)
    lambda <- exp(Lambda[i0])
  }
  
  theta <- (Sigi.sqrt %*% U)%*%(1/(1+lambda*s)*AStWY)
  
  res <- list(fitted.values = as.vector(B%*%theta), B = B,theta=theta,s = s,
              knots=knots,p=p,m=m,
              lambda=lambda,argvals.new = tnew, mu.new = as.vector(Bnew%*%theta))
  class(res) <- "pspline.face"
  return(res)
  
}


### About:
## The function is to estimate the mean and covariance function
## from a cluster of functions/longitudinal observations. The output
## will also predict each curve.

## Function arguments:
### "data" is a data frame with three arguments:
# (1) "argvals": observation times;
# (2) "subj": subject indices;
# (3) "y": values of observations;
# Note that: we only handle complete data, so missing values are not allowed at this moment

## "newdata" is of the same strucutre as "data". 
## To predict at new values of "argvals", make the corresponding "y" NA

### "center" means if we want to compute population mean

### "argvals.new" if we want the estimated covariance function at "argvals.new"; if NULL,
### then 100 equidistant points in the range of "argvals" in "data"

### "knots" is the number of interior knots for B-spline basis functions to be used; 

face.sparse.inner <- function(data, newdata = NULL, W = NULL,
                              center=TRUE,argvals.new=NULL,
                              knots=7, knots.option="equally-spaced",
                              p=3,m=2,lambda=NULL,lambda_mean=NULL,
                              search.length=14,
                              lower=-3,upper=10, 
                              calculate.scores=FALSE,pve=0.99){
  
  #########################
  ####step 0: read in data
  #########################
  face:::check.data(data)
  if(!is.null(newdata)){ face:::check.data(newdata,type="predict")}
  
  y <- data$y
  t <- data$argvals
  subj <- data$subj
  tnew <- argvals.new
  if(is.null(tnew)) tnew <- seq(min(t),max(t),length=100)
  
  fit_mean <- NULL
  
  knots.initial <- knots
  #########################
  ####step 1: demean
  #########################
  r <- y
  mu.new <- rep(0,length(tnew))
  if(center){
    fit_mean <- face:::pspline(data,argvals.new=tnew,knots=knots.initial,lambda=lambda_mean)
    mu.new <- fit_mean$mu.new
    r <- y - fit_mean$fitted.values 
  }
  #########################
  ####step 2:raw estimates
  #########################
  indW <- F # whether identity W
  if(is.null(W)) indW <- T
  
  raw <- face:::raw.construct(data.frame("argvals" = t, "subj" = subj, "y" = as.vector(r)))
  C <- raw$C
  st <- raw$st
  N <- raw$st
  N2 <- raw$N2
  if(indW) W <- raw$W  
  n0 <- raw$n0
  
  delta <- Matrix((st[,1]==st[,2]) * 1) # sparse
  
  #########################
  ####step 3: smooth
  #########################
  knots <- face:::construct.knots(t,knots,knots.option,p)
  
  List <- face:::pspline.setting(st[,1],knots=knots,p,m,type="simple",knots.option=knots.option)
  B1 <- List$B
  B1 <- Matrix(B1)
  DtD <- List$P
  
  B2 = spline.des(knots=knots, x=st[,2], ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  c = dim(B1)[2]
  c2 = c*(c+1)/2
  B = Matrix(t(KhatriRao(Matrix(t(B2)),Matrix(t(B1)))))
  G = Matrix(duplication.matrix(c))
  
  BtWB = matrix(0,nrow=c^2,ncol=c^2)
  Wdelta = c()
  WC = c()
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    B3 = Matrix(matrix(B[seq,],nrow=length(seq)))
    W3 = W[[i]] # don't form a large W
    BtWB = BtWB + crossprod(B3, W3%*%B3)
    Wdelta <- c(Wdelta,as.matrix(W3 %*% delta[seq]))
    WC <- c(WC,as.matrix(W3 %*% C[seq]))
  }
  
  GtBtWBG = crossprod(G,BtWB%*%G)
  
  BG = B%*%G # sparse
  detWde <- crossprod(delta,Wdelta) # detWde = sum(delta)
  GtBtWdelta <- crossprod(BG,Wdelta)
  XtWX <- rbind(cbind(GtBtWBG,GtBtWdelta), cbind(t(GtBtWdelta),detWde))
  
  eSig = eigen(XtWX,symmetric=TRUE)
  V = eSig$vectors
  E = eSig$values
  E = E + 0.000001*max(E)
  Sigi_sqrt = face:::matrix.multiply(V,1/sqrt(E))%*%t(V)
  
  P = crossprod(G,Matrix(suppressMessages(kronecker(diag(c),DtD))))%*%G
  
  Q = bdiag(P,0)
  tUQU = crossprod(Sigi_sqrt,(Q%*%Sigi_sqrt))
  Esig = eigen(tUQU,symmetric=TRUE)
  
  U = Esig$vectors
  s = Esig$values
  A0 <- Sigi_sqrt%*%U
  X <- cbind(BG,delta)
  A = as.matrix(X%*%A0) # F=XA dense
  
  AtA = crossprod(A) # diff
  f = crossprod(A,C) # diff
  ftilde = crossprod(A,WC) # diff
  
  c2 <- c2 + 1
  g <- rep(0, c2)
  G1 <- matrix(0,c2,c2)
  mat_list <- list()
  
  for(i in 1:n0){
    seq = (sum(N2[1:i])-N2[i]+1):(sum(N2[1:i]))
    Ai = matrix(A[seq,],nrow=length(seq))
    AitAi = crossprod(Ai) #t(Ai)%*%Ai
    Wi = W[[i]]
    
    fi = crossprod(Ai,C[seq]) # t(Fi)Ci
    Ji = crossprod(Ai,Wi%*%C[seq])
    Li = crossprod(Ai,Wi%*%Ai)
    g = g + Ji*fi
    G1 = G1 + AitAi*(Ji%*%t(ftilde))
    
    LList <- list()
    LList[[1]] = AitAi
    LList[[2]] = Li
    mat_list[[i]] = LList
    
  }
  
  
  Lambda <- seq(lower,upper,length=search.length)
  Gcv <- 0*Lambda
  gcv <- function(x){
    lambda <- exp(x)
    d <- 1/(1+lambda*s)
    ftilde_d <- ftilde*d
    cv0 <- -2*sum(ftilde_d*f)
    cv1 <-  sum(ftilde_d*(AtA%*%ftilde_d))
    cv2 <-  2*sum(d*g)
    cv3 <-  -4*sum(d*(G1%*%d))
    cv4 <- sum(unlist(sapply(mat_list,function(x){
      a <- x[[1]]%*%ftilde_d
      b <- x[[2]]%*%ftilde_d
      2*sum(a*b*d)
    })))
    cv <- cv0 + cv1 + cv2 + cv3 + cv4
    return(cv)
  }
  if(is.null(lambda)){
    Lambda <- seq(lower,upper,length=search.length)
    Length <- length(Lambda)
    Gcv <- rep(0,Length)
    for(i in 1:Length) 
      Gcv[i] <- gcv(Lambda[i])
    i0 <- which.min(Gcv)
    lambda <- exp(Lambda[i0])
  }
  
  alpha <- face:::matrix.multiply(A0,1/(1+lambda*s))%*%ftilde
  Theta <- G %*% alpha[1:c2-1]
  Theta <- matrix(Theta,c,c)         # parameter estimated (sym)
  sigma2 <- alpha[c2]
  if(sigma2 <= 0.000001) {                                               
    warning("error variance cannot be non-positive, reset to 1e-6!")    
    sigma2 <- 0.000001                                                  
  }
  
  Eigen <- eigen(Theta,symmetric=TRUE)
  Eigen$values[Eigen$values<0] <- 0
  npc <- sum(Eigen$values>0) #which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1]
  if(npc >1){
    Theta <- face:::matrix.multiply(Eigen$vectors[,1:npc],Eigen$values[1:npc])%*%t(Eigen$vectors[,1:npc])
    Theta_half <- face:::matrix.multiply(Eigen$vectors[,1:npc],sqrt(Eigen$values[1:npc]))
  }
  if(npc==1){
    Theta <- Eigen$values[1]*suppressMessages(kronecker(Eigen$vectors[,1],t(Eigen$vectors[,1])))
    Theta_half <- sqrt(Eigen$values[1])*Eigen$vectors[,1]
  }
  Eigen <- eigen(Theta,symmetric=TRUE)
  
  #########################
  ####step 4: calculate estimated covariance function
  #########################
  Bnew = spline.des(knots=knots, x=tnew, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Gmat <- crossprod(Bnew) / nrow(Bnew)
  eig_G <- eigen(Gmat, symmetric = T)
  G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
  G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
  
  Chat.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) 
  Chat.diag.new = as.vector(diag(Chat.new))  
  Cor.new = diag(1/sqrt(Chat.diag.new))%*%Chat.new%*%diag(1/sqrt(Chat.diag.new))
  
  Eigen.new = eigen(as.matrix(G_half%*%Matrix(Theta)%*%G_half),symmetric=TRUE)
  # Eigen.new = eigen(Chat.new,symmetric=TRUE)
  npc = which.max(cumsum(Eigen$values)/sum(Eigen$values)>pve)[1] #determine number of PCs
  
  eigenfunctions = matrix(Bnew%*%G_invhalf%*%Eigen.new$vectors[,1:min(npc,length(tnew))],ncol=min(npc,length(tnew)))
  eigenvalues = Eigen.new$values[1:min(npc,length(tnew))]
  # eigenfunctions = eigenfunctions*sqrt(length(tnew))/sqrt(max(tnew)-min(tnew))
  # eigenvalues = eigenvalues/length(tnew)*(max(tnew)-min(tnew))
  
  
  #########################
  ####step 5: calculate variance
  #########################
  var.error.hat <- rep(sigma2,length(t))
  var.error.new <- rep(sigma2,length(tnew))
  
  
  
  Chat.raw.new = as.matrix(tcrossprod(Bnew%*%Matrix(Theta),Bnew)) + diag(var.error.new) 
  Chat.raw.diag.new = as.vector(diag(Chat.raw.new)) 
  Cor.raw.new = diag(1/sqrt(Chat.raw.diag.new))%*%Chat.raw.new%*%diag(1/sqrt(Chat.raw.diag.new))
  #########################
  ####step 6: prediction
  #########################
  if(is.null(newdata) && calculate.scores==T){
    newdata = data
  }
  if(!is.null(newdata)){
    
    mu.pred <- rep(0,length(newdata$argvals))
    var.error.pred <- rep(sigma2,length(newdata$argvals))
    if(center){
      mu.pred <- face:::predict.pspline.face(fit_mean,newdata$argvals)
    }
    
    subj.pred = newdata$subj
    subj_unique.pred = unique(subj.pred)
    y.pred = newdata$y
    Chat.diag.pred = 0*y.pred
    se.pred = 0*y.pred
    
    scores = list(subj=subj_unique.pred,
                  scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                  u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
    )
    
    Ulist = list()
    for(i in 1:length(subj_unique.pred)){
      sel.pred = which(subj.pred==subj_unique.pred[i])
      lengthi = length(sel.pred)
      
      pred.points <- newdata$argvals[sel.pred]
      mu.predi <- mu.pred[sel.pred]
      var.error.predi <- var.error.pred[sel.pred]
      
      y.predi = y.pred[sel.pred] - mu.predi
      sel.pred.obs = which(!is.na(y.predi))
      obs.points <- pred.points[sel.pred.obs]
      if(!is.null(obs.points)){
        var <- mean(var.error.predi[sel.pred.obs])
        if(var==0&length(sel.pred.obs) < npc)
          stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
        B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
        Chati = tcrossprod(B3i%*%Theta,B3i)
        Chat.diag.pred[sel.pred] = diag(Chati)
        if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
        if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
        Vi.inv = as.matrix(solve(Chati + Ri))
        Vi.pred = tcrossprod(B3i.pred%*%Theta,B3i.pred)
        Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
        ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
        scores$u[i,] = as.vector(ui)
        y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
        temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
        if(length(sel.pred.obs) >1){
          se.pred[sel.pred] = sqrt(diag(Vi.pred - temp%*%Vi.inv%*%t(temp)))
        }
        if(length(sel.pred.obs) ==1){
          se.pred[sel.pred] = sqrt(Vi.pred[1,1] - Vi.inv[1,1]*temp%*%t(temp))
        }
        
        Ulist[[i]] = ui
        
        ## predict scores
        if(calculate.scores==TRUE){ 
          temp = matrix(t(eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(eigenfunctions[,1]^2)
          temp = as.matrix(temp)
          scores$scores[i,1:npc] = temp[,1]
        }
      }
    }
  }## if(is.null(newdata))
  if(is.null(newdata)){
    y.pred=NULL
    mu.pred = NULL
    var.error.pred = NULL
    Chat.diag.pred = NULL
    se.pred = NULL
    scores=NULL
    
  }
  
  
  
  res <- list(newdata=newdata, W = W, y.pred = y.pred, Theta=Theta,argvals.new=tnew, 
              mu.new = mu.new, Chat.new=Chat.new, var.error.new = var.error.new,
              Cor.new = Cor.new, eigenfunctions = eigenfunctions, eigenvalues = eigenvalues,
              Cor.raw.new = Cor.raw.new, Chat.raw.diag.new = Chat.raw.diag.new,
              rand_eff = scores, calculate.scores=calculate.scores,
              mu.hat = fit_mean$fitted.values,var.error.hat = var.error.hat,
              mu.pred = mu.pred, var.error.pred = var.error.pred, Chat.diag.pred = Chat.diag.pred,
              se.pred = se.pred,
              fit_mean = fit_mean, lambda_mean=fit_mean$lambda,
              lambda=lambda,Gcv=Gcv,Lambda=Lambda,knots=knots,knots.option=knots.option,s=s,npc=npc, p = p, m=m,
              center=center,pve=pve,sigma2=sigma2, r = r, DtD = DtD,
              U = Eigen.new$vectors[,1:npc],G_invhalf = G_invhalf, Ulist = Ulist, Bnew = Bnew)
  
  class(res) <- "face.sparse"
  return(res)
}




#####################################################################
face.sparse <- function(data, newdata = NULL,
                        center=TRUE,argvals.new=NULL,
                        knots=7, 
                        # knots.option="equally-spaced",
                        p=3,m=2,lambda=NULL,lambda_mean=NULL,
                        search.length=14,
                        lower=-3,upper=10, lower2 = -3, upper2=5,
                        calculate.scores=FALSE,pve=0.998,two_step=FALSE){
  
  if(!two_step){
    fit <- face.sparse.inner(data = data, newdata=newdata, W = NULL,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower,upper=upper,search.length=search.length,
                             calculate.scores=calculate.scores, pve = pve)
    
  }else{
    fit <- face.sparse.inner(data = data, newdata=NULL, W = NULL,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower,upper=upper,search.length=search.length,
                             pve = pve)
    
    W <- face:::weight.construct(fit,data)
    
    
    fit <- face.sparse.inner(data = data, newdata=newdata, W = W,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower2,upper=upper2,search.length=length(lower2:upper2),
                             calculate.scores=calculate.scores,pve=pve)
  }
  
  
  return(fit)
}








####################### Irregular to Matrix Function #################################

irreg2mat <- function(ydata, binning=FALSE, maxbins=1e3){
  
  ## drop any missings:
  ydata <- ydata[complete.cases(ydata), ]
  
  ## turn into row/column indices for new matrix
  nobs <- length(unique(ydata$.id))
  # make sure newid takes values 1:nobs
  newid <- as.numeric(as.factor(ydata$.id))
  
  ## bin y-index, if necessary
  bins <- sort(unique(ydata$.index))
  if(binning && (length(bins) > maxbins)){
    # linear binning;
    # bin-borders go from just below min to just above max
    # TODO: quantile-based binning?
    binvalues <- seq((1-.001*sign(bins[1]))*bins[1], 
                     (1+.001*sign(bins[length(bins)]))*bins[length(bins)],
                     l=maxbins+1)
    bins <- binvalues
    binvalues <- head(filter(binvalues, c(.5, .5)), -1)
  } else {
    binvalues <- bins
    bins <- c((1-.001*sign(bins[1]))*bins[1], 
              bins[-length(bins)],
              (1+.001*sign(bins[length(bins)]))*bins[length(bins)])
    # take care of edge cases:
    if(bins[1] == 0) bins[1] <- -.001
    if(bins[length(bins)] == 0) bins[length(bins)] <- .001
  }
  newindex <- cut(ydata$.index, breaks=bins, include.lowest=TRUE)
  
  Y <- matrix(NA, nrow=nobs, ncol=nlevels(newindex))
  colnames(Y) <- binvalues
  attr(Y, "index") <- binvalues
  Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
  return(Y)
}




