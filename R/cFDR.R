##' Compute estimated cFDR values from a set of pairs of p values, or adjusted p-values if controls are shared.
##'
##' For a set of p-values or adjusted p values \eqn{P_{j}}{P_i} for the principal phenotype, and a set of p values \eqn{P_{j}}{P_j} for the conditional phenotype, the cFDR for SNP \eqn{x} is estimated as
##' 
##' \deqn{
##'   \widehat{cFDR}(x)=P_{i}(x) \frac{\textrm{(# SNPs with $P_{j} \leq P_{j}(x)$}{\textrm{(# SNPs with $P_{i} \leq P_{i}(x), P_{j} \leq P_{j}(x))}}}
##' }{
##'   est. cFDR(x) = P_i(x) * (# P_j \le P_j(x)) / (# P_i \le P_i(x) AND P_j \le P_j(x))
##' }
##'
##' The procedure can be sped up by specifying a subset of SNPs at which to estimate the cFDR. This is most effective if a boundary of the subset can be written as \eqn{P_{j} = f(P_{i})}{P_j = f(P_i)}, where \eqn{f} is non-increasing. Equivalently, if \eqn{S} is the set of all p-value pairs, the subset \eqn{S'} should satisfy the condition
##' 
##' \deqn{
##'   \forall \left( p_{i}',p_{j}' \in S', p_{i},p_{j} \in S \right) \: p_{i} \leq  p_{i}', p_{j} \leq p_{j}' \implies p_{i},p_{j} \in S'
##' }{
##'   for all (p_i',p_j') in S', (p_i,p_j) in S:   p_i \le p_i', q_j \le p_j   =>    (p_i,p_j) in S'.
##' }
##'
##' An example of such a subset is the set of \eqn{p_{i},p_{j}}{p_i,p_j} with \eqn{p_{i}^2 + p_{j}^2 < \alpha}{p_i^2 + p_j^2 < \alpha}.
##' 
##' The estimated cFDR is a conservative estimate (upper bound) on the quantity 
##' 
##' \deqn{
##'   Pr(H_{0}^{i} | P_{i}<P_{i}(x),P_{j}<P_{j}(x))
##'  }{
##'   Pr(H0 | P_i<P_i(x),P_j<P_j(x))
##'  }
##' 
##' where \eqn{H_{0}^{i}}{H0} is the null hypothesis at SNP \eqn{x} for the principal phenotype (ie, SNP \eqn{x} is not associated with the principal phenotype)
##' 
##' @title cfdr
##' @param P_i vector of p-values for the principal phenotype.
##' @param P_j vector of p-values or adjusted p-values for the conditional phenotype. If controls are shared between GWAS, the p-values must be adjusted using the function \code{\link{exp_quantile}}.
##' @param sub set of indices; only compute cFDR at this subset of p-value pairs. Typically, only SNPs with low values of P_i or P_j are of interest; it is possible to increase the speed of the procedure by only calculating cFDR at such SNPs.
##' @export
##' @author James Liley
##' @examples
##' require(mnormt)
##' P = 2*pnorm(-abs(rmnorm(10000,varcov=diag(c(0.8,0.8))+0.2)))
##' P_i = P[,1]; P_j=P[,2]
##' W = which((qnorm(P_i/2)^2) + (qnorm(P_j/2)^2) > 4)
##'
##' # Using function
##' C = cfdr(P_i,P_j,sub=W); C[W[1]]
##' 
##' # Explicit computation
##' P_i[W[1]]/ (length(which (P_i <= P_i[1] & P_j <= P_j[1])) / length(which (P_j <= P_j[1])) )
##'

cfdr = function(P_i,P_j,sub=1:length(P_i)) {
  
  # various error handlers
  if (!is.numeric(P_i) | !is.numeric(P_j) || (max(P_i,na.rm=TRUE)>1) || (max(P_j,na.rm=TRUE)>1) || (min(P_i,na.rm=TRUE)<0) || (min(P_j, na.rm=TRUE)<0)) stop("Arguments P_i and P_j must both be vectors of p values or adjusted p values")
  if (!(all(sub %in% 1:length(P_i)) | all(sub %in% names(P_i)))) stop("Parameter sub must be a valid list of indices of P_i")
  
  # check sub, add warning
  
  
  C=rep(NA,length(P_i)); names(C)=names(P_i)
  wx=rep(0,length(P_i)); wx[sub]=1;
  
  ww=which(!is.na(P_i+P_j))
  P_i=P_i[ww]; P_j=P_j[ww]
  wx=wx[ww]; subx=which(wx==1)
  
  # For each P_i,P_j pair in 'sub', count the number of SNPs with smaller P_i and P_j values
  R_i=rank(P_i); R_jx=rank(P_j)
  o=order(R_i)
  wx=wx[o]; suby=which(wx==1)
  R_j=R_jx[o]
  N=length(P_i)
  maxl=rep(1,N)
  out=rep(1,N)
  for (i in 1:length(suby)) {
    out[suby[i]]=maxl[R_j[suby[i]]]
    maxl[R_j[suby[i]]:N]=1+maxl[R_j[suby[i]]:N]
  }
  
  # Divide by number of SNPs with smaller P_j values, then divide P_j values by the result to get cFDR.
  Cx=rep(1,length(ww))
  Cx[subx]=P_i[subx]/(out[order(o)][subx]/R_jx[subx])
  
  C[ww]=Cx
  C
}





##' Transform principal p-values for estimating cFDR in a shared control design. Computes the 'expected quantile' of a p value.
##'
##' Computes the probability that a p value at some SNP for a phenotype \eqn{i} is less than some cutoff \eqn{p_{i}}{p_i} given that the p value at that SNP for a second phenotype \eqn{j} is less than a cutoff \eqn{p_{j}}{p_j}, under the null hypothesis that the SNP is not associated with phenotype \eqn{i}.
##'
##' If the GWAS for phenotypes \eqn{i} and \eqn{j} share a number of controls, effect sizes for the two phenotypes will be correlated, even at null SNPs. This leads to dependence of the expected quantile both on the number of shared controls and study sizes (through paramter \code{\link{rho}}), and the distribution of Z values for phenotype \eqn{j} across all SNPs, not all of which will necessarily be null for phenotype \eqn{j}.
##'
##' The distribution of Z scores for phenotype \eqn{j} can be specified in two ways, depending on the parameter method. In the first case, we assume that the distribution of 'true' Z scores (Z scores which we would observe if the observed allele frequencies exactly matched the population allele frequencies) has a mixture distribution of 0 with probability \code{\link{pi0}}, and N(0,\code{\link{sigma}}^2) with probability 1-\code{\link{pi0}}. These parameters can be estimated from an observed distribution of Z scores using the function \code{\link{fit.em}} .
##'
##' In the second case, the distribution of true Z scores can be specified as an array of Z scores and corresponding probabilities (or PDF heights).
##'
##' @title exp_quantile
##' @param p Observed p value or set of p values for principal phenotype
##' @param pc Observed p value or set of p values for conditional phenotype. Must be of the same length as \code{\link{p}}.
##' @param rho Correlation between z values induced by shared controls
##' @param  method Set to 1 to assume Z scores for the conditional phenotype follow a mixed distribution characterised by \code{\link{pi0}}, \code{\link{sigma}}; 2 to specify assume they follow empirical distribution characterised by \code{\link{dist}}. 
##' @param pi0 Relevant only if \code{\link{method}} = 1. Parameter for distribution of Z scores for the conditional phenotype (proportion of SNPs null for conditional phenotype). 
##' @param sigma Relevant only if \code{\link{method}} = 1. Parameter for distribution of Z scores for the conditional phenotype (variance of true effect sizes for SNPs non-null for conditional phenotype).
##' @param dist Relevant only if \code{\link{method}} = 2. Matrix with two rows; first is a list of potential values for Z scores for conditional phenotype, second is a list of probabilities associated with those values.
##' @return A list of transformed p values (expected quantiles) of the same length as \code{\link{p}}.
##' @export
##' @author James Liley
##' @examples
##' p <- 0.2; pc <- 0.1;
##' rho <- 0.25
##' sigma <- 5; pi0 <- 0.95
##' exp_quantile(p,pc,rho,method=1,pi0=pi0,sigma=sigma)
##' # Generally the expected quantile is close to p. In this case, the SNP is 'probably' null for phenotype j, so the MAF in controls is somewhat aberrant from expected (as pc<0.1). This means that the expected quantile is higher than the p value.
##'
##' p<-0.2; pc <-1e-8;
##' exp_quantile(p,pc,rho,method=1,pi0=pi0,sigma=sigma)
##' # If a low value of pc is observed, it is 'likely' to have arisen from a non-null SNP for phenotype j. This p value is quite reasonable for a non-null SNP, given the distribution of z values for such SNPs (sigma=5), so the aberration from expected MAF amongst controls is minimal.
##'
##' pi0=1
##' exp_quantile(p,pc,rho,method=1,pi0=pi0,sigma=sigma)
##' # If, on the other hand, pi0=1 (no non-null SNPs for phenotype j), and we happen to observe a p value of 1e-8 for the conditional phenotype, we conclude that MAF in controls is likely to be very aberrant from expected, so the expected quantile is markedly higher than the p value.
##'
##' rho <- 0
##' exp_quantile(p,pc,rho,method=1,pi0=pi0,sigma=sigma)
##' # If rho=0 (no shared controls) z values for the two phenotypes are independent, and the expected quantile equals parameter p.
##'
##' p3 = rbind((-500:500)/20, dnorm((-500:500)/20,sd=3))
##' exp_quantile(5e-4,5e-3,0.3,method=2,pi0=0.9,dist=p3)
##' exp_quantile(5e-4,5e-3,0.3,method=1,pi0=0.9,sigma=3)
##' # Demonstrates specification of an empirical distribution for Z scores for the conditional phenotype.
##' 
##' 

exp_quantile = function(p,pc,rho,method=1,pi0=0.5,sigma=1,dist = rbind((1:10)/10,dnorm((1:10)/10,sd=3))) {
  
  require(pbivnorm)
  
  if (rho==0) return(p) else {
    
    z = -qnorm(p/2)
    zc = -qnorm(pc/2);
    
    if (method==1) {
      num = (pi0*pbivnorm(-z,-zc,rho=rho)) + (pi0*pbivnorm(-z,-zc,rho=-rho)) + ((1-pi0)*pbivnorm(-z,-zc/sqrt(1+(sigma^2)),rho=rho/sqrt(1+(sigma^2)))) +  ((1-pi0)*pbivnorm(-z,-zc/sqrt(1+(sigma^2)),rho=-rho/sqrt(1+(sigma^2))))
      denom = (pi0*pnorm(-zc)) + (1-pi0)*pnorm(-zc/sqrt(1+(sigma^2)))
      return(num/denom)
    }
    if (method==2) {
      dist[,2]=dist[,2]/(trap(dist[,1],dist[,2]))
      num = matrix(0,length(z), dim(dist)[2]); denom=num;
      for (i in 1:dim(dist)[2]) {
        num[,i] = pbivnorm(-z,-zc-dist[1,i],rho=rho) + pbivnorm(-z,-zc+dist[1,i],rho=rho) + pbivnorm(-z,-zc-dist[1,i],rho=-rho) + pbivnorm(-z,-zc+dist[1,i],rho=-rho)
        denom[,i] = pnorm(-zc-dist[1,i]) + pnorm(-zc+dist[1,i])
      }
      num_null = 2*((pbivnorm(-z,-zc,rho=rho) + pbivnorm(-z,-zc,rho=-rho))); denom_null = 2*pnorm(-zc); 
      num_nonnull=num_null; denom_nonnull=denom_null; # initialise
      for (i in 1:length(z)) {
        num_nonnull[i] = trap(dist[1,],num[i,]*dist[2,]); 
        denom_nonnull[i] = trap(dist[1,],denom[i,]*dist[2,])
      }
      return(((pi0*num_null) + ((1-pi0)*(num_nonnull)))/((pi0*denom_null) + ((1-pi0)*(denom_nonnull))))
    }
  }
}

# Trapezoid rule
trap = function(x,y) {
  0.5*((c(x,0)-c(0,x))[2:(length(x))] %*% (c(y,0)+c(0,y))[2:(length(y))])
}




##' Evaluate the conditional expected quantile of a p value in the shared control design conditional on a p value for a second phenotype
##'
##' Computes the probability that a p value for a phenotype \eqn{i} is less than some cutoff \code{\link{p}} given that the p value for a second phenotype \eqn{j} is equal to a value \code{\link{pc}}, under the null hypothesis for phenotype \eqn{i}.
##'
##' If the GWAS for phenotypes \eqn{i} and \eqn{j} share a number of controls, effect sizes for the two phenotypes will be correlated, even at null SNPs. This leads to dependence of the expected quantile both on the number of shared controls and study sizes (through paramter \code{\link{rho}}), and the distribution of Z scores for phenotype \eqn{j} across all SNPs, not all of which will necessarily be null for phenotype \eqn{j}.
##'
##' The 'true' Z scores (Z scores which we would observe if the observed allele frequencies exactly matched the population allele frequencies) are assumed to follow a mixture distribution of 0 with probability \code{\link{pi0}}, and N(0,\code{\link{sigma}}^2) with probability 1-\code{\link{pi0}}. These parameters can be estimated from an observed distribution of Z scores using the function \code{\link{fit.em}} .
##'
##' The value of the function is the first derivative of the function \code{\link{exp_quantile}} with respect to \eqn{P_[j]}{P_j}.
##'
##' @title exp_quantile_point
##' @param p Observed p value or set of p values for principal phenotype
##' @param pc Observed p value or set of p values for conditional phenotype. Must be of the same length as \code{\link{p}}.
##' @param rho Correlation between Z scores induced by shared controls. Output from \code{\link{cor_shared}}.
##' @param pi0 Parameter for distribution of Z scores for the conditional phenotype (proportion of SNPs null for conditional phenotype).
##' @param sigma Parameter for distribution of Z scores for the conditional phenotype (variance of true effect sizes for SNPs non-null for conditional phenotype).
##' @return A list of probabilities of the same length as \code{\link{p}}.
##' @export
##' @author James Liley
##' @examples
##' p <- 0.2; pc <- 0.1;
##' rho <- 0.25
##' sigma <- 5; pi0 <- 0.95
##' exp_quantile_point(p,pc,rho,pi0=pi0,sigma=sigma)
##' # Generally the expected quantile is close to p. 

exp_quantile_point = function(p,pc,rho,pi0=1,sigma=1) {
  
  require(pbivnorm)
  
  z = -qnorm(p/2)
  zc = -qnorm(pc/2);
  
  num = (dnorm(-zc)*pi0*
           (pnorm(-z,mean = -rho*zc,sd = sqrt(1-(rho^2)))+pnorm(-z,mean = rho*zc,sd = sqrt(1-(rho^2))))) + 
    (dnorm(-zc,sd = sqrt(1+(sigma^2)))*(1-pi0)*
       (pnorm(-z,mean = -rho*zc/sqrt(1+(sigma^2)),sd = sqrt(1-(rho^2)+(sigma^2))/sqrt(1+(sigma^2))) +
          pnorm(-z,mean = rho*zc/sqrt(1+(sigma^2)),sd = sqrt(1-(rho^2)+(sigma^2))/sqrt(1+(sigma^2)))))
  denom = ((pi0*dnorm(-zc)) + (1-pi0)*dnorm(-zc,sd = sqrt(1+(sigma^2))))
  
  return(num/denom)
}




##' Compute correlation between z values arising from shared controls
##'
##' Uses a formula from Zaykin et al (2010)
##'
##' @title cor_shared
##' @param n0i Number of controls in GWAS for principal phenotype
##' @param n0j Number of controls in GWAS for conditional phenotype
##' @param n1i Number of cases in GWAS for principal phenotype
##' @param n1j Number of cases in GWAS for conditional phenotype
##' @param overlap Number of controls shared between studies
##' @return Asymptotic correlation between z scores.
##' @export
##' @author James Liley
##' @examples
##'
##' n0i = 1500; n0j = 1200; n1i = 800; n1j = 700;
##' overlap = 1000;
##' cor_shared(n0i,n0j,n1i,n1j,overlap)
##'
##'
##' # Simulation
##' n_snps = 10000; sim_maf = runif(n_snps,min=0.1,max=0.5) # simulated MAF
##' sim_con_i = matrix(runif(n_snps*n0i),n_snps,n0i)<sim_maf # simulated genotypes for controls for phenotype i
##' sim_con_j = matrix(runif(n_snps*n0j),n_snps,n0j)<sim_maf; sim_con_j[,1:overlap]=sim_con_i[,1:overlap] # simulated genotypes for controls for phenotype j, with shared controls
##' sim_case_i = matrix(runif(n_snps*n1i),n_snps,n1i)<sim_maf; sim_case_j = matrix(runif(n_snps*n1j),n_snps,n1j)<sim_maf # simulated genotypes for cases
##' om0i=rowMeans(sim_con_i); om0j=rowMeans(sim_con_j); om1i = rowMeans(sim_case_i); om1j = rowMeans(sim_case_j) # observed minor allele frequencies
##' or_i = om0i*(1-om1i)/(om1i*(1-om0i)); or_j = om0j*(1-om1j)/(om1j*(1-om0j)) # Odds ratios
##' se_i = sqrt( (1/(2*om0i*n0i)) + (1/(2*(1-om0i)*n0i)) + (1/(2*om1i*n1i)) + (1/(2*(1-om1i)*n1i)) )
##' se_j = sqrt( (1/(2*om0j*n0j)) + (1/(2*(1-om0j)*n0j)) + (1/(2*om1j*n1j)) + (1/(2*(1-om1j)*n1j)) ) # Standard errors
##' effect_i = log(or_i)/se_i; effect_j = log(or_j)/se_j # Effect sizes
##' 
##' cor(effect_i,effect_j)
##' cor_shared(n0i,n0j,n1i,n1j,overlap)
##' 
##' cor(abs(effect_i),abs(effect_j))
##' cor_shared(n0i,n0j,n1i,n1j,overlap)^2

cor_shared = function(n0i, n0j, overlap, n1i, n1j) {
  
  N01 = n0i-overlap; N02=n0j-overlap
  sqrt(1/( (1+(N02/overlap))*(1+(N01/overlap))*(1+(N01/n1i)+(overlap/n1i))*(1+(N02/n1j)+(overlap/n1j)) ));
  
}










##' Compute an upper bound on the false discovery rate amongst SNPs with cFDR less than some cutoff \eqn{\alpha} .
##' 
##' Bound is based on estimating the area of the region of the unit square containing all potential p-value pairs \eqn{p_{i},p_{j}}{p_i,p_j} such that \eqn{\widehat{cFDR}(p_{i},p_{j}) \leq \alpha}{est. cFDR(p_i,p_j) \le \alpha}. It is typically conservative.
##' 
##' Computation requires parametrisation of the joint distribution of Z scores for the conditional phenotype at SNPs null for the principal phenotype. This is assumed to be mixture-Gaussian, consisting of N(0,I_2) with probability \code{\link{pi0}} and N(0,(1,rho; rho,\code{\link{sigma}}^2)) with probability 1-\code{\link{pi0}}. Values of \code{\link{pi0}} and \code{\link{sigma}} can be obtained from the fitted parameters of the null model usign the function \code{\link{fit.em}}.
##' 
##' The probability is computed using a numerical integral over the (+/+) quadrant and the range and resolution of the integral can be set.
##' @title c2a
##' @param Z n x 2 matrix of Z scores; Z[,1] corresponding to the principal phenotype, Z[,2] to the conditional phenotype
##' @param cutoffs vector of cFDR cutoffs at which to compute overall FDR
##' @param pi0 proportion of SNPs not associated with conditional phenotype
##' @param sigma standard deviation of observed conditional Z scores in SNPs associated with the conditional phenotype
##' @param rho covariance between principal and conditional Z scores arising due to shared controls; output from \code{\link{cor_shared}}.
##' @param xmax compute integral over [0,\code{\link{xmax}}] x [0,\code{\link{ymax}}] as approximation to [0,inf] x [0,inf]
##' @param ymax compute integral over [0,\code{\link{xmax}}] x [0,\code{\link{ymax}}] as approximation to [0,inf] x [0,inf]
##' @param res compute integral at gridpoints with this spacing.
##' @param vector of FDR values corresponding to \code{\link{cutoffs}}
##' @author James Liley
##' @export
##' @return list of FDR values
##' @examples
##' nn=100000
##' 
##' Z=abs(rbind(rmnorm(0.8*nn,varcov=diag(2)), rmnorm(0.15*nn,varcov=rbind(c(1,0),c(0,2^2))), rmnorm(0.05*nn,varcov=rbind(c(3^2,2),c(2,4^2)))));
##' P=2*pnorm(-abs(Z))
##' 
##' X=cfdr(P[,1],P[,2],sub=which(Z[,1]^2 + Z[,2]^2 > 6))
##' Xm=which(X<0.05); Xsub=Xm[order(runif(length(Xm)))[1:100]] # sample of cFDR values 
##' 
##' true_fdr=rep(0,100); for (i in 1:100) true_fdr[i]=length(which(X[1:(0.95*nn)] <= X[Xsub[i]]))/length(which(X<=X[Xsub[i]])) # true FDR values (empirical)
##' fdr=c2a(Z,X[Xsub],pi0=0.95,sigma=2) # estimated FDR using area method
##' 
##' plot(true_fdr,fdr,xlab="True FDR",ylab="Estimated",col="red"); points(true_fdr,X[Xsub],col="blue"); abline(0,1); legend(0.1*max(true_fdr),0.7*max(fdr),c("Area method", "cFDR"),col=c("red", "blue"),pch=c(1,1)) # cFDR underestimates true FDR; area method gives good approximation.
##'
c2a=function(Z,cutoffs,pi0=0.5,sigma=1,rho=0,xmax=12,ymax=12,res=0.01) {
  
  if (!is.numeric(pi0) | !is.numeric(sigma) || (pi0>1) | (pi0<0) | (sigma<0)) stop("Parameters pi0 and sigma must be set; pi0 must be in [0,1] and sigma must be non-negative. ")
  
  gx=seq(0,xmax,res); gy=seq(0,ymax,res); gridx=matrix(gx,length(gx),length(gy)); gridy=t(matrix(gy,length(gy),length(gx))); 
  Zx=cbind(as.numeric(gridx),as.numeric(gridy))
  
  ww=which(Z[,1]^2 + Z[,2]^2 > 3.5^2)
  qn=matrix(0,length(gx),length(gy)) # qn[i,j] is set to the number of elements of Z with Z[,1]>gx[i] and Z[,2]>gx[j]
  e1=round(ecdf(gx)(Z[ww,1])*length(gx)); e2=round(ecdf(gy)(Z[ww,2])*length(gy)); nx=length(gx); ny=length(gy)
  for (ii in 1:length(ww)) {
    qn[1:e1[ii],1:e2[ii]]=qn[1:e1[ii],1:e2[ii]]+1 # block out lower-left rectangle for every point
  }
  pj=exp_quantile(2*pnorm(-Zx[,1]),2*pnorm(-Zx[,2]),rho=rho,pi0=pi0, sigma=sigma) # transformed p value
  qd=round(dim(Z)[1]*(1-ecdf(Z[,2])(Zx[,2]))) # denominator of observed quantile
  cf=pj*(1+qd)/(1+as.vector(qn)); cf[which(Zx[,1]^2 + Zx[,2]^2 <3.6^2)]=1; cf=matrix(cf,length(gx),length(gy))
  
  kk=4*((pi0*dmnorm(Zx,varcov=diag(2)))+ ((1-pi0)*dmnorm(Zx,varcov=rbind(c(1,rho),c(rho,sigma^2))))) # probability distribution of null SNPs
  
  ycut=round(ecdf(cf)(cutoffs)*length(cf)) # number of cf values less than each value of cutoffs
  null_L=cumsum(kk[order(cf)])[ycut]/sum(kk) # integrate pdf of null SNPs over region with cfdr less than cutoffs. Expected proportion of null SNPs lying in region with cfdr < cut.
  
  # xij is defined such that xij[i,j]=probability of null SNP falling in region x>gx[i],y>gy[j]
  xij=t(apply(apply(matrix(kk,length(gx),length(gy))[length(gx):1,length(gy):1], 2, cumsum), 1, cumsum))[length(gx):1,length(gy):1]
  maxm=  cummax(as.numeric(xij)[order(cf)])[ycut]/sum(kk) # maxm[j] is the maximum expected number of null SNPs in a region (x>xx,y>yy) with cfdr(xx,yy)<cutoffs[j]
  
  cutoffs*null_L/maxm # final FDR is bounded above by the ratio of expected number of SNPs in L divided by the expected number of SNPs in M
  
}







##' Fit a specific two Guassian mixture distribution to a set of Z values.
##'
##' Assumes 
##' Z ~ N(0,1) with probability  \code{\link{pi0}}, Z ~ N(0,1 + \code{\link{sigma}}^2) with probability 1-\code{\link{pi0}}
##'
##' We define 'true' Z scores as the Z scores that would be obtained if MAF for both groups exactly matched the corresponding MAFs in the population, or equivalently the expected values of Z scores. If 'true' Z scores are distributed following a 'spike and tail' model of 0 with probability \code{\link{pi0}} and N(0,\code{\link{sigma}}^2) with probability 1-\code{\link{pi0}}, then observed Z scores follow the above distribution.
##' @title fit.em
##' @param Z numeric vector of observed data
##' @param sigma_init initial value for \code{\link{sigma}}
##' @param pi0_init initial value for \code{\link{pi0}}
##' @param tol how small a change lhood prompts continued optimization
##' @param maxit maximum number of iterations
##' @return a list containing fitted \code{\link{pi0}}, fitted \code{\link{sigma}}, and a record of fitted parameters at eachs stage of the E-M algorithm.
##' @export
##' @author Chris Wallace, James Liley
##' @examples
##' sigma <- 2
##' pi0 <- 0.8
##' n <- 10000; n0=round(pi0*n); n1=n-n0
##' Z <- c(rnorm(n0,0,1),rnorm(n1,0,sqrt(1+ (sigma^2))))
##' fit<-fit.em(Z)
##' fit$pi0
##' fit$sigma
##' fit$history
fit.em <- function(Z,pi0_init=0.9,sigma_init=1,tol=1e-4,verbose=TRUE,maxit=1e4) {
  
  ## probabilities of group0, group1 membership
  p <- c(pi0_init,1-pi0_init)
  px <- matrix(p,length(Z),2,byrow=TRUE)
  
  ## parameter vector
  pars <- c(pi0_init,sigma_init)
  

  pars.fail <- function(pars) if ((pars[1]<0) | (pars[1]>1) | (pars[2]<0)) return(TRUE) else return(FALSE)
  
  ## likelihood for first group
  lh1=function(pi) pi*dnorm(Z)
  
  ## likelihood for second group
  lh2=function(pi,sigma) pi*dnorm(Z,sd=sigma)
  
  ## likelihood function to be maximized
  lhood <- function(pars, sumlog=TRUE) {
   e=lh1(pars[1]) + lh2(1-pars[1],pars[2]) 
   e[which(e==0)] <- 1e-64
   e[which(is.infinite(e))] <- 1e64
   if(sumlog) return(-sum(log(e))) else return(-e)
  }
  
  nit <- 1
  df <- 1
  value <- matrix(NA,maxit,3,dimnames=list(NULL,c("pi0","sigma","lhood")))
  value[nit,] <- c(pars, lhood(pars))
  while(df>tol & nit<maxit) {
    
    nit <- nit+1
    
    ## E step
    px[,1]=lh1(pars[1]);
    px[,2]=lh2(1-pars[1],pars[2])

    px <- px/rowSums(px) ## normalise
    px[is.nan(px)] <- 0
    
    ## M step
    pars[1] = mean(px[,1])
    pars[2] <- sqrt( max(1,sum( px[,2] * Z^2 )) / sum(px[,2]) ) # enforce sigma >= 1
    value[nit,] <- c(pars, lhood(pars))
    df <- abs(value[nit,3] - value[nit-1,3])
    
  }
  return(list(pi0=pars[1], sigma=sqrt(pars[2]^2 - 1),
              history=value[1:nit,]))
}



