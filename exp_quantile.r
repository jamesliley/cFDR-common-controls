##' Evaluate the conditional expected quantile of a p value in the shared control design.
##'
##' Computes the probability that a p value for a phenotype i is less than some cutoff p_i given that the p value for a second phenotype j is less than a cutoff p_j, under the null hypothesis for phenotype i.
##'
##' If the GWAS for phenotypes i and j share a number of controls, effect sizes for the two phenotypes will be correlated, even at null SNPs. This entails that the expected quantile depends both on the number of shared controls, through a paramter 'rho', and the distribution of z values for phenotype j across all SNPs, not all of which will be null for j.
##'
##' The distribution of z values for phenotype j can be specified in two ways, depending on the parameter 'method'. In the first case, we assume that the distribution of 'true' z values (z values which we would observe if the observed allele frequencies exactly matched the population allele frequencies) has a mixture distribution of 0 with probability pi0, and N(0,sigma^2) with probability (1-pi0). These parameters can be estimated from an observed distribution of z values using the function fit.em.
##'
##' In the second case, the distribution of 'true' z values can be specified as an n x 2 array of z values and corresponding probabilities.
##'
##' @title exp_quantile
##' @param p Observed p value or set of p values for principal phenotype
##' @param pc Observed p value or set of p values for conditional phenotype. Must be of the same length as p.
##' @param rho Correlation between z values induced by shared controls
##' @param  method; Set to 1 to assume Zj has a mixed distribution characterised by pi0j, sigmaj; 2 to specify an empirical distribution. 
##' @param pi0j Relevant only if method == 1. Parameter for distribution of Zj. 
##' @param sigmaj Relevant only if method == 1. Parameter for distribution of Zj
##' @param prior Relevant only if method == 2. Matrix with two rows; first is a list of potential values for Zj, second is a list of probabilities associated with those values.
##' @return A list of probabilities of the same length as p.
##' @export
##' @author James Liley
##' @examples
##' p <- 0.2; pc <- 0.1;
##' rho <- 0.25
##' sigma <- 5; pi0 <- 0.95
##' exp_quantile(p,pc,rho,method=1,pi0j=pi0,sigmaj=sigma)
##' # Generally the expected quantile is close to p. In this case, the SNP is 'probably' null for phenotype j, so the MAF in controls is somewhat aberrant from expected (as pc<0.1). This means that the expected quantile is higher than the p value.
##'
##' p<-0.2; pc <-1e-8;
##' exp_quantile(p,pc,rho,method=1,pi0j=pi0,sigmaj=sigma)
##' # If a low value of pc is observed, it is 'likely' to have arisen from a non-null SNP for phenotype j. This p value is quite reasonable for a non-null SNP, given the distribution of z values for such SNPs (sigma=5), so the aberration from expected MAF amongst controls is minimal.
##'
##' pi0=1
##' exp_quantile(p,pc,rho,method=1,pi0j=pi0,sigmaj=sigma)
##' # If, on the other hand, pi0=1 (no non-null SNPs for phenotype j), and we happen to observe a p value of 1e-8 for the conditional phenotype, we conclude that MAF in controls is likely to be very aberrant from expected, so the expected quantile is markedly higher than the p value.
##'
##' rho <- 0
##' exp_quantile(p,pc,rho,method=1,pi0j=pi0,sigmaj=sigma)
##' # If rho=0 (no shared controls) z values for the two phenotypes are independent, and the expected quantile equals pc.
##'
##' p3 = rbind((-500:500)/20, dnorm((-500:500)/20,sd=3))
##' exp_quantile(5e-4,5e-3,0.3,method=2,pi0j=0.9,prior=p3)
##' exp_quantile(5e-4,5e-3,0.3,method=1,pi0j=0.9,sigmaj=3)
##' # Demonstrates specification of a prior.

exp_quantile = function(p,pc,rho,method=1,pi0j=0.9,sigmaj=3,prior = rbind((1:10)/10,dnorm((1:10)/10,sd=3))) {

require(pbivnorm)

z = -qnorm(p/2)
zc = -qnorm(pc/2);

if (method==1) {
 num = (pi0j*pbivnorm(-z,-zc,rho=rho)) + (pi0j*pbivnorm(-z,-zc,rho=-rho)) + ((1-pi0j)*pbivnorm(-z,-zc/sqrt(1+(sigmaj^2)),rho=rho/sqrt(1+(sigmaj^2)))) +  ((1-pi0j)*pbivnorm(-z,-zc/sqrt(1+(sigmaj^2)),rho=-rho/sqrt(1+(sigmaj^2))))
 denom = (pi0j*pnorm(-zc)) + (1-pi0j)*pnorm(-zc/sqrt(1+(sigmaj^2)))
 return(num/denom)
}
if (method==2) {
 prior[,2]=prior[,2]/(trap(prior[,1],prior[,2]))
 num = matrix(0,length(z), dim(prior)[2]); denom=num;
 for (i in 1:dim(prior)[2]) {
  num[,i] = pbivnorm(-z,-zc-prior[1,i],rho=rho) + pbivnorm(-z,-zc+prior[1,i],rho=rho) + pbivnorm(-z,-zc-prior[1,i],rho=-rho) + pbivnorm(-z,-zc+prior[1,i],rho=-rho)
  denom[,i] = pnorm(-zc-prior[1,i]) + pnorm(-zc+prior[1,i])
  }
 num_null = 2*((pbivnorm(-z,-zc,rho=rho) + pbivnorm(-z,-zc,rho=-rho))); denom_null = 2*pnorm(-zc); 
 num_nonnull=num_null; denom_nonnull=denom_null; # initialise
 for (i in 1:length(z)) {
  num_nonnull[i] = trap(prior[1,],num[i,]*prior[2,]); 
  denom_nonnull[i] = trap(prior[1,],denom[i,]*prior[2,])
 }
 return(((pi0j*num_null) + ((1-pi0j)*(num_nonnull)))/((pi0j*denom_null) + ((1-pi0j)*(denom_nonnull))))
}
}


# Trapezoid rule

trap = function(x,y) {
0.5*((c(x,0)-c(0,x))[2:(length(x))] %*% (c(y,0)+c(0,y))[2:(length(y))])
}
