##' Evaluate the conditional expected quantile of a p value in the shared control design conditional on a p value bound for a second phenotype
##'
##' Computes the probability that a p value for a phenotype i is less than some cutoff 'p' given that the p value for a second phenotype j is equal to a value 'pc', under the null hypothesis for phenotype i.
##'
##' If the GWAS for phenotypes i and j share a number of controls, effect sizes for the two phenotypes will be correlated, even at null SNPs. This entails that the expected quantile depends both on the number of shared controls, through a paramter 'rho', and the distribution of z values for phenotype j across all SNPs, not all of which will be null for j.
##'
##' The 'true' z values (z values which we would observe if the observed allele frequencies exactly matched the population allele frequencies) is assumed to follow a mixture distribution of 0 with probability pi0, and N(0,sigma^2) with probability (1-pi0). These parameters can be estimated from an observed distribution of z values using the function fit.em.
##'
##' The value of the function is the first derivative of the exp_quantile function with respect to P_j.
##'
##' @title exp_quantile_point
##' @param p Observed p value or set of p values for principal phenotype
##' @param pc Observed p value or set of p values for conditional phenotype. Must be of the same length as p.
##' @param rho Correlation between z values induced by shared controls
##' @param pi0j Parameter for distribution of Zj. Proportion of SNPs null for phenotype j.
##' @param sigmaj Parameter for distribution of Zj. Variance of true effect sizes for SNPs non-null for phenotype j.
##' @return A list of probabilities of the same length as p.
##' @export
##' @author James Liley
##' @examples
##' p <- 0.2; pc <- 0.1;
##' rho <- 0.25
##' sigma <- 5; pi0 <- 0.95
##' exp_quantile_point(p,pc,rho,pi0j=pi0,sigmaj=sigma)
##' # Generally the expected quantile is close to p. 

exp_quantile_point = function(p,pc,rho,pi0j=1,sigmaj=0) {

require(pbivnorm)

z = -qnorm(p/2)
zc = -qnorm(pc/2);

 num = (dnorm(-zc)*pi0j*
   (pnorm(-z,mean = -rho*zc,sd = sqrt(1-(rho^2)))+pnorm(-z,mean = rho*zc,sd = sqrt(1-(rho^2))))) + 
  (dnorm(-zc,sd = sqrt(1+(sigmaj^2)))*(1-pi0j)*
   (pnorm(-z,mean = -rho*zc/sqrt(1+(sigmaj^2)),sd = sqrt(1-(rho^2)+(sigmaj^2))/sqrt(1+(sigmaj^2))) +
    pnorm(-z,mean = rho*zc/sqrt(1+(sigmaj^2)),sd = sqrt(1-(rho^2)+(sigmaj^2))/sqrt(1+(sigmaj^2)))))
 denom = ((pi0j*dnorm(-zc)) + (1-pi0j)*dnorm(-zc,sd = sqrt(1+(sigmaj^2))))

return(num/denom)
}
