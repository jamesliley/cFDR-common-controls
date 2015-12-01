##' Compute cFDR values from a set of pairs of p values
##'
##' Assumes a distribution of z values for conditional phenotype distributed as 0 with probability pi0j, N(0,sigmaj^2) with probability 1-pi0j. 
##'
##' @title cfdr
##' @param p_vals_ij an array of dimension n x 2; entry [i,1] is the p value for the ith SNP at the principal phenotype, [i,2] is the p value for the ith SNP at the jth phenotype.
##' @param rho correlation between z values arising from shared controls; output from cor_shared
##' @param method Passed to exp_quantile. Set to 1 to assume Zj has a mixed distribution characterised by pi0j, sigmaj; 2 to specify an empirical distribution. 
##' @param pi0j Passed to exp_quantile. Relevant only if method == 1. Parameter for distribution of Zj. 
##' @param sigmaj Passed to exp_quantile. Relevant only if method == 1. Parameter for distribution of Zj
##' @param prior Passed to exp_quantile. Relevant only if method == 2. Matrix with two rows; first is a list of potential values for Zj, second is a list of probabilities associated with those values.
##' @param all set to 1 to compute cfdr at all possible values; 0 false to only compute at values for which it will potentially be less than 0.1 (faster).
##' @return Array of dimensions n x 2; entry [i,1] is cFDR(p_vals[i,1]|p_vals[i,2]); entry [i,2] is cFDR(p_vals[i,2]|p_vals[i,1]) 
##' @export
##' @author James Liley
##' @examples
##' require(mnormt)
##' rho = cor_shared(1500,1200,1000,700,800)
##' z_vals = rmnorm(10000,mean=c(0,0),varcov=rbind(c(1,rho),c(rho,1)))
##' p_vals = 2*pnorm(-abs(z_vals))
##'
##' rec_z = -qnorm(p_vals/2) # Recovered z values
##' f_n = fit.em(c(rec_z,-rec_z))
##' pi0 = f_n$pars[1]; sigma=sqrt(f_n$pars[2]) 
##'
##'
##' cond_fdr = cfdr(p_vals,rho, method=1, pi0j=pi0,sigmaj=sigma)
##' 
##' # Explicit computation
##' p_vals[1,]
##' (a=length(which(p_vals[,2]<=p_vals[1,2])))
##' (b=length(which((p_vals[,1]<=p_vals[1,1]) & (p_vals[,2]<=p_vals[1,2]))))
##' (c=exp_quantile(p_vals[1,1],p_vals[1,2],rho,pi0j=pi0,sigmaj=sigma))
##' c*a/b
##' cond_fdr[1]

cfdr = function(p_vals_ij,rho,method=1, pi0j = 1, sigmaj = 1, prior = rbind((1:10)/10,dnorm((1:10)/10,sd=3)),all=FALSE) {

f_i = rep(0,dim(p_vals_ij)[1])

if (! all) {
# Compute cutoff c_i
p_i = p_vals_ij[,1]; p_j = p_vals_ij[,2];
ord=order(p_i); cdf=(1:length(ord))/length(ord) # Quantiles of ordered list of p values
p_i = p_i[ord];
fdr1 = p_i/cdf;
c_i = p_i[min(which(fdr1>0.1))]

r1 = which(log(p_vals_ij[,1])+log(p_vals_ij[,2])<0.5*log(c_i)) # values not in triangle (p_i,p_j<1; -log10(p_i) + -log10(p_j) <-log10(c_i))
r2 = which(p_vals_ij[,1]<(10^-0.75)) # 
#r3 = which(p_vals_ij[,1]>(10^-12)) # SNPs with values lower than this 
r = intersect(r1,r2) # intersect(intersect(r1,r2),r3) # Values with ambiguous FDR

total = 1:length(f_i);
f_i[setdiff(total,r1)]=1; # approximately
f_i[setdiff(total,r2)]=1; # approximately
#f_i[setdiff(total,r3)]=0;

} else {
r = 1:dim(p_vals_ij)[1]
}

p_i = p_vals_ij[,1]; # Remaining points
p_j = p_vals_ij[,2];


# Compute empirical FDR at remaining points

cdf = rep(0,length(p_i))

oj = order(p_j)
p_i1 = p_i[oj]; p_j1 = p_j[oj] #
or = order(oj)[r] # index of snps in R in p_i1/p_j1

cr = rep(0,length(r))
for (i in 1:length(r)) {
 cdf[r[i]] = length(which(p_i1[1:or[i]]<=p_i1[or[i]]))/or[i]
}

# Equivalent technique

#for (i in 1:length(r)) {
# cdf[r[i]] = length(which(p_i<=p_i[r[i]] & p_j<=p_j[r[i]]))/length(which(p_j<=p_j[r[i]]))
#}

f_i[r]=exp_quantile(p_vals_ij[r,1],p_vals_ij[r,2],rho,method=method,pi0j=pi0j,sigmaj=sigmaj,prior=prior)/cdf[r] # Computation of FDR

return(f_i)
}






##' Compute correlation between z values arising from shared controls
##'
##' Utilises a formula from Zaykin et al (2010)
##'
##' @title cor_shared
##' @param n0i Number of controls in GWAS for phenotype i
##' @param n0j Number of controls in GWAS for phenotype j
##' @param n1i Number of cases in GWAS for phenotype i
##' @param n1j Number of cases in GWAS for phenotype j
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
##' # Simulation - heavy memory usage
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
