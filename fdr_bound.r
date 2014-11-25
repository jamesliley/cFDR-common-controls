##' Computation of overall false discovery rate given a cFDR cutoff, given distribution parameters for z values for conditional phenotype at SNPs null for principal phenotype
##'
##' If, for some SNP with principal and conditional p values equal to p_i and p_j respectively, we have cFDR(p_i|p_j) = a, then amongst the set of SNPs for which P_i <= p_i and P_j <= p_j, the expected false discovery rate (proportion of false positives) is less than a.
##'
##' However, amongst the set of SNPs for which cFDR(p_i|p_j) <= a, the expected false discovery rate is not controlled at a. A bound on the expected false discovery rate can be computed using the area of the region in which cFDR<a, and an assumed distribution of z-values of null SNPs.
##'
##' This function provides an implementation of that calculation.
##'
##' @title c2a
##' @param p_vals an array of dimension n x 2; entry [i,1] is the p value for the ith SNP at the principal phenotype, [i,2] is the p value for the ith SNP at the jth phenotype.
##' @param cfdr the cFDR of p_vals[,1]|p_vals[,2]; output from function cfdr.
##' @param cutoff the significance cutoff used for cFDR.
##' @param rho correlation between z values arising from shared controls; output from cor_shared
##' @param pi0_n proportion of SNPs null for principal phenotype which are also null for conditional phenotype
##' @param sigma_n standard deviation of z-scores for SNPs null for conditional phenotype but non-null for conditional phenotype
##' @return upper bound on the expected false discovery rate amongst all SNPs for which cFDR<cutoff.
##' @export
##' @author James Liley
##' @examples
##'
 
c2a=function(p_vals,cfdr,cutoff,rho,pi0_n,sigma_n) {
w = which(!is.na(cfdr))
L = which(cfdr[w]<cutoff)
ufdr = p_vals[w,1]*length(w)/rank(p_vals[w,1]) 

ul = which.max(p_vals[w[L],1]); lr = which.min(p_vals[w[-L],1]);
left = -qnorm(p_vals[w[L[ul]],1]/2); upper = -qnorm(p_vals[w[L[ul]],2]/2)
right = -qnorm(p_vals[w[-L][lr],1]/2)+1; lower = 0 #-qnorm(p_vals[w[-L][lr],2]/2)

if (left >= right | length(L)<50 | max(p_vals[w[which(ufdr<cutoff)],1])> max(p_vals[which(cfdr<cutoff),1]) ) {
 return(cutoff)
} else {

# Values left,right,upper and lower bound the curved region of L. L also contains M* and the region V_3 bounded by
#  left<zi<right,zj>upper.



V_1 = (pi0_n*pmnorm(c(-left, -upper),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1)))) +
 ((1-pi0_n)*pmnorm(c(-left, -upper),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1+(sigma_n^2)))))
# integral over area bounded by minz<zi, a_z<zj; where a_z = z_j(which(zi=minz)
V_2 = (pi0_n*pmnorm(c(-right,-upper),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1)))) +
 ((1-pi0_n)*pmnorm(c(-right,-upper),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1+(sigma_n^2)))))
# integral over area bounded by maxz<zi, a_z<zj

V_3 = V_1 - V_2
# integral over area bounded by minp<zi<maxp; a_z<zj


V_4 = (pi0_n*pmnorm(c(-left, -upper),mean = c(0,0),varcov = rbind(c(1,-rho),c(-rho,1)))) +
 ((1-pi0_n)*pmnorm(c(-left, -upper),mean = c(0,0),varcov = rbind(c(1,-rho),c(-rho,1+(sigma_n^2)))))
# integral over area bounded by minz<zi, a_z<zj; where a_z = z_j(which(zi=minz)
V_5 = (pi0_n*pmnorm(c(-right,-upper),mean = c(0,0),varcov = rbind(c(1,-rho),c(-rho,1)))) +
 ((1-pi0_n)*pmnorm(c(-right,-upper),mean = c(0,0),varcov = rbind(c(1,-rho),c(-rho,1+(sigma_n^2)))))
# integral over area bounded by maxz<zi, a_z<zj

V_6 = V_4 - V_5
# integral over area bounded by -minp<-zi<-maxp; a_z<zj



xx = seq(left, right,length.out = 100); yy = seq(lower,upper,length.out = 100);
x0 = matrix(xx,length(xx),length(xx)); y0 = t(matrix(yy,length(yy),length(yy)))

# Messy numerical necessities and interpolation over region bounded by lower, upper, left and right
potential = which(p_vals[w,1]<2.2*pnorm(-abs(left))) 
k = 20; pot_s= potential[(1:floor(length(potential)/k))*k] 

 zz = interp(c(-qnorm(p_vals[w[potential],1]/2),-qnorm(p_vals[w[pot_s],1]/2)),c(-qnorm(p_vals[w[potential],2]/2),rep(0,length(pot_s))),c(cfdr[w[potential]],ufdr[pot_s]),xo = xx,yo=yy,duplicate="median") # Interpolation to cFDR values for a grid
# zz$z[i,j] is now the interpolated value of cFDR(i,j)

density_1 = ff(c(x0),c(y0),rho,pi0_n,sigma_n);
density_2 = ff(c(x0),c(y0),-rho,pi0_n,sigma_n);
int = (sum(density_1[which(zz$z<cutoff)]) + sum(density_2[which(zz$z<cutoff)]))*(upper-lower)*(right-left)/(length(xx)*length(yy))
# Integral over rectangle [left, right] x [lower,upper] of function zz.


if (max(zz$z[which(y0 <3)],na.rm=TRUE)>cutoff) {V_M = 2*pnorm(-max(x0[which(y0<3)][which(zz$z[which(y0<3)]>cutoff)]))} else {V_M=min(max(p_vals[w[which(ufdr<cutoff)],1]),p_vals[w[-L][lr],1])}
# integral over M

V_L = 2*int + 2*(V_3 + V_6) + 2*pnorm(-right)

c_r = approx(p_vals[w,1], ufdr,V_M)$y # c_r is the cFDR of the upper right border of M*, equivalent to the UFDR at that point.

if (FALSE) {
plot(-qnorm(p_vals[,1]/2),-qnorm(p_vals[,2]/2),pch=".",xlim = c(0,10),ylim = c(0,10))
points(-qnorm(p_vals[w[L],1]/2),-qnorm(p_vals[w[L],2]/2),col="blue",cex= 0.3)
points(-qnorm(p_vals[w[which(ufdr<cutoff)],1]/2),-qnorm(p_vals[w[which(ufdr<cutoff)],2]/2),col="green",cex= 1)
points(-qnorm(p_vals[w[L],1]/2),-qnorm(p_vals[w[L],2]/2),col="blue",cex= 0.3)
lines(c(left, left,right,right),c(lower,upper,upper,lower),col="red");
wx = which(zz$z<cutoff); points(x0[wx],y0[wx],col="red",pch=".")
lines(c(-qnorm(2.5e-8),-qnorm(2.5e-8)),c(0,10),col="green")
lines(c(-qnorm(V_M/2),-qnorm(V_M/2)),c(0,10),col="gray")
}

return(V_L* c_r /V_M)

}
}



##' PDF of z-value pairs for SNPs which are null for principal phenotype
##'
##' Form of distribution is normal with mean (0,0), variance (1,rho; rho,1) with probability pi0_n
##'
##' Parameters pi0_n and sigma_n can be estimated from distribution of z_j; see fit.em function and implementation of fdr_bound function
##'
##' @title ff
##' @param z_i z value for principal phenotype at which to evaluate pdf
##' @param z_j z value for conditional phenotype at which to evaluate pdf
##' @param rho correlation between z values arising from shared controls; output from cor_shared
##' @param pi0 proportion of SNPs null for principal phenotype which are also null for conditional phenotype
##' @param sigma standard deviation of z-scores for SNPs null for conditional phenotype but non-null for conditional phenotype
##' @return value of pdf at z_i,z_j
##' @export
##' @author James Liley
##' @examples
##'
 
ff = function(z_i,z_j,rho=0,pi0=0.9,sigma=3) {
require(mnormt)
(pi0*dmnorm(cbind(z_i,z_j),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1)))) +
 ((1-pi0)*dmnorm(cbind(z_i,z_j),mean = c(0,0),varcov = rbind(c(1,rho),c(rho,1+(sigma^2)))))
}





##' PDF of p-value pairs for SNPs which are null for principal phenotype
##'
##' Form of distribution of z values is normal with mean (0,0), variance (1,rho; rho,1) with probability pi0_n. Distribution transformed to p-value space.
##'
##' Parameters pi0_n and sigma_n can be estimated from distribution of z_j; see fit.em function and implementation of fdr_bound function
##'
##' @title pdf_p
##' @param p_i p value for principal phenotype at which to evaluate pdf
##' @param p_j p value for conditional phenotype at which to evaluate pdf
##' @param rho correlation between z values arising from shared controls; output from cor_shared
##' @param pi0 proportion of SNPs null for principal phenotype which are also null for conditional phenotype
##' @param sigma standard deviation of z-scores for SNPs null for conditional phenotype but non-null for conditional phenotype
##' @return value of pdf at p_i,p_j
##' @export
##' @author James Liley
##' @examples
##'
 
pdf_p = function(p_i,p_j,rho=0,pi0=0.9,sigma=3) {
 (ff(-qnorm(p_i/2),-qnorm(p_j/2),rho,pi0,sigma) + 
 ff(qnorm(p_i/2),-qnorm(p_j/2),rho,pi0,sigma) + 
 ff(-qnorm(p_i/2),qnorm(p_j/2),rho,pi0,sigma) + 
 ff(qnorm(p_i/2),qnorm(p_j/2),rho,pi0,sigma)) *
 (pi/2)*exp(0.5*(qnorm(p_i/2)^2 + qnorm(p_j/2)^2))
 }
