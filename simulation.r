# Pleiotropy-informed detection of genetic variants using a shared control group
#
# Simulation to demonstrate adjustment of p values
# Assumes haploid organisms for simplicity
#

pi0 = 0.7 #0.95
sigma = 1.5

n_snps = 20000;
n_cases_p = 1200; # Number of cases in principal phenotype
n_cases_c = 1300; # Number of cases in conditional phenotype

n_controls_p = 1000;
n_controls_c = 1100;
overlap = 1000;



# Set up population effect sizes for conditional phenotype

pop_z_c = rnorm(n_snps,mean=0,sd=sigma)

#pop_z_c = runif(n_snps,min=-sigma,max=sigma)   # Alternative: uniform distribution
#pop_z_c = rt(n_snps,df=2) # t distribution
#pop_z_c = sample(c(sigma,-sigma),n_snps,replace=TRUE) # discrete distribution

pop_z_c[which(runif(n_snps)<pi0)]=0 # Effect size is 0 for pi0 % of SNPs



# Set up minor allele population frequencies to match population effect sizes; all effect sizes positive in this simulation

maf_control=runif(n_snps,min=0.1,max=0.5)
maf_case_p = maf_control
maf_case_c = maf_control+ ((pop_z_c*sqrt(maf_control*(1-maf_control)*((1/(n_controls_c))+(1/(n_cases_c)))/2)));



# Check population effect sizes match what they should

s_or = maf_case_c*(1-maf_control)/(maf_control*(1-maf_case_c))
se = sqrt((1/(maf_control*n_controls_c)) + (1/((1-maf_control)*n_controls_c)) + (1/(maf_case_c*n_cases_c)) + (1/((1-maf_case_c)*n_cases_c)))

z = log(s_or)/se

#plot(abs(z),abs(pop_z_c))



# Set up sample genotypes

g_p_control=matrix(runif(n_snps*n_controls_p)<rep(maf_control,n_controls_p),n_snps,n_controls_p)
g_c_control=matrix(runif(n_snps*n_controls_c)<rep(maf_control,n_controls_c),n_snps,n_controls_c)
g_c_control[,1:overlap]=g_p_control[,1:overlap] # Overlap in control groups

g_p_case = matrix(runif(n_snps*n_cases_p)<rep(maf_case_p,n_cases_p),n_snps,n_cases_p)
g_c_case = matrix(runif(n_snps*n_cases_c)<rep(maf_case_c,n_cases_c),n_snps,n_cases_c)



# Observed MAFs

obs_maf_control_p = rowMeans(g_p_control)
obs_maf_control_c = rowMeans(g_c_control)

obs_maf_case_p = rowMeans(g_p_case)
obs_maf_case_c = rowMeans(g_c_case)



# Observed effect sizes

log_or_p = log( obs_maf_control_p*(1-obs_maf_case_p)/(obs_maf_case_p*(1-obs_maf_control_p)) )
se_p = sqrt((1/(obs_maf_control_p*n_controls_p)) + (1/((1-obs_maf_control_p)*n_controls_p)) + (1/(obs_maf_case_p*n_cases_p)) + (1/((1-obs_maf_case_p)*n_cases_p)))#*(1/sqrt(2))*

z_p = log_or_p/se_p
p_p = 2*pnorm(-abs(z_p))

log_or_c = log( obs_maf_control_c*(1-obs_maf_case_c)/(obs_maf_case_c*(1-obs_maf_control_c)) )
se_c = sqrt((1/(obs_maf_control_c*n_controls_c)) + (1/((1-obs_maf_control_c)*n_controls_c)) + (1/(obs_maf_case_c*n_cases_c)) + (1/((1-obs_maf_case_c)*n_cases_c)))#*(1/sqrt(2))

z_c = log_or_c/se_c
p_c = 2*pnorm(-abs(z_c))



# Recovery of pi0, sigma, and rho

z_z = abs(z_c);
z_z = z_z[which(!is.na(z_z))]
z_z = z_z*sample(c(-1,1),length(z_z),replace=TRUE)

ps = fit.em(z_z,pi0=0.9,s2=1,tol=1e-4,verbose=TRUE,maxit=1000) 
pi0_obs = ps$pars[1]
sigma_obs = sqrt(ps$pars[2])





# Computation of expected correlation between effect sizes due to shared controls, as per Zaykin et al.

N1=n_cases_p; N2=n_cases_c; N0=overlap; N01 = n_controls_p-N0; N02=n_controls_c-N0
rho = sqrt(1/( (1+(N02/N0))*(1+(N01/N0))*(1+(N01/N1)+(N0/N1))*(1+(N02/N2)+(N0/N2)) ));





# Our function for computing expected quantile

exp_quantile = function(p,pc,rho,pi0j=0.9,sigmaj=3) {

z = -qnorm(p/2)
zc = -qnorm(pc/2);

 num = (pi0j*pbivnorm(-z,-zc,rho=rho)) + (pi0j*pbivnorm(-z,-zc,rho=-rho)) + ((1-pi0j)*pbivnorm(-z,-zc/sqrt(1+(sigmaj^2)),rho=rho/sqrt(1+(sigmaj^2)))) +  ((1-pi0j)*pbivnorm(-z,-zc/sqrt(1+(sigmaj^2)),rho=-rho/sqrt(1+(sigmaj^2))))
 denom = (pi0j*pnorm(-zc) + (1-pi0j)*pnorm(-zc/sqrt(1+(sigmaj^2))))
 return(num/denom)
}





# Plot showing inflation

#pdf("inflation.pdf")

c_cut = 0.05
subp = which(p_c<c_cut)

plot(sort(p_p[subp]),main="Principal p values with conditional p < 0.05",ylab="p value",xlab="Quantile")
lines(c(0,length(subp)),c(0,1),col="red")

#dev.off()





# Plot showing correction using our algorithm

#pdf("correction1.pdf")

c_cut = 0.05
subp = which(p_c<c_cut)

plot(sort(p_p[subp]),main="Principal p values and adjusted p values",ylab="p value",xlab="Quantile")
points(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=pi0_obs,sigmaj=sigma_obs)[order(p_p[subp])],col="blue",cex=0.5)
lines(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=pi0_obs,sigmaj=sigma_obs)[order(p_p[subp])],col="blue")
lines(c(0,length(subp)),c(0,1),col="red")

l_lab = c("Raw p","Estimated quantile")
#legend(0.6,1,l_lab,c("black","blue"))

#dev.off()




# Plot showing correction using our algorithm and correction assuming all SNPs are null for conditional phenotype

pdf("/Users/James/University/Cambridge/Year 1/Projects/Shared Controls/Figures/Simulation/correction2.pdf")

c_cut = 0.05
subp = which(p_c<c_cut)

plot(sort(p_p[subp]),main="Principal p values and adjusted p values",ylab="p value/adjusted p value",xlab="Quantile")
points(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=pi0_obs,sigmaj=sigma_obs)[order(p_p[subp])],col="blue",cex=0.5,pch=16)
lines(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=pi0_obs,sigmaj=sigma_obs)[order(p_p[subp])],col="blue",lwd=2)

points(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=1,sigmaj=0)[order(p_p[subp])],col="red",cex=0.5,pch=16)
lines(exp_quantile(p_p[subp],rep(c_cut,length(subp)),rho, pi0j=1,sigmaj=0)[order(p_p[subp])],col="red",lwd=2)

lines(c(0,length(subp)),c(0,1),col="red")

l_lab = c("Raw p","Assumed null","Estimated quantile")
legend(0.6,1,l_lab,c("black","red","blue"))

dev.off()



