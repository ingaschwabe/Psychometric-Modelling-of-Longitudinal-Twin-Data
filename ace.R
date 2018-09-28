#============================================================================#
# Following JAGS script fits a longitudinal ACE latent trait-state AE twin model
# T = total number of time points (equal to 2 here)
# n_dz = total number of DZ twin families; n_dz = total number of MZ families 
# n_items = total number of questionnaire items
#
# Note that this script works for 2 time points but can easily be adjusted 
# to work for more time points. 
# Any questions: Contact Inga Schwabe (I.Schwabe@uvt.nl)
#============================================================================#

model{
    # Twin model part for MZ twins
    for (fam in 1:n_mz){
        
        a1_mz[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        c_mz[fam, 1:T] ~ dmnorm(mu_c[1:T], tau_c[1:T,1:T])
        
        # Modelling phenotypes:
        pheno_mz_twin1[fam,1:T] ~ dmnorm(mu[1:T] + a1_mz[fam,1:T] + 
                                         c_mz[fam,1:T], tau_e[1:T,1:T])
        pheno_mz_twin2[fam,1:T] ~ dmnorm(mu[1:T] + a1_mz[fam,1:T] + 
                                         c_mz[fam,1:T], tau_e[1:T,1:T])
        
        # IRT part: Rasch model for both time points 
        for (k in 1:n_items){
            logit(p_mz_twin1_t1[fam,k]) <- pheno_mz_twin1[fam,1] - beta[k]
            data_mz_twin1_t1[fam,k] ~ dbern(p_mz_twin1_t1[fam,k])
            
            logit(p_mz_twin2_t1[fam,k]) <- pheno_mz_twin2[fam,1] - beta[k]
            data_mz_twin2_t1[fam,k] ~ dbern(p_mz_twin2_t1[fam,k])
        }
        
        for (k in 1:n_items){
            logit(p_mz_twin1_t2[fam,k]) <- pheno_mz_twin1[fam,2] - beta[k]
            data_mz_twin1_t2[fam,k] ~ dbern(p_mz_twin1_t2[fam,k])
            
            logit(p_mz_twin2_t2[fam,k]) <- pheno_mz_twin2[fam,2] - beta[k]
            data_mz_twin2_t2[fam,k] ~ dbern(p_mz_twin2_t2[fam,k])
        }
    }
    
    # Twin model part for DZ twins 
    for (fam in 1:n_dz){
        a1_dz[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        a2_dz[fam,1:T] <- a1_dz[fam,1:T]/sqrt(2)
        
        a3_dz_twin1[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        a3_dz_twin2[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        
        a4_dz_twin1[fam,1:T] <- a3_dz_twin1[fam,1:T]/sqrt(2)
        a4_dz_twin2[fam,1:T] <- a3_dz_twin2[fam,1:T]/sqrt(2)
        
        c_dz[fam,1:T] ~ dmnorm(mu_c[1:T], tau_c[1:T,1:T])
        
        # Modelling phenotypes:
        pheno_dz_twin1[fam,1:T] ~ dmnorm(mu[1:T] + a2_dz[fam,1:T] + 
                                         a4_dz_twin1[fam,1:T] + c_dz[fam,1:T], 
                                         tau_e[1:T,1:T])
        pheno_dz_twin2[fam,1:T] ~ dmnorm(mu[1:T] + a2_dz[fam,1:T] + 
                                         a4_dz_twin2[fam,1:T] + c_dz[fam,1:T], 
                                         tau_e[1:T,1:T])
        
        # IRT part (Rasch model for both time points): 
        for (k in 1:n_items){
            logit(p_dz_twin1_t1[fam,k]) <- pheno_dz_twin1[fam,1] - beta[k]
            data_dz_twin1_t1[fam,k] ~ dbern(p_dz_twin1_t1[fam,k])
            
            logit(p_dz_twin2_t1[fam,k]) <- pheno_dz_twin2[fam,1] - beta[k]
            data_dz_twin2_t1[fam,k] ~ dbern(p_dz_twin2_t1[fam,k])
        }
        
        for (k in 1:n_items){
            logit(p_dz_twin1_t2[fam,k]) <- pheno_dz_twin1[fam,2] - beta[k]
            data_dz_twin1_t2[fam,k] ~ dbern(p_dz_twin1_t2[fam,k])
            
            logit(p_dz_twin2_t2[fam,k]) <- pheno_dz_twin2[fam,2] - beta[k]
            data_dz_twin2_t2[fam,k] ~ dbern(p_dz_twin2_t2[fam,k])
        }
    }
    
    # Prior distributions  					  							     
    # Normal prior for beta item parameters: 
    for (k in 1:n_items){
        beta[k] ~ dnorm(0, .1)
    } 
    
    # The expected value of the first time point is set to 0 to identify the 
    # measurement scale. A normal prior is set on the remaining time point averages
    mu[1] <- 0 
    mu[2] ~ dnorm(0, .1)
    
    # Inverse Wishart prior distribution for covariance matrices 
    tau_a[1:T,1:T] ~ dwish(omega_a[,], T+1) 
    tau_e[1:T,1:T] ~ dwish(omega_e[,], T+1)
    tau_c[1:T,1:T] ~ dwish(omega_c[,], T+1)
    
    # Calculate inverse: 
    Sigma_a[1:T,1:T] <- inverse(tau_a[,])
    Sigma_e[1:T,1:T] <- inverse(tau_e[,])
    Sigma_c[1:T,1:T] <- inverse(tau_c[,])
}