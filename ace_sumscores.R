#============================================================================#
# Following JAGS script fits a longitudinal ACE twin model using sum-scores
# T = total number of time points 
# n_dz = total number of DZ twin families; n_dz = total number of MZ families 
#
# Any questions: Contact Inga Schwabe (I.Schwabe@uvt.nl)
#============================================================================##

model{
    
    # Twin model part for MZ twins
    for (fam in 1:n_mz){
        a1_mz[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        c_mz[fam,1:T] ~ dmnorm(mu_c[1:T], tau_c[1:T,1:T])
        
        # Modelling phenotypes: 
        data_mz_twin1[fam,1:T] ~ dmnorm(mu[1:T] + a1_mz[fam,1:T] + 
                                        c_mz[fam,1:T], tau_e[1:T,1:T])
        data_mz_twin2[fam,1:T] ~ dmnorm(mu[1:T] + a1_mz[fam,1:T] + 
                                        c_mz[fam,1:T], tau_e[1:T,1:T])
    }
    
    # Twin model part for MZ twins
    for (fam in 1:n_dz){
        a1_dz[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        a2_dz[fam,1:T] <- a1_dz[fam,1:T]/sqrt(2)
        
        a3_dz_twin1[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        a3_dz_twin2[fam,1:T] ~ dmnorm(mu_a[1:T], tau_a[1:T,1:T])
        
        a4_dz_twin1[fam,1:T] <- a3_dz_twin1[fam,1:T]/sqrt(2)
        a4_dz_twin2[fam,1:T] <- a3_dz_twin2[fam,1:T]/sqrt(2)
        
        c_dz[fam,1:T] ~ dmnorm(mu_c[1:T], tau_c[1:T,1:T])
        
        # Modelling phenotypes: 
        data_dz_twin1[fam,1:T] ~ dmnorm(mu[1:T] + a2_dz[fam,1:T] + 
                                        a4_dz_twin1[fam,1:T] + c_dz[fam,1:T], 
                                        tau_e[1:T,1:T])
        data_dz_twin2[fam,1:T] ~ dmnorm(mu[1:T] + a2_dz[fam,1:T] + 
                                        a4_dz_twin2[fam,1:T] + c_dz[fam,1:T], 
                                        tau_e[1:T,1:T])        
    }
    
    # Prior distributions
    
    # Normal prior for expected values of time points 
    for (timepoint in 1:T){
        mu[timepoint] ~ dnorm(0, .1)
    }

    # Inverse Wishart prior distribution for covariance matrices 
    tau_a[1:T,1:T] ~ dwish(omega_a[,], T+1) 
    tau_e[1:T,1:T] ~ dwish(omega_e[,], T+1)
    tau_c[1:T,1:T] ~ dwish(omega_c[,], T+1)
    
    # Calculate inverse 
    Sigma_a[1:T,1:T] <- inverse(tau_a[,])
    Sigma_e[1:T,1:T] <- inverse(tau_e[,])
    Sigma_c[1:T,1:T] <- inverse(tau_c[,])
}