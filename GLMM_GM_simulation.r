# combined first, then split as different files


library("lme4")
# source("ul_ci.r")




# data generating process of unbalanced Binomial Data with one Gaussian random effect
sim_ul_group_time_paper <- function(
  x1_type = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma = 0.5,  # standard deviation for random effects
  beta0 = -0.3,  # intercept of linear part
  beta1 = -3,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2 = 2,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3 = 0.2,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern = 0.5  # probability for Bernoulli distribution in baseline variable
){

  n_group = 4  # number of groups, 
  N=sum(ni)  # Total observations for all objects

  dat <- matrix(0, ncol = 7, nrow = N)  # initializ data matrix
  
  # loop to generate id, x1, x2, x3
  k=1  # object indicator
  for(i in 1:n_group) {   # loop for each group
    # input id
    if(x3_type==1){  # x3 is gender
      id = k:(k + ni[i] - 1)
    }
    if(x3_type==2){  # x3 is time      
      if(i==1) {id <- c(1:ni[1]) }    
      if(i==3) {id <- c((ni[1]+1):(ni[1]+ni[3])) }  # id are the same for group 1,2 ; group 3,4; as group 1,2 have same patients, and group 3,4 have same patients.
    }
    # input x1i
    if(x1_type==1 && x3_type==1 ){  # x1~Bern(p), x3 is gender
      x1i <- rbinom(ni[i], 1, p_bern)  # x1~Bern(p)
    }
    if(x1_type==1 && x3_type==2){  # x1~Bern(p), x3 is time
      if(i==1||i==3) {x1i <- rbinom(ni[i], 1, p_bern)}  # randomly generate x1i~ Bernoulli(p_bern), generate new x1i only for group 1,3.
    }
    if(x1_type==2 && x3_type==1){  # x1~Unif(0,1), x3 is gender
      x1i <- runif(ni[i], min=0, max=1)  # x1~Bern(p)
    }
    if(x1_type==2 && x3_type==2){  # x1~Unif(0,1), x3 is time
      if(i==1||i==3) {x1i <- runif(ni[i], min=0, max=1)}  # randomly generate x1i~ Bernoulli(p_bern), generate new x1i only for group 1,3.
    }
    # input x2 and x3
    if(i==1||i==2) {x2i <- 1 } else {x2i = 0}   # x2i=1 for group 1,2, i.e. treat=1;x2i=0 for group 3,4, i.e. treat=0
    if(i==1||i==3) {x3i <- 0 } else {x3i = 1}   # x3i=0 for group 1,3, i.e. time=0;x3i=1 for group 2,4, time=1
    
    dat[k:(k + ni[i] - 1), 1] = id[1:ni[i]]         # input id into data matrix;  
    dat[k:(k + ni[i] - 1), 3] = x1i[1:ni[i]]        # input baseline variable into data matrix;
    dat[k:(k + ni[i] - 1), 4] = x2i        # input treat indicator into data matrix; 
    dat[k:(k + ni[i] - 1), 5] = x3i        # input time indicator into data matrix;
    
    k <- k + ni[i]       # object indictor moves to next group
  }
  
  id = dat[,1]; x1 = dat[,3]; x2 = dat[,4]; x3 = dat[,5];  # get id, x1, x2, x3 for all of 4 groups from data matrix     
  
  # Loop to generate y and ui(bi)
  k <- 1   # object indicator
  # generate ui for RE
  for(i in 1:n_group) {   # loop for each group
    x1i=dat[k:(k + ni[i] - 1), 3]      # get x1i from data matrix
    x2i=dat[k:(k + ni[i] - 1), 4]      # get x2i from data matrix
    x3i=dat[k:(k + ni[i] - 1), 5]      # get x3i from data matrix
    if(x3_type==1){  # x3 is gender  
      ui <- sigma * rnorm(ni[i])
    }
    if(x3_type==2){  # x3 is time
      if(i==1||i==3) {ui <- sigma * rnorm(ni[i])};  # random effects are same for group 1,2 and same for  group 3,4   
    }  
    eta <- beta0+beta1*x1i+beta2*x2i+beta3*x3i+ui[1:ni[i]];     # generate link or eta for negative binomial 
    mu <- exp(eta)/(1 + exp(eta))        # generate mu(i.e. prob) by logit link
    yi <- rbinom(ni[i],1,mu)  # generate y by logistic distribution
    
    dat[k:(k + ni[i] - 1), 2] <- yi    # input yi into 2nd column of data matrix
    dat[k:(k + ni[i] - 1), 6] <- ui[1:ni[i]]    # input uiinto 5th column of data matrix
    dat[k:(k + ni[i] - 1), 7] <- mu    # input mu ie. prob into 7th column of data matrix      
    
    k <- k + ni[i]          # indicator move to next object
  }
  
  y=dat[,2];b <- dat[,6];mu <- dat[,7]  # get y, b and mu for all of 4 groups from data matrix     
  
  mydata = data.frame(id, y, x1, x2, x3, b, mu)   # generate data.frame 
  
  return(mydata)   # return the Negative Binomial Data by data.frame
}


# data generating process of unbalanced Negative Binomial Data with one Gaussian random effect
sim_unb_group_time_paper <- function(
  x1_type = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma = 0.1,  # standard deviation for random effects
  size = 50,  # size for negative binomial data; for logistic data, no need to input this value
  beta0 = 0.3,  # intercept of linear part
  beta1 = -0.2,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2 = 0.3,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3 = 0.4,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern = 0.5  # probability for Bernoulli distribution in baseline variable
){  

  n_group = 4   # number of groups, 
  N=sum(ni)  # Total observations for all objects
  
  dat <- matrix(0, ncol = 7, nrow = N)  #initializ data matrix
  
  #loop to generate id, x1, x2, x3
  k=1  # object indicator
  for(i in 1:n_group) {   # loop for each group
    # input id
    if(x3_type==1){  # x3 is gender
      id = k:(k + ni[i] - 1)
    }
    if(x3_type==2){  # x3 is time      
      if(i==1) {id <- c(1:ni[1]) }    
      if(i==3) {id <- c((ni[1]+1):(ni[1]+ni[3])) }  # id are the same for group 1,2 ; group 3,4; as group 1,2 have same patients, and group 3,4 have same patients.
    }
    #input x1i
    if(x1_type==1 && x3_type==1 ){  # x1~Bern(p), x3 is gender
      x1i <- rbinom(ni[i], 1, p_bern)  #x1~Bern(p)
    }
    if(x1_type==1 && x3_type==2){  # x1~Bern(p), x3 is time
      if(i==1||i==3) {x1i <- rbinom(ni[i], 1, p_bern)}  # randomly generate x1i~ Bernoulli(p_bern), generate new x1i only for group 1,3. In paper prob=0.2
    }
    if(x1_type==2 && x3_type==1){  # x1~Unif(0,1), x3 is gender
      x1i <- runif(ni[i], min=0, max=1)  #x1~Bern(p)
    }
    if(x1_type==2 && x3_type==2){  # x1~Unif(0,1), x3 is time
      if(i==1||i==3) {x1i <- runif(ni[i], min=0, max=1)}  # randomly generate x1i~ Bernoulli(p_bern), generate new x1i only for group 1,3. In paper prob=0.2
    }
    #input x2 and x3
    if(i==1||i==2) {x2i <- 1 } else {x2i = 0}   # x2i=1 for group 1,2, i.e. Treat=1;x2i=0 for group 3,4, i.e. Treat=0
    if(i==1||i==3) {x3i <- 0 } else {x3i = 1}   # x3i=0 for group 1,3, i.e. time=0;x3i=1 for group 2,4, time=1
    
    dat[k:(k + ni[i] - 1), 1] = id[1:ni[i]]         # input id into data matrix;  
    dat[k:(k + ni[i] - 1), 3] = x1i[1:ni[i]]        # input baseline variable into data matrix;
    dat[k:(k + ni[i] - 1), 4] = x2i        # input treat indicator into data matrix; 
    dat[k:(k + ni[i] - 1), 5] = x3i        # input time indicator into data matrix;
    
    k <- k + ni[i]       # object indictor moves to next group
  }
  
  id = dat[,1]; x1 = dat[,3]; x2 = dat[,4]; x3 = dat[,5];  # get id, x1, x2, x3 for all of 4 groups from data matrix     
  
  # Loop to generate y and ui(bi)
  k <- 1   # object indicator
  # generate ui for RE
  for(i in 1:n_group) {   # loop for each group
    x1i=dat[k:(k + ni[i] - 1), 3]      # get x1i from data matrix
    x2i=dat[k:(k + ni[i] - 1), 4]      # get x2i from data matrix
    x3i=dat[k:(k + ni[i] - 1), 5]      # get x3i from data matrix
    if(x3_type==1){  # x3 is gender  
      ui <- sigma * rnorm(ni[i])
    }
    if(x3_type==2){  # x3 is time
      if(i==1||i==3) {ui <- sigma * rnorm(ni[i])};  # random effects are same for group 1,2 and same for  group 3,4   
    }  
    eta <- beta0+beta1*x1i+beta2*x2i+beta3*x3i+ui[1:ni[i]];     # generate link or eta for negative binomial 
    mu <- exp(eta)        # generate mu by NB link
    yi <- rnbinom(ni[i],mu=mu,size=size) # generate y by NB distribution
    
    dat[k:(k + ni[i] - 1), 2] <- yi    # input yi into 2nd column of data matrix
    dat[k:(k + ni[i] - 1), 6] <- ui[1:ni[i]]    # input uiinto 5th column of data matrix
    dat[k:(k + ni[i] - 1), 7] <- mu    # input mu ie. prob into 7th column of data matrix      
    
    k <- k + ni[i]          # indicator move to next object
  }
  
  y=dat[,2];b <- dat[,6];mu <- dat[,7]  # get y, b and mu for all of 4 groups from data matrix     
  
  mydata = data.frame(id, y, x1, x2, x3, b, mu)   # generate data.frame 
  
  return(mydata)   # return the Negative Binomial Data by data.frame
}



#############################################################################








ul_ci <- function(
  x1_type_sim = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type_sim = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni_sim = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma_sim = 0.5,  # standard deviation for random effects
  size_sim = 50,  # size for negative binomial data; for logistic data, no need to input this value
  beta0_sim = -0.3,  # intercept of linear part
  beta1_sim = -3,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2_sim = 2,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3_sim = 0.2,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern_sim = 0.5,  # probability for Bernoulli distribution in baseline variable
  alpha_sim = 0.05  # significance level for the confidence/prediction interval
){
  

  # generate logistic data with one random variable for simulation
  mydata = sim_ul_group_time_paper(x1_type_sim, x3_type_sim, ni_sim, sigma_sim, beta0_sim, beta1_sim, beta2_sim, beta3_sim, p_bern_sim)
    
  # fit the logistic data by glmm methods in "lme4" package, suppress the warnings
  glmm.laplace = glmer(y ~ x1 + x2 + x3 + (1|id), data=mydata, family="binomial", control=glmerControl(check.conv.singular=.makeCC(action="ignore", tol=1e-4), check.conv.grad=.makeCC("ignore", tol=1e-4))) 

  result <- matrix(ncol = 17, nrow = 4)  # initialize result matrix
  
  
  # get X,Z covariates matrix, i.e. fixed and random design matrix
  Z_raw <- getME(glmm.laplace, "Z")
  X_raw <- getME(glmm.laplace, "X")
  Z <- as.matrix(Z_raw)
  X <- as.matrix(X_raw)
  colnames(X) <- NULL

  # get number of groups from input data 
  n_group = length(unique(X[,3]))*length(unique(X[,4]))  

  # get number of obseravations for each group
  ni <- c(0,0,0,0)
  ni[1] <- sum(X[,3]==1 & X[,4]==0)
  ni[2] <- sum(X[,3]==1 & X[,4]==1)
  ni[3] <- sum(X[,3]==0 & X[,4]==0)
  ni[4] <- sum(X[,3]==0 & X[,4]==1)

  
  # Calculate C matrix for variance of predicted group mean
  weights_working_vector <- weights(glmm.laplace,type="working")
  W_wave <- diag(weights_working_vector)
  nrow(W_wave);ncol(W_wave)
  
  sigma_hat <- getME(glmm.laplace,"theta")[[1]]
  if(sigma_hat == 0) next
  
  G <- diag(rep(sigma_hat,ncol(Z)))
  G_inv = diag(rep(1/sigma_hat,ncol(Z)))
  nrow(G);ncol(G)
  
  C_wave <- cbind(rbind(t(X)%*%W_wave%*%X,t(Z)%*%W_wave%*%X),rbind(t(X)%*%W_wave%*%Z,t(Z)%*%W_wave%*%Z+G_inv))
  nrow(C_wave);ncol(C_wave)

  
  
  
  
  # vector to store coverage prability (whether population group mean inside/outside the interval)
  cp_logitlink = rep(0,n_group)    
  cp_clt = rep(0,n_group)
  cp_lognorm = rep(0,n_group)
  
  # input result table for each group i
  k <- 1  # Initialize pointer to find id / object indicator
  for(i in 1:n_group){
    
    result[i,1] <- ni[i];      # input number of observations  n_obs
    result[i,2] <- i;   # input i
    if(i==1||i==2) {result[i,3] <- 1 } else {result[i,3] <- 0}    # input treat
    if(i==1||i==3) {result[i,4] <- 0 } else {result[i,4] <- 1}    # input time
    
    # get observed y for each group
    y_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),2];  
    result[i,5] <- mean(y_i)   # input Y_bar, mean of true value of Y; observed mean for each group
    
    # get adjusted z_i and z_i for group i
    x_i <- X[k:(k+ni[i]-1),];
    z_i <- Z[k:(k+ni[i]-1),];  
    
    sigma_hat <- getME(glmm.laplace,"theta")[[1]]  # get estimated sigma
    beta_hat_raw <- getME(glmm.laplace,"beta")
    beta_hat <- as.numeric(beta_hat_raw)  # get esimated betas for FE
    b_hat <- getME(glmm.laplace,"b")   # get estimated bi same for all patients in group i
    
    # use estimated beta, b and adjusted X, Z to estimate eta_hat  
    FE_hat_i <- x_i%*%beta_hat   # FE part (beta_hat use the same)
    RE_hat_i <- z_i%*%b_hat    # RE part
    
    # estimate eta(link)
    eta_hat_i <- as.numeric(FE_hat_i + RE_hat_i);  
    
    # get mu_hat for each patient in group i, mu=Eg^{-1}(eta)
    c_hat = (1+0.346*sigma_hat^2)^(-0.5)
    mu_hat_i <- exp(c_hat*x_i %*% beta_hat)/(1+exp(c_hat*x_i %*% beta_hat)) 
    # mu_hat_i <- exp(eta_hat_i) 
    # estimated adjusted group mean, the esimators for our paper
    mu_z_hat_i <- mean(mu_hat_i)
    result[i,6] <- mu_z_hat_i  
    
    # get derivatives of mu_hat by beta 
    partial_mu_hat_beta_hat_i <-  diag(as.vector(exp(c_hat*x_i %*% beta_hat)/(1+exp(c_hat*x_i %*% beta_hat))^2)) %*% x_i * c_hat; 
    # partial_mu_hat_beta_hat_i 
    
    # get derivatives of mu_hat by sigma
    partial_mu_hat_sigma_hat_i <- -0.346 * sigma_hat * c_hat^3 * diag(as.vector(exp(c_hat*x_i %*% beta_hat)/(1+exp(c_hat*x_i %*% beta_hat))^2))  %*% (x_i %*% beta_hat);
    # partial_mu_hat_sigma_hat_i
    
    # get derivatives of mu_hat by psi(parameters)  
    psi_hat <- c(sigma_hat,beta_hat);psi_hat
    partial_mu_hat_psi_hat_i <- cbind(partial_mu_hat_sigma_hat_i,partial_mu_hat_beta_hat_i); partial_mu_hat_psi_hat_i
    partial_mu_z_hat_psi_hat_i <- apply(partial_mu_hat_psi_hat_i,2,mean)  #derivate of adjusted group mean on psi
    
    
    ## extract var-cov matrix of sigma and beta from glmer.nb function
    ## Now pull out the full variance-covariance matrix:
    ## extract Hessian
    hh <- glmm.laplace@optinfo$derivs$Hessian
    ## invert
    vv <- solve(hh)
    ## double/symmetrize (move from deviance to log-likelihood scale)
    Cov_sigma_hat_beta_hat_i <- unname(as.matrix(forceSymmetric(vv + t(vv)))) # var-cov matrix of sigma and beta: 1st col is sigma
    
    n <- ni[i]  # get n for each group
    
    # Use delta method to calculate Var of estimated adjusted group mean
    Var_mu_z_hat_i <- as.numeric(partial_mu_z_hat_psi_hat_i%*% Cov_sigma_hat_beta_hat_i %*%partial_mu_z_hat_psi_hat_i) ; Var_mu_z_hat_i
    
    SD_mu_z_hat_i <- Var_mu_z_hat_i^0.5   # Calculated Standard deviation of mu_z_hat
    
    result[i,8] <- Var_mu_z_hat_i  # input variance of mu_z_hat 
    result[i,9] <- SD_mu_z_hat_i   # input SD of mu_z_hat 
    
    alpha_sig <- alpha_sim   # significant level of CLT method
    Z_alpha <- qnorm(1-alpha_sig/2);Z_alpha 
    
    #######Calculation of CI###############  
    
    ##generate confidence interval of mu_z_hat for logit link function: by formulas in paper
    CI_mu_z_hat_lower_i_logitlink <- (mu_z_hat_i/(1-mu_z_hat_i))*exp(-(Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i)))/(1+(mu_z_hat_i/(1-mu_z_hat_i))*exp(-(Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i))))     #lower bound of CI
    CI_mu_z_hat_upper_i_logitlink <- (mu_z_hat_i/(1-mu_z_hat_i))*exp((Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i)))/(1+(mu_z_hat_i/(1-mu_z_hat_i))*exp((Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i))))     #upper bound of CI
    
    
    # CLT method for CI
    CI_mu_z_hat_lower_i_clt <- mu_z_hat_i - Z_alpha * SD_mu_z_hat_i
    CI_mu_z_hat_upper_i_clt <- mu_z_hat_i + Z_alpha * SD_mu_z_hat_i
    
    ## no lognormal for this, set 0 instead
    CI_mu_z_hat_lower_i_lognorm <- 0
    CI_mu_z_hat_upper_i_lognorm <- 0
    
    result[i,10] <- CI_mu_z_hat_lower_i_clt   # input lower bound of CI 
    result[i,11] <- CI_mu_z_hat_upper_i_clt   # input upper bound of CI
    
    # generate mu_z_star_hat, i.e estimator in formula (2) of Qu and Luo's paper
    x1_i <- mean(X[k:(k+ni[i]-1),2])    # mean of x1(baseline) for group_i
    x2_i <- result[i,3]  # group  for i
    x3_i <- result[i,4]  # time for i
    
    mu_z_star_hat <- exp(c_hat*c(1,x1_i,x2_i,x3_i) %*% beta_hat)/(1+exp(c_hat*c(1,x1_i,x2_i,x3_i) %*% beta_hat)) # eta->mu for logit link 
    result[i,7] <- mu_z_star_hat   # input mu_z_star_hat
    
    # Calculate the population value of mu_l (i.e expectation)
    eta_ij <- x_i %*%  c(beta0_sim, beta1_sim, beta2_sim, beta3_sim)
    
    c = (1+0.346*sigma_sim^2)^(-0.5)
    mu_E_ij  <- exp(c*eta_ij)/(1+exp(c*eta_ij))   #population mu
    mu_l <- mean(mu_E_ij)    #population group mean for mu
    result[i,13] <- mu_l  # input mu_l
    
    # Calculate the population value of mu_bar (i.e average of each mu)
    mu_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),7];  # observed y for each group
    mu_bar_i  <- mean(mu_i)     # mu_bar : mean of each mu_i
    result[i,12] <- mu_bar_i   # input mu_bar
    
    if((CI_mu_z_hat_lower_i_logitlink <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_logitlink)){ 
      cp_logitlink[i] <- cp_logitlink[i]+1      # count number of Y_bar inCI for group+i
    } 
    
    if((CI_mu_z_hat_lower_i_clt <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_clt)){ 
      cp_clt[i] <- cp_clt[i]+1      # count number of Y_bar inCI for group+i
    }
    
    if((CI_mu_z_hat_lower_i_lognorm <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_lognorm)){ 
      cp_lognorm[i] <- cp_lognorm[i]+1      # count number of Y_bar inCI for group+i
    }
    
    k <- k + ni[i] 
  }
  
  result[, 14]=mean(result[, 6])
  result[, 15]=mean(result[, 13])
  result[, 16]=mean(result[, 5])
  result[, 17]=mean(result[, 12])
  
  # input names for each column of the table
  colnames(result) <- c("n_obs","group","treat","time", "Y_bar","mu_z_hat","mu_z_star_hat", "Var_mu_z_hat","SD_mu_z_hat","CI_lower","CI_upper", "mu_bar", "mu_l", "mu_z_hat_total", "mu_l_total", "y_bar_total","mu_bar_total")
  
  result = data.frame(result)
  
  #print(result,digits=5)
  
  cp_table = cbind(result$n_obs, result$group, result$treat, result$time)
  
  # generate output for estimation of paramter
  sigma_hat = getME(glmm.laplace,"theta")[[1]]  # get estimated sigma for RE
  sigma_hat = round(sigma_hat,digits=5)  # only output 5 digits
  
  beta_hat = as.numeric(getME(glmm.laplace, "beta"))  # get esimated betas for FE
  beta_hat = round(beta_hat, digits=5)
  
  beta_hat_output = paste(beta_hat, collapse=" ")# change vector into paste string form for output  
  
  group = result$group
  treat = result$treat
  time = result$time
  mu = result$mu_l   # true group mean
  mu = round(mu, digits=5)
  mu_z_hat = result$mu_z_hat   # group mean/ marginal mean by adjustment 
  mu_z_hat = round(mu_z_hat, digits=5)
  y_bar = result$Y_bar   # group mean /marginal mean
  y_bar = round(y_bar, digits=5)
  mu_z_star_hat = result$mu_z_star_hat  # repsonse at mean value
  mu_z_star_hat = round(mu_z_star_hat, digits=5)
  
  # save all results 
  z = list()
  #Z$model_fit = glmm.laplace
  z$result = result
  z$x1_type_sim = x1_type_sim
  z$x3_type_sim = x3_type_sim
  z$ni_sim = ni_sim
  z$sigma_sim = ni_sim
  z$beta0_sim = beta0_sim
  z$beta1_sim = beta1_sim
  z$beta2_sim = beta2_sim
  z$beta3_sim = beta3_sim
  z$p_bern_sim = p_bern_sim
  z$sigma_hat = sigma_hat
  z$beta_hat_output = beta_hat_output
  z$group = group
  z$treat = treat
  z$time = time
  z$mu = mu
  z$mu_z_hat = mu_z_hat
  z$y_bar = y_bar
  z$mu_z_star_hat = mu_z_star_hat
  z$cp_logitlink = cp_logitlink
  z$cp_clt = cp_clt
  z$cp_lognorm = cp_lognorm
  z$CI_mu_z_hat_lower_i_logitlink = CI_mu_z_hat_lower_i_logitlink 
  z$CI_mu_z_hat_upper_i_logitlink = CI_mu_z_hat_upper_i_logitlink 
  z$CI_mu_z_hat_lower_i_clt = CI_mu_z_hat_lower_i_clt
  z$CI_mu_z_hat_upper_i_clt = CI_mu_z_hat_upper_i_clt
  z$CI_mu_z_hat_lower_i_lognorm = CI_mu_z_hat_lower_i_lognorm
  z$CI_mu_z_hat_upper_i_lognorm = CI_mu_z_hat_upper_i_lognorm
  
  return(z)
}



########################################################################################




ul_pi <- function(
  x1_type_sim = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type_sim = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni_sim = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma_sim = 0.5,  # standard deviation for random effects
  size_sim = 50,  # size for negative binomial data; for logistic data, no need to input this value
  beta0_sim = -0.3,  # intercept of linear part
  beta1_sim = -3,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2_sim = 2,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3_sim = 0.2,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern_sim = 0.5,  # probability for Bernoulli distribution in baseline variable
  alpha_sim = 0.05  # significance level for the confidence/prediction interval
){

  # generate logistic data with one random variable for simulation
  mydata = sim_ul_group_time_paper(x1_type_sim, x3_type_sim, ni_sim, sigma_sim, beta0_sim, beta1_sim, beta2_sim, beta3_sim, p_bern_sim)
  
  # fit the logistic data by glmm methods in "lme4" package, suppress the warnings
  glmm.laplace = glmer(y ~ x1 + x2 + x3 + (1|id), data=mydata, family="binomial", control=glmerControl(check.conv.singular=.makeCC(action="ignore", tol=1e-4), check.conv.grad=.makeCC("ignore", tol=1e-4))) 

  result <- matrix(ncol = 17, nrow = 4)  # initialize result matrix

  # generate X,Z matrix, i.e. fixed and random design matrix
  Z_raw <- getME(glmm.laplace, "Z")
  X_raw <- getME(glmm.laplace, "X")
  Z <- as.matrix(Z_raw)
  X <- as.matrix(X_raw)
  colnames(X) <- NULL
  #nrow(X);ncol(X);nrow(Z);ncol(Z)
  
  n_group=length(unique(X[,3]))*length(unique(X[,4]))  
  
  #special for unbalanced case########
  ni <- c(0,0,0,0)
  ni[1] <- sum(X[,3]==1 & X[,4]==0)
  ni[2] <- sum(X[,3]==1 & X[,4]==1)
  ni[3] <- sum(X[,3]==0 & X[,4]==0)
  ni[4] <- sum(X[,3]==0 & X[,4]==1)
  ni
  ####################################
  
  
  
  ###################pi######################################  
  
  weights_working_vector <- weights(glmm.laplace,type="working")
  W_wave <- diag(weights_working_vector)
  nrow(W_wave);ncol(W_wave)
  
  sigma_hat <- getME(glmm.laplace,"theta")[[1]]
  if(sigma_hat == 0) next
  
  G <- diag(rep(sigma_hat,ncol(Z)))
  G_inv = diag(rep(1/sigma_hat,ncol(Z)))
  nrow(G);ncol(G)
  
  C_wave <- cbind(rbind(t(X)%*%W_wave%*%X,t(Z)%*%W_wave%*%X),rbind(t(X)%*%W_wave%*%Z,t(Z)%*%W_wave%*%Z+G_inv))
  nrow(C_wave);ncol(C_wave)
  
  #########pi###############################################  
  
  
  
  cp_logitlink = rep(0,n_group)
  cp_clt = rep(0,n_group)
  cp_lognorm = rep(0,n_group)
  
  # input result table for each group i
  k <- 1  # Initialize pointer to find id / object indicator
  for(i in 1:4){
    
    result[i,1] <- ni[i];      # input number of observations  n_obs
    result[i,2] <- i;   # input i
    if(i==1||i==2) {result[i,3] <- 1 } else {result[i,3] <- 0}    # input treat
    if(i==1||i==3) {result[i,4] <- 0 } else {result[i,4] <- 1}    # input time
    
    # get observed y for each group
    y_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),2];  
    result[i,5] <- mean(y_i)   # input Y_bar, mean of true value of Y; observed mean for each group
    
    # get adjusted z_i and z_i for group i
    x_i <- X[k:(k+ni[i]-1),];
    z_i <- Z[k:(k+ni[i]-1),];  
    
    
    #####################pi############
    
    
    L_i <- rbind(t(x_i),t(z_i));   #generate L
    #nrow(L_i);ncol(L_i)
    
    CMSEP_eta_hat_i <- t(L_i)%*%solve(C_wave)%*%L_i   #generate variance for eta
    #nrow(CMSEP_eta_hat_i);ncol(CMSEP_eta_hat_i)
    
    mu_hat_i <- getME(glmm.laplace,"mu")[k:(k+ni[i]-1)]  #y_hat or mu_hat: estimation of mu
    eta_hat_i <- log(mu_hat_i/(1-mu_hat_i))
    mu_z_hat_i <- mean(mu_hat_i)
    g_inv_deriv_eta_hat_i <- exp(eta_hat_i)/(1+exp(eta_hat_i))^2  #g-1()
    
    result[i,6] <- mu_z_hat_i   #input mu_z_hat 
    
    n <- ni[i]  #get n for each group
    
    Var_mu_z_hat_i <- 1/(ni[i]^2)*g_inv_deriv_eta_hat_i%*%CMSEP_eta_hat_i%*%g_inv_deriv_eta_hat_i
    SD_mu_z_hat_i <- Var_mu_z_hat_i[1,1]^0.5
    
    result[i,8] <- Var_mu_z_hat_i   #input variance of mu_z_hat 
    result[i,9] <- SD_mu_z_hat_i   #input SD of mu_z_hat 
    
    beta_hat_raw <- getME(glmm.laplace,"beta")
    beta_hat <- as.numeric(beta_hat_raw)  # get esimated betas for FE

    alpha_sig <- alpha_sim   # significant level of CLT method
    Z_alpha <- qnorm(1-alpha_sig/2);Z_alpha 
    
    #######Calculation of CI###############  
    
    ##generate confidence interval of mu_z_hat for logit link function: by formulas in paper
    CI_mu_z_hat_lower_i_logitlink <- (mu_z_hat_i/(1-mu_z_hat_i))*exp(-(Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i)))/(1+(mu_z_hat_i/(1-mu_z_hat_i))*exp(-(Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i))))     #lower bound of CI
    CI_mu_z_hat_upper_i_logitlink <- (mu_z_hat_i/(1-mu_z_hat_i))*exp((Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i)))/(1+(mu_z_hat_i/(1-mu_z_hat_i))*exp((Z_alpha*SD_mu_z_hat_i)/(mu_z_hat_i*(1-mu_z_hat_i))))     #upper bound of CI
    
    # CLT method for CI
    CI_mu_z_hat_lower_i_clt <- mu_z_hat_i - Z_alpha * SD_mu_z_hat_i
    CI_mu_z_hat_upper_i_clt <- mu_z_hat_i + Z_alpha * SD_mu_z_hat_i
    
    # no lognormal for this, set 0 instead
    CI_mu_z_hat_lower_i_lognorm <- 0
    CI_mu_z_hat_upper_i_lognorm <- 0
    
    result[i,10] <- CI_mu_z_hat_lower_i_clt   # input lower bound of CI 
    result[i,11] <- CI_mu_z_hat_upper_i_clt   # input upper bound of CI
    
    # generate mu_z_star_hat, i.e estimator in formula (2) of Qu and Luo's paper
    x1_i <- mean(X[k:(k+ni[i]-1),2])    # mean of x1(baseline) for group_i
    x2_i <- result[i,3]  # group  for i
    x3_i <- result[i,4]  # time for i
    
    if(x3_type_sim==1){
      b_hat <- as.numeric(getME(glmm.laplace,"b"))[k:(k+ni[i]-1)]
    }
    if(x3_type_sim==2){
      if((k==1)|(k==2)){
        b_hat <- as.numeric(getME(glmm.laplace,"b"))[1:ni[i]]
      }
      if((k==3)|(k==4)){
        b_hat <- as.numeric(getME(glmm.laplace,"b"))[(ni[1]+1):(ni[1]+ni[i])]
      }
    }
    mu_z_star_hat <- mean(exp(as.numeric(c(1,x1_i,x2_i,x3_i)%*%beta_hat)+b_hat)/(1+exp(as.numeric(c(1,x1_i,x2_i,x3_i)%*%beta_hat)+b_hat)))
    result[i,7] <- mu_z_star_hat   # input mu_z_star_hat
    
    ##############pi###################    
    
    mu_l  <- mean(mydata[k:(k+ni[i]-1),7])     # mu_bar : mean of each mu_i
    result[i,13] <- mu_l  # input mu_l
    # Y_bar <- mean(y_i)
    
    ###############pi##############   
    
    # Calculate the population value of mu_bar (i.e average of each mu)
    mu_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),7];  # observed y for each group
    mu_bar_i  <- mean(mu_i)     # mu_bar : mean of each mu_i
    result[i,12] <- mu_bar_i   # input mu_bar
    
    if((CI_mu_z_hat_lower_i_logitlink <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_logitlink)){ 
      cp_logitlink[i] <- cp_logitlink[i]+1      # count number of Y_bar inCI for group+i
    } 
    
    if((CI_mu_z_hat_lower_i_clt <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_clt)){ 
      cp_clt[i] <- cp_clt[i]+1      # count number of Y_bar inCI for group+i
    }
    
    if((CI_mu_z_hat_lower_i_lognorm <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_lognorm)){ 
      cp_lognorm[i] <- cp_lognorm[i]+1      # count number of Y_bar inCI for group+i
    }
    
    k <- k + ni[i] 
  }
  
  result[, 14]=mean(result[, 6])
  result[, 15]=mean(result[, 13])
  result[, 16]=mean(result[, 5])
  result[, 17]=mean(result[, 12])
  
  # input names for each column of the table
  colnames(result) <- c("n_obs","group","treat","time", "Y_bar","mu_z_hat","mu_z_star_hat", "Var_mu_z_hat","SD_mu_z_hat","CI_lower","CI_upper", "mu_bar", "mu_l", "mu_z_hat_total", "mu_l_total", "y_bar_total","mu_bar_total")
  
  result = data.frame(result)
  
  #print(result,digits=5)
  
  cp_table = cbind(result$n_obs, result$group, result$treat, result$time)
  
  # generate output for estimation of paramter
  sigma_hat = getME(glmm.laplace,"theta")[[1]]  # get estimated sigma for RE
  sigma_hat = round(sigma_hat,digits=5)  # only output 5 digits
  
  beta_hat = as.numeric(getME(glmm.laplace, "beta"))  # get esimated betas for FE
  beta_hat = round(beta_hat, digits=5)
  
  beta_hat_output = paste(beta_hat, collapse=" ")# change vector into paste string form for output  
  
  
  group = result$group
  treat = result$treat
  time = result$time
  mu = result$mu_l   # true group mean
  mu = round(mu, digits=5)
  mu_z_hat = result$mu_z_hat   # group mean/ marginal mean by adjustment 
  mu_z_hat = round(mu_z_hat, digits=5)
  y_bar = result$Y_bar   # group mean /marginal mean
  y_bar = round(y_bar, digits=5)
  mu_z_star_hat = result$mu_z_star_hat  # repsonse at mean value
  mu_z_star_hat = round(mu_z_star_hat, digits=5)

  # save all results 
  z = list()
  #Z$model_fit = glmm.laplace
  z$result = result
  z$x1_type_sim = x1_type_sim
  z$x3_type_sim = x3_type_sim
  z$ni_sim = ni_sim
  z$sigma_sim = ni_sim
  z$beta0_sim = beta0_sim
  z$beta1_sim = beta1_sim
  z$beta2_sim = beta2_sim
  z$beta3_sim = beta3_sim
  z$p_bern_sim = p_bern_sim
  z$sigma_hat = sigma_hat
  z$beta_hat_output = beta_hat_output
  z$group = group
  z$treat = treat
  z$time = time
  z$mu = mu
  z$mu_z_hat = mu_z_hat
  z$y_bar = y_bar
  z$mu_z_star_hat = mu_z_star_hat
  z$cp_logitlink = cp_logitlink
  z$cp_clt = cp_clt
  z$cp_lognorm = cp_lognorm
  z$CI_mu_z_hat_lower_i_logitlink = CI_mu_z_hat_lower_i_logitlink 
  z$CI_mu_z_hat_upper_i_logitlink = CI_mu_z_hat_upper_i_logitlink 
  z$CI_mu_z_hat_lower_i_clt = CI_mu_z_hat_lower_i_clt
  z$CI_mu_z_hat_upper_i_clt = CI_mu_z_hat_upper_i_clt
  z$CI_mu_z_hat_lower_i_lognorm = CI_mu_z_hat_lower_i_lognorm
  z$CI_mu_z_hat_upper_i_lognorm = CI_mu_z_hat_upper_i_lognorm
  
  return(z)
}

  

##############################################################################################




unb_ci <- function(
  x1_type_sim = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type_sim = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni_sim = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma_sim = 0.1,  # standard deviation for random effects
  size_sim = 50,  # size for negative binomial data; for logistic data, no need to input this value
  beta0_sim = 0.3,  # intercept of linear part
  beta1_sim = -0.2,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2_sim = 0.3,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3_sim = 0.4,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern_sim = 0.5,  # probability for Bernoulli distribution in baseline variable
  alpha_sim = 0.05  # significance level for the confidence/prediction interval
){
  
  # generate logistic data with one random variable for simulation
  mydata = sim_unb_group_time_paper(x1_type_sim, x3_type_sim, ni_sim, sigma_sim, size_sim, beta0_sim, beta1_sim, beta2_sim, beta3_sim, p_bern_sim)
  
  # fit the logistic data by glmm methods in "lme4" package, suppress the warnings
  glmm.laplace = glmer.nb(y ~ x1 + x2 + x3 + (1|id), data=mydata, verbose=FALSE, control=glmerControl(check.conv.singular=.makeCC(action="ignore", tol=1e-4), check.conv.grad=.makeCC("ignore", tol=1e-4)))   #,verbose=TRUE is used to show the iteration, can be closed for CP simulation

  result <- matrix(ncol = 17, nrow = 4)  # initialize result matrix
  
  
  # generate X,Z matrix, i.e. fixed and random design matrix
  Z_raw <- getME(glmm.laplace, "Z")
  X_raw <- getME(glmm.laplace, "X")
  Z <- as.matrix(Z_raw)
  X <- as.matrix(X_raw)
  colnames(X) <- NULL
  #nrow(X);ncol(X);nrow(Z);ncol(Z)
  
  n_group=length(unique(X[,3]))*length(unique(X[,4]))  #In case num of groups>4
  
  #special for unbalanced case
  ni <- c(0,0,0,0)
  ni[1] <- sum(X[,3]==1 & X[,4]==0)
  ni[2] <- sum(X[,3]==1 & X[,4]==1)
  ni[3] <- sum(X[,3]==0 & X[,4]==0)
  ni[4] <- sum(X[,3]==0 & X[,4]==1)
  ni
  
  ####################################
  
  sigma_hat <- getME(glmm.laplace,"theta")[[1]]  # get estimated sigma
  if(sigma_hat == 0) next   #remove the bug case of sigma_hat=0
  
  cp_loglink = rep(0,n_group)
  cp_clt = rep(0,n_group)
  cp_lognorm = rep(0,n_group)
  
  # input result table for each group i
  k <- 1  # Initialize pointer to find id / object indicator
  for(i in 1:4){
    
    result[i,1] <- ni[i];      # input number of observations  n_obs
    result[i,2] <- i;   # input i
    if(i==1||i==2) {result[i,3] <- 1 } else {result[i,3] <- 0}    # input treat
    if(i==1||i==3) {result[i,4] <- 0 } else {result[i,4] <- 1}    # input time
    
    # get observed y for each group
    y_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),2];  
    result[i,5] <- mean(y_i)   # input Y_bar, mean of true value of Y; observed mean for each group
    
    # get adjusted z_i and z_i for group i
    x_i <- X[k:(k+ni[i]-1),];
    z_i <- Z[k:(k+ni[i]-1),];  
    
    beta_hat_raw <- getME(glmm.laplace,"beta")
    beta_hat <- as.numeric(beta_hat_raw)  # get esimated betas for FE
    b_hat <- getME(glmm.laplace,"b")   # get estimated bi same for all patients in group i
    
    # use estimated beta, b and adjusted X, Z to estimate eta_hat  
    FE_hat_i <- x_i%*%beta_hat   # FE part (beta_hat use the same)
    RE_hat_i <- z_i%*%b_hat    # RE part
    
    # estimate eta(link)
    eta_hat_i <- as.numeric(FE_hat_i + RE_hat_i);  
    
    # get mu_hat for each patient in group i, mu=Eg^{-1}(eta)
    mu_hat_i <- exp(x_i %*% beta_hat + 0.5*sigma_hat^2) 
    # mu_hat_i <- exp(eta_hat_i) 
    # estimated adjusted group mean, the esimators for our paper
    mu_z_hat_i <- mean(mu_hat_i)
    result[i,6] <- mu_z_hat_i  
    
    # get derivatives of mu_hat by beta_hat 
    partial_mu_hat_beta_hat_i <- diag(as.vector(mu_hat_i)) %*% x_i; 
    
    # get derivatives of mu_hat by sigma_hat
    partial_mu_hat_sigma_hat_i <- sigma_hat * mu_hat_i;
    
    # get derivatives of mu_hat by psi(parameters)  
    psi_hat <- c(sigma_hat,beta_hat);psi_hat
    partial_mu_hat_psi_hat_i <- cbind(partial_mu_hat_sigma_hat_i,partial_mu_hat_beta_hat_i); partial_mu_hat_psi_hat_i
    partial_mu_z_hat_psi_hat_i <- apply(partial_mu_hat_psi_hat_i,2,mean)  #derivate of adjusted group mean on psi
    
    
    ## extract var-cov matrix of sigma and beta from glmer.nb function
    ## Now pull out the full variance-covariance matrix:
    ## extract Hessian
    hh <- glmm.laplace@optinfo$derivs$Hessian
    ## invert
    vv <- solve(hh)
    ## double/symmetrize (move from deviance to log-likelihood scale)
    Cov_sigma_hat_beta_hat_i <- unname(as.matrix(forceSymmetric(vv + t(vv)))) # var-cov matrix of sigma and beta: 1st col is sigma
    
    n <- ni[i]  #get n for each group
    
    # Sum of lognormal approximation to calculate Variance of mu_z_hat
    mu_x_sn_i <- x_i %*% beta_hat + 0.5*sigma_hat^2 
    partial_mu_x_sn_sigma_beta <- cbind(sigma_hat,x_i)
    M_x_sn_i <- as.matrix(partial_mu_x_sn_sigma_beta  %*%   Cov_sigma_hat_beta_hat_i  %*%  t(partial_mu_x_sn_sigma_beta))
    Var_mu_z_hat_i <- (1/ni[i]^2) * as.numeric( t(exp(mu_x_sn_i+ 0.5 *diag(M_x_sn_i )))  %*% (exp(M_x_sn_i)-1) %*% exp(mu_x_sn_i+ 0.5 *diag(M_x_sn_i)))
    
    
    SD_mu_z_hat_i <- Var_mu_z_hat_i^0.5   # Calculated Standard deviation of mu_z_hat
    
    result[i,8] <- Var_mu_z_hat_i  # input variance of mu_z_hat 
    result[i,9] <- SD_mu_z_hat_i   # input SD of mu_z_hat 
    
    alpha_sig <- alpha_sim   # significant level of CLT method
    Z_alpha <- qnorm(1-alpha_sig/2);Z_alpha 
    
    #######Calculation of CI###############  
    
    ## generate confidence interval of mu_z_hat for log link function: by formulas in paper
    CI_mu_z_hat_lower_i_loglink <- mu_z_hat_i*exp(-Z_alpha*SD_mu_z_hat_i/mu_z_hat_i) #lower bound of CI
    CI_mu_z_hat_upper_i_loglink <- mu_z_hat_i*exp(Z_alpha*SD_mu_z_hat_i/mu_z_hat_i)  #upper bound of CI
    
    # CLT method for CI
    CI_mu_z_hat_lower_i_clt <- mu_z_hat_i - Z_alpha * SD_mu_z_hat_i
    CI_mu_z_hat_upper_i_clt <- mu_z_hat_i + Z_alpha * SD_mu_z_hat_i
    
    # lognormal approximation method for CI(sum of depended lognormals)
    mu_lognormal_i <- log(ni[i] * mu_z_hat_i) - 0.5* log (Var_mu_z_hat_i/mu_z_hat_i^2+1)   #calculated mu for logmormal distribution
    sd_lognormal_i <- (log (Var_mu_z_hat_i/mu_z_hat_i^2+1))^0.5   #calculated variance for logmormal distribution
    
    CI_mu_z_hat_lower_i_lognorm <-(1/ni[i])* qlnorm(0.025, meanlog = mu_lognormal_i, sdlog = sd_lognormal_i, lower.tail = TRUE, log.p = FALSE)
    CI_mu_z_hat_upper_i_lognorm <-(1/ni[i])* qlnorm(0.975, meanlog = mu_lognormal_i, sdlog = sd_lognormal_i, lower.tail = TRUE, log.p = FALSE)  
    
    result[i,10] <- CI_mu_z_hat_lower_i_clt   # input lower bound of CI 
    result[i,11] <- CI_mu_z_hat_upper_i_clt   # input upper bound of CI
    
    # generate mu_z_star_hat, i.e estimator in formula (2) of Qu and Luo's paper
    x1_i <- mean(X[k:(k+ni[i]-1),2])    # mean of x1(baseline) for group_i
    x2_i <- result[i,3]  # group  for i
    x3_i <- result[i,4]  # time for i
    
    mu_z_star_hat <- exp(c(1,x1_i,x2_i,x3_i) %*% beta_hat+0.5*sigma_hat^2) 
    result[i,7] <- mu_z_star_hat   # input mu_z_star_hat
    
    # Calculate the population value of mu_l (i.e expectation)
    eta_ij <- x_i %*%  c(beta0_sim, beta1_sim, beta2_sim, beta3_sim)
    
    mu_E_ij <- exp(eta_ij + 0.5*sigma_sim^2)   #population mu
    mu_l <- mean(mu_E_ij)  #population group mean for mu
    result[i,13] <- mu_l  # input mu_l
    
    # Calculate the population value of mu_bar (i.e average of each mu)
    mu_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),7];  # observed y for each group
    mu_bar_i  <- mean(mu_i)     # mu_bar : mean of each mu_i
    result[i,12] <- mu_bar_i   # input mu_bar
    
    if((CI_mu_z_hat_lower_i_loglink <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_loglink)){ 
      cp_loglink[i] <- cp_loglink[i]+1      # count number of Y_bar inCI for group+i
    } 
    
    if((CI_mu_z_hat_lower_i_clt <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_clt)){ 
      cp_clt[i] <- cp_clt[i]+1      # count number of Y_bar inCI for group+i
    }
    
    if((CI_mu_z_hat_lower_i_lognorm <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_lognorm)){ 
      cp_lognorm[i] <- cp_lognorm[i]+1      # count number of Y_bar inCI for group+i
    }
    
    k <- k + ni[i] 
  }
  
  result[, 14]=mean(result[, 6])
  result[, 15]=mean(result[, 13])
  result[, 16]=mean(result[, 5])
  result[, 17]=mean(result[, 12])
  
  # input names for each column of the table
  colnames(result) <- c("n_obs","group","treat","time", "Y_bar","mu_z_hat","mu_z_star_hat", "Var_mu_z_hat","SD_mu_z_hat","CI_lower","CI_upper", "mu_bar", "mu_l", "mu_z_hat_total", "mu_l_total", "y_bar_total","mu_bar_total")
  
  result = data.frame(result)
  
  cp_table = cbind(result$n_obs, result$group, result$treat, result$time)
  
  # generate output for estimation of paramter
  sigma_hat = getME(glmm.laplace,"theta")[[1]]  # get estimated sigma for RE
  sigma_hat = round(sigma_hat,digits=5)  # only output 5 digits
  
  size_hat = getME(glmm.laplace,"glmer.nb.theta") # estimated size parameter of NB
  size_hat = round(size_hat,digits=5)
  
  beta_hat = as.numeric(getME(glmm.laplace, "beta"))  # get esimated betas for FE
  beta_hat = round(beta_hat, digits=5)
  
  beta_hat_output = paste(beta_hat, collapse=" ")# change vector into paste string form for output  
  
  
  group = result$group
  treat = result$treat
  time = result$time
  mu = result$mu_l   # true group mean
  mu = round(mu, digits=5)
  mu_z_hat = result$mu_z_hat   # group mean/ marginal mean by adjustment 
  mu_z_hat = round(mu_z_hat, digits=5)
  y_bar = result$Y_bar   # group mean /marginal mean
  y_bar = round(y_bar, digits=5)
  mu_z_star_hat = result$mu_z_star_hat  # repsonse at mean value
  mu_z_star_hat = round(mu_z_star_hat, digits=5)
  
  # save all results 
  z = list()
  #Z$model_fit = glmm.laplace
  z$result = result
  z$x1_type_sim = x1_type_sim
  z$x3_type_sim = x3_type_sim
  z$ni_sim = ni_sim
  z$sigma_sim = ni_sim
  z$beta0_sim = beta0_sim
  z$beta1_sim = beta1_sim
  z$beta2_sim = beta2_sim
  z$beta3_sim = beta3_sim
  z$p_bern_sim = p_bern_sim
  z$sigma_hat = sigma_hat
  z$beta_hat_output = beta_hat_output
  z$group = group
  z$treat = treat
  z$time = time
  z$mu = mu
  z$mu_z_hat = mu_z_hat
  z$y_bar = y_bar
  z$mu_z_star_hat = mu_z_star_hat
  z$cp_loglink = cp_loglink
  z$cp_clt = cp_clt
  z$cp_lognorm = cp_lognorm
  z$CI_mu_z_hat_lower_i_loglink = CI_mu_z_hat_lower_i_loglink 
  z$CI_mu_z_hat_upper_i_loglink = CI_mu_z_hat_upper_i_loglink 
  z$CI_mu_z_hat_lower_i_clt = CI_mu_z_hat_lower_i_clt
  z$CI_mu_z_hat_upper_i_clt = CI_mu_z_hat_upper_i_clt
  z$CI_mu_z_hat_lower_i_lognorm = CI_mu_z_hat_lower_i_lognorm
  z$CI_mu_z_hat_upper_i_lognorm = CI_mu_z_hat_upper_i_lognorm
  
  return(z)
}



#############################################################################

unb_pi <- function(
  x1_type_sim = 1,  # =1, x1~Bernoulli(p); =2, x1~Unif[0,1]
  x3_type_sim = 2,  # =1, x3 is gender, different observation at different X3=0 or 1; =2, x3 is time, same observation at different time
  ni_sim = c(200,180,200,160),  # a vector to store number of observation in each group
  sigma_sim = 0.1,  # standard deviation for random effects
  size_sim = 50,  # size for negative binomial data; for logistic data, no need to input this value
  beta0_sim = 0.3,  # intercept of linear part
  beta1_sim = -0.2,  # slope for baseline variable: x1_type=1: x1~Bernoulli(p_bern); x1_type=2: x1~Uniform(0,1)
  beta2_sim = 0.3,  # slope for treat indicator: x2=0: control group; x2=1: treatment group; 
  beta3_sim = 0.4,  # slope for time indicator: x1_type=1: time 0 or 1; x1_type=2: gender 0 or 1
  p_bern_sim = 0.5,  # probability for Bernoulli distribution in baseline variable
  alpha_sim = 0.05  # significance level for the confidence/prediction interval
){
  
  # generate logistic data with one random variable for simulation
  mydata = sim_unb_group_time_paper(x1_type_sim, x3_type_sim, ni_sim, sigma_sim, size_sim, beta0_sim, beta1_sim, beta2_sim, beta3_sim, p_bern_sim)
  
  # fit the logistic data by glmm methods in "lme4" package, suppress the warnings
  glmm.laplace = glmer.nb(y ~ x1 + x2 + x3 + (1|id), data=mydata, verbose=FALSE, control=glmerControl(check.conv.singular=.makeCC(action="ignore", tol=1e-1), check.conv.grad=.makeCC("ignore", tol=1e-1)))  
  
  result <- matrix(ncol = 17, nrow = 4)  # initialize result matrix
  
  
  # generate X,Z matrix, i.e. fixed and random design matrix
  Z_raw <- getME(glmm.laplace, "Z")
  X_raw <- getME(glmm.laplace, "X")
  Z <- as.matrix(Z_raw)
  X <- as.matrix(X_raw)
  colnames(X) <- NULL
  #nrow(X);ncol(X);nrow(Z);ncol(Z)
  
  n_group=length(unique(X[,3]))*length(unique(X[,4]))  #In case num of groups>4
  
  #special for unbalanced case
  ni <- c(0,0,0,0)
  ni[1] <- sum(X[,3]==1 & X[,4]==0)
  ni[2] <- sum(X[,3]==1 & X[,4]==1)
  ni[3] <- sum(X[,3]==0 & X[,4]==0)
  ni[4] <- sum(X[,3]==0 & X[,4]==1)
  ni
  
  ####################################
  
  ###################pi######################################    
  
  weights_working_vector <- weights(glmm.laplace,type="working")
  W_wave <- diag(weights_working_vector)
  nrow(W_wave);ncol(W_wave)
  
  sigma_hat <- getME(glmm.laplace,"theta")[[1]]
  #if(sigma_hat <= sigma_sim/10) next
  if(sigma_hat == 0) next
  
  #size_hat = getME(glmm.laplace,"glmer.nb.theta")  #new
  #if(size_hat <= size_sim*20) next   #new
  
  G <- diag(rep(sigma_hat,ncol(Z)))
  G_inv = diag(rep(1/sigma_hat,ncol(Z)))
  nrow(G);ncol(G)
  
  C_wave <- cbind(rbind(t(X)%*%W_wave%*%X,t(Z)%*%W_wave%*%X),rbind(t(X)%*%W_wave%*%Z,t(Z)%*%W_wave%*%Z+G_inv))
  nrow(C_wave);ncol(C_wave)
  
  ###################pi######################################     
  
  cp_loglink = rep(0,n_group)
  cp_clt = rep(0,n_group)
  cp_lognorm = rep(0,n_group)
  
  # input result table for each group i
  k <- 1  # Initialize pointer to find id / object indicator
  for(i in 1:4){
    
    result[i,1] <- ni[i];      # input number of observations  n_obs
    result[i,2] <- i;   # input i
    if(i==1||i==2) {result[i,3] <- 1 } else {result[i,3] <- 0}    # input treat
    if(i==1||i==3) {result[i,4] <- 0 } else {result[i,4] <- 1}    # input time
    
    # get observed y for each group
    y_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),2];  
    result[i,5] <- mean(y_i)   # input Y_bar, mean of true value of Y; observed mean for each group
    
    # get adjusted z_i and z_i for group i
    x_i <- X[k:(k+ni[i]-1),];
    z_i <- Z[k:(k+ni[i]-1),];  
    
    
    #####################pi############
    
    
    L_i <- rbind(t(x_i),t(z_i));   #generate L
    #nrow(L_i);ncol(L_i)
    
    CMSEP_eta_hat_i <- t(L_i)%*%solve(C_wave)%*%L_i   #generate variance for eta
    #nrow(CMSEP_eta_hat_i);ncol(CMSEP_eta_hat_i)
    
    mu_hat_i <- getME(glmm.laplace,"mu")[k:(k+ni[i]-1)]  #y_hat or mu_hat: estimation of mu
    eta_hat_i <- log(mu_hat_i)   #change!!
    mu_z_hat_i <- mean(mu_hat_i)
    g_inv_deriv_eta_hat_i <- exp(eta_hat_i) #g-1()   #change!!
    
    result[i,6] <- mu_z_hat_i   #input mu_z_hat 
    
    n <- ni[i]  #get n for each group
    
    Var_mu_z_hat_i <- 1/(ni[i]^2)*g_inv_deriv_eta_hat_i%*%CMSEP_eta_hat_i%*%g_inv_deriv_eta_hat_i
    Var_mu_z_hat_i = Var_mu_z_hat_i[1,1] 
    SD_mu_z_hat_i <- Var_mu_z_hat_i^0.5
    
    result[i,8] <- Var_mu_z_hat_i   #input variance of mu_z_hat 
    result[i,9] <- SD_mu_z_hat_i   #input SD of mu_z_hat 

    beta_hat_raw <- getME(glmm.laplace,"beta")
    beta_hat <- as.numeric(beta_hat_raw)  # get esimated betas for FE

    result[i,8] <- Var_mu_z_hat_i  # input variance of mu_z_hat 
    result[i,9] <- SD_mu_z_hat_i   # input SD of mu_z_hat 
    
    alpha_sig <- alpha_sim   # significant level of CLT method
    Z_alpha <- qnorm(1-alpha_sig/2);Z_alpha 
    
    #######Calculation of CI###############  
    
    ## generate confidence interval of mu_z_hat for log link function: by formulas in paper
    CI_mu_z_hat_lower_i_loglink <- mu_z_hat_i*exp(-Z_alpha*SD_mu_z_hat_i/mu_z_hat_i) #lower bound of CI
    CI_mu_z_hat_upper_i_loglink <- mu_z_hat_i*exp(Z_alpha*SD_mu_z_hat_i/mu_z_hat_i)  #upper bound of CI
    
    # CLT method for CI
    CI_mu_z_hat_lower_i_clt <- mu_z_hat_i - Z_alpha * SD_mu_z_hat_i
    CI_mu_z_hat_upper_i_clt <- mu_z_hat_i + Z_alpha * SD_mu_z_hat_i
    
    # no lognormal for this, set 0 instead
    CI_mu_z_hat_lower_i_lognorm = 0
    CI_mu_z_hat_upper_i_lognorm = 0
    
    result[i,10] <- CI_mu_z_hat_lower_i_clt   # input lower bound of CI 
    result[i,11] <- CI_mu_z_hat_upper_i_clt   # input upper bound of CI
    
    # generate mu_z_star_hat, i.e estimator in formula (2) of Qu and Luo's paper
    x1_i <- mean(X[k:(k+ni[i]-1),2])    # mean of x1(baseline) for group_i
    x2_i <- result[i,3]  # group  for i
    x3_i <- result[i,4]  # time for i
    
    if(x3_type_sim==1){
      b_hat <- as.numeric(getME(glmm.laplace,"b"))[k:(k+ni[i]-1)]
    }
    if(x3_type_sim==2){
      if((k==1)|(k==2)){
        b_hat <- as.numeric(getME(glmm.laplace,"b"))[1:ni[i]]
      }
      if((k==3)|(k==4)){
        b_hat <- as.numeric(getME(glmm.laplace,"b"))[(ni[1]+1):(ni[1]+ni[i])]
      }
    }
    mu_z_star_hat <- mean(exp(as.numeric(c(1,x1_i,x2_i,x3_i)%*%beta_hat)+b_hat))
    result[i,7] <- mu_z_star_hat   # input mu_z_star_hat
    
    ##############pi###################    
    
    mu_l  <- mean(mydata[k:(k+ni[i]-1),7])     # mu_bar : mean of each mu_i
    result[i,13] <- mu_l  # input mu_l
    # Y_bar <- mean(y_i)
    
    ###############pi##############   
    
    # Calculate the population value of mu_bar (i.e average of each mu)
    mu_i <- mydata[(mydata[,4]==result[i,3] & mydata[,5]==result[i,4]),7];  # observed y for each group
    mu_bar_i  <- mean(mu_i)     # mu_bar : mean of each mu_i
    result[i,12] <- mu_bar_i   # input mu_bar
    
    if((CI_mu_z_hat_lower_i_loglink <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_loglink)){ 
      cp_loglink[i] <- cp_loglink[i]+1      # count number of Y_bar inCI for group+i
    } 
    
    if((CI_mu_z_hat_lower_i_clt <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_clt)){ 
      cp_clt[i] <- cp_clt[i]+1      # count number of Y_bar inCI for group+i
    }
    
    if((CI_mu_z_hat_lower_i_lognorm <=  mu_l) && (mu_l <= CI_mu_z_hat_upper_i_lognorm)){ 
      cp_lognorm[i] <- cp_lognorm[i]+1      # count number of Y_bar inCI for group+i
    }
    
    k <- k + ni[i] 
  }
  
  result[, 14]=mean(result[, 6])
  result[, 15]=mean(result[, 13])
  result[, 16]=mean(result[, 5])
  result[, 17]=mean(result[, 12])
  
  # input names for each column of the table
  colnames(result) <- c("n_obs","group","treat","time", "Y_bar","mu_z_hat","mu_z_star_hat", "Var_mu_z_hat","SD_mu_z_hat","CI_lower","CI_upper", "mu_bar", "mu_l", "mu_z_hat_total", "mu_l_total", "y_bar_total","mu_bar_total")
  
  result = data.frame(result)
  

  cp_table = cbind(result$n_obs, result$group, result$treat, result$time)
  
  # generate output for estimation of paramter
  sigma_hat = getME(glmm.laplace,"theta")[[1]]  # get estimated sigma for RE
  sigma_hat = round(sigma_hat,digits=5)  # only output 5 digits
  
  size_hat = getME(glmm.laplace,"glmer.nb.theta") # estimated size parameter of NB
  size_hat = round(size_hat,digits=5)
  
  beta_hat = as.numeric(getME(glmm.laplace, "beta"))  # get esimated betas for FE
  beta_hat = round(beta_hat, digits=5)
  
  beta_hat_output = paste(beta_hat, collapse=" ")# change vector into paste string form for output  
  
  
  group = result$group
  treat = result$treat
  time = result$time
  mu = result$mu_l   # true group mean
  mu = round(mu, digits=5)
  mu_z_hat = result$mu_z_hat   # group mean/ marginal mean by adjustment 
  mu_z_hat = round(mu_z_hat, digits=5)
  y_bar = result$Y_bar   # group mean /marginal mean
  y_bar = round(y_bar, digits=5)
  mu_z_star_hat = result$mu_z_star_hat  # repsonse at mean value
  mu_z_star_hat = round(mu_z_star_hat, digits=5)
  
  
  # save all results 
  z = list()
  #Z$model_fit = glmm.laplace
  z$result = result
  z$x1_type_sim = x1_type_sim
  z$x3_type_sim = x3_type_sim
  z$ni_sim = ni_sim
  z$sigma_sim = ni_sim
  z$beta0_sim = beta0_sim
  z$beta1_sim = beta1_sim
  z$beta2_sim = beta2_sim
  z$beta3_sim = beta3_sim
  z$p_bern_sim = p_bern_sim
  z$sigma_hat = sigma_hat
  z$beta_hat_output = beta_hat_output
  z$group = group
  z$treat = treat
  z$time = time
  z$mu = mu
  z$mu_z_hat = mu_z_hat
  z$y_bar = y_bar
  z$mu_z_star_hat = mu_z_star_hat
  z$cp_loglink = cp_loglink
  z$cp_clt = cp_clt
  z$cp_lognorm = cp_lognorm
  z$CI_mu_z_hat_lower_i_loglink = CI_mu_z_hat_lower_i_loglink 
  z$CI_mu_z_hat_upper_i_loglink = CI_mu_z_hat_upper_i_loglink 
  z$CI_mu_z_hat_lower_i_clt = CI_mu_z_hat_lower_i_clt
  z$CI_mu_z_hat_upper_i_clt = CI_mu_z_hat_upper_i_clt
  z$CI_mu_z_hat_lower_i_lognorm = CI_mu_z_hat_lower_i_lognorm
  z$CI_mu_z_hat_upper_i_lognorm = CI_mu_z_hat_upper_i_lognorm
  
  return(z)
  
}  


###########################################################################



# result for simulation
set.seed(100)

# logistic data with marginal group mean and confidence interval
my_sim = ul_ci()
my_sim

# logistic data with conditional group mean and prediction interval
my_sim = ul_pi()
my_sim

# negativ binomial data with marginal group mean and confidence interval
my_sim = unb_ci()
my_sim

# negativ binomial data with conditional group mean and prediction interval
my_sim = unb_pi()
my_sim




# setup of parameters
alpha = 0.05  # signifcant level for CI: e.g. 1-0.05=95% confidence inverval
rep_sim = 1  # replication times for simulation, but in cluster, each job we usually only simulate once 
#set.seed(100)



# simulation for coverage probabiliy in several repetitions
ul_ci_cp1_m = matrix(NA,ncol=4,nrow=rep_sim)
ul_ci_cp2_m = matrix(NA,ncol=4,nrow=rep_sim)
ul_ci_cp3_m = matrix(NA,ncol=4,nrow=rep_sim)
for(i in 1:rep_sim){
  my_sim = ul_ci()
  ul_ci_cp1_m[i,] = my_sim$cp_logitlink
  ul_ci_cp2_m[i,] = my_sim$cp_clt
  ul_ci_cp3_m[i,] = my_sim$cp_lognorm
}
my_sim$y_bar
# display cp in the simulations
apply(ul_ci_cp1_m, 2, mean)
apply(ul_ci_cp2_m, 2, mean)
apply(ul_ci_cp3_m, 2, mean)



# simulation for coverage probabiliy in several repetitions
ul_pi_cp1_m = matrix(0,ncol=4,nrow=rep_sim)
ul_pi_cp2_m = matrix(0,ncol=4,nrow=rep_sim)
ul_pi_cp3_m = matrix(0,ncol=4,nrow=rep_sim)
for(i in 1:rep_sim){
  my_sim = ul_pi()
  ul_pi_cp1_m[i,] = my_sim$cp_logitlink
  ul_pi_cp2_m[i,] = my_sim$cp_clt
  ul_pi_cp3_m[i,] = my_sim$cp_lognorm
}

# display cp in the simulations
apply(ul_pi_cp1_m, 2, mean)
apply(ul_pi_cp2_m, 2, mean)
apply(ul_pi_cp3_m, 2, mean)



# simulation for coverage probabiliy in several repetitions
unb_ci_cp1_m = matrix(0,ncol=4,nrow=rep_sim)
unb_ci_cp2_m = matrix(0,ncol=4,nrow=rep_sim)
unb_ci_cp3_m = matrix(0,ncol=4,nrow=rep_sim)
for(i in 1:rep_sim){
  my_sim = unb_ci()
  unb_ci_cp1_m[i,] = my_sim$cp_loglink
  unb_ci_cp2_m[i,] = my_sim$cp_clt
  unb_ci_cp3_m[i,] = my_sim$cp_lognorm
}

# display cp in the simulations
apply(unb_ci_cp1_m, 2, mean)
apply(unb_ci_cp2_m, 2, mean)
apply(unb_ci_cp3_m, 2, mean)



# simulation for coverage probabiliy in several repetitions
unb_pi_cp1_m = matrix(0,ncol=4,nrow=rep_sim)
unb_pi_cp2_m = matrix(0,ncol=4,nrow=rep_sim)
unb_pi_cp3_m = matrix(0,ncol=4,nrow=rep_sim)
for(i in 1:rep_sim){
  my_sim = unb_pi()
  unb_pi_cp1_m[i,] = my_sim$cp_loglink
  unb_pi_cp2_m[i,] = my_sim$cp_clt
  unb_pi_cp3_m[i,] = my_sim$cp_lognorm
}

# display cp in the simulations
apply(unb_pi_cp1_m, 2, mean)
apply(unb_pi_cp2_m, 2, mean)
apply(unb_pi_cp3_m, 2, mean)




