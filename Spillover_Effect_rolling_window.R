return <- read_excel("return_G7.xlsx")
retmat = return
T = length(retmat$Date)
sdmat = matrix(0,1199,6)
ermat = matrix(0,1199,6)
for(k in 1:6){
  r_vec = as.numeric(unlist(retmat[,(k+1)]))
  # transformation function to restrict the parameter values
  trans <- function(theta){
    theta_trans = theta
    theta_trans[1] = -abs(theta[1])
    theta_trans[2] = abs(theta[2])
    theta_trans[3] = (0.8+exp(theta[3]))/(1+exp(theta[3])) # persistent probability of regime 0 between 0 and 1
    theta_trans[4] = (0.8+exp(theta[4]))/(1+exp(theta[4])) # persistent probability of regime 1 between 0 and 1
    theta_trans[5] = abs(theta[5])
    theta_trans[6] = abs(theta[6])
    theta_trans[7] = exp(theta[7])/(1+exp(theta[7])+exp(theta[8]))
    theta_trans[8] = exp(theta[8])/(1+exp(theta[7])+exp(theta[8]))
    return(theta_trans)
  }
  
  # Likelihood function
  fn <- function(theta){
    
    theta = trans(theta) # this transformation guarantee the input parameter values are always reasonable
    
    b0_e = theta[1]
    b1_e = theta[2]
    p00_e = theta[3]
    p11_e = theta[4]
    c0_e = theta[5]
    c1_e = theta[6]
    alpha_e = theta[7]
    beta_e = theta[8]
    
    p01_e = 1-p00_e
    p10_e = 1-p11_e
    
    p_0_ss = (1-p11_e)/(2-p00_e-p11_e)
    p_1_ss = 1-p_0_ss
    
    pr0_L = p_0_ss
    pr1_L = p_1_ss
    
    h0_e = var(r_vec)
    r0_e = 0
    h_e_1_vec = matrix(0,T,1)
    h_e_1_vec[1] = h0_e
    
    LnL = 0
    
    for(i in 1:T){
      if (i==1){
        rt_1 = 0
      } else{
        rt_1 = r_vec[i-1]
      }
      rt = r_vec[i]
      
      #Distribution of y given st and st-1
      
      #st-1=0 st=0
      h0_e = c0_e + alpha_e*(rt_1-b0_e)^2 + beta_e*h_e_1_vec[i]
      fy_00 = (1/sqrt(h0_e))*exp(-(rt-b0_e)^2/(2*h0_e))*p00_e*pr0_L
      
      #st-1=0 st=1
      h1_e = c1_e + alpha_e*(rt_1-b0_e)^2 + beta_e*h_e_1_vec[i]
      fy_01 = (1/sqrt(h1_e))*exp(-(rt-b1_e)^2/(2*h1_e))*p01_e*pr0_L
      
      #st-1=1 st=0
      h0_e = c0_e + alpha_e*(rt_1-b1_e)^2 + beta_e*h_e_1_vec[i]
      fy_10 = (1/sqrt(h0_e))*exp(-(rt-b0_e)^2/(2*h0_e))*p10_e*pr1_L
      
      #st-1=1 st=1
      h1_e = c1_e + alpha_e*(rt_1-b1_e)^2 + beta_e*h_e_1_vec[i]
      fy_11 = (1/sqrt(h1_e))*exp(-(rt-b1_e)^2/(2*h1_e))*p11_e*pr1_L
      
      fy = fy_00 + fy_01 + fy_10 + fy_11
      
      LnL = LnL + log(fy)
      
      #updating process
      
      pr0_u = (fy_00 + fy_10)/fy
      pr1_u = (fy_01 + fy_11)/fy
      
      h_e_1_vec[i+1] = h0_e*pr0_u + h1_e*pr1_u
      
      pr0_L = pr0_u
      pr1_L = pr1_u
    }
    return(-LnL)
  }
  
  if(k<2){
    garchms <- optim(c(0.04,0.06,0.85,0.95,0.1,0.05,0.1,0.1),fn,hessian=TRUE)
  } else if(k>1 & k<3){
    garchms <- optim(c(0.05,0.05,0.85,0.95,0.1,0.05,0.1,0.1),fn,hessian=TRUE)
  } else if(k>2 & k<5) {
    garchms <- optim(c(0.03,0.05,0.85,0.95,0.1,0.05,0.1,0.1),fn,hessian=TRUE)
  } else {
    garchms <- optim(c(0.02,0.05,0.85,0.95,0.1,0.05,0.1,0.1),fn,hessian=TRUE)
  }
  
  #print(c(trans(garchms$par)))
  
  #AICc = 2*garchms$value+2*length(garchms$par)+2*(garchms$par)*(garchms$par+1)/(T-garchms$par-1) # information criteria
  
  b0_e = trans(garchms$par)[1]
  b1_e = trans(garchms$par)[2]
  p00_e = trans(garchms$par)[3]
  p11_e = trans(garchms$par)[4]
  c0_e = trans(garchms$par)[5]
  c1_e = trans(garchms$par)[6]
  alpha_e = trans(garchms$par)[7]
  beta_e = trans(garchms$par)[8]
  
  p01_e = 1-p00_e
  p10_e = 1-p11_e
  
  p_0_ss = (1-p11_e)/(2-p00_e-p11_e)
  p_1_ss = 1-p_0_ss
  
  pr0_L = p_0_ss
  pr1_L = p_1_ss
  
  h0_e = var(r_vec)
  r0_e = 0
  h_e_1_vec = matrix(0,T,1)
  h_e_1_vec[1] = h0_e
  
  pr0 = matrix(0,T,1)
  pr1 = matrix(0,T,1)
  
  for(i in 1:T){
    if (i==1){
      rt_1 = 0
    } else{
      rt_1 = r_vec[i-1]
    }
    rt = r_vec[i]
    
    #Distribution of y given st and st-1
    
    #st-1=0 st=0
    h0_e = c0_e + alpha_e*(rt_1-b0_e)^2 + beta_e*h_e_1_vec[i]
    fy_00 = (1/sqrt(h0_e))*exp(-(rt-b0_e)^2/(2*h0_e))*p00_e*pr0_L
    
    #st-1=0 st=1
    h1_e = c1_e + alpha_e*(rt_1-b0_e)^2 + beta_e*h_e_1_vec[i]
    fy_01 = (1/sqrt(h1_e))*exp(-(rt-b1_e)^2/(2*h1_e))*p01_e*pr0_L
    
    #st-1=1 st=0
    h0_e = c0_e + alpha_e*(rt_1-b1_e)^2 + beta_e*h_e_1_vec[i]
    fy_10 = (1/sqrt(h0_e))*exp(-(rt-b0_e)^2/(2*h0_e))*p10_e*pr1_L
    
    #st-1=1 st=1
    h1_e = c1_e + alpha_e*(rt_1-b1_e)^2 + beta_e*h_e_1_vec[i]
    fy_11 = (1/sqrt(h1_e))*exp(-(rt-b1_e)^2/(2*h1_e))*p11_e*pr1_L
    
    fy = fy_00 + fy_01 + fy_10 + fy_11
    
    #updating process
    
    pr0_u = (fy_00 + fy_10)/fy
    pr1_u = (fy_01 + fy_11)/fy
    
    h_e_1_vec[i+1] = h0_e*pr0_u + h1_e*pr1_u
    
    pr0_L = pr0_u
    pr1_L = pr1_u
    
    pr0[i] = pr0_u
    pr1[i] = pr1_u
  }
  
  sdmat[,k] = as.vector((h_e_1_vec[2:(T+1)]))
  ermat[,k] = as.vector(pr0*b0_e+pr1*b1_e)
}

theta1mat = matrix(0,1148,10)
theta2mat = matrix(0,1148,10)
theta3mat = matrix(0,1148,10)
theta4mat = matrix(0,1148,10)
theta5mat = matrix(0,1148,10)
theta6mat = matrix(0,1148,10)

#t = 1148

for(t in 1:1148){
sp500vol = sdmat[2:(52+t-1),1]
sp500vol_1 = sdmat[1:((52+t-1)-1),1]
TSX60vol = sdmat[2:(52+t-1),2]
TSX60vol_1 = sdmat[1:((52+t-1)-1),2]
Nikkeivol = sdmat[2:(52+t-1),3]
Nikkeivol_1 = sdmat[1:((52+t-1)-1),3]
MIB30vol = sdmat[2:(52+t-1),4]
MIB30vol_1 = sdmat[1:((52+t-1)-1),4]
DAX40vol = sdmat[2:(52+t-1),5]
DAX40vol_1 = sdmat[1:((52+t-1)-1),5]
CAC40vol = sdmat[2:(52+t-1),6]
CAC40vol_1 = sdmat[1:((52+t-1)-1),6]

# Vector Autoregressive Part

var1 <- lm(sp500vol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

var2 <- lm(TSX60vol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

var3 <- lm(Nikkeivol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

var4 <- lm(MIB30vol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

var5 <- lm(DAX40vol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

var6 <- lm(CAC40vol~sp500vol_1+TSX60vol_1+Nikkeivol_1
           +MIB30vol_1+DAX40vol_1+CAC40vol_1)

PHI1 = matrix(c(summary(var1)$coefficients[2:7,1],
                summary(var2)$coefficients[2:7,1],
                summary(var3)$coefficients[2:7,1],
                summary(var4)$coefficients[2:7,1],
                summary(var5)$coefficients[2:7,1],
                summary(var6)$coefficients[2:7,1]
),6,6
)

PHI1 = t(PHI1)

residual = data.frame(var1$residuals,var2$residuals,var3$residuals,var4$residuals,var5$residuals,var6$residuals)

Sigma = cov(residual)

A0 = diag(6)

A1 = PHI1 %*% A0
A2 = PHI1 %*% A1
A3 = PHI1 %*% A2

e1 = matrix(c(1,0,0,0,0,0),6,1)
e2 = matrix(c(0,1,0,0,0,0),6,1)
e3 = matrix(c(0,0,1,0,0,0),6,1)
e4 = matrix(c(0,0,0,1,0,0),6,1)
e5 = matrix(c(0,0,0,0,1,0),6,1)
e6 = matrix(c(0,0,0,0,0,1),6,1)

# SP500

theta_11 = sqrt(Sigma[1,1])^(-1)*((t(e1)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e1)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_12 = sqrt(Sigma[2,2])^(-1)*((t(e1)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e2)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_13 = sqrt(Sigma[3,3])^(-1)*((t(e1)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e3)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_14 = sqrt(Sigma[4,4])^(-1)*((t(e1)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e4)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_15 = sqrt(Sigma[5,5])^(-1)*((t(e1)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e5)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_16 = sqrt(Sigma[6,6])^(-1)*((t(e1)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e1)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e1)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e1)%*%A3%*%Sigma%*%e6)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_1st = c(theta_11/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16),
              theta_12/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16),
              theta_13/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16),
              theta_14/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16),
              theta_15/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16),
              theta_16/(theta_11+theta_12+theta_13+theta_14+theta_15+theta_16))

# TSX60

theta_21 = sqrt(Sigma[1,1])^(-1)*((t(e2)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e1)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_22 = sqrt(Sigma[2,2])^(-1)*((t(e2)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e2)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_23 = sqrt(Sigma[3,3])^(-1)*((t(e2)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e3)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_24 = sqrt(Sigma[4,4])^(-1)*((t(e2)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e4)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_25 = sqrt(Sigma[5,5])^(-1)*((t(e2)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e5)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_26 = sqrt(Sigma[6,6])^(-1)*((t(e2)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e2)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e2)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e2)%*%A3%*%Sigma%*%e6)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_2st = c(theta_21/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26),
              theta_22/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26),
              theta_23/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26),
              theta_24/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26),
              theta_25/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26),
              theta_26/(theta_21+theta_22+theta_23+theta_24+theta_25+theta_26))

# Nikkei

theta_31 = sqrt(Sigma[1,1])^(-1)*((t(e3)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e1)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_32 = sqrt(Sigma[2,2])^(-1)*((t(e3)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e2)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_33 = sqrt(Sigma[3,3])^(-1)*((t(e3)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e3)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_34 = sqrt(Sigma[4,4])^(-1)*((t(e3)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e4)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_35 = sqrt(Sigma[5,5])^(-1)*((t(e3)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e5)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_36 = sqrt(Sigma[6,6])^(-1)*((t(e3)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e3)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e3)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e3)%*%A3%*%Sigma%*%e6)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_3st = c(theta_31/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36),
              theta_32/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36),
              theta_33/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36),
              theta_34/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36),
              theta_35/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36),
              theta_36/(theta_31+theta_32+theta_33+theta_34+theta_35+theta_36))

# MIB30

theta_41 = sqrt(Sigma[1,1])^(-1)*((t(e4)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e1)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_42 = sqrt(Sigma[2,2])^(-1)*((t(e4)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e2)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_43 = sqrt(Sigma[3,3])^(-1)*((t(e4)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e3)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_44 = sqrt(Sigma[4,4])^(-1)*((t(e4)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e4)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_45 = sqrt(Sigma[5,5])^(-1)*((t(e4)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e5)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_46 = sqrt(Sigma[6,6])^(-1)*((t(e4)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e4)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e4)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e4)%*%A3%*%Sigma%*%e6)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_4st = c(theta_41/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46),
              theta_42/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46),
              theta_43/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46),
              theta_44/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46),
              theta_45/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46),
              theta_46/(theta_41+theta_42+theta_43+theta_44+theta_45+theta_46))

# DAX40

theta_51 = sqrt(Sigma[1,1])^(-1)*((t(e5)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e1)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_52 = sqrt(Sigma[2,2])^(-1)*((t(e5)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e2)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_53 = sqrt(Sigma[3,3])^(-1)*((t(e5)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e3)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_54 = sqrt(Sigma[4,4])^(-1)*((t(e5)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e4)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_55 = sqrt(Sigma[5,5])^(-1)*((t(e5)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e5)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_56 = sqrt(Sigma[6,6])^(-1)*((t(e5)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e5)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e5)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e5)%*%A3%*%Sigma%*%e6)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_5st = c(theta_51/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56),
              theta_52/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56),
              theta_53/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56),
              theta_54/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56),
              theta_55/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56),
              theta_56/(theta_51+theta_52+theta_53+theta_54+theta_55+theta_56))

# CAC40

theta_61 = sqrt(Sigma[1,1])^(-1)*((t(e6)%*%A0%*%Sigma%*%e1)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e1)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e1)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e1)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_62 = sqrt(Sigma[2,2])^(-1)*((t(e6)%*%A0%*%Sigma%*%e2)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e2)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e2)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e2)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_63 = sqrt(Sigma[3,3])^(-1)*((t(e6)%*%A0%*%Sigma%*%e3)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e3)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e3)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e3)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_64 = sqrt(Sigma[4,4])^(-1)*((t(e6)%*%A0%*%Sigma%*%e4)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e4)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e4)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e4)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_65 = sqrt(Sigma[5,5])^(-1)*((t(e6)%*%A0%*%Sigma%*%e5)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e5)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e5)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e5)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_66 = sqrt(Sigma[6,6])^(-1)*((t(e6)%*%A0%*%Sigma%*%e6)^2
                                  +(t(e6)%*%A1%*%Sigma%*%e6)^2
                                  +(t(e6)%*%A2%*%Sigma%*%e6)^2
                                  +(t(e6)%*%A3%*%Sigma%*%e6)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_6st = c(theta_61/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66),
              theta_62/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66),
              theta_63/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66),
              theta_64/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66),
              theta_65/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66),
              theta_66/(theta_61+theta_62+theta_63+theta_64+theta_65+theta_66))

# Spillover to all other 

theta1all = sum(theta_1st[2:6])
theta2all = sum(theta_2st[1:6])-theta_2st[2]
theta3all = sum(theta_3st[1:6])-theta_3st[3]
theta4all = sum(theta_4st[1:6])-theta_4st[4]
theta5all = sum(theta_5st[1:6])-theta_5st[5]
theta6all = sum(theta_6st[1:6])-theta_6st[6]

# Spillover from all other

thetaall1 = theta_2st[1]+theta_3st[1]+theta_4st[1]+theta_5st[1]+theta_6st[1]
thetaall2 = theta_1st[2]+theta_3st[2]+theta_4st[2]+theta_5st[2]+theta_6st[2]
thetaall3 = theta_1st[3]+theta_2st[3]+theta_4st[3]+theta_5st[3]+theta_6st[3]
thetaall4 = theta_1st[4]+theta_2st[4]+theta_3st[4]+theta_5st[4]+theta_6st[4]
thetaall5 = theta_1st[5]+theta_2st[5]+theta_3st[5]+theta_4st[5]+theta_6st[5]
thetaall6 = theta_1st[6]+theta_2st[6]+theta_3st[6]+theta_4st[6]+theta_5st[6]

# Good Volatility

sp500gdvol = sdmat[2:(52+t-1),1]*(retmat[2:(52+t-1),2]>ermat[2:(52+t-1),1])
sp500gdvol_1 = sdmat[1:((52+t-1)-1),1]*(retmat[1:((52+t-1)-1),2]>ermat[1:((52+t-1)-1),1])
TSX60gdvol = sdmat[2:(52+t-1),2]*(retmat[2:(52+t-1),3]>ermat[2:(52+t-1),2])
TSX60gdvol_1 = sdmat[1:((52+t-1)-1),2]*(retmat[1:((52+t-1)-1),3]>ermat[1:((52+t-1)-1),2])
Nikkeigdvol = sdmat[2:(52+t-1),3]*(retmat[2:(52+t-1),4]>ermat[2:(52+t-1),3])
Nikkeigdvol_1 = sdmat[1:((52+t-1)-1),3]*(retmat[1:((52+t-1)-1),4]>ermat[1:((52+t-1)-1),3])
MIB30gdvol = sdmat[2:(52+t-1),4]*(retmat[2:(52+t-1),5]>ermat[2:(52+t-1),4])
MIB30gdvol_1 = sdmat[1:((52+t-1)-1),4]*(retmat[1:((52+t-1)-1),5]>ermat[1:((52+t-1)-1),4])
DAX40gdvol = sdmat[2:(52+t-1),5]*(retmat[2:(52+t-1),6]>ermat[2:(52+t-1),5])
DAX40gdvol_1 = sdmat[1:((52+t-1)-1),5]*(retmat[1:((52+t-1)-1),6]>ermat[1:((52+t-1)-1),5])
CAC40gdvol = sdmat[2:(52+t-1),6]*(retmat[2:(52+t-1),7]>ermat[2:(52+t-1),6])
CAC40gdvol_1 = sdmat[1:((52+t-1)-1),6]*(retmat[1:((52+t-1)-1),7]>ermat[1:((52+t-1)-1),6])

# Vector Autoregressive Part

var1 <- lm(sp500gdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

var2 <- lm(TSX60gdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

var3 <- lm(Nikkeigdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

var4 <- lm(MIB30gdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

var5 <- lm(DAX40gdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

var6 <- lm(CAC40gdvol~sp500gdvol_1+TSX60gdvol_1+Nikkeigdvol_1
           +MIB30gdvol_1+DAX40gdvol_1+CAC40gdvol_1)

PHI1 = matrix(c(summary(var1)$coefficients[2:7,1],
                summary(var2)$coefficients[2:7,1],
                summary(var3)$coefficients[2:7,1],
                summary(var4)$coefficients[2:7,1],
                summary(var5)$coefficients[2:7,1],
                summary(var6)$coefficients[2:7,1]
),6,6
)

PHI1 = t(PHI1)

residual = data.frame(var1$residuals,var2$residuals,var3$residuals,var4$residuals,var5$residuals,var6$residuals)

Sigma = cov(residual)

A0 = diag(6)

A1 = PHI1 %*% A0
A2 = PHI1 %*% A1
A3 = PHI1 %*% A2

e1 = matrix(c(1,0,0,0,0,0),6,1)
e2 = matrix(c(0,1,0,0,0,0),6,1)
e3 = matrix(c(0,0,1,0,0,0),6,1)
e4 = matrix(c(0,0,0,1,0,0),6,1)
e5 = matrix(c(0,0,0,0,1,0),6,1)
e6 = matrix(c(0,0,0,0,0,1),6,1)

# SP500

theta_gd_11 = sqrt(Sigma[1,1])^(-1)*((t(e1)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e1)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_12 = sqrt(Sigma[2,2])^(-1)*((t(e1)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e2)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_13 = sqrt(Sigma[3,3])^(-1)*((t(e1)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e3)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_14 = sqrt(Sigma[4,4])^(-1)*((t(e1)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e4)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_15 = sqrt(Sigma[5,5])^(-1)*((t(e1)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e5)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_16 = sqrt(Sigma[6,6])^(-1)*((t(e1)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e6)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd_1st = c(theta_gd_11/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16),
                 theta_gd_12/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16),
                 theta_gd_13/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16),
                 theta_gd_14/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16),
                 theta_gd_15/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16),
                 theta_gd_16/(theta_gd_11+theta_gd_12+theta_gd_13+theta_gd_14+theta_gd_15+theta_gd_16))

# TSX60

theta_gd_21 = sqrt(Sigma[1,1])^(-1)*((t(e2)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e1)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_22 = sqrt(Sigma[2,2])^(-1)*((t(e2)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e2)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_23 = sqrt(Sigma[3,3])^(-1)*((t(e2)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e3)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_24 = sqrt(Sigma[4,4])^(-1)*((t(e2)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e4)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_25 = sqrt(Sigma[5,5])^(-1)*((t(e2)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e5)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_26 = sqrt(Sigma[6,6])^(-1)*((t(e2)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e6)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd_2st = c(theta_gd_21/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26),
                 theta_gd_22/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26),
                 theta_gd_23/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26),
                 theta_gd_24/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26),
                 theta_gd_25/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26),
                 theta_gd_26/(theta_gd_21+theta_gd_22+theta_gd_23+theta_gd_24+theta_gd_25+theta_gd_26))

# Nikkei

theta_gd_31 = sqrt(Sigma[1,1])^(-1)*((t(e3)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e1)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_32 = sqrt(Sigma[2,2])^(-1)*((t(e3)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e2)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_33 = sqrt(Sigma[3,3])^(-1)*((t(e3)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e3)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_34 = sqrt(Sigma[4,4])^(-1)*((t(e3)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e4)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_35 = sqrt(Sigma[5,5])^(-1)*((t(e3)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e5)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_36 = sqrt(Sigma[6,6])^(-1)*((t(e3)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e6)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd_3st = c(theta_gd_31/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36),
                 theta_gd_32/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36),
                 theta_gd_33/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36),
                 theta_gd_34/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36),
                 theta_gd_35/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36),
                 theta_gd_36/(theta_gd_31+theta_gd_32+theta_gd_33+theta_gd_34+theta_gd_35+theta_gd_36))

# MIB30

theta_gd_41 = sqrt(Sigma[1,1])^(-1)*((t(e4)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e1)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_42 = sqrt(Sigma[2,2])^(-1)*((t(e4)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e2)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_43 = sqrt(Sigma[3,3])^(-1)*((t(e4)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e3)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_44 = sqrt(Sigma[4,4])^(-1)*((t(e4)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e4)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_45 = sqrt(Sigma[5,5])^(-1)*((t(e4)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e5)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_46 = sqrt(Sigma[6,6])^(-1)*((t(e4)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e6)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd_4st = c(theta_gd_41/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46),
                 theta_gd_42/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46),
                 theta_gd_43/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46),
                 theta_gd_44/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46),
                 theta_gd_45/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46),
                 theta_gd_46/(theta_gd_41+theta_gd_42+theta_gd_43+theta_gd_44+theta_gd_45+theta_gd_46))

# DAX40

theta_gd_51 = sqrt(Sigma[1,1])^(-1)*((t(e5)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e1)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_52 = sqrt(Sigma[2,2])^(-1)*((t(e5)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e2)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_53 = sqrt(Sigma[3,3])^(-1)*((t(e5)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e3)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_54 = sqrt(Sigma[4,4])^(-1)*((t(e5)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e4)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_55 = sqrt(Sigma[5,5])^(-1)*((t(e5)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e5)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_56 = sqrt(Sigma[6,6])^(-1)*((t(e5)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e6)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd_5st = c(theta_gd_51/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56),
                 theta_gd_52/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56),
                 theta_gd_53/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56),
                 theta_gd_54/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56),
                 theta_gd_55/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56),
                 theta_gd_56/(theta_gd_51+theta_gd_52+theta_gd_53+theta_gd_54+theta_gd_55+theta_gd_56))

# CAC40

theta_gd_61 = sqrt(Sigma[1,1])^(-1)*((t(e6)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e1)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_62 = sqrt(Sigma[2,2])^(-1)*((t(e6)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e2)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_63 = sqrt(Sigma[3,3])^(-1)*((t(e6)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e3)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_64 = sqrt(Sigma[4,4])^(-1)*((t(e6)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e4)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_65 = sqrt(Sigma[5,5])^(-1)*((t(e6)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e5)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_66 = sqrt(Sigma[6,6])^(-1)*((t(e6)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e6)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd_6st = c(theta_gd_61/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66),
                 theta_gd_62/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66),
                 theta_gd_63/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66),
                 theta_gd_64/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66),
                 theta_gd_65/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66),
                 theta_gd_66/(theta_gd_61+theta_gd_62+theta_gd_63+theta_gd_64+theta_gd_65+theta_gd_66))

# Spillover to all other 

theta1all_gd = sum(theta_gd_1st[2:6])
theta2all_gd = sum(theta_gd_2st[1:6])-theta_gd_2st[2]
theta3all_gd = sum(theta_gd_3st[1:6])-theta_gd_3st[3]
theta4all_gd = sum(theta_gd_4st[1:6])-theta_gd_4st[4]
theta5all_gd = sum(theta_gd_5st[1:6])-theta_gd_5st[5]
theta6all_gd = sum(theta_gd_6st[1:6])-theta_gd_6st[6]

# Spillover from all other

thetaall1_gd = theta_gd_2st[1]+theta_gd_3st[1]+theta_gd_4st[1]+theta_gd_5st[1]+theta_gd_6st[1]
thetaall2_gd = theta_gd_1st[2]+theta_gd_3st[2]+theta_gd_4st[2]+theta_gd_5st[2]+theta_gd_6st[2]
thetaall3_gd = theta_gd_1st[3]+theta_gd_2st[3]+theta_gd_4st[3]+theta_gd_5st[3]+theta_gd_6st[3]
thetaall4_gd = theta_gd_1st[4]+theta_gd_2st[4]+theta_gd_3st[4]+theta_gd_5st[4]+theta_gd_6st[4]
thetaall5_gd = theta_gd_1st[5]+theta_gd_2st[5]+theta_gd_3st[5]+theta_gd_4st[5]+theta_gd_6st[5]
thetaall6_gd = theta_gd_1st[6]+theta_gd_2st[6]+theta_gd_3st[6]+theta_gd_4st[6]+theta_gd_5st[6]

# Good Volatility with 0 threshold

sp500gd0vol = sdmat[2:(52+t-1),1]*(retmat[2:(52+t-1),2]>0)
sp500gd0vol_1 = sdmat[1:((52+t-1)-1),1]*(retmat[1:((52+t-1)-1),2]>0)
TSX60gd0vol = sdmat[2:(52+t-1),2]*(retmat[2:(52+t-1),3]>0)
TSX60gd0vol_1 = sdmat[1:((52+t-1)-1),2]*(retmat[1:((52+t-1)-1),3]>0)
Nikkeigd0vol = sdmat[2:(52+t-1),3]*(retmat[2:(52+t-1),4]>0)
Nikkeigd0vol_1 = sdmat[1:((52+t-1)-1),3]*(retmat[1:((52+t-1)-1),4]>0)
MIB30gd0vol = sdmat[2:(52+t-1),4]*(retmat[2:(52+t-1),5]>0)
MIB30gd0vol_1 = sdmat[1:((52+t-1)-1),4]*(retmat[1:((52+t-1)-1),5]>0)
DAX40gd0vol = sdmat[2:(52+t-1),5]*(retmat[2:(52+t-1),6]>0)
DAX40gd0vol_1 = sdmat[1:((52+t-1)-1),5]*(retmat[1:((52+t-1)-1),6]>0)
CAC40gd0vol = sdmat[2:(52+t-1),6]*(retmat[2:(52+t-1),7]>0)
CAC40gd0vol_1 = sdmat[1:((52+t-1)-1),6]*(retmat[1:((52+t-1)-1),7]>0)

var1 <- lm(sp500gd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

var2 <- lm(TSX60gd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

var3 <- lm(Nikkeigd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

var4 <- lm(MIB30gd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

var5 <- lm(DAX40gd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

var6 <- lm(CAC40gd0vol~sp500gd0vol_1+TSX60gd0vol_1+Nikkeigd0vol_1
           +MIB30gd0vol_1+DAX40gd0vol_1+CAC40gd0vol_1)

PHI1 = matrix(c(summary(var1)$coefficients[2:7,1],
                summary(var2)$coefficients[2:7,1],
                summary(var3)$coefficients[2:7,1],
                summary(var4)$coefficients[2:7,1],
                summary(var5)$coefficients[2:7,1],
                summary(var6)$coefficients[2:7,1]
),6,6
)

PHI1 = t(PHI1)

residual = data.frame(var1$residuals,var2$residuals,var3$residuals,var4$residuals,var5$residuals,var6$residuals)

Sigma = cov(residual)

A0 = diag(6)

A1 = PHI1 %*% A0
A2 = PHI1 %*% A1
A3 = PHI1 %*% A2

e1 = matrix(c(1,0,0,0,0,0),6,1)
e2 = matrix(c(0,1,0,0,0,0),6,1)
e3 = matrix(c(0,0,1,0,0,0),6,1)
e4 = matrix(c(0,0,0,1,0,0),6,1)
e5 = matrix(c(0,0,0,0,1,0),6,1)
e6 = matrix(c(0,0,0,0,0,1),6,1)

# SP500

theta_gd0_11 = sqrt(Sigma[1,1])^(-1)*((t(e1)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e1)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_12 = sqrt(Sigma[2,2])^(-1)*((t(e1)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e2)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_13 = sqrt(Sigma[3,3])^(-1)*((t(e1)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e3)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_14 = sqrt(Sigma[4,4])^(-1)*((t(e1)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e4)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_15 = sqrt(Sigma[5,5])^(-1)*((t(e1)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e5)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_16 = sqrt(Sigma[6,6])^(-1)*((t(e1)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e6)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_gd0_1st = c(theta_gd0_11/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16),
                  theta_gd0_12/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16),
                  theta_gd0_13/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16),
                  theta_gd0_14/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16),
                  theta_gd0_15/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16),
                  theta_gd0_16/(theta_gd0_11+theta_gd0_12+theta_gd0_13+theta_gd0_14+theta_gd0_15+theta_gd0_16))

# TSX60

theta_gd0_21 = sqrt(Sigma[1,1])^(-1)*((t(e2)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e1)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_22 = sqrt(Sigma[2,2])^(-1)*((t(e2)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e2)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_23 = sqrt(Sigma[3,3])^(-1)*((t(e2)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e3)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_24 = sqrt(Sigma[4,4])^(-1)*((t(e2)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e4)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_25 = sqrt(Sigma[5,5])^(-1)*((t(e2)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e5)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_26 = sqrt(Sigma[6,6])^(-1)*((t(e2)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e6)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_gd0_2st = c(theta_gd0_21/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26),
                  theta_gd0_22/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26),
                  theta_gd0_23/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26),
                  theta_gd0_24/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26),
                  theta_gd0_25/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26),
                  theta_gd0_26/(theta_gd0_21+theta_gd0_22+theta_gd0_23+theta_gd0_24+theta_gd0_25+theta_gd0_26))

# Nikkei

theta_gd0_31 = sqrt(Sigma[1,1])^(-1)*((t(e3)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e1)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_32 = sqrt(Sigma[2,2])^(-1)*((t(e3)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e2)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_33 = sqrt(Sigma[3,3])^(-1)*((t(e3)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e3)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_34 = sqrt(Sigma[4,4])^(-1)*((t(e3)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e4)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_35 = sqrt(Sigma[5,5])^(-1)*((t(e3)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e5)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_36 = sqrt(Sigma[6,6])^(-1)*((t(e3)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e6)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_gd0_3st = c(theta_gd0_31/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36),
                  theta_gd0_32/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36),
                  theta_gd0_33/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36),
                  theta_gd0_34/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36),
                  theta_gd0_35/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36),
                  theta_gd0_36/(theta_gd0_31+theta_gd0_32+theta_gd0_33+theta_gd0_34+theta_gd0_35+theta_gd0_36))

# MIB30

theta_gd0_41 = sqrt(Sigma[1,1])^(-1)*((t(e4)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e1)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_42 = sqrt(Sigma[2,2])^(-1)*((t(e4)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e2)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_43 = sqrt(Sigma[3,3])^(-1)*((t(e4)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e3)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_44 = sqrt(Sigma[4,4])^(-1)*((t(e4)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e4)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_45 = sqrt(Sigma[5,5])^(-1)*((t(e4)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e5)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_46 = sqrt(Sigma[6,6])^(-1)*((t(e4)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e6)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_gd0_4st = c(theta_gd0_41/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46), 
                  theta_gd0_42/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46),
                  theta_gd0_43/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46),
                  theta_gd0_44/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46),
                  theta_gd0_45/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46),
                  theta_gd0_46/(theta_gd0_41+theta_gd0_42+theta_gd0_43+theta_gd0_44+theta_gd0_45+theta_gd0_46))

# DAX40

theta_gd0_51 = sqrt(Sigma[1,1])^(-1)*((t(e5)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e1)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_52 = sqrt(Sigma[2,2])^(-1)*((t(e5)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e2)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_53 = sqrt(Sigma[3,3])^(-1)*((t(e5)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e3)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_54 = sqrt(Sigma[4,4])^(-1)*((t(e5)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e4)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_55 = sqrt(Sigma[5,5])^(-1)*((t(e5)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e5)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_56 = sqrt(Sigma[6,6])^(-1)*((t(e5)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e6)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_gd0_5st = c(theta_gd0_51/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56),
                  theta_gd0_52/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56),
                  theta_gd0_53/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56),
                  theta_gd0_54/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56),
                  theta_gd0_55/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56),
                  theta_gd0_56/(theta_gd0_51+theta_gd0_52+theta_gd0_53+theta_gd0_54+theta_gd0_55+theta_gd0_56))

# CAC40

theta_gd0_61 = sqrt(Sigma[1,1])^(-1)*((t(e6)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e1)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_62 = sqrt(Sigma[2,2])^(-1)*((t(e6)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e2)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_63 = sqrt(Sigma[3,3])^(-1)*((t(e6)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e3)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_64 = sqrt(Sigma[4,4])^(-1)*((t(e6)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e4)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_65 = sqrt(Sigma[5,5])^(-1)*((t(e6)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e5)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_66 = sqrt(Sigma[6,6])^(-1)*((t(e6)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e6)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_gd0_6st = c(theta_gd0_61/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66),
                  theta_gd0_62/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66),
                  theta_gd0_63/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66),
                  theta_gd0_64/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66),
                  theta_gd0_65/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66),
                  theta_gd0_66/(theta_gd0_61+theta_gd0_62+theta_gd0_63+theta_gd0_64+theta_gd0_65+theta_gd0_66))

# Spillover to all other 

theta1all_gd0 = sum(theta_gd0_1st[2:6])
theta2all_gd0 = sum(theta_gd0_2st[1:6])-theta_gd0_2st[2]
theta3all_gd0 = sum(theta_gd0_3st[1:6])-theta_gd0_3st[3]
theta4all_gd0 = sum(theta_gd0_4st[1:6])-theta_gd0_4st[4]
theta5all_gd0 = sum(theta_gd0_5st[1:6])-theta_gd0_5st[5]
theta6all_gd0 = sum(theta_gd0_6st[1:6])-theta_gd0_6st[6]

# Spillover from all other

thetaall1_gd0 = theta_gd0_2st[1]+theta_gd0_3st[1]+theta_gd0_4st[1]+theta_gd0_5st[1]+theta_gd0_6st[1]
thetaall2_gd0 = theta_gd0_1st[2]+theta_gd0_3st[2]+theta_gd0_4st[2]+theta_gd0_5st[2]+theta_gd0_6st[2]
thetaall3_gd0 = theta_gd0_1st[3]+theta_gd0_2st[3]+theta_gd0_4st[3]+theta_gd0_5st[3]+theta_gd0_6st[3]
thetaall4_gd0 = theta_gd0_1st[4]+theta_gd0_2st[4]+theta_gd0_3st[4]+theta_gd0_5st[4]+theta_gd0_6st[4]
thetaall5_gd0 = theta_gd0_1st[5]+theta_gd0_2st[5]+theta_gd0_3st[5]+theta_gd0_4st[5]+theta_gd0_6st[5]
thetaall6_gd0 = theta_gd0_1st[6]+theta_gd0_2st[6]+theta_gd0_3st[6]+theta_gd0_4st[6]+theta_gd0_5st[6]

#-----------------------Bad Spillover---------------------

# Bad Volatility

sp500bdvol = sdmat[2:(52+t-1),1]*(retmat[2:(52+t-1),2]<ermat[2:(52+t-1),1])
sp500bdvol_1 = sdmat[1:((52+t-1)-1),1]*(retmat[1:((52+t-1)-1),2]<ermat[1:((52+t-1)-1),1])
TSX60bdvol = sdmat[2:(52+t-1),2]*(retmat[2:(52+t-1),3]<ermat[2:(52+t-1),2])
TSX60bdvol_1 = sdmat[1:((52+t-1)-1),2]*(retmat[1:((52+t-1)-1),3]<ermat[1:((52+t-1)-1),2])
Nikkeibdvol = sdmat[2:(52+t-1),3]*(retmat[2:(52+t-1),4]<ermat[2:(52+t-1),3])
Nikkeibdvol_1 = sdmat[1:((52+t-1)-1),3]*(retmat[1:((52+t-1)-1),4]<ermat[1:((52+t-1)-1),3])
MIB30bdvol = sdmat[2:(52+t-1),4]*(retmat[2:(52+t-1),5]<ermat[2:(52+t-1),4])
MIB30bdvol_1 = sdmat[1:((52+t-1)-1),4]*(retmat[1:((52+t-1)-1),5]<ermat[1:((52+t-1)-1),4])
DAX40bdvol = sdmat[2:(52+t-1),5]*(retmat[2:(52+t-1),6]<ermat[2:(52+t-1),5])
DAX40bdvol_1 = sdmat[1:((52+t-1)-1),5]*(retmat[1:((52+t-1)-1),6]<ermat[1:((52+t-1)-1),5])
CAC40bdvol = sdmat[2:(52+t-1),6]*(retmat[2:(52+t-1),7]<ermat[2:(52+t-1),6])
CAC40bdvol_1 = sdmat[1:((52+t-1)-1),6]*(retmat[1:((52+t-1)-1),7]<ermat[1:((52+t-1)-1),6])

# Vector Autoregressive Part

var1 <- lm(sp500bdvol~sp500bdvol_1+TSX60bdvol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40bdvol_1+CAC40bdvol_1)

var2 <- lm(TSX60bdvol~sp500bdvol_1+TSX60bdvol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40bdvol_1+CAC40bdvol_1)

var3 <- lm(Nikkeibdvol~sp500bdvol_1+TSX60bdvol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40bdvol_1+CAC40bdvol_1)

var4 <- lm(MIB30bdvol~sp500bdvol_1+TSX60bdvol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40vol_1+CAC40bdvol_1)

var5 <- lm(DAX40bdvol~sp500bdvol_1+TSX60vol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40bdvol_1+CAC40bdvol_1)

var6 <- lm(CAC40bdvol~sp500bdvol_1+TSX60bdvol_1+Nikkeibdvol_1
           +MIB30bdvol_1+DAX40bdvol_1+CAC40bdvol_1)

PHI1 = matrix(c(summary(var1)$coefficients[2:7,1],
                summary(var2)$coefficients[2:7,1],
                summary(var3)$coefficients[2:7,1],
                summary(var4)$coefficients[2:7,1],
                summary(var5)$coefficients[2:7,1],
                summary(var6)$coefficients[2:7,1]
),6,6
)

PHI1 = t(PHI1)

residual = data.frame(var1$residuals,var2$residuals,var3$residuals,var4$residuals,var5$residuals,var6$residuals)

Sigma = cov(residual)

A0 = diag(6)

A1 = PHI1 %*% A0
A2 = PHI1 %*% A1
A3 = PHI1 %*% A2

e1 = matrix(c(1,0,0,0,0,0),6,1)
e2 = matrix(c(0,1,0,0,0,0),6,1)
e3 = matrix(c(0,0,1,0,0,0),6,1)
e4 = matrix(c(0,0,0,1,0,0),6,1)
e5 = matrix(c(0,0,0,0,1,0),6,1)
e6 = matrix(c(0,0,0,0,0,1),6,1)

# SP500

theta_bd_11 = sqrt(Sigma[1,1])^(-1)*((t(e1)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e1)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_12 = sqrt(Sigma[2,2])^(-1)*((t(e1)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e2)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_13 = sqrt(Sigma[3,3])^(-1)*((t(e1)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e3)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_14 = sqrt(Sigma[4,4])^(-1)*((t(e1)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e4)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_15 = sqrt(Sigma[5,5])^(-1)*((t(e1)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e5)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_16 = sqrt(Sigma[6,6])^(-1)*((t(e1)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e1)%*%A3%*%Sigma%*%e6)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd_1st = c(theta_bd_11/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16),
                 theta_bd_12/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16),
                 theta_bd_13/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16),
                 theta_bd_14/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16),
                 theta_bd_15/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16),
                 theta_bd_16/(theta_bd_11+theta_bd_12+theta_bd_13+theta_bd_14+theta_bd_15+theta_bd_16))

# TSX60

theta_bd_21 = sqrt(Sigma[1,1])^(-1)*((t(e2)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e1)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_22 = sqrt(Sigma[2,2])^(-1)*((t(e2)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e2)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_23 = sqrt(Sigma[3,3])^(-1)*((t(e2)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e3)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_24 = sqrt(Sigma[4,4])^(-1)*((t(e2)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e4)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_25 = sqrt(Sigma[5,5])^(-1)*((t(e2)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e5)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_26 = sqrt(Sigma[6,6])^(-1)*((t(e2)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e2)%*%A3%*%Sigma%*%e6)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd_2st = c(theta_bd_21/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26),
                 theta_bd_22/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26),
                 theta_bd_23/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26),
                 theta_bd_24/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26),
                 theta_bd_25/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26),
                 theta_bd_26/(theta_bd_21+theta_bd_22+theta_bd_23+theta_bd_24+theta_bd_25+theta_bd_26))

# Nikkei

theta_bd_31 = sqrt(Sigma[1,1])^(-1)*((t(e3)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e1)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_32 = sqrt(Sigma[2,2])^(-1)*((t(e3)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e2)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_33 = sqrt(Sigma[3,3])^(-1)*((t(e3)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e3)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_34 = sqrt(Sigma[4,4])^(-1)*((t(e3)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e4)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_35 = sqrt(Sigma[5,5])^(-1)*((t(e3)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e5)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_36 = sqrt(Sigma[6,6])^(-1)*((t(e3)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e3)%*%A3%*%Sigma%*%e6)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd_3st = c(theta_bd_31/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36),
                 theta_bd_32/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36),
                 theta_bd_33/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36),
                 theta_bd_34/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36),
                 theta_bd_35/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36),
                 theta_bd_36/(theta_bd_31+theta_bd_32+theta_bd_33+theta_bd_34+theta_bd_35+theta_bd_36))

# MIB30

theta_bd_41 = sqrt(Sigma[1,1])^(-1)*((t(e4)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e1)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_42 = sqrt(Sigma[2,2])^(-1)*((t(e4)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e2)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_43 = sqrt(Sigma[3,3])^(-1)*((t(e4)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e3)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_44 = sqrt(Sigma[4,4])^(-1)*((t(e4)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e4)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_45 = sqrt(Sigma[5,5])^(-1)*((t(e4)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e5)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_46 = sqrt(Sigma[6,6])^(-1)*((t(e4)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e4)%*%A3%*%Sigma%*%e6)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd_4st = c(theta_bd_41/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46),
                 theta_bd_42/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46),
                 theta_bd_43/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46),
                 theta_bd_44/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46),
                 theta_bd_45/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46),
                 theta_bd_46/(theta_bd_41+theta_bd_42+theta_bd_43+theta_bd_44+theta_bd_45+theta_bd_46))

# DAX40

theta_bd_51 = sqrt(Sigma[1,1])^(-1)*((t(e5)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e1)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_52 = sqrt(Sigma[2,2])^(-1)*((t(e5)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e2)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_53 = sqrt(Sigma[3,3])^(-1)*((t(e5)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e3)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_54 = sqrt(Sigma[4,4])^(-1)*((t(e5)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e4)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_55 = sqrt(Sigma[5,5])^(-1)*((t(e5)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e5)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_56 = sqrt(Sigma[6,6])^(-1)*((t(e5)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e5)%*%A3%*%Sigma%*%e6)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd_5st = c(theta_bd_51/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56),
                 theta_bd_52/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56),
                 theta_bd_53/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56),
                 theta_bd_54/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56),
                 theta_bd_55/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56),
                 theta_bd_56/(theta_bd_51+theta_bd_52+theta_bd_53+theta_bd_54+theta_bd_55+theta_bd_56))

# CAC40

theta_bd_61 = sqrt(Sigma[1,1])^(-1)*((t(e6)%*%A0%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e1)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e1)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_62 = sqrt(Sigma[2,2])^(-1)*((t(e6)%*%A0%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e2)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e2)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_63 = sqrt(Sigma[3,3])^(-1)*((t(e6)%*%A0%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e3)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e3)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_64 = sqrt(Sigma[4,4])^(-1)*((t(e6)%*%A0%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e4)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e4)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_65 = sqrt(Sigma[5,5])^(-1)*((t(e6)%*%A0%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e5)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e5)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_66 = sqrt(Sigma[6,6])^(-1)*((t(e6)%*%A0%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A1%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A2%*%Sigma%*%e6)^2
                                     +(t(e6)%*%A3%*%Sigma%*%e6)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd_6st = c(theta_bd_61/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66),
                 theta_bd_62/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66),
                 theta_bd_63/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66),
                 theta_bd_64/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66),
                 theta_bd_65/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66),
                 theta_bd_66/(theta_bd_61+theta_bd_62+theta_bd_63+theta_bd_64+theta_bd_65+theta_bd_66))

# Spillover to all other 

theta1all_bd = sum(theta_bd_1st[2:6])
theta2all_bd = sum(theta_bd_2st[1:6])-theta_bd_2st[2]
theta3all_bd = sum(theta_bd_3st[1:6])-theta_bd_3st[3]
theta4all_bd = sum(theta_bd_4st[1:6])-theta_bd_4st[4]
theta5all_bd = sum(theta_bd_5st[1:6])-theta_bd_5st[5]
theta6all_bd = sum(theta_bd_6st[1:6])-theta_bd_6st[6]

# Spillover from all other

thetaall1_bd = theta_bd_2st[1]+theta_bd_3st[1]+theta_bd_4st[1]+theta_bd_5st[1]+theta_bd_6st[1]
thetaall2_bd = theta_bd_1st[2]+theta_bd_3st[2]+theta_bd_4st[2]+theta_bd_5st[2]+theta_bd_6st[2]
thetaall3_bd = theta_bd_1st[3]+theta_bd_2st[3]+theta_bd_4st[3]+theta_bd_5st[3]+theta_bd_6st[3]
thetaall4_bd = theta_bd_1st[4]+theta_bd_2st[4]+theta_bd_3st[4]+theta_bd_5st[4]+theta_bd_6st[4]
thetaall5_bd = theta_bd_1st[5]+theta_bd_2st[5]+theta_bd_3st[5]+theta_bd_4st[5]+theta_bd_6st[5]
thetaall6_bd = theta_bd_1st[6]+theta_bd_2st[6]+theta_bd_3st[6]+theta_bd_4st[6]+theta_bd_5st[6]

# Bad Volatility with 0 threshold


sp500bd0vol = sdmat[2:(52+t-1),1]*(retmat[2:(52+t-1),2]<0)
sp500bd0vol_1 = sdmat[1:((52+t-1)-1),1]*(retmat[1:((52+t-1)-1),2]<0)
TSX60bd0vol = sdmat[2:(52+t-1),2]*(retmat[2:(52+t-1),3]<0)
TSX60bd0vol_1 = sdmat[1:((52+t-1)-1),2]*(retmat[1:((52+t-1)-1),3]<0)
Nikkeibd0vol = sdmat[2:(52+t-1),3]*(retmat[2:(52+t-1),4]<0)
Nikkeibd0vol_1 = sdmat[1:((52+t-1)-1),3]*(retmat[1:((52+t-1)-1),4]<0)
MIB30bd0vol = sdmat[2:(52+t-1),4]*(retmat[2:(52+t-1),5]<0)
MIB30bd0vol_1 = sdmat[1:((52+t-1)-1),4]*(retmat[1:((52+t-1)-1),5]<0)
DAX40bd0vol = sdmat[2:(52+t-1),5]*(retmat[2:(52+t-1),6]<0)
DAX40bd0vol_1 = sdmat[1:((52+t-1)-1),5]*(retmat[1:((52+t-1)-1),6]<0)
CAC40bd0vol = sdmat[2:(52+t-1),6]*(retmat[2:(52+t-1),7]<0)
CAC40bd0vol_1 = sdmat[1:((52+t-1)-1),6]*(retmat[1:((52+t-1)-1),7]<0)

var1 <- lm(sp500bd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

var2 <- lm(TSX60bd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

var3 <- lm(Nikkeibd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

var4 <- lm(MIB30bd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

var5 <- lm(DAX40bd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

var6 <- lm(CAC40bd0vol~sp500bd0vol_1+TSX60bd0vol_1+Nikkeibd0vol_1
           +MIB30bd0vol_1+DAX40bd0vol_1+CAC40bd0vol_1)

PHI1 = matrix(c(summary(var1)$coefficients[2:7,1],
                summary(var2)$coefficients[2:7,1],
                summary(var3)$coefficients[2:7,1],
                summary(var4)$coefficients[2:7,1],
                summary(var5)$coefficients[2:7,1],
                summary(var6)$coefficients[2:7,1]
),6,6)

PHI1 = t(PHI1)

residual = data.frame(var1$residuals,var2$residuals,var3$residuals,var4$residuals,var5$residuals,var6$residuals)

Sigma = cov(residual)

A0 = diag(6)

A1 = PHI1 %*% A0
A2 = PHI1 %*% A1
A3 = PHI1 %*% A2

e1 = matrix(c(1,0,0,0,0,0),6,1)
e2 = matrix(c(0,1,0,0,0,0),6,1)
e3 = matrix(c(0,0,1,0,0,0),6,1)
e4 = matrix(c(0,0,0,1,0,0),6,1)
e5 = matrix(c(0,0,0,0,1,0),6,1)
e6 = matrix(c(0,0,0,0,0,1),6,1)

# SP500

theta_bd0_11 = sqrt(Sigma[1,1])^(-1)*((t(e1)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e1)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd0_12 = sqrt(Sigma[2,2])^(-1)*((t(e1)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e2)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd0_13 = sqrt(Sigma[3,3])^(-1)*((t(e1)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e3)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd0_14 = sqrt(Sigma[4,4])^(-1)*((t(e1)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e4)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))

theta_bd0_15 = sqrt(Sigma[5,5])^(-1)*((t(e1)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e5)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))


theta_bd0_16 = sqrt(Sigma[6,6])^(-1)*((t(e1)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e1)%*%A3%*%Sigma%*%e6)^2)/
  (t(e1)%*%A0%*%Sigma%*%A0%*%e1+(t(e1)%*%A1%*%Sigma%*%A0%*%e1)
   +(t(e1)%*%A2%*%Sigma%*%A0%*%e1)+(t(e1)%*%A3%*%Sigma%*%A0%*%e1))


theta_bd0_1st = c(theta_bd0_11/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16),
                  theta_bd0_12/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16),
                  theta_bd0_13/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16),
                  theta_bd0_14/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16),
                  theta_bd0_15/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16),
                  theta_bd0_16/(theta_bd0_11+theta_bd0_12+theta_bd0_13+theta_bd0_14+theta_bd0_15+theta_bd0_16))

# TSX60

theta_bd0_21 = sqrt(Sigma[1,1])^(-1)*((t(e2)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e1)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_22 = sqrt(Sigma[2,2])^(-1)*((t(e2)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e2)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_23 = sqrt(Sigma[3,3])^(-1)*((t(e2)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e3)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_24 = sqrt(Sigma[4,4])^(-1)*((t(e2)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e4)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_25 = sqrt(Sigma[5,5])^(-1)*((t(e2)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e5)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_26 = sqrt(Sigma[6,6])^(-1)*((t(e2)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e2)%*%A3%*%Sigma%*%e6)^2)/
  (t(e2)%*%A0%*%Sigma%*%A0%*%e2+(t(e2)%*%A1%*%Sigma%*%A0%*%e2)
   +(t(e2)%*%A2%*%Sigma%*%A0%*%e2)+(t(e2)%*%A3%*%Sigma%*%A0%*%e2))

theta_bd0_2st = c(theta_bd0_21/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26),
                  theta_bd0_22/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26),
                  theta_bd0_23/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26),
                  theta_bd0_24/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26),
                  theta_bd0_25/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26),
                  theta_bd0_26/(theta_bd0_21+theta_bd0_22+theta_bd0_23+theta_bd0_24+theta_bd0_25+theta_bd0_26))

# Nikkei

theta_bd0_31 = sqrt(Sigma[1,1])^(-1)*((t(e3)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e1)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_32 = sqrt(Sigma[2,2])^(-1)*((t(e3)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e2)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_33 = sqrt(Sigma[3,3])^(-1)*((t(e3)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e3)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_34 = sqrt(Sigma[4,4])^(-1)*((t(e3)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e4)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_35 = sqrt(Sigma[5,5])^(-1)*((t(e3)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e5)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_36 = sqrt(Sigma[6,6])^(-1)*((t(e3)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e3)%*%A3%*%Sigma%*%e6)^2)/
  (t(e3)%*%A0%*%Sigma%*%A0%*%e3+(t(e3)%*%A1%*%Sigma%*%A0%*%e3)
   +(t(e3)%*%A2%*%Sigma%*%A0%*%e3)+(t(e3)%*%A3%*%Sigma%*%A0%*%e3))

theta_bd0_3st = c(theta_bd0_31/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36),
                  theta_bd0_32/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36),
                  theta_bd0_33/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36),
                  theta_bd0_34/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36),
                  theta_bd0_35/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36),
                  theta_bd0_36/(theta_bd0_31+theta_bd0_32+theta_bd0_33+theta_bd0_34+theta_bd0_35+theta_bd0_36))

# MIB30

theta_bd0_41 = sqrt(Sigma[1,1])^(-1)*((t(e4)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e1)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_42 = sqrt(Sigma[2,2])^(-1)*((t(e4)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e2)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_43 = sqrt(Sigma[3,3])^(-1)*((t(e4)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e3)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_44 = sqrt(Sigma[4,4])^(-1)*((t(e4)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e4)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_45 = sqrt(Sigma[5,5])^(-1)*((t(e4)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e5)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_46 = sqrt(Sigma[6,6])^(-1)*((t(e4)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e4)%*%A3%*%Sigma%*%e6)^2)/
  (t(e4)%*%A0%*%Sigma%*%A0%*%e4+(t(e4)%*%A1%*%Sigma%*%A0%*%e4)
   +(t(e4)%*%A2%*%Sigma%*%A0%*%e4)+(t(e4)%*%A3%*%Sigma%*%A0%*%e4))

theta_bd0_4st = c(theta_bd0_41/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46),
                  theta_bd0_42/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46),
                  theta_bd0_43/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46),
                  theta_bd0_44/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46),
                  theta_bd0_45/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46),
                  theta_bd0_46/(theta_bd0_41+theta_bd0_42+theta_bd0_43+theta_bd0_44+theta_bd0_45+theta_bd0_46))

# DAX40

theta_bd0_51 = sqrt(Sigma[1,1])^(-1)*((t(e5)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e1)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_52 = sqrt(Sigma[2,2])^(-1)*((t(e5)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e2)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_53 = sqrt(Sigma[3,3])^(-1)*((t(e5)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e3)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_54 = sqrt(Sigma[4,4])^(-1)*((t(e5)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e4)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_55 = sqrt(Sigma[5,5])^(-1)*((t(e5)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e5)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_56 = sqrt(Sigma[6,6])^(-1)*((t(e5)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e5)%*%A3%*%Sigma%*%e6)^2)/
  (t(e5)%*%A0%*%Sigma%*%A0%*%e4+(t(e5)%*%A1%*%Sigma%*%A0%*%e5)
   +(t(e5)%*%A2%*%Sigma%*%A0%*%e4)+(t(e5)%*%A3%*%Sigma%*%A0%*%e5))

theta_bd0_5st = c(theta_bd0_51/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56),
                  theta_bd0_52/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56),
                  theta_bd0_53/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56),
                  theta_bd0_54/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56),
                  theta_bd0_55/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56),
                  theta_bd0_56/(theta_bd0_51+theta_bd0_52+theta_bd0_53+theta_bd0_54+theta_bd0_55+theta_bd0_56))

# CAC40

theta_bd0_61 = sqrt(Sigma[1,1])^(-1)*((t(e6)%*%A0%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e1)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e1)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_62 = sqrt(Sigma[2,2])^(-1)*((t(e6)%*%A0%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e2)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e2)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_63 = sqrt(Sigma[3,3])^(-1)*((t(e6)%*%A0%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e3)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e3)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_64 = sqrt(Sigma[4,4])^(-1)*((t(e6)%*%A0%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e4)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e4)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_65 = sqrt(Sigma[5,5])^(-1)*((t(e6)%*%A0%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e5)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e5)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_66 = sqrt(Sigma[6,6])^(-1)*((t(e6)%*%A0%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A1%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A2%*%Sigma%*%e6)^2
                                      +(t(e6)%*%A3%*%Sigma%*%e6)^2)/
  (t(e6)%*%A0%*%Sigma%*%A0%*%e6+(t(e6)%*%A1%*%Sigma%*%A0%*%e6)
   +(t(e6)%*%A2%*%Sigma%*%A0%*%e6)+(t(e6)%*%A3%*%Sigma%*%A0%*%e6))

theta_bd0_6st = c(theta_bd0_61/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66),
                  theta_bd0_62/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66),
                  theta_bd0_63/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66),
                  theta_bd0_64/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66),
                  theta_bd0_65/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66),
                  theta_bd0_66/(theta_bd0_61+theta_bd0_62+theta_bd0_63+theta_bd0_64+theta_bd0_65+theta_bd0_66))

# Spillover to all other 

theta1all_bd0 = sum(theta_bd0_1st[2:6])
theta2all_bd0 = sum(theta_bd0_2st[1:6])-theta_bd0_2st[2]
theta3all_bd0 = sum(theta_bd0_3st[1:6])-theta_bd0_3st[3]
theta4all_bd0 = sum(theta_bd0_4st[1:6])-theta_bd0_4st[4]
theta5all_bd0 = sum(theta_bd0_5st[1:6])-theta_bd0_5st[5]
theta6all_bd0 = sum(theta_bd0_6st[1:6])-theta_bd0_6st[6]

# Spillover from all other

thetaall1_bd0 = theta_bd0_2st[1]+theta_bd0_3st[1]+theta_bd0_4st[1]+theta_bd0_5st[1]+theta_bd0_6st[1]
thetaall2_bd0 = theta_bd0_1st[2]+theta_bd0_3st[2]+theta_bd0_4st[2]+theta_bd0_5st[2]+theta_bd0_6st[2]
thetaall3_bd0 = theta_bd0_1st[3]+theta_bd0_2st[3]+theta_bd0_4st[3]+theta_bd0_5st[3]+theta_bd0_6st[3]
thetaall4_bd0 = theta_bd0_1st[4]+theta_bd0_2st[4]+theta_bd0_3st[4]+theta_bd0_5st[4]+theta_bd0_6st[4]
thetaall5_bd0 = theta_bd0_1st[5]+theta_bd0_2st[5]+theta_bd0_3st[5]+theta_bd0_4st[5]+theta_bd0_6st[5]
thetaall6_bd0 = theta_bd0_1st[6]+theta_bd0_2st[6]+theta_bd0_3st[6]+theta_bd0_4st[6]+theta_bd0_5st[6]

# Differences in two different methods

theta1mat[t,]=matrix(c(theta1all,thetaall1,theta1all_gd,thetaall1_gd,theta1all_gd0,thetaall1_gd0,theta1all_bd,thetaall1_bd,theta1all_bd0,thetaall1_bd0),1,10)
theta2mat[t,]=matrix(c(theta2all,thetaall2,theta2all_gd,thetaall2_gd,theta2all_gd0,thetaall2_gd0,theta2all_bd,thetaall2_bd,theta2all_bd0,thetaall2_bd0),1,10)
theta3mat[t,]=matrix(c(theta3all,thetaall3,theta3all_gd,thetaall3_gd,theta3all_gd0,thetaall3_gd0,theta3all_bd,thetaall3_bd,theta3all_bd0,thetaall3_bd0),1,10)
theta4mat[t,]=matrix(c(theta4all,thetaall4,theta4all_gd,thetaall4_gd,theta4all_gd0,thetaall4_gd0,theta4all_bd,thetaall4_bd,theta4all_bd0,thetaall4_bd0),1,10)
theta5mat[t,]=matrix(c(theta5all,thetaall5,theta5all_gd,thetaall5_gd,theta5all_gd0,thetaall5_gd0,theta5all_bd,thetaall5_bd,theta5all_bd0,thetaall5_bd0),1,10)
theta6mat[t,]=matrix(c(theta6all,thetaall6,theta6all_gd,thetaall6_gd,theta6all_gd0,thetaall6_gd0,theta6all_bd,thetaall6_bd,theta6all_bd0,thetaall6_bd0),1,10)
}