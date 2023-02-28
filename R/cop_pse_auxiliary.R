cal_ryz <- function(ryxGz, R2){
    R2yz = (ryxGz^2 - R2)/(ryxGz^2 - 1)
    if (R2yz >= 0) {
        ryz = sqrt(R2yz)
    } else {
        stop("Error! R2yz < 0!")
    }
    return(ryz)
}

cal_rxz <- function(var_x, var_y, R2, df, std_err){
    R2xz = 1 - ((var_y * (1 - R2))/(var_x * df * std_err^2))
    if (R2xz <= 0) {stop("Error! R2xz < 0!")} 
    ## Note output the number of R2xz
    rxz = sqrt(R2xz)
    return(rxz)
}

cal_rxy <- function(ryxGz, rxz, ryz){
    rxy = ryxGz * sqrt((1 - rxz^2)*(1 - ryz^2)) + rxz * ryz
    return(rxy)
}

cal_delta_star <- function(FR2max, R2, R2_uncond, est_eff, eff_thr, var_x, var_y, est_uncond, rxz, n_obs){
    if (FR2max > .99) {FR2max = .99}
    # if (FR2max < R2 + inci) {FR2max = R2 + inci} check with Ken what this means
    if (FR2max > R2) {D = sqrt(FR2max - R2)}
    
    #elements for computing Oster's delta_star
    bt_m_b = est_eff - eff_thr
    rt_m_ro_t_syy = (R2 - R2_uncond) * var_y
    b0_m_b1 = est_uncond - est_eff
    rm_m_rt_t_syy = (FR2max - R2) * var_y
    
    t_x = var_x * (n_obs / (n_obs - 1)) * (1 - rxz^2)
    ## adjust df for var_x 
    ## var_x is population variance, need sample variance from x
    ## this adjustment is to get closer to what robomit generates as they run regression using the sample data 
    num1 = bt_m_b * rt_m_ro_t_syy * t_x
    num2 = bt_m_b * var_x * t_x * b0_m_b1^2
    num3 = 2 * bt_m_b^2 * (t_x * b0_m_b1 * var_x)
    num4 = bt_m_b^3 * (t_x * var_x - t_x^2)
    num = num1 + num2 + num3 + num4
    den1 = rm_m_rt_t_syy * b0_m_b1 * var_x
    den2 = bt_m_b * rm_m_rt_t_syy * (var_x - t_x)
    den3 = bt_m_b^2 * (t_x * b0_m_b1 * var_x)
    den4 = bt_m_b^3 * (t_x * var_x - t_x^2)
    den = den1 + den2 + den3 + den4
    #obtain delta_star which is Oster's delta
    delta_star = num / den
    return(delta_star)
}

# this is the original approach
# cal_delta_exact <- function(ryx, ryz, rxz, beta_thr, FR2max, R2, sdx, sdz){
#     # setting up simple values to align with Mathematica
#     y = ryx
#     z = ryz
#     w = rxz
#     b = beta_thr ## check with Ken, is beta_thr here or eff_thr here? 
#     v = 1 - rxz^2 # this is to simplify calculation later
#     D = sqrt(FR2max - R2) # same above
#     
#     kden1 = sqrt(b^2 * v^2 + 
#                      2 * b * v * w * z - 
#                      2 * b * v * y + 
#                      D^2 * v + 
#                      w^2 * z^2 -
#                      2 * w * y * z +
#                      y^2)
#     knum1 = sqrt(v) * (b * v + w * z - y)
#     rcvx = rxcv = - knum1 / kden1
# 
#     # If regular conditions hold, calculate rycv from rxcv 
#     yescalc = 0
#     if (abs(rcvx) < 1 && (rcvx^2 / v) < 1){
#         yescalc = 1
#         rcvy = rycv = 
#         D * sqrt(1 - (rcvx^2 / v)) +
#         (ryx * rcvx) / (v) -
#         (ryz * rcvx * rxz) / (v)
#         delta_exact = rcvx / rxz
#         result = list(rxcv, rycv, delta_exact)
#     }
#     
#     if (yescalc == 1) {
#         return(result)
#     } else {
#         stop("Error!")
#     }
# }

# see test_cop for updated approach to calculate delta exact 

verify_reg_Gzcv = function(n_obs, sdx, sdy, sdz, sdcv, 
                           rxy, rxz, rzy, rcvy, rcvx, rcvz){
    
    model <- 'Y ~ beta1 * X + beta2 * Z + beta3 * CV'
    
    ccvy = rcvy * sdcv * sdy # cov(cv, y)
    ccvx = rcvx * sdcv * sdx # cov(cv, x)
    ccvz = rcvz * sdcv * sdz
    cxy = rxy * sdx * sdy
    czy = rzy * sdz * sdy
    cxz = rxz * sdx * sdz
    
    # set up the covariance matrix
    cov.matrix <- matrix(c(sdy^2, cxy, czy, ccvy,
                           cxy, sdx^2, cxz, ccvx,
                           czy, cxz, sdz^2, ccvz,
                           ccvy, ccvx, ccvz, sdcv^2), 4, 4)
    rownames(cov.matrix) <- colnames(cov.matrix) <- c("Y", "X", "Z", "CV")
    
    # Check if model can be run
    flag_cov <- tryCatch(
        expr = {
            lavaan::sem(model, 
                sample.cov = cov.matrix, 
                sample.nobs = n_obs)
        },
        error = function(e){
            flag_cov = F
            return(flag_cov)
        },
        warning = function(w){
            flag_cov = F
            return(flag_cov)
        }
    )
    #if model can be run to verify true delta, then run it can save results
    if (class(flag_cov) == "lavaan") {
        fit <- lavaan::sem(model,
                      sample.cov = cov.matrix,
                      sample.nobs = n_obs)
        ## the R2 extracted from summary is NOT right, do the calculation below
        R2 <- (sdy^2 - lavaan::parameterEstimates(fit)[4,]$est) / sdy^2
        betaX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$est
        seX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$se
        betaZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$est
        seZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$se
        betaCV <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta3',]$est
        seCV <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta3',]$se
    }

    #get regression based on true delta in terms of standardized coefficent
    cor.matrix <- matrix(c(1,rxy, rzy, rcvy,
                           rxy, 1, rxz, rcvx,
                           rzy, rxz, 1, rcvz,
                           rcvy, rcvx, rcvz, 1), 4, 4)
    rownames(cor.matrix) <- colnames(cor.matrix) <- c("Y", "X", "Z", "CV")
    
    # check to see if model can be run
    flag_cor <- tryCatch(
        expr = {
            lavaan::sem(model, 
                sample.cov = cor.matrix, 
                sample.nobs = n_obs)
        },
        error = function(e){
            flag_cor = F
            return(flag_cor)
        },
        warning = function(w){
            flag_cor = F
            return(flag_cor)
        }
    )
    
    # if model can be run, then run it
    if (class(flag_cor) == "lavaan") {
        fit <- lavaan::sem(model,
                   sample.cov = cor.matrix,
                   sample.nobs = n_obs)
        std_R2 <- 1 - lavaan::parameterEstimates(fit)[4,]$est
        std_betaX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$est
        std_seX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$se
        std_betaZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$est
        std_seZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$se
        std_betaCV <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta3',]$est
        std_seCV <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta3',]$se
    }
    
    if (class(flag_cor) == "lavaan" && class(flag_cov) == "lavaan") {
        result = list(R2, betaX, seX, betaZ, seZ, betaCV, seCV,
                      std_R2, std_betaX, std_seX,
                      cov.matrix, cor.matrix)
        return(result)
    } else {
        stop("Error!")
    }
}

cal_pse1 <- function(thr, kryx){
    # calculations for preserving standard error
    i1 <- complex(real = 1, imaginary = -sqrt(3))
    i2 <- complex(real = 1, imaginary = sqrt(3))
    temp <- -((2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)/(3 * (-1 + kryx^2))) - 
        ((2^(1/3)) * (-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + kryx^2)))/
        (3 * (-1 + kryx^2) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                  12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 + 
                                  6 * thr * kryx^5 + 7 * kryx^6 + sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                                       39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                                       30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 
                                                                                      4 * (-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + 
                                                                                                                                                                   kryx^2))^3)))^(1/3)) + 
        1/(3 * 2^(1/3) * (-1 + kryx^2)) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                               12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 +
                                               6 * thr * kryx^5 + 7 * kryx^6 + 
                                               sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                    39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                    30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 4 * (-(2 - thr^2 + 2 * thr * kryx - 
                                                                                                                                                      2 * kryx^2)^2 + 
                                                                                                                                                    6 * (-1 + kryx^2) * 
                                                                                                                                                    (thr^2 - 2 * thr * kryx + kryx^2))^3)))^(1/3)
    # calculations generate correlations conditional on Z: rxcvGz and rycvGz

    if (Re(temp) < 0) {
        output <- "Not valid."} 
    else {
        rxcvGz_sepreserve <- Re(sqrt(temp))
        if (abs(rxcvGz_sepreserve) <= 1){
            rycvGz_sepreserve <- (kryx - thr * (1 - rxcvGz_sepreserve^2))/rxcvGz_sepreserve
            output <- list(rxcvGz_sepreserve, rycvGz_sepreserve)
        } else {
            output <- "Not valid."
        }
    }
    
    return(output)
}

cal_pse2 <- function(thr, kryx){
    # calculations for preserving standard error
    i1 <- complex(real = 1, imaginary = -sqrt(3))
    i2 <- complex(real = 1, imaginary = sqrt(3))
    temp <- -((2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)/(3 * (-1 + kryx^2))) + 
        ((i2)*(-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + kryx^2)))/
        (3 * 2^(2/3) * (-1 + kryx^2) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                            12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 + 
                                            6 * thr * kryx^5 + 7 * kryx^6 + sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                                                 39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                                                 30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 
                                                                                                4 * (-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + 
                                                                                                                                                                             kryx^2))^3)))^(1/3)) - 
        1/(6 * 2^(1/3) * (-1 + kryx^2)) * (i1) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                                      12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 +
                                                      6 * thr * kryx^5 + 7 * kryx^6 + 
                                                      sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                           39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                           30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 4 * (-(2 - thr^2 + 2 * thr * kryx - 
                                                                                                                                                             2 * kryx^2)^2 + 
                                                                                                                                                           6 * (-1 + kryx^2) * 
                                                                                                                                                           (thr^2 - 2 * thr * kryx + kryx^2))^3)))^(1/3)
    # calculations generate correlations conditional on Z: rxcvGz and rycvGz

    if (Re(temp) < 0) {
        output <- "Not valid."} 
    else {
        rxcvGz_sepreserve <- Re(sqrt(temp))
        if (abs(rxcvGz_sepreserve) <= 1){
            rycvGz_sepreserve <- (kryx - thr * (1 - rxcvGz_sepreserve^2))/rxcvGz_sepreserve
            output <- list(rxcvGz_sepreserve, rycvGz_sepreserve)
        } else {
            output <- "Not valid."
        }
    }
    
    return(output)
    
}

cal_pse3 <- function(thr, kryx){
    # calculations for preserving standard error
    i1 <- complex(real = 1, imaginary = -sqrt(3))
    i2 <- complex(real = 1, imaginary = sqrt(3))
    temp <- -((2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)/(3 * (-1 + kryx^2))) + 
        ((i1)*(-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + kryx^2)))/
        (3 * 2^(2/3) * (-1 + kryx^2) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                            12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 + 
                                            6 * thr * kryx^5 + 7 * kryx^6 + sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                                                 39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                                                 30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 
                                                                                                4 * (-(2 - thr^2 + 2 * thr * kryx - 2 * kryx^2)^2 + 6 * (-1 + kryx^2) * (thr^2 - 2 * thr * kryx + 
                                                                                                                                                                             kryx^2))^3)))^(1/3)) - 
        1/(6 * 2^(1/3) * (-1 + kryx^2)) * (i2) * (-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 39 * kryx^2 + 
                                                      12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 30 * kryx^4 - 27 * thr^2 * kryx^4 +
                                                      6 * thr * kryx^5 + 7 * kryx^6 + 
                                                      sqrt(as.complex((-16 + 15 * thr^2 + 6 * thr^4 + 2 * thr^6 - 30 * thr * kryx - 24 * thr^3 * kryx - 12 * thr^5 * kryx + 
                                                                           39 * kryx^2 + 12 * thr^2 * kryx^2 + 18 * thr^4 * kryx^2 + 24 * thr * kryx^3 + 8 * thr^3 * kryx^3 - 
                                                                           30 * kryx^4 - 27 * thr^2 * kryx^4 + 6 * thr * kryx^5 + 7 * kryx^6)^2 + 4 * (-(2 - thr^2 + 2 * thr * kryx - 
                                                                                                                                                             2 * kryx^2)^2 + 
                                                                                                                                                           6 * (-1 + kryx^2) * 
                                                                                                                                                           (thr^2 - 2 * thr * kryx + kryx^2))^3)))^(1/3)
    # calculations generate correlations conditional on Z: rxcvGz and rycvGz
    if (Re(temp) < 0) {
        output <- "Not valid."} 
    else {
        rxcvGz_sepreserve <- Re(sqrt(temp))
        if (abs(rxcvGz_sepreserve) <= 1){
            rycvGz_sepreserve <- (kryx - thr * (1 - rxcvGz_sepreserve^2))/rxcvGz_sepreserve
            output <- list(rxcvGz_sepreserve, rycvGz_sepreserve)
        } else {
            output <- "Not valid."
        }
    }
    
    return(output)
    
}

verify_manual <- function(rxy, rxz, rxcv, ryz, rycv, rzcv, sdy, sdx){
    beta <- (rxy + rycv * rxz * rzcv + ryz * rxcv * rzcv - rxy * rzcv^2 - rycv * rxcv - ryz * rxz) /
        (1 + 2 * rxcv * rzcv * rxz - rxcv^2 - rzcv^2 - rxz^2)
    eff <- beta * sdy / sdx
    return(beta)
}

verify_reg_Gz = function(n_obs, sdx, sdy, sdz, rxy, rxz, rzy){
    
    model <- 'Y ~ beta1 * X + beta2 * Z'
    
    cxy = rxy * sdx * sdy
    czy = rzy * sdz * sdy
    cxz = rxz * sdx * sdz
    
    # set up the covariance matrix
    cov.matrix <- matrix(c(sdy^2, cxy, czy, 
                           cxy, sdx^2, cxz,
                           czy, cxz, sdz^2), 3, 3)
    rownames(cov.matrix) <- colnames(cov.matrix) <- c("Y", "X", "Z")
    
    # Check if model can be run
    flag_cov <- tryCatch(
        expr = {
            lavaan::sem(model, 
                sample.cov = cov.matrix, 
                sample.nobs = n_obs)
        },
        error = function(e){
            flag_cov = F
            return(flag_cov)
        },
        warning = function(w){
            flag_cov = F
            return(flag_cov)
        }
    )
    #if model can be run to verify true delta, then run it can save results
    if (class(flag_cov) == "lavaan") {
        fit <- lavaan::sem(model,
                   sample.cov = cov.matrix,
                   sample.nobs = n_obs)
        ## the R2 extracted from summary is NOT right, do the calculation below
        R2 <- (sdy^2 - lavaan::parameterEstimates(fit)[3,]$est) / sdy^2
        betaX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$est
        seX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$se
        betaZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$est
        seZ <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta2',]$se
   }

    if (class(flag_cov) == "lavaan") {
        result = list(R2, betaX, seX, betaZ, seZ)
        return(result)
    } else {
        stop("Error!")
    }
}

verify_reg_uncond = function(n_obs, sdx, sdy, rxy){
    
    model <- 'Y ~ beta1 * X'
    cxy = rxy * sdx * sdy
    
    # set up the covariance matrix
    cov.matrix <- matrix(c(sdy^2, cxy,  
                           cxy, sdx^2), 2, 2)
    rownames(cov.matrix) <- colnames(cov.matrix) <- c("Y", "X")
    
    # Check if model can be run
    flag_cov <- tryCatch(
        expr = {
            lavaan::sem(model, 
                sample.cov = cov.matrix, 
                sample.nobs = n_obs)
        },
        error = function(e){
            flag_cov = F
            return(flag_cov)
        },
        warning = function(w){
            flag_cov = F
            return(flag_cov)
        }
    )
    #if model can be run to verify true delta, then run it can save results
    if (class(flag_cov) == "lavaan") {
        fit <- lavaan::sem(model,
                   sample.cov = cov.matrix,
                   sample.nobs = n_obs)
        ## the R2 extracted from summary is NOT right, do the calculation below
        R2 <- (sdy^2 - lavaan::parameterEstimates(fit)[2,]$est) / sdy^2
        betaX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$est
        seX <- lavaan::parameterEstimates(fit)[lavaan::parameterEstimates(fit)$label == 'beta1',]$se
    }
    
    if (class(flag_cov) == "lavaan") {
        result = list(R2, betaX, seX)
        return(result)
    } else {
        stop("Error!")
    }
}

gen_fTable <- function(rxcvGz, rycvGz, rcvz, rxy, rxz, rzy, 
                       n_obs, sdx, sdy, sdz, sdcv, 
                       ryxGz){
    # convert conditional correlations to unconditional correlations to be used in new regression
    rxcv = rxcvGz * sqrt((1 - rcvz^2) * (1 - rxz^2)) + rxz * rcvz
    rycv = rycvGz * sqrt((1 - rcvz^2) * (1 - rzy^2)) + rzy * rcvz
    
    verify_pse_reg_M3 = verify_reg_Gzcv(n_obs, sdx, sdy, sdz, sdcv, rxy, rxz, rzy, rycv, rxcv, rcvz)
    verfiy_pse_manual_thr = verify_manual(rxy, rxz, rxcv, rzy, rycv, rcvz, sdy, sdx)
    cov_pse = verify_pse_reg_M3[[11]]
    
    # prepare some other values in the final Table (long output)
    R2_M3 = as.numeric(verify_pse_reg_M3[1])
    eff_x_M3 = as.numeric(verify_pse_reg_M3[2]) # should be equivalent or very close to eff_thr
    se_x_M3 = as.numeric(verify_pse_reg_M3[3])
    beta_x_M3 = as.numeric(verify_pse_reg_M3[9]) # should be equivalent or very close to thr
    t_x_M3 = eff_x_M3 / se_x_M3 
    eff_z_M3 = as.numeric(verify_pse_reg_M3[4])
    se_z_M3 = as.numeric(verify_pse_reg_M3[5])
    eff_cv_M3 = as.numeric(verify_pse_reg_M3[6])
    se_cv_M3 = as.numeric(verify_pse_reg_M3[7])
    
    verify_pse_reg_M2 = verify_reg_Gz(n_obs, sdx, sdy, sdz, rxy, rxz, rzy)
    R2_M2 = as.numeric(verify_pse_reg_M2[1])
    eff_x_M2 = as.numeric(verify_pse_reg_M2[2]) # should be equivalent or very close to est_eff
    se_x_M2 = as.numeric(verify_pse_reg_M2[3]) # should be equivalent or very close to std_err
    eff_z_M2 = as.numeric(verify_pse_reg_M2[4]) 
    se_z_M2 = as.numeric(verify_pse_reg_M2[5]) 
    t_x_M2 = eff_x_M2 / se_x_M2 
    
    verify_pse_reg_M1 = verify_reg_uncond(n_obs, sdx, sdy, rxy)
    R2_M1 = as.numeric(verify_pse_reg_M1[1]) # should be equivalent or very close to rxy^2
    eff_x_M1 = as.numeric(verify_pse_reg_M1[2]) # should be equivalent or very close to rxy*sdy/sdx
    se_x_M1 = as.numeric(verify_pse_reg_M1[3]) 
    t_x_M1 = eff_x_M1 / se_x_M1 
    
    fTable <- matrix(c(R2_M1, R2_M2, R2_M3, # R2 for three reg models
                       eff_x_M1, eff_x_M2, eff_x_M3, # unstd reg coef for X in three reg  models
                       se_x_M1, se_x_M2, se_x_M3, # unstd reg se for X in three reg models
                       rxy, ryxGz, beta_x_M3, # std reg coef for X in three reg models
                       t_x_M1, t_x_M2, t_x_M3, # t values for X in three reg models
                       NA, eff_z_M2, eff_z_M3, # reg coef for Z in three reg models
                       NA, se_z_M2, se_z_M3, # se for Z in three reg models
                       NA, eff_z_M2 / se_z_M2, eff_z_M3 / se_z_M3, # t for Z in three reg models,
                       NA, NA, eff_cv_M3, # reg coef for CV in three reg models
                       NA, NA, se_cv_M3, # se for CV in three reg models
                       NA, NA, eff_cv_M3 / se_cv_M3), # t for CV in three reg models
                     nrow = 11, ncol = 3, byrow = T) 
    
    rownames(fTable) <- c("R2", "coef_X", "SE_X", "std_coef_X", "t_X",
                          "coef_Z", "SE_Z", "t_Z",
                          "coef_CV", "SE_CV", "t_CV")
    
    colnames(fTable) <- c("M1:X", "M2:X,Z", "M3:X,Z,CV")
    
    output <- list("correlation between X and CV conditional on Z" = rxcvGz, 
                   "correlation between Y and CV conditional on Z" = rycvGz, 
                   "correlation between X and CV" = rxcv, 
                   "correlation between Y and CV" = rycv,
                   "covariance matrix" = cov_pse, 
                   "Table" = fTable)
    
    return(output)
}


