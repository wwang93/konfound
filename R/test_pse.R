test_pse <- function(est_eff,
                     std_err,
                     n_obs, # input df for now
                     eff_thr,
                     sdx,
                     sdy,
                     R2){
    
    ## test_pse(est_eff = .5, std_err = .056, n_obs = 6174, 
    ##         eff_thr = .1, sdx = 0.22, sdy = 1, R2 = .3,to_return = "full")
    ## prepare input
    var_x = sdx^2
    var_y = sdy^2
    var_z = sdz = 1
    
    ## now standardize 
    beta_thr = eff_thr * sdx / sdy
    beta = est_eff * sdx / sdy
    SE = std_err * sdx / sdy
    
    ## observed regression, reg y on x Given z
    tyxGz = beta / SE  ### this should be equal to est_eff / std_err
    ryxGz = tyxGz / sqrt(n_obs + tyxGz^2) 
    
    ## make sure R2 due to x alone is not larger than overall or observed R2
    if (ryxGz^2 > R2) {stop("Error! ryxGz^2 > R2")}
    
    ## calculate ryz, rxz, rxy
    ryz = rzy = cal_ryz(ryxGz, R2)
    rxz = rzx = cal_rxz(var_x, var_y, R2, n_obs, std_err)
    rxy = ryx = cal_rxy(ryxGz, rxz, ryz)
    
    thr = eff_thr * sdx / sdy
    sdz = sdcv = 1
    rcvz = rzcv = 0
    
    # to record how many & which solutions we have are valid
    solvalid = c(F, F, F)
    finaloutput = list()
    
    Gz_pse1 <- cal_pse1(thr, ryxGz)
    if (is.list(Gz_pse1)) {
        rxcvGz1 = as.numeric(Gz_pse1[1])
        rycvGz1 = as.numeric(Gz_pse1[2])
        solvalid[1] = T
        rawoutput1 = gen_fTable(rxcvGz1, rycvGz1, rcvz, rxy, rxz, rzy, 
                             n_obs, sdx, sdy, sdz, sdcv, 
                             ryxGz)
        finaloutput = c(finaloutput,rawoutput1)
    }
    
    Gz_pse2 <- cal_pse2(thr, ryxGz)
    if (is.list(Gz_pse2)) {
        rxcvGz2 = as.numeric(Gz_pse2[1])
        rycvGz2 = as.numeric(Gz_pse2[2])
        solvalid[2] = T
        rawoutput2 = gen_fTable(rxcvGz2, rycvGz2, rcvz, rxy, rxz, rzy, 
                             n_obs, sdx, sdy, sdz, sdcv, 
                             ryxGz)
        finaloutput = c(finaloutput,rawoutput2)
    }
    
    Gz_pse3 <- cal_pse3(thr, ryxGz)
    if (is.list(Gz_pse3)) {
        rxcvGz3 = as.numeric(Gz_pse3[1])
        rycvGz3 = as.numeric(Gz_pse3[2])
        solvalid[3] = T
        rawoutput3 = gen_fTable(rxcvGz3, rycvGz3, rcvz, rxy, rxz, rzy, 
                             n_obs, sdx, sdy, sdz, sdcv, 
                             ryxGz)
        finaloutput = c(finaloutput,rawoutput3)
    }
    
    if (sum(solvalid) > 0) {
        cat("This function calculates the conditions that set the estimated effect approximately equal to the threshold while preserving the standard error.")
        cat("\n")
        cat(sprintf("There are %i possible solutions, as specified below.", sum(solvalid)))
        cat("\n")
        return(finaloutput)
    }
    
    
    # if (to_return == "print") {
    #    cat("This function calculates the conditions that set the estimated effect approximately equal to the threshold while preserving the standard error.")
    #    cat("\n")
    #    cat(sprintf("The correlation between X and CV is %.3f, and the correlation between Y and CV is %.3f.", rxcv, rycv))
    #    cat("\n")
    #    cat(sprintf("Conditional on the covariates, the correlation between X and CV is %.3f, and the correlation between Y and CV is %.3f.", rxcvGz, rycvGz))
    #    cat("\n")
    #    cat(sprintf("Including such CV, the coefficient changes to %.3f, and standard error is %.3f.", eff_x_M3, se_x_M3))
    #    cat("\n")
    #    cat("Use to_return = raw_ouput to see more specific results.")
    #}
        
}
    