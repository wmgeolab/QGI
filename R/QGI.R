library(doParallel)
library(MatchIt)
library(devtools)
library(roxygen2)

#' Plot Propensity Model
#'
#' Matches samples of groups of treated and control observations that have a spatial distribution (i.e., latitude and longitude coordinates), so that matched groups have (1) similar covariate distributions [via a propensity score matching], and (2) do not result in paired matches being geographically proximate.  The package also provides a distance-decay estimate of the treatment effect (i.e., the impact a treatment had on a defined outcome) following the procedure outlined in Runfola and Batra et al. (2020) <https://doi.org/10.3390/su12083225>.
#' @param rDF Takes a dataset with a 'distance' column...
#' @param outcome_var Outcome variable (i.e., under 5 child mortality)
#' @param controlVars List of all independent variables to be used to model propensity and outcome
#' @param sample_density Total number of runs
#' @param lower_dist_bound Smallest distance (in decimal degrees) to calculate results for. Defaults to 0.01.
#' @param upper_dist_bound Largest distance (in decimal degrees) to calculate results for. Defaults to 0.5.
#' @param maximum_match_diff From 0 to 1, the largest difference in match quality allowed. Defaults to 0.9.
#' @keywords
#' @export
#' @examples
#' plotPropensityModel(read.csv("../rDF.csv"), independent_vars, x,x,x, 30*16)

#Allow user to change the name of the distance column.  Right now it is hard coded as "Distance".

QGI <- function(rDF, 
                                controlVars, 
                                outcome_var, 
                                outcome_label, 
                                sample_density, 
                                lower_dist_bound=0.01, 
                                upper_dist_bound=0.5, 
                                maximum_match_diff=0.9, 
                                untreated_thresh=0.3,
                                task = 'gef') {
  
  cores=detectCores()
  cl <- makeForkCluster(cores[1]-1) 
  registerDoParallel(cl)
  
  mc = 0
  #define sampling frame, with as much coverage as possible
  #given a defined density
  
  if(sample_density <= 1){
    sample_density = 2
    print("Sample Density set to minimum allowed (2)")
  }

  #trials = seq(upper_dist_bound, lower_dist_bound, (lower_dist_bound-upper_dist_bound)/(sample_density-1))
  trials = seq(lower_dist_bound, upper_dist_bound, (upper_dist_bound-lower_dist_bound)/(sample_density-1))
  
  mc <- foreach(i = 1:length(trials), .combine=rbind) %dopar% {
    set.seed(i)
    
    rDF$distanceTrue <- rDF$distance
    
    rDF$Treatment <- 999
    
    #Change the distance for the threshold.
    #We actually need to do this a few times, to construct the distance bands.
    dist_thresh <- trials[i] 
    if (task == 'fao'){
      #rDF$Treatment[rDF$distanceTrue > 0.1 & rDF$distanceTrue < dist_thresh] <- 0
      #rDF$Treatment[rDF$distanceTrue > dist_thresh] <- 0
      #rDF$Treatment[rDF$distanceTrue == 0] <- 1
      rDF$Treatment[rDF$distanceTrue <= 0 + dist_thresh] <- 1
      rDF$Treatment[rDF$distanceTrue >= untreated_thresh + dist_thresh] <- 0
    }else {
      rDF$Treatment[rDF$distanceTrue <= dist_thresh] <- 1
      rDF$Treatment[rDF$distanceTrue > (dist_thresh + upper_dist_bound)] <- 0
    }
      
    rDF <- rDF[!rDF$Treatment == 999,]
    
    rDF[, !(names(rDF) == "distance")]
    
    #Choose pscore variables
    pVars <- c(controlVars,"Treatment")
    #Choose analysis variables
    aVars <- c(controlVars,"Treatment","outcome")
    
    #construct propensity score formula
    #the 3 debug lines are required for our data, the else case is the simpler generalized version
    if (task == 'uganda') {
      exclude_some_var <- rep(controlVars)
      exclude_some_var <- exclude_some_var[!exclude_some_var %in% c("toJoin_LSMS_2011_geogAsset_h14q5.y","extracted_data_gdp_grid.na.mean.y")]
      f1 <- as.formula(paste("Treatment", paste(exclude_some_var, collapse = " + "), sep = " ~ "))
    } else {
      f1 <- as.formula(paste("Treatment", paste(controlVars, collapse = " + "), sep = " ~ "))
    }
    
    pscore.Calc <- matchit(f1, data=na.omit(rDF[pVars], cols=pVars), method="nearest", distance="logit")
    
    matched.dta <- match.data(pscore.Calc)
    
    
    #construct outcome modeling formula
    #the last 3 debug lines are required for our data, the else case is the simpler generalized version
    if (task == 'uganda') {
      #we can assume user already has a reasonable $outcome measurement.
      matched.dta$outcome <- (matched.dta$toJoin_LSMS_2011_geogAsset_h14q5.y - matched.dta$toJoin_LSMS_2009_geogAsset_h14q5.y) * 0.00027
      
      exclude_some_var2 <- c("Treatment", controlVars)
      exclude_some_var2 <- exclude_some_var2[!exclude_some_var2 %in% c("toJoin_LSMS_2011_geogAsset_h14q5.y","extracted_data_gdp_grid.na.mean.y")]
      f2 <- as.formula(paste("outcome", paste(exclude_some_var2, collapse = " + "), sep = " ~ "))
    } else {
      idx = as.numeric(rownames(matched.dta))
      matched.dta$outcome <- rDF[idx,][[outcome_var]]
      f2 <- as.formula(paste("outcome", paste(c("Treatment", controlVars), collapse = " + "), sep = " ~ "))
    }
    
    trtModel <- lm(f2, data= na.omit(matched.dta[aVars], cols=pVars))
    
    
    match_diff = abs(summary(pscore.Calc)$sum.matched[[1]][1] - summary(pscore.Calc)$sum.matched[[2]][1])
    
    return( list(i, dist_thresh, trtModel$coef["Treatment"][[1]], nrow(rDF), match_diff, summary(trtModel)$r.squared, coef(summary(trtModel))[2,4], coef(summary(trtModel))[2,2], nrow(matched.dta)))
    
  }

  mc_df <- data.frame(t(matrix(unlist(mc), nrow=length(mc[1,]), byrow=T)))
  names(mc_df) <- c("ItId", "thresh", "coef", "obs", "match_diff", "R2", "TreatSig", "StdError", "SampleSize")
  
  #Drop cases where not enough data was available to conduct the matches
  #according to the user set parameter
  mc_df <- mc_df[!mc_df["thresh"] == 1,]
  
  #Drop cases with a really,really high match_diff
  #Generally indicates no decent matches were made in matching
  mc_df <- mc_df[!mc_df["match_diff"] > maximum_match_diff,]
  
  #Good (Runs with above-average matches, for now):
  #Ignore this I think:
  #good_matches <- mc_df[!mc_df$match_diff > mean(mc_df$match_diff),]
  
  
  #Scale for viz and weighting
  mc_df$matchQual_scale = 1 - ((mc_df$match_diff - min(mc_df$match_diff)) / (max(mc_df$match_diff) - min(mc_df$match_diff)))
  
  mc_df$size = (mc_df$matchQual_scale) # + mc_df$sampleSize_scale) / 2
  
  #Change Threshold from Degrees to KM for interpretation
  mc_df$thresh_km = mc_df$thresh * 110.567
  #mc_df$thresh_km = mc_df$thresh*0.001
  
  mc_df$b = mc_df$thresh_km**2
  mc_df$c = mc_df$thresh_km**3
  mc_df$d = mc_df$thresh_km**4
  
  mc_df <- mc_df[order(mc_df$thresh_km),]
  
  #Weighted function based on match quality
  mean_mdl = lm(coef ~ thresh_km + b + c, data = mc_df, weights =size)
  std_mdl = lm(StdError ~ thresh_km + b + c, data = mc_df, weights =size)
  #mean_mdl = lm(coef ~ thresh_km + b + c, data = mc_df)
  #std_mdl = lm(StdError ~ thresh_km + b + c, data = mc_df)
 
  par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = FALSE)        
  
  ylim_upper <- max(mc_df$coef + 1.96*mc_df$StdError, na.rm = TRUE)
  ylim_lower <- min(mc_df$coef - 1.96*mc_df$StdError, na.rm = TRUE)
  plot(mc_df$thresh_km, mc_df$coef, cex = mc_df$size, xlab="Distance (km)", ylab=outcome_label, ylim=c(min(0,ylim_lower), max(ylim_upper,0)))
  #plot(mc_df$thresh_km, mc_df$coef, cex = mc_df$size, xlab="Distance (km)", ylab=outcome_label, ylim=c(-0.5,0.5))
  
  newx <- seq(min(mc_df$thresh_km), max(mc_df$thresh_km), length.out=1000)
  
  if (task == 'fao'){
    newx <- rev(newx)
  }
  
  preds <- predict(mean_mdl, newdata=data.frame(thresh_km=newx, b=newx**2, c=newx**3), interval='confidence')
  preds_std <- predict(std_mdl, newdata=data.frame(thresh_km=newx, b=newx**2, c=newx**3), interval='confidence')
  polygon(c(rev(newx), newx), c(rev( preds[ ,3] + 1.96*preds_std[1:1000]),preds[ ,2] - 1.96*preds_std[1:1000]), col=rgb(0.2, 0.2, 0.25,0.25), border = NA)
  
  #polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col=rgb(0, 1, 0,0.25), border = NA)
  lines(newx, preds[ ,3] + 1.96*preds_std[1:1000], lty = 'dashed', col = 'black')
  lines(newx, preds[ ,2] - 1.96*preds_std[1:1000], lty = 'dashed', col = 'black')
  lines(newx, preds[ ,3], lty = 'dashed', col = 'yellow')
  lines(newx, preds[ ,2], lty = 'dashed', col = 'yellow')
  legend(5, 400, legend=c("Weighted Best Fit", "Lower Match Quality", "Higher Match Quality"), col=c("blue","black", "black"), lty=c(1,NA, NA), pch=c(NA,1,1), cex=0.8, pt.cex=c(NA, 0.5, 1))
  lines(mc_df$thresh_km, fitted(mean_mdl), col="blue")
  abline(h = 0, lty = 2)
  
  #Overall impact 
  print(mean(unlist(mc_df[mc_df$thresh_km >=2.66 & mc_df$thresh_km <=6,]["coef"][1])))
  upper_std = preds[ ,3] + 1.96*preds_std[1:1000]
  print(mean(unlist(mc_df[mc_df$thresh_km >=min(newx[which(upper_std<0)]) & mc_df$thresh_km <=max(newx[which(upper_std<0)]),]["coef"][1])))
  print(paste('significant distance intervel: ',min(newx[which(upper_std<0)]),max(newx[which(upper_std<0)])))
  lower_std = preds[ ,2] - 1.96*preds_std[1:1000]
  print(mean(unlist(mc_df[mc_df$thresh_km >=min(newx[which(lower_std>0)]) & mc_df$thresh_km <=max(newx[which(lower_std>0)]),]["coef"][1])))
  print(paste('significant distance intervel: ',min(newx[which(lower_std>0)]),max(newx[which(lower_std>0)])))
  print(paste('plotted distance range: ',min(mc_df$thresh_km),' ',max(mc_df$thresh_km)))
  
  return(c(mc_df, min(newx[which(lower_std>0)]),max(newx[which(lower_std>0)])))
  
}