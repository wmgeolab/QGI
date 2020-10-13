library(foreach)
library(doParallel)
library(MatchIt)
library(devtools)

#' Plot Propensity Model
#'
#' Matches samples of groups of treated and control observations that have a spatial distribution (i.e., latitude and longitude coordinates), so that matched groups have (1) similar covariate distributions [via a propensity score matching], and (2) do not result in paired matches being geographically proximate.  The package also provides a distance-decay estimate of the treatment effect (i.e., the impact a treatment had on a defined outcome) following the procedure outlined in Runfola and Batra et al. (2020) <https://doi.org/10.3390/su12083225>.
#' @param df A R dataframe containing your observations
#' @param distanceCol The column in your dataframe containing the distance in kilometers between each observation and the intervention locations
#' @param controlVars List of all independent variables to be used to model propensity and outcome
#' @param outcomeVar The column in your dataframe containing the outcome variable 
#' @param density Total number of estimations to be made, equally spaced between the lower and upper distance bound. Defaults to 50.
#' @param lowerDistBound Smallest distance (in kilometers) to calculate results for. Defaults to 0km.
#' @param upperDistBound Largest distance (in kilometers) to calculate results for. Defaults to the median of distance in the data.
#' @param enforcedMinimumDistance Minimum distance (in kilometers) between a location considered treated and untreated, across the full dataset.  Defaults to 5km.
#' @param maximum_match_diff From 0 to 1, the largest difference in match quality allowed. 
#' @param cores The number of computer cores to use to fit the QGI model.  Defaults to number of cores on machine - 1.
#' @keywords
#' @export
#' @examples
#' QGI(read.csv("../df.csv"), independent_vars, x,x,x, 30*16)

QGI <- function(df, 
                distanceCol = "distance",
                controlVars, 
                outcomeVar = "outcome", 
                density = 50, 
                lowerDistBound=0.0, 
                upperDistBound="Default", 
                maximum_match_diff=0.9, 
                enforcedMinimumDistance=0.5,
                strategyNA = "Omit",
                pScoreMethod = "linear",
                pScoreDistance = "logit",
                matchDiscardStrategy = "both",
                pScoreReEstimate = TRUE,
                minN = 10,
                matchQualityWeighting = TRUE,
                figFile = "Default",
                cores = "Default") 
  {
  
  if(cores == "Default")
    {
      cores=parallel::detectCores()
    }

  if(strategyNA == "Omit")
    {
      df <- na.omit(df)
    }

  cl <- parallel::makeForkCluster(cores[1]-1) 
  doParallel::registerDoParallel(cl)

if(enforcedMinimumDistance <= 0.1)
{
  print("WARNING: Small enforced minimum distance values can result in duplicate cases in treatment and control sets, depending on the resolution of your data.  Use with extreme caution.")
}  


  if(upperDistBound == "Default")
    {
    upperDistBound = as.numeric(summary(df[,distanceCol])[3])
    }
  
  mc = 0
  #define sampling frame, with as much coverage as possible
  #given a defined density
  
  if(density <= 1)
    {
    density = 2
    print("Sample Density set to minimum allowed (2)")
    }

  #trials = seq(upperDistBound, lowerDistBound, (lowerDistBound-upperDistBound)/(density-1))
  trials = seq(lowerDistBound, upperDistBound, (upperDistBound-lowerDistBound)/(density-1))
  
  mc <- foreach(i = 1:length(trials), .combine=rbind) %dopar% {    
    
  #Rename the distance column to avoid later confusion with p-score distances
  df$geogDistance <- df[,distanceCol]
  
  #Set initial treatment value to 999 (representative of a null)
  df$Treatment <- 999
  
  #Identify the distance threshold for this iteration
  dist_thresh <- trials[i] 

  #Create a binary indicating if the location is considered "treated" this iteration
  #Or control.
  df$Treatment[df$geogDistance <= (0 + dist_thresh)] <- 1
  df$Treatment[df$geogDistance >= (enforcedMinimumDistance + dist_thresh)] <- 0
  
  #Count the number of treatment and control cases for output.
  df$Control = 0
  df$Control[df$Treatment == 0] <- 1

  treatCount = sum(df$Treatment)
  controlCount = sum(df$Control)



  #Remove any cases that were not assigned to treatment or control
  #Can happen if there are data errors (i.e., NAs in distance)
  df <- df[!df$Treatment == 999,]
  
  #Remove the distance column in prep for propensity matching
  df[, !(names(df) == "distance")]
  
  #Pscore variables
  pVars <- c(controlVars,"Treatment")

  #Analysis variables
  aVars <- c(controlVars,"Treatment","outcome")
  
  #construct propensity score formula
  f1 <- as.formula(paste("Treatment", paste(controlVars, collapse = " + "), sep = " ~ "))
  
  #Calculate our propensity scores using matchit
  pscoreCalc <- matchit(f1, 
                          data=df[pVars], 
                          method=pScoreMethod, 
                          distance=pScoreDistance,
                          discard=matchDiscardStrategy,
                          reestimate=pScoreReEstimate)
  
  matched.df <- match.data(pscoreCalc)
    
  if(nrow(matched.df) > minN)
  {
    #construct outcome modeling formula
    idx = as.numeric(rownames(matched.df))
    matched.df$outcome <- df[idx,][[outcomeVar]]
    f2 <- as.formula(paste("outcome", paste(c("Treatment", controlVars), collapse = " + "), sep = " ~ "))
    
    trtModel <- lm(f2, data=matched.df[aVars])
    
    
    match_diff = abs(summary(pscoreCalc)$sum.matched[[1]][1] - summary(pscoreCalc)$sum.matched[[2]][1])
    
    return( list(i, dist_thresh, trtModel$coef["Treatment"][[1]], nrow(df), match_diff, summary(trtModel)$r.squared, coef(summary(trtModel))[2,4], coef(summary(trtModel))[2,2], nrow(matched.df), treatCount, controlCount))
  } 

  }


  #Create visualizations and outputs
  tDF <- data.frame(t(matrix(unlist(mc), nrow=11, byrow=T)))
  colnames(tDF) <- c("ItId", "thresh", "coef", "obs", "match_diff", "R2", "TreatSig", "StdError", "ItSampleSize", "ItTreatmentCount", "ItControlCount")
  
  
  # #Scale for data visualization and (optionally) weighting
  tDF$matchWeight = 1 - ((tDF$match_diff - min(tDF$match_diff)) / (max(tDF$match_diff) - min(tDF$match_diff)))

  # #Third-order polynomial for the distance-decay
  tDF$b = tDF$thresh**2
  tDF$c = tDF$thresh**3
  tDF$d = tDF$thresh**4
  
  tDF <- tDF[order(tDF$thresh),]

  if(matchQualityWeighting == TRUE)
  {
    mean_mdl = lm(coef ~ thresh + b + c, data = tDF, weights =matchWeight)
    std_mdl = lm(StdError ~ thresh + b + c, data = tDF, weights =matchWeight)
  } else {
    mean_mdl = lm(coef ~ thresh + b + c, data = tDF)
    std_mdl = lm(StdError ~ thresh + b + c, data = tDF)
  }

  
  


  

  newx <- seq(min(tDF$thresh), max(tDF$thresh), length.out=1000)

  preds <- predict(mean_mdl, newdata=data.frame(thresh=newx, b=newx**2, c=newx**3), interval='confidence')
  preds_std <- predict(std_mdl, newdata=data.frame(thresh=newx, b=newx**2, c=newx**3), interval='confidence')

#================= Viz
if(figFile != "Default"){
png(figFile)
}
par(mfrow = c(1, 1),     # 2x2 layout
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(3, 3, 0, 0), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = FALSE)   
  ylim_upper <- max(tDF$coef + 1.96*tDF$StdError, na.rm = TRUE)
  ylim_lower <- min(tDF$coef - 1.96*tDF$StdError, na.rm = TRUE)

  #Calculate mean outcome for later use
  meanOutcome = mean(df[[outcomeVar]], na.rm=TRUE)

  plot(tDF$thresh, tDF$coef, cex = tDF$matchWeight, xlab="Distance (km)", ylab=paste(outcomeVar, " relative to mean of ", round(meanOutcome,3)), ylim=c(min(0,ylim_lower), max(ylim_upper,0)))
  
  polygon(c(rev(newx), newx), c(rev( preds[ ,3] + 1.96*preds_std[1:1000]),preds[ ,2] - 1.96*preds_std[1:1000]), col=rgb(0.2, 0.2, 0.25,0.25), border = NA)
  
  #polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col=rgb(0, 1, 0,0.25), border = NA)
  lines(newx, preds[ ,3] + 1.96*preds_std[1:1000], lty = 'dashed', col = 'black')
  lines(newx, preds[ ,2] - 1.96*preds_std[1:1000], lty = 'dashed', col = 'black')
  lines(newx, preds[ ,3], lty = 'dashed', col = 'yellow')
  lines(newx, preds[ ,2], lty = 'dashed', col = 'yellow')
  legend(5, 
         400, 
         legend=c("Best Fit", "Lower Match Quality", "Higher Match Quality"), 
         col=c("blue","black", "black"), 
         lty=c(1,NA, NA), 
         pch=c(NA,1,1), 
         cex=0.8, 
         pt.cex=c(NA, 0.5, 1))
  lines(tDF$thresh, fitted(mean_mdl), col="blue")
  abline(h = 0, lty = 2)

dev.off()
#================= End Viz


  retClass = list(figure = figFile, 
                  distanceModels = tDF)
  class(retClass) <- c("figure", 
                       "distanceModels")

  return(retClass)

}

  #plot(tDF$thresh, tDF$coef, cex = tDF$size, xlab="Distance (km)", ylab=outcome_label, ylim=c(-0.5,0.5))
  #return(c(tDF, min(newx[which(lower_std>0)]),max(newx[which(lower_std>0)])))
  
