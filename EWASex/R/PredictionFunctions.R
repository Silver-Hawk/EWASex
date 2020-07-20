#' Precomputed means and SD from 49 CpGs.
#'
#' Means and SD from 49 CpGs, computed from MADT, BWT, and LSADT2 cohorts. See @source.
#'
#' @docType data
#'
#' @usage data(MeansAndSD49)
#'
#' @format An object of class \code{"data.frame"}; The keys represents genders, here 1 represent a "male" and 2 represent a "female".
#'
#' @keywords datasets
#'
#' @references Not published yet.
#'
#' @source Combined means and SDs computed from Middle Aged Danish Twins (MADT),
#' Longitudinal Study of Aging Danish Twins 2 (LSADT2), and, Birth-weight Discordant Twins (BWT).
#' (See https://doi.org/10.1017/thg.2012.77 and https://doi.org/10.1017/thg.2012.77 for dataset descriptions)
#'
#' @examples
#' data(MeansAndSD49)
#'
"MeansAndSD49"

#' Get the 49 goldset CpGs.
#'
#' Allows for easy subsetting of a data.frame.
#'
#' @return A character vector with the names of the 49 CpGs that best predict sex.
#' @examples
#' CpGs <- getGoldCpGNames()
#'
#' @export
getGoldCpGNames <- function() {
  c("cg05257947", "cg25940844", "cg06615444", "cg05958126", "cg22969661", "cg00114625", "cg00666173", "cg14295915",
    "cg00026186", "cg24741068", "cg13574945", "cg15132216", "cg10723556", "cg23696472", "cg22604777", "cg25225807",
    "cg15977272", "cg16221895", "cg24186901", "cg21201934", "cg04317640", "cg13203541", "cg00832270", "cg24790801",
    "cg27501723", "cg20208613", "cg18689730", "cg15565409", "cg17363084", "cg10981178", "cg26505478", "cg18102950",
    "cg14098973", "cg11673471", "cg07861180", "cg03773146", "cg13244998", "cg01120894", "cg10717149", "cg25206026",
    "cg01742836", "cg20662859", "cg06136002", "cg00963467", "cg03670113", "cg13766601", "cg05872808", "cg06042004",
    "cg00813156")
}

#' Computes the vectors for the means and SDs that are used for the prediction part.
#'
#' The vectors can be stored to predict at a later time or used on another dataset that doesn't contain gender information.
#' If you are only intrested in predicting gender, you can use \code{getPredictions} with \code{means} = "none".
#'
#' @param genders Gender information indicating the genders. This can be a vector consiting of 1's and 2's,
#' or characters, like, "male" and "female".
#' @param df The data.frame consiting of the CpG information. This can be any CpGs, but it is recommended to use the gold set of 49 CpGs.
#' @param margin Whether to compute the means and SDs column or row-wise. Default is 2, corresponding to one sample per row and a CpG for each column.
#' @return The means and SDs of \code{df} for each gender in \code{genders}
#' @examples
#' genders <- pheno$genders # get genders from your phenotype file
#' df <- wholeDataFrame[,getGoldCpGNames()] # subset to only use the 49 CpGs
#'
#' MeansAndSD <- getMeansAndSD(genders, df, margin=2)
#' @export
getMeansAndSD <- function(genders, df, margin=2) {
  keys_ = unique(genders)

  if(margin==1)
    df = t(df)

  ret_list = list()
  for (k_ in keys_) {
    ret_list[[k_]] = list()
    ret_list[[k_]][['mean']] = apply(df[genders == k_,,drop=F], 2, mean, na.rm=T)
    ret_list[[k_]][['sd']] = apply(df[genders == k_,,drop=F], 2, sd, na.rm=T)
  }
  names(ret_list) <- keys_

  return(ret_list)
}



#' Computes the errors, normalized errors and performs the sex predictions.
#'
#' The vectors can be stored to predict at a later time or used on another dataset that doesn't contain gender information.
#' If you are only intrested in predicting gender, you can use \code{getPredictions} with \code{means} = "none".
#'
#' @param genders Gender information indicating the genders. This can be a vector consiting of 1's and 2's,
#' or characters, like, "male" and "female".
#' @param df The data.frame consiting of the CpG information. This can be any CpGs, but it is recommended to use the gold set of 49 CpGs.
#' @param margin Whether to compute the means and SDs column or row-wise. Default is 2, corresponding to one sample per row and a CpG for each column.
#' @param means Means and SDs already computed using the \code{getMeansAndSD} function.
#' @return The means and SDs of \code{df} for each gender in \code{genders}
#' @examples
#' # Example 1:
#' # if you want to predict using the gender information of the dataset
#' genders <- pheno$genders # get genders from your phenotype file
#' df <- wholeDataFrame[,getGoldCpGNames()] # subset to only use the 49 CpGs
#'
#' predictions <- getPredictions(genders, df, margin=2, means="none")
#'
#' # Example 2:
#' # if you already have means and SD from the current of another dataset
#' df <- wholeDataFrame[,getGoldCpGNames()] # subset to only use the 49 CpGs
#'
#' predictions <- getPredictions(genders="none", df, margin=2, means=MeansAndSD)
#' @export
getPredictions <- function(genders, df, margin=2, means="none") {
  # convert to usuable data.frame
  if (margin == 1) {
    df = t(df)
  }

  # if means are not given, calculate them for the current set and genders
  if (is.character(means)) {
    means = getMeansAndSD(genders, df, margin=2)
  }

  # so far, the script only supports two genders
  keys_ = names(means)

  E1 <- data.frame(gender=keys_[1], cpg=colnames(df), mean=means[[keys_[1]]][['mean']], sd=means[[keys_[1]]][['sd']])
  E2 <- data.frame(gender=keys_[2], cpg=colnames(df), mean=means[[keys_[2]]][['mean']], sd=means[[keys_[2]]][['sd']])

  E1testVals <- abs(df - E1$mean) / E1$sd
  E1test <- rowSums(E1testVals) / ncol(df)
  E2testVals <- abs(df - E2$mean) / E2$sd
  E2test <- rowSums(E2testVals) / ncol(df)

  returnFrame <- data.frame(Error1=E1test,
                            Error2=E2test,
                            NormError1=E1test/max(E1test),
                            NormError2=E2test/max(E2test))

  returnFrame$predictedGender <- ifelse(returnFrame$NormError1 < returnFrame$NormError2, keys_[1], keys_[2])

  return (returnFrame)
}
