# Copyright (C) 2025 Yuhan Tian
#
# This program is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software Foundation; either version
# 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Yuhan Tian, <tyh9293@gmail.com>


#' Generate VAR Data
#'
#' This function generates Vector Auto-Regressive (VAR) data based on the provided parameters.
#'
#' @param n Integer. The length of the VAR data to be generated.
#' @param As List. A list containing the transition matrices for the VAR process.
#' @param Sig Matrix. The covariance matrix of errors.
#' @param h Integer. The order of the VAR process.
#' @param isOldProvided Logical. If TRUE, the VAR data will be generated based on the
#'        last observations from the previous segment of data. Defaults to FALSE.
#' @param oldxs Matrix. A \code{p} by \code{h} matrix containing the last observations from
#'        the previous segment of data. Required if \code{isOldProvided = TRUE}.
#'
#' @return A data matrix of dimensions \code{p} by \code{n}.
#' @importFrom MASS mvrnorm
#' @examples
#' # Example usage
#' As <- list(matrix(c(0.5, 0.2, 0.1, 0.4), 2, 2))
#' Sig <- diag(2)
#' data <- generateVAR(n = 100, As = As, Sig = Sig, h = 1)
#' @export
generateVAR <- function(n, As, Sig, h, isOldProvided = FALSE, oldxs = NULL) {
  p <- ncol(As[[1]])

  if (!isOldProvided) {
    oldxs <- matrix(0, p, h)
  }

  data <- cbind(oldxs, matrix(0, p, n + 100))

  for (t in (1 + h):(n + h + 100)) {
    temp <- 0
    for (i in 1:h) {
      temp <- temp + As[[i]] %*% data[, t - i]
    }
    data[, t] <- temp + mvrnorm(n = 1, mu = rep(0, p), Sigma = Sig)
  }

  if (isOldProvided) {
    return(data[, (1 + h):(n + h)])
  } else {
    return(data[, (1 + h + 100):(n + h + 100)])
  }
}


#' VAR_cpDetect_Online: Sequential change point Detection for Vector Auto-Regressive Models
#'
#' This function performs sequential change point detection in high-dimensional time series data modeled as a Vector Auto-Regressive (VAR) process, targeting changes in the transition matrices that encode temporal and cross-correlations.
#'
#' @param data A matrix where rows represent different dimensions (features) and columns represent observations. The first \code{n0} columns are treated as historical data.
#' @param n0 Integer. The size of the historical data (number of columns in \code{data} treated as historical).
#' @param w Integer. The size of the sliding window used for calculating test statistics; referred to as the pre-specified detection delay.
#' @param alpha Numeric. The desired false alarm rate, where 1/alpha represents the targeted average run length (ARL), which should exceed the length of the data to be monitored.
#' @param h Integer. The order of the VAR process.
#' @param RLmode Logical. If \code{TRUE}, the algorithm terminates when the first alarm is issued.
#' @param needRefine Logical. If \code{TRUE}, a refinement process is conducted to pinpoint the change point location.
#' @param refineSize Numeric. The proportion of the new window size to the original window size, used during refinement.
#' @param needRefineCorrection Logical. If \code{TRUE}, a confirmation step is performed during the refinement process to verify the detected change point.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{RL}}{The index (ignoring historical data) of the last observation read by the algorithm when the first alarm is issued. This is returned only if \code{RLmode = TRUE}.}
#'   \item{\code{cp_refined}}{The refined estimate for the location (ignoring historical data) of the change point. This is returned only if \code{RLmode = TRUE} and \code{needRefine = TRUE}.}
#'   \item{\code{alarm_locations}}{A vector of indices (ignoring historical data) where alarms were raised. This is returned only if \code{RLmode = FALSE}.}
#'   \item{\code{cp_locations}}{A vector of refined change point locations (ignoring historical data), corresponding 1-to-1 with the \code{alarm_locations}. This is returned only if \code{RLmode = FALSE} and \code{needRefine = TRUE}.}
#' }
#'
#' @details This function fits a VAR model to the historical data using the l1 penalty and calculates test statistics for the sliding window to detect change points. If refinement is enabled, a second step narrows down the change point location. Optionally, a correction step can verify the detected change points.
#' @importFrom stats qnorm
#' @examples
#' library(MASS)
#' set.seed(2024)
#' As <- list(matrix(c(0.5, 0.2, 0.1, 0.4), 2, 2))
#' As_new <- list(matrix(c(-0.5, 0.2, 0.1, -0.4), 2, 2))
#' Sig <- diag(2)
#' data_IC <- generateVAR(n = 400, As = As, Sig = Sig, h = 1)
#' data_OC <- generateVAR(n = 100, As = As_new, Sig = Sig, h = 1,
#'                        isOldProvided = TRUE, oldxs = data_IC[, ncol(data_IC)])
#' data <- cbind(data_IC, data_OC)
#' result <- VAR_cpDetect_Online(data, n0 = 300, w = 20, alpha = 1/200, h = 1,
#'                               RLmode = TRUE, needRefine = TRUE, refineSize = 1/5,
#'                               needRefineCorrection = TRUE)
#' print(result)
#'
#' @export
VAR_cpDetect_Online <- function(data, n0, w, alpha, h, RLmode = TRUE, needRefine = TRUE, refineSize = 1/5, needRefineCorrection = TRUE){
  lastindex <- 0
  alarm_locations <- NULL
  cp_locations <- NULL
  p <- nrow(data)

  data_estimation <- data[, (lastindex + 1):(lastindex + n0)]
  fit <- fitVAR(t(data_estimation), p = h, penalty = "ENET")
  SSE <- 0
  SSE_2 <- 0
  for (t in (h + 1):n0) {
    temp_sum <- 0
    for (l in 1:h) {
      temp_sum <- temp_sum + fit$A[[l]]%*%data[, (t - l)]
    }
    SSE <- SSE + sum((temp_sum - data[, t])^2)
    SSE_2 <- SSE_2 + sum((temp_sum - data[, t])^4)
  }
  sigma2hat <- SSE/((n0 - h)*p)
  varhat <- SSE_2/((n0 - h)*p) - sigma2hat^2
  lastindex <- lastindex + n0

  data_in_use <- data[, (lastindex + 1 - h):(lastindex + w)]
  ts_ori <- 0
  head <- 0
  for (t in (1 + h):(w + h)) {
    temp_sum1 <- 0
    for (l in 1:h) {
      temp_sum1 <- temp_sum1 + fit$A[[l]]%*%data_in_use[, (t - l)]
    }
    temp_sum2 <- sum((data_in_use[, t] - temp_sum1)^2)
    ts_ori <- ts_ori + temp_sum2
    if (t == (1 + h)) {
      head <- temp_sum2
    }
  }
  ts_nor <- sqrt((p*w)/varhat)*(ts_ori/(p*w) - sigma2hat)
  if (abs(ts_nor) > qnorm(1 - (alpha/2))) {
    alarm_locations <- append(alarm_locations, (lastindex - n0))
    if (needRefine) {
      w_new <- ceiling(w*refineSize)
      i <- h
      while (i <= (w + h - w_new)) {
        ts_oriprime <- 0
        data_in_use_in_use <- data_in_use[, (i + 1 - h):(i + w_new)]
        for (tprime in (h + 1):(h + w_new)) {
          tempsum <- 0
          for (l in 1:h) {
            tempsum <- tempsum + fit$A[[l]]%*%data_in_use_in_use[, (tprime - l)]
          }
          ts_oriprime <- ts_oriprime + sum((tempsum - data_in_use_in_use[, tprime])^2)
        }
        ts_norprime <- sqrt((p*w_new)/varhat)*(ts_oriprime/(p*w_new) - sigma2hat)
        if (abs(ts_norprime) > qnorm(1 - (alpha/2))) {
          if (RLmode) {
            return(list(RL = (lastindex + w - n0), cp_refined = (lastindex + (i - h) - n0)))
          }
          else {
            cp_locations <- append(cp_locations, (lastindex + (i - h) - n0))
            break
          }
        }
        i <- i + 1
      }
      if(i == (w + h - w_new + 1)){
        if (needRefineCorrection) {
          alarm_locations <- alarm_locations[-length(alarm_locations)]
        }
        else {
          if (RLmode) {
            return(list(RL = (lastindex + w - n0), cp_refined = (lastindex + w - n0)))
          }
          else {
            cp_locations <- append(cp_locations, (lastindex + w - n0))
          }
        }
      }
    }
    if(RLmode&(!needRefine)){
      return(list(RL = (lastindex + w - n0)))
    }
  }
  lastindex <- lastindex + 1

  while(lastindex <= (ncol(data) - w)){
    data_in_use <- data[, (lastindex + 1 - h):(lastindex + w)]
    ts_ori <- ts_ori - head
    temp_sum1 <- 0
    for (l in 1:h) {
      temp_sum1 <- temp_sum1 + fit$A[[l]]%*%data_in_use[, (h + 1 - l)]
    }
    head <- sum((data_in_use[, (h + 1)] - temp_sum1)^2)
    temp_sum1 <- 0
    for (l in 1:h) {
      temp_sum1 <- temp_sum1 + fit$A[[l]]%*%data_in_use[, (w + h - l)]
    }
    ts_ori <- ts_ori + sum((data_in_use[, (w + h)] - temp_sum1)^2)
    ts_nor <- sqrt((p*w)/varhat)*(ts_ori/(p*w) - sigma2hat)
    if (abs(ts_nor) > qnorm(1 - (alpha/2))) {
      alarm_locations <- append(alarm_locations, (lastindex - n0))
      if (needRefine) {
        w_new <- ceiling(w*refineSize)
        i <- h
        while (i <= (w + h - w_new)) {
          ts_oriprime <- 0
          data_in_use_in_use <- data_in_use[, (i + 1 - h):(i + w_new)]
          for (tprime in (h + 1):(h + w_new)) {
            tempsum <- 0
            for (l in 1:h) {
              tempsum <- tempsum + fit$A[[l]]%*%data_in_use_in_use[, (tprime - l)]
            }
            ts_oriprime <- ts_oriprime + sum((tempsum - data_in_use_in_use[, tprime])^2)
          }
          ts_norprime <- sqrt((p*w_new)/varhat)*(ts_oriprime/(p*w_new) - sigma2hat)
          if (abs(ts_norprime) > qnorm(1 - (alpha/2))) {
            if (RLmode) {
              return(list(RL = (lastindex + w - n0), cp_refined = (lastindex + (i - h) - n0)))
            }
            else {
              cp_locations <- append(cp_locations, (lastindex + (i - h) - n0))
              break
            }
          }
          i <- i + 1
        }
        if(i == (w + h - w_new + 1)){
          if (needRefineCorrection) {
            alarm_locations <- alarm_locations[-length(alarm_locations)]
          }
          else {
            if (RLmode) {
              return(list(RL = (lastindex + w - n0), cp_refined = (lastindex + w - n0)))
            }
            else {
              cp_locations <- append(cp_locations, (lastindex + w - n0))
            }
          }
        }
      }
      if(RLmode&(!needRefine)){
        return(list(RL = (lastindex + w - n0)))
      }

    }
    lastindex <- lastindex + 1
  }
  if (RLmode) {
    return(list(RL = (ncol(data) - n0), cp_refined = (ncol(data) - n0)))
  }
  return(list(alarm_locations = alarm_locations, cp_locations = cp_locations))
}

#' Identify the Beginning of the Alarm Clusters
#'
#' This function clusters alarms into groups and identifies the starting points of the alarm clusters.
#' If the next alarm occurs within a specified window size (\code{w}) from the current alarm,
#' it will be considered part of the current cluster. Otherwise, a new cluster will be formed.
#'
#' @param alarms A numeric vector. The alarms raised during the monitoring process.
#' @param w An integer. The window size used to group alarms into clusters.
#'
#' @return A numeric vector containing the starting points of the alarm clusters.
#'         If the next alarm is within \code{w} observations of the current alarm,
#'         the next alarm will be considered part of the current alarm cluster.
#'         Otherwise, a new cluster is formed and the next alarm is considered the beginning
#'         of a new alarm cluster.
#'
#' @examples
#' # Example usage:
#' alarms <- c(10, 15, 30, 35, 60)
#' change_points <- get_cps(alarms, w = 10)
#'
#' @export
get_cps <- function(alarms, w) {
  ans <- NULL
  prev <- -w - 1
  for (num in alarms) {
    if (num > (prev + w)) {
      ans <- append(ans, num)
    }
    prev <- num
  }
  return(ans)
}



#' S&P 500 Daily Log Returns and Corresponding Dates
#'
#' This dataset contains daily log returns for 186 stocks in the S&P 500
#' index from February 6, 2004, to March 2, 2016. The daily log returns
#' are calculated using the adjusted daily closing prices. The dataset also
#' contains the corresponding dates for each log return.
#'
#' The dataset is provided as an `.RData` file containing:
#' \itemize{
#'   \item \code{sp500$log_daily_return}: A matrix of daily log returns with 3037 rows (trading days) and 186 columns (stocks).
#'   \item \code{sp500$date}: A vector of length 3037 containing the dates for each daily log return.
#' }
#'
#' @name sp500
#' @docType data
#' @usage data(sp500)
#' @format A list with two elements:
#' \describe{
#'   \item{\code{sp500$log_daily_return}}{A matrix with dimensions 3037 (rows, trading days) by 186 (columns, stocks).}
#'   \item{\code{sp500$date}}{A vector of length 3037, containing the dates for each trading day.}
#' }
#' @source Data from the S&P 500 stock index (2004-2016).
#' @examples
#' # Example Usage: Applying Change Point Detection to S&P 500 Data
#' # This is an example of how to apply the change point detection method
#' # (using the VAR_cpDetect_Online function) on the daily log return
#' # dataset from the S&P 500 (stored in the sp500 dataset). The code
#' # below calculates the average return volatility for all stocks, applies
#' # the change point detection algorithm, and plots the results with detected
#' # change points shown as vertical red and black lines.
#' \donttest{
#' # Load the dataset
#' data(sp500)
#'
#' # Set parameters
#' library(ggplot2)
#' set.seed(2024)
#' n_sp <- nrow(sp500$log_daily_return)
#' p_sp <- ncol(sp500$log_daily_return)
#'
#' # Calculate average return volatility for all data points
#' volatility_sum <- rep(0, (n_sp - 21))
#' for(col in 1:p_sp){
#'   temp <- as.numeric(sp500$log_daily_return[, col])
#'   temp1 <- rep(0, (n_sp - 21))
#'   for(row in 1:(n_sp - 21)){
#'     temp1[row] <- sd(temp[(row):(row + 21)])
#'   }
#'   volatility_sum <- volatility_sum + temp1
#' }
#' volatility_ave <- volatility_sum / p_sp
#'
#' # Apply change point detection method
#' n0 <- 200
#' w <- 22
#' alpha <- 1 / 5000
#'
#' res <- VAR_cpDetect_Online(t(sp500$log_daily_return), n0, w, alpha, 1, FALSE, TRUE, 5 / w, TRUE)
#' res_sp <- res$alarm_locations + n0
#' res_sp_cps <- res$cp_locations + n0
#' # Get the estimated starting points of each alarm cluster
#' cps_est_sp <- unique(res_sp_cps[which(res_sp %in% get_cps(res_sp, w))])
#'
#' # Prepare data for plotting
#' y_values <- c(volatility_ave)
#' x_values <- sp500$date[1:(n_sp - 21)]
#' df <- data.frame(y_values, x_values)
#' plot_sp <- ggplot(df, aes(y = y_values, x = x_values)) +
#'   geom_line() +
#'   theme(legend.position = "none") +
#'   labs(title = "", x = "", y = "") +
#'   scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
#'   geom_vline(xintercept = sp500$date[res_sp], linetype = "solid", color = "red", alpha = .1) +
#'   geom_vline(xintercept = sp500$date[cps_est_sp], linetype = "solid", color = "black")
#'
#' # Print the detected change points
#' sp500$date[cps_est_sp] # The dates for the starting of the alarm clusters
#' plot_sp
#' }
#' @seealso \code{\link{VAR_cpDetect_Online}}, \code{\link{get_cps}}
NULL


#' Fit VAR Model with Elastic Net via Cross Validation
#'
#' Estimates a (possibly high-dimensional) VAR model using penalized least squares
#' with an elastic net penalty and cross validation.
#' This function is adapted from the \emph{sparsevar} package
#' (<https://github.com/svazzole/sparsevar/tree/master>), which is distributed under
#' the GNU General Public License v2. The code has been modified to specifically implement
#' the elastic net penalty (penalty = "ENET") and cross validation (method = "cv").
#'
#' @param data A numeric matrix or data frame with time series data (observations in rows,
#'   variables in columns).
#' @param p Integer. The order of the VAR model.
#' @param ... Additional options for estimation. Global options include:
#'   \itemize{
#'     \item \code{threshold}: Logical. If \code{TRUE}, all entries smaller than the oracle
#'           threshold are set to zero.
#'     \item \code{scale}: Logical. Whether to scale the data (default is \code{FALSE}).
#'     \item \code{nfolds}: Integer. The number of folds used for cross validation (default is 10).
#'     \item \code{parallel}: Logical. If \code{TRUE}, use multicore backend (default is \code{FALSE}).
#'     \item \code{ncores}: Integer. If \code{parallel = TRUE}, specify the number of cores to use.
#'     \item \code{alpha}: Numeric. The elastic net mixing parameter (default is 1, i.e. LASSO).
#'     \item \code{type.measure}: Character. The error measure for CV (e.g., \code{"mse"} or \code{"mae"}).
#'     \item \code{nlambda}: Integer. The number of lambda values to use in cross validation (default is 100).
#'     \item \code{leaveOut}: Integer. In time slice validation, leave out the last observations (default is 15).
#'     \item \code{horizon}: Integer. The forecast horizon to use for estimating error (default is 1).
#'     \item \code{lambda}: Either a numeric vector of lambda values or a string indicating which
#'           lambda to use (default is \code{"lambda.min"}).
#'     \item \code{return_fit}: Logical. If \code{TRUE}, return the complete fit object.
#'   }
#'
#' @return A list with the following components:
#'   \item{mu}{A vector of means for each variable.}
#'   \item{A}{A list (of length \code{p}) of the estimated coefficient matrices for the VAR process.}
#'   \item{fit}{(Optional) The complete results of the penalized least squares estimation.}
#'   \item{lambda}{The chosen lambda value (by cross validation).}
#'   \item{mse}{The minimum mean squared error from cross validation.}
#'   \item{mse_sd}{The standard deviation of the mean squared error.}
#'   \item{time}{Elapsed time for the estimation.}
#'   \item{series}{The (possibly transformed) input time series.}
#'   \item{residuals}{The residuals of the VAR model.}
#'   \item{sigma}{The estimated variance/covariance matrix of the residuals.}
#'
#' @references The original source code is adapted from the
#'   \href{https://github.com/svazzole/sparsevar/tree/master}{sparsevar package},
#'   which is distributed under the GNU General Public License v2.
#'
#' @export
fitVAR <- function(data, p = 1, ...) {
  opt <- list(...)

  # convert data to matrix if necessary
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }

  cnames <- colnames(data)

  # Use cross validation to find lambda and fit the model
  out <- cvVAR(data, p, opt)

  # Add variable names to the estimated matrices if available
  if (!is.null(cnames)) {
    for (k in 1:length(out$A)) {
      colnames(out$A[[k]]) <- cnames
      rownames(out$A[[k]]) <- cnames
    }
  }

  return(out)
}


#' Cross-Validated VAR Estimation using Elastic Net
#'
#' This internal function performs cross validation for VAR estimation using the elastic net
#' penalty. It prepares the data, calls the elastic net CV routine, reshapes the estimated coefficients,
#' applies optional thresholding, computes residuals, and estimates the error covariance.
#'
#' @param data A numeric matrix with time series data (observations in rows, variables in columns).
#' @param p Integer. The order of the VAR model.
#' @param opt List. A list of options (see \code{fitVAR} for details).
#'
#' @return A list with components:
#'   \item{mu}{Vector of means of the original series.}
#'   \item{A}{List of VAR coefficient matrices (one for each lag).}
#'   \item{fit}{The complete elastic net CV fit (if requested).}
#'   \item{lambda}{The optimal lambda value chosen by CV.}
#'   \item{mse}{The minimum mean squared error from CV.}
#'   \item{mse_sd}{Standard deviation of the MSE.}
#'   \item{time}{Elapsed time for the ENET estimation.}
#'   \item{series}{The transformed series (after centering/scaling).}
#'   \item{residuals}{Residuals from the VAR model.}
#'   \item{sigma}{Estimated covariance matrix of the residuals.}
#'
#' @keywords internal
cvVAR <- function(data, p, opt = NULL) {
  nc <- ncol(data)
  nr <- nrow(data)

  threshold <- ifelse(!is.null(opt$threshold), opt$threshold, FALSE)
  threshold_type <- ifelse(!is.null(opt$threshold_type), opt$threshold_type, "soft")
  return_fit <- ifelse(!is.null(opt$return_fit), opt$return_fit, FALSE)

  # Transform the dataset into design matrices for VAR estimation
  tr_dt <- transformData(data, p, opt)

  # Fit the elastic net model via cross validation
  t <- Sys.time()
  fit <- cvVAR_ENET(tr_dt$X, tr_dt$y, nvar = nc, opt)
  elapsed <- Sys.time() - t

  # Extract the lambda option
  lambda <- ifelse(is.null(opt$lambda), "lambda.min", opt$lambda)

  # Extract coefficients (ignoring the intercept) and reshape into a matrix
  Avector <- stats::coef(fit, s = lambda)
  A <- matrix(Avector[2:length(Avector)],
              nrow = nc, ncol = nc * p,
              byrow = TRUE)

  mse <- min(fit$cvm)

  # Apply thresholding if requested
  if (threshold == TRUE) {
    A <- applyThreshold(A, nr, nc, p, type = threshold_type)
  }

  # Split the coefficient matrix into a list (one matrix per lag)
  A <- splitMatrix(A, p)

  # Compute residuals from the VAR model
  res <- computeResiduals(tr_dt$series, A)

  # Extract the standard deviation of the CV error
  ix <- which(fit$cvm == min(fit$cvm))
  mse_sd <- fit$cvsd[ix]

  # Create and return the output list
  output <- list()
  output$mu <- tr_dt$mu
  output$A <- A

  if (return_fit == TRUE) {
    output$fit <- fit
  }

  output$lambda <- fit$lambda.min
  output$mse <- mse
  output$mse_sd <- mse_sd
  output$time <- elapsed
  output$series <- tr_dt$series
  output$residuals <- res

  # Estimate the error covariance matrix
  output$sigma <- estimateCovariance(res)

  attr(output, "class") <- "var"
  attr(output, "type") <- "fit"
  return(output)
}


#' Cross Validation for Elastic Net VAR Estimation
#'
#' This internal function performs cross validation using elastic net (ENET)
#' estimation via the \code{glmnet} package. It supports parallel processing if requested.
#'
#' @param X A numeric matrix of predictors.
#' @param y Numeric vector of responses.
#' @param nvar Integer. The number of variables in the original VAR (number of columns in data).
#' @param opt List. A list of options including:
#'   \itemize{
#'     \item \code{alpha}: The elastic net mixing parameter (default = 1).
#'     \item \code{nlambda}: Number of lambda values (default = 100).
#'     \item \code{type.measure}: Error measure for CV (default = "mse").
#'     \item \code{nfolds}: Number of folds for CV (default = 10).
#'     \item \code{parallel}: Logical. Whether to use parallel processing (default = FALSE).
#'     \item \code{ncores}: Number of cores for parallel processing (default = 1).
#'     \item \code{lambdas_list}: Optionally, a user-specified list of lambdas.
#'     \item \code{folds_ids}: Optionally, user-specified fold IDs for CV.
#'   }
#'
#' @return An object of class \code{cv.glmnet} as returned by \code{glmnet::cv.glmnet}.
#'
#' @keywords internal
cvVAR_ENET <- function(X, y, nvar, opt) {
  a <- ifelse(is.null(opt$alpha), 1, opt$alpha)
  nl <- ifelse(is.null(opt$nlambda), 100, opt$nlambda)
  tm <- ifelse(is.null(opt$type.measure), "mse", opt$type.measure)
  nf <- ifelse(is.null(opt$nfolds), 10, opt$nfolds)
  parall <- ifelse(is.null(opt$parallel), FALSE, opt$parallel)
  ncores <- ifelse(is.null(opt$ncores), 1, opt$ncores)

  # Define lambda values to use
  if (!is.null(opt$lambdas_list)) {
    lambdas_list <- opt$lambdas_list
  } else {
    lambdas_list <- c(0)
  }

  # Assign fold IDs if provided
  if (is.null(opt$folds_ids)) {
    folds_ids <- numeric(0)
  } else {
    nr <- nrow(X)
    folds_ids <- rep(sort(rep(seq(nf), length.out = nr / nvar)), nvar)
  }

  # Call cv.glmnet with or without parallel processing
  if (parall == TRUE) {
    if (ncores < 1) {
      stop("The number of cores must be > 1")
    } else {
      cl <- doParallel::registerDoParallel(cores = ncores)
      if (length(folds_ids) == 0) {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, nlambda = nl,
                                     type.measure = tm, nfolds = nf,
                                     parallel = TRUE, standardize = FALSE)
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, lambda = lambdas_list,
                                     type.measure = tm, nfolds = nf,
                                     parallel = TRUE, standardize = FALSE)
        }
      } else {
        if (length(lambdas_list) < 2) {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, nlambda = nl,
                                     type.measure = tm, foldid = folds_ids,
                                     parallel = TRUE, standardize = FALSE)
        } else {
          cvfit <- glmnet::cv.glmnet(X, y,
                                     alpha = a, lambda = lambdas_list,
                                     type.measure = tm, foldid = folds_ids,
                                     parallel = TRUE, standardize = FALSE)
        }
      }
    }
  } else {
    if (length(folds_ids) == 0) {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, nlambda = nl,
                                   type.measure = tm, nfolds = nf,
                                   parallel = FALSE, standardize = FALSE)
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, lambda = lambdas_list,
                                   type.measure = tm, nfolds = nf,
                                   parallel = FALSE, standardize = FALSE)
      }
    } else {
      if (length(lambdas_list) < 2) {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, nlambda = nl,
                                   type.measure = tm, foldid = folds_ids,
                                   parallel = FALSE, standardize = FALSE)
      } else {
        cvfit <- glmnet::cv.glmnet(X, y,
                                   alpha = a, lambda = lambdas_list,
                                   type.measure = tm, foldid = folds_ids,
                                   parallel = FALSE, standardize = FALSE)
      }
    }
  }

  return(cvfit)
}


#' Transform Data for VAR Estimation
#'
#' Transforms the input time series data into the design matrices required for VAR estimation.
#' This includes centering, optional scaling, and constructing the lagged predictor matrix.
#'
#' @param data A numeric matrix or data frame with time series data (observations in rows,
#'   variables in columns).
#' @param p Integer. The order of the VAR model (number of lags).
#' @param opt List. Options for data transformation. Supported options include:
#'   \itemize{
#'     \item \code{scale}: Logical. Whether to scale the data columns (default is \code{FALSE}).
#'     \item \code{center}: Logical. Whether to center the data columns (default is \code{TRUE}).
#'   }
#'
#' @return A list with the following components:
#'   \item{X}{The design matrix (via the Kronecker product) for lagged predictors.}
#'   \item{y}{A vectorized response corresponding to the lagged data.}
#'   \item{series}{The (centered and possibly scaled) original time series matrix.}
#'   \item{mu}{A row vector of the column means used for centering.}
#'
#' @keywords internal
transformData <- function(data, p, opt) {
  nr <- nrow(data)
  nc <- ncol(data)
  data <- as.matrix(data)

  # Determine whether to scale and/or center the data
  scale_flag <- ifelse(is.null(opt$scale), FALSE, opt$scale)
  center_flag <- ifelse(is.null(opt$center), TRUE, opt$center)

  if (center_flag == TRUE) {
    m <- colMeans(data)
    cm <- matrix(rep(m, nrow(data)), nrow = nrow(data), byrow = TRUE)
    data <- data - cm
  } else {
    m <- rep(0, nc)
  }

  if (scale_flag == TRUE) {
    data <- apply(X = data, MARGIN = 2, FUN = scale)
  }

  # Construct lagged matrices
  tmpX <- data[1:(nr - 1), ]
  tmpY <- data[2:nr, ]
  tmpX <- duplicateMatrix(tmpX, p)
  tmpY <- tmpY[p:nrow(tmpY), ]
  y <- as.vector(tmpY)

  # Create the design matrix using the Kronecker product
  I <- Matrix::Diagonal(nc)
  X <- kronecker(I, tmpX)

  output <- list()
  output$X <- X
  output$y <- y
  output$series <- data
  output$mu <- t(m)
  return(output)
}


#' Apply Thresholding to VAR Coefficients
#'
#' Applies a thresholding rule to a coefficient matrix by setting entries below a
#' certain threshold to zero. Two types of thresholding are available: "soft" and "hard".
#'
#' @param a_mat Numeric matrix. The coefficient matrix to be thresholded.
#' @param nr Integer. The number of rows in the original data.
#' @param nc Integer. The number of variables (columns) in the original data.
#' @param p Integer. The order of the VAR model.
#' @param type Character. The type of threshold to apply; either \code{"soft"} (default)
#'   or \code{"hard"}.
#'
#' @return The thresholded coefficient matrix.
#'
#' @keywords internal
applyThreshold <- function(a_mat, nr, nc, p, type = "soft") {
  if (type == "soft") {
    tr <- 1 / sqrt(p * nc * log(nr))
  } else if (type == "hard") {
    tr <- (nc) ^ (-0.49)
  } else {
    stop("Unknown threshold type. Possible values are: \"soft\" or \"hard\"")
  }
  l_mat <- abs(a_mat) >= tr
  a_mat <- a_mat * l_mat
  return(a_mat)
}


#' Compute VAR Model Residuals
#'
#' Computes the residuals from a VAR model by subtracting the fitted values (obtained
#' from the estimated coefficient matrices) from the original time series data.
#'
#' @param data A numeric matrix of the original time series (observations in rows).
#' @param A List. A list of VAR coefficient matrices (one for each lag).
#'
#' @return A numeric matrix of residuals.
#'
#' @keywords internal
computeResiduals <- function(data, A) {
  nr <- nrow(data)
  nc <- ncol(data)
  p <- length(A)
  res <- matrix(0, nrow = nr, ncol = nc)
  f <- matrix(0, nrow = nr, ncol = nc)

  for (i in 1:p) {
    tmpD <- rbind(matrix(0, nrow = i, ncol = nc), data[1:(nrow(data) - i), ])
    tmpF <- t(A[[i]] %*% t(tmpD))
    f <- f + tmpF
  }
  res <- data - f
  return(res)
}


#' Estimate Covariance Matrix from Residuals
#'
#' Estimates the covariance (or variance) matrix of the residuals using shrinkage estimation.
#' This function utilizes \code{corpcor::cov.shrink} for covariance estimation.
#'
#' @param res A numeric matrix of residuals from the VAR model.
#' @param ... Additional arguments passed to \code{corpcor::cov.shrink} (if any).
#'
#' @return A numeric covariance matrix.
#'
#' @importFrom corpcor cov.shrink
#' @keywords internal
estimateCovariance <- function(res, ...) {
  nc <- ncol(res)
  s <- corpcor::cov.shrink(res, verbose = FALSE)
  sigma <- matrix(0, nrow = nc, ncol = nc)
  for (i in 1:nc) {
    for (j in 1:nc) {
      sigma[i, j] <- s[i, j]
    }
  }
  return(sigma)
}


#' Construct Lagged Design Matrix for VAR
#'
#' Duplicates the original data matrix to create a lagged predictor matrix for VAR estimation.
#'
#' @param data A numeric matrix with time series data (observations in rows).
#' @param p Integer. The order of the VAR model (number of lags).
#'
#' @return A numeric matrix with duplicated columns corresponding to lagged observations.
#'
#' @keywords internal
duplicateMatrix <- function(data, p) {
  nr <- nrow(data)
  nc <- ncol(data)
  outputData <- data
  if (p > 1) {
    for (i in 1:(p - 1)) {
      tmpData <- matrix(0, nrow = nr, ncol = nc)
      tmpData[(i + 1):nr, ] <- data[1:(nr - i), ]
      outputData <- cbind(outputData, tmpData)
    }
  }
  outputData <- outputData[p:nr, ]
  return(outputData)
}


#' Split Coefficient Matrix into VAR Lags
#'
#' Splits a matrix of estimated coefficients into a list of matrices,
#' each corresponding to one lag of the VAR model.
#'
#' @param M A numeric matrix of coefficients.
#' @param p Integer. The order of the VAR model (number of lags).
#'
#' @return A list of \code{p} matrices, each of dimension (number of variables) x (number of variables).
#'
#' @keywords internal
splitMatrix <- function(M, p) {
  nr <- nrow(M)
  A <- list()
  for (i in 1:p) {
    ix <- ((i - 1) * nr) + (1:nr)
    A[[i]] <- M[ , ix]
  }
  return(A)
}

