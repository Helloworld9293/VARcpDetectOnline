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
#' @importFrom sparsevar fitVAR
#' @importFrom stats qnorm
#' @examples
#' library(sparsevar)
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
#'
#' @examples
#' # Example Usage: Applying Change Point Detection to S&P 500 Data
#' # This is an example of how to apply the change point detection method
#' # (using the VAR_cpDetect_Online function) on the daily log return
#' # dataset from the S&P 500 (stored in the sp500 dataset). The code
#' # below calculates the average return volatility for all stocks, applies
#' # the change point detection algorithm, and plots the results with detected
#' # change points shown as vertical red and black lines.
#' \dontrun{
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
