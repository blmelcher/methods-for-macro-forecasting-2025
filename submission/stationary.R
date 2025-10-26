# Check stationarity and transform time series to stationarity if needed.

check_stationarity <- function(ts_data, min_length = 10) {
  # Return TRUE if series appears stationary by ADF test, FALSE otherwise.
  # NA values are removed before testing. If the series is too short, return FALSE.
  vec <- na.omit(as.numeric(ts_data))
  if (length(vec) < min_length) {
    warning("Series too short for ADF test; returning FALSE.")
    return(FALSE)
  }
  adf <- tryCatch(
    tseries::adf.test(vec),
    error = function(e) {
      warning("ADF test failed: ", conditionMessage(e))
      # Return non-stationary on error
      return(list(p.value = 1))
    }
  )
  is_stat <- adf$p.value < 0.05
  if (is_stat) {
    message("Series stationary (ADF p-value = ", signif(adf$p.value, 3), ").")
  } else {
    message("Series non-stationary (ADF p-value = ", signif(adf$p.value, 3), ").")
  }
  return(is_stat)
}

make_stationary <- function(ts_data, use_log = TRUE, align = TRUE, tiny = 1e-6) {
  # Make a series approximately stationary.
  # - If use_log = TRUE, apply log (adding a small offset if needed) and then first difference.
  # - If use_log = FALSE, apply first difference.
  # - If align = TRUE (default), return a vector with the same length as input
  #   (leading NA(s) inserted to align differenced output with original timestamps).
  orig_len <- length(ts_data)
  orig_names <- names(ts_data)
  vec <- as.numeric(ts_data)

  # If already stationary, return original (aligned as requested)
  if (check_stationarity(vec)) {
    message("No transformation applied; series already stationary.")
    if (align) {
      out <- vec
      if (!is.null(orig_names)) names(out) <- orig_names
      return(out)
    } else {
      return(vec)
    }
  }

  # Apply transformation
  if (use_log) {
    # handle non-positive values by adding a small offset
    if (any(vec <= 0, na.rm = TRUE)) {
      offset <- abs(min(vec, na.rm = TRUE)) + tiny
      warning("Non-positive values found; adding offset = ", signif(offset, 6), " before log.")
      vec_adj <- vec + offset
    } else {
      vec_adj <- vec
    }
    trans <- diff(log(vec_adj))
    message("Applied log + first difference.")
  } else {
    trans <- diff(vec)
    message("Applied first difference.")
  }

  # Return aligned or raw differenced series
  if (align) {
    out <- rep(NA_real_, orig_len)
    if (orig_len > 1) out[-1] <- as.numeric(trans)
    if (!is.null(orig_names)) names(out) <- orig_names
    return(out)
  } else {
    if (!is.null(orig_names)) names(trans) <- orig_names[-1]
    return(as.numeric(trans))
  }
}