.PBFitDynr <- function(i,
                       path,
                       prefix,
                       dynr_initial,
                       dynr_measurement,
                       dynr_noise,
                       dynr_dynamics,
                       optimization_flag,
                       hessian_flag,
                       verbose,
                       weight_flag,
                       debug_flag,
                       perturb_flag,
                       xtol_rel,
                       stopval,
                       ftol_rel,
                       ftol_abs,
                       maxeval,
                       maxtime) {
  fn <- file.path(
    path,
    paste0(
      prefix,
      "_",
      "fit",
      "_",
      i,
      ".Rds"
    )
  )
  run <- function(path,
                  prefix,
                  dynr_initial,
                  dynr_measurement,
                  dynr_noise,
                  dynr_dynamics,
                  optimization_flag,
                  hessian_flag,
                  verbose,
                  weight_flag,
                  debug_flag,
                  perturb_flag,
                  xtol_rel,
                  stopval,
                  ftol_rel,
                  ftol_abs,
                  maxeval,
                  maxtime) {
    temp <- tempdir()
    on.exit(
      unlink(temp)
    )
    fn_data <- file.path(
      path,
      paste0(
        prefix,
        "_",
        "data",
        "_",
        i,
        ".Rds"
      )
    )
    fn_fit <- file.path(
      path,
      paste0(
        prefix,
        "_",
        "fit",
        "_",
        i,
        ".Rds"
      )
    )
    dynr_model <- dynr::dynr.model(
      data = .DynrData(
        object = readRDS(
          file = fn_data
        )
      ),
      initial = dynr_initial,
      measurement = dynr_measurement,
      dynamics = dynr_dynamics,
      noise = dynr_noise,
      outfile = tempfile(
        pattern = "dynr_",
        tmpdir = temp,
        fileext = ".c"
      )
    )
    dynr_model@options$xtol_rel <- xtol_rel
    dynr_model@options$stopval <- stopval
    dynr_model@options$ftol_rel <- ftol_rel
    dynr_model@options$ftol_abs <- ftol_abs
    dynr_model@options$maxeval <- maxeval
    dynr_model@options$maxtime <- maxtime
    if (verbose) {
      fit <- dynr::dynr.cook(
        dynr_model,
        optimization_flag = optimization_flag,
        hessian_flag = hessian_flag,
        verbose = verbose,
        weight_flag = weight_flag,
        debug_flag = debug_flag,
        perturb_flag = perturb_flag
      )
    } else {
      utils::capture.output(
        fit <- dynr::dynr.cook(
          dynr_model,
          optimization_flag = optimization_flag,
          hessian_flag = hessian_flag,
          verbose = verbose,
          weight_flag = weight_flag,
          debug_flag = debug_flag,
          perturb_flag = perturb_flag
        )
      )
    }
    saveRDS(
      object = fit,
      file = fn_fit
    )
    invisible()
  }
  tryCatch(
    {
      if (file.exists(fn)) {
        readRDS(
          file = fn
        )
      } else {
        run(
          path = path,
          prefix = prefix,
          dynr_initial = dynr_initial,
          dynr_measurement = dynr_measurement,
          dynr_noise = dynr_noise,
          dynr_dynamics = dynr_dynamics,
          optimization_flag = optimization_flag,
          hessian_flag = hessian_flag,
          verbose = verbose,
          weight_flag = weight_flag,
          debug_flag = debug_flag,
          perturb_flag = perturb_flag,
          xtol_rel = xtol_rel,
          stopval = stopval,
          ftol_rel = ftol_rel,
          ftol_abs = ftol_abs,
          maxeval = maxeval,
          maxtime = maxtime
        )
      }
    },
    error = function(e) {
      run(
        path = path,
        prefix = prefix,
        dynr_initial = dynr_initial,
        dynr_measurement = dynr_measurement,
        dynr_noise = dynr_noise,
        dynr_dynamics = dynr_dynamics,
        optimization_flag = optimization_flag,
        hessian_flag = hessian_flag,
        verbose = verbose,
        weight_flag = weight_flag,
        debug_flag = debug_flag,
        perturb_flag = perturb_flag,
        xtol_rel = xtol_rel,
        stopval = stopval,
        ftol_rel = ftol_rel,
        ftol_abs = ftol_abs,
        maxeval = maxeval,
        maxtime = maxtime
      )
    }
  )
  invisible()
}
