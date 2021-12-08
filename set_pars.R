global_pars <- c(
  lod = 40							# Limit of detection
)
run_pars_list <- list(
  # UNINFORMATIVE PRIORS: ---------------------------------------------------
  list(tpsd=2,
       dpmin=0,
       dpmean_prior=global_pars[["lod"]]/2,
       dpsd_prior=global_pars[["lod"]]/6,
       wpmin=0.25,
       wpmax=14,
       wpmean_prior=14/2,
       wpsd_prior=14/6,
       wrmin=2,
       wrmax=30,
       wrmean_prior=30/2,
       wrsd_prior=30/6,
       sigma_max=10,
       sigma_prior_scale=5,
       lambda=0.01,
       fpmean=1/log(10),
       trapfit=FALSE
  ),
  # PRIORS INFORMED BY PREVIOUS ANALYSIS: -----------------------------------
  list(tpsd=2,
       dpmin=0,
       dpmean_prior=global_pars[["lod"]]/2,
       dpsd_prior=global_pars[["lod"]]/6,
       wpmin=0.25,
       wpmax=14,
       wpmean_prior=2.7,
       wpsd_prior=14/6,
       wrmin=2,
       wrmax=30,
       wrmean_prior=7.4,
       wrsd_prior=30/6,
       sigma_max=10,
       sigma_prior_scale=5,
       lambda=0.01,
       fpmean=1/log(10),
       trapfit=FALSE
  ),
  # LOW PRIORS: -------------------------------------------------------------
  list(tpsd=2,
       dpmin=0,
       dpmean_prior=global_pars[["lod"]]/2,
       dpsd_prior=global_pars[["lod"]]/6,
       wpmin=0.25,
       wpmax=14,
       wpmean_prior=0,
       wpsd_prior=14/6,
       wrmin=2,
       wrmax=30,
       wrmean_prior=0,
       wrsd_prior=30/6,
       sigma_max=10,
       sigma_prior_scale=5,
       lambda=0.01,
       fpmean=1/log(10),
       trapfit=FALSE
  ))

