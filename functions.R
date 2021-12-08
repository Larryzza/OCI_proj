make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
  # Use "reduce" here
  out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
    as_tibble %>% 
    mutate(iteration=1:n()) %>% 
    pivot_longer(-iteration) %>% 
    separate(name, c("param","id"), sep="_") %>%
    pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}

# For parsing the 'parameters' output from Stan: 
parseparam <- function(extracted_params, parname, n_indiv){
  as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}


launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(
    launch_shinystan(fit) #You can replace this with any other Shiny app
  )
}

make_shared_params_df <- function(extracted_params, parnames){
  # Use "reduce" here
  out <- reduce(lapply(parnames, function(x) 
    as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
  ), cbind) %>%
    as_tibble() %>%
    mutate(iteration=1:n())
}

# For generating names for a matrix-turned-data frame: 
makenames <- function(parname, n_indiv){
  unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}

truncnormmean <- function(mu, sigma, T0, T1){
  alpha <- (T0 - mu)/sigma
  beta <- (T1 - mu)/sigma
  out <- mu + 
    (dnorm(alpha, 0, 1) - dnorm(beta, 0, 1))/
    (pnorm(beta, 0, 1) - pnorm(alpha, 0, 1))*sigma
  return(out)
}
