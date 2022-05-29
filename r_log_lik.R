## no need for mixl...

data("Train", package="mlogit")
head(Train, 3)
Train$ID <- Train$id
Train$CHOICE <- as.numeric(Train$choice)

## for each individual return log likelihood (or sum over all individuals)
## betas is named vector!
## see Molloy equation (3)
PnAt <- function(betas, data)
{
  nom <- exp(betas["B_price"] * data["price_A"] / 1000 + betas["B_time"] * data["time_A"] / 60 + betas["B_change"] * data["change_A"])
  denom <- nom + exp(betas["ASC_B"] + betas["B_price"] * data["price_B"] / 1000 + betas["B_timeB"] * data["time_B"] / 60 + betas["B_change"] * data["change_B"])
  
  return(as.numeric(nom / denom))
}


PnBt <- function(betas, data)
{
  nom <- exp(betas["ASC_B"] + betas["B_price"] * data["price_B"] / 1000 + betas["B_timeB"] * data["time_B"] / 60 + betas["B_change"] * data["change_B"])
  denom <- nom + exp(betas["B_price"] * data["price_A"] / 1000 + betas["B_time"] * data["time_A"] / 60 + betas["B_change"] * data["change_A"])
  
  return(as.numeric(nom / denom))
}


r_log_lik <- function(betas, data, e, verbose = TRUE)
{
  p <-
    apply(data, 1, function(choice_task) {
      nm <- names(choice_task)
      ct <- stringr::str_remove_all(choice_task, "[A-Za-z]")
      choice_task <- as.numeric(ct)
      names(choice_task) <- nm
      if (choice_task["CHOICE"] == 1) {
        p <- PnAt(betas, choice_task)
      } else {
        p <- PnBt(betas, choice_task)
      }
    })
  
  out <- sum(log(p))
  
  if (verbose)
  {
    cat("=")
    info <- sprintf("Call: %3d  -->  LL: %+4.2f", e$i, out)
    if (e$i %% 30 == 0) cat("\n", info, "\n")
    e$i <- e$i + 1
  }
  
  return(out)
}

e <- new.env()
e$i <- 1

LL <- function(betas) r_log_lik(betas, data = Train, e = e)
start <- c("B_price" = 0, "B_time" = 0, "B_change" = 0, "ASC_B" = 0, "B_timeB" = 0)
mL <- maxLik::maxLik(LL, start = start, method = "BFGS")






## comparison

devtools::load_all()


## mnl without rand
mnl_test <- "
    U_A =          @B_price * $price_A / 1000 + @B_time * $time_A / 60 + @B_change * $change_A; 
    U_B = @ASC_B + @B_price * $price_B / 1000 + @B_timeB * $time_B / 60;
  "

model_spec <- mixl::specify_model(mnl_test, Train)

#only take starting values that are needed
est <- stats::setNames(c(0,0,0,0,0), c("B_price", "B_time", "B_timeB", 
                                       "B_change", "ASC_B"))

availabilities <- mixl::generate_default_availabilities(Train, model_spec$num_utility_functions)

model <- mixl::estimate(model_spec, est, Train, availabilities = availabilities, nDraws = 20)  

summary(model)





## alternatively we can pass r_log_lik to mixl (kind of unnecessary though)

r_model_spec <- mixl::specify_model(mnl_test, Train, logLik = r_log_lik)
model <- mixl::estimate(r_model_spec, est, Train, availabilities = availabilities)
