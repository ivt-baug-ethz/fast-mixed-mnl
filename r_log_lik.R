devtools::load_all()


data("Train", package="mlogit")
head(Train, 3)
Train$ID <- Train$id
Train$CHOICE <- as.numeric(Train$choice)

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




## for each individual return log likelihood
## betas is named vector!
## see Molloy equation (3)
PnAt <- function(betas, data)
{
  nom <- exp(betas["B_price"] * data["price_A"] / 1000 + betas["B_time"] * data["time_A"] / 60 + betas["B_change"] * data["change_A"])
  denom <- nom + exp(betas["ASC_B"] + betas["B_price"] * data["price_B"] / 1000 + betas["B_time"] * data["time_B"] / 60 + betas["B_change"] * data["change_B"])
  
  return(as.numeric(nom / denom))
}


PnBt <- function(betas, data)
{
  nom <- exp(betas["ASC_B"] + betas["B_price"] * data["price_B"] / 1000 + betas["B_time"] * data["time_B"] / 60 + betas["B_change"] * data["change_B"])
  denom <- nom + exp(betas["B_price"] * data["price_A"] / 1000 + betas["B_time"] * data["time_A"] / 60 + betas["B_change"] * data["change_A"])
  
  return(as.numeric(nom / denom))
}


r_log_lik <- function(betas, data, Nindividuals = NULL, availabilities = NULL, nullableDraws = NULL, nDraws = NULL, P = NULL, weights = NULL, num_threads = 1L, p_indices = FALSE)
{
  p <-
    apply(data, 1, function(choice_task) {
      nm <- names(choice_task)
      choice_task <- as.numeric(choice_task)
      names(choice_task) <- nm
      if (choice_task["CHOICE"] == 1) {
        p <- PnAt(betas, choice_task)
      } else {
        p <- PnBt(betas, choice_task)
      }
    })
  
  return(sum(log(p)))
}


r_log_lik(betas, data)

r_model_spec <- mixl::specify_model(mnl_test, Train, logLik = r_log_lik)

debugonce(mixl::estimate)
model <- mixl::estimate(r_model_spec, est, Train, availabilities = availabilities)
