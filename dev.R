devtools::load_all()

## draws are created in estimate (passed as matrix). If draws matrix is not specified
## the function mixl::create_halton_draws is used. Alternatively, you can specify a
## draws matrix (following any distribution you'd like) and then this will be used
## in the utility specification...

## mixl::convert_to_valid_cpp (see preprocessor.R) actually also parses the utility!

# simple model ------------------------------------------------------------

data("Train", package="mlogit")
head(Train, 3)
Train$ID <- Train$id
Train$CHOICE <- as.numeric(Train$choice)

mnl_test <- "
    ASC_B_RND 	= @ASC_B 	+ draw_2 * @SIGMA_B;

    U_A =             @B_price * $price_A / 1000 + @B_time * $time_A / 60 + @B_change * $change_A; 
    U_B = ASC_B_RND + @B_price * $price_B / 1000 + @B_timeB * $time_B / 60;
  "

model_spec <- mixl::specify_model(mnl_test, Train)

#only take starting values that are needed
est <- stats::setNames(c(0,0,0,0,0,0), c("B_price", "B_time", "B_timeB", 
                                         "B_change", "ASC_B", "SIGMA_B"))

availabilities <- mixl::generate_default_availabilities(Train, model_spec$num_utility_functions)

model <- mixl::estimate(model_spec, est, Train, availabilities = availabilities, nDraws = 20)  

summary(model)



# check where cpp log lik function ----------------------------------------

debugonce(mixl::specify_model)
model_spec <- mixl::specify_model(mnl_test, Train)    ## write e1 to global env

## debug specify_model -> e1$cpp_code -> which compiles -> logLik -> returns cpp_container (env) with logLik function!
ls(model_spec)
model_spec$logLik

## how to get there (extracted and slightly adjusted from specify_model)
template <- "loglik.cpp"
template_location <- system.file("include", "mixl", template, package = "mixl")
cpp_template <- readr::read_file(template_location)

# e1$cpp_code <-
cpp_code <- convert_to_valid_cpp(cpp_template, e1=e1)

#### some stuff left out ####

# Rcpp::sourceCpp(code = e1$cpp_code, env = cpp_container, ...)
Rcpp::sourceCpp(code = cpp_code)

## and tada: logLik is a callable function in global env!

## can adjust like so: LL[i] = std::log(s) - lognDraws;



# downstream implications -------------------------------------------------

## see ll2 line 81: ll2 is a function(betas) which gets passed to maxLik
## write ll2 to global env -> debugonce(ll2) -> write betas to global env
debugonce(mixl::estimate)
model <- mixl::estimate(model_spec, est, Train, availabilities = availabilities, nDraws = 20)

ll2(betas)

## ll2 is simply a wrapper around logLik passing all the arguments required to
## logLik - with the exemption of betas, which is the vector of coefficients to
## be estimated!



# extend specify_model ----------------------------------------------------

## TODO: extend function to pass alternative c++ likelihood function (passing R functions will work as well!)
## TODO: write R logLik equivalent and test

## for each individual return log likelihood
## betas is named vector!
r_log_lik <- function(betas, data)
{
  #### -> see script r_log_lik.R
  NULL
}




# test --------------------------------------------------------------------

Rcpp::sourceCpp("logLik_dani.cpp")

model_spec_lik <- mixl::specify_model(mnl_test, Train, logLik = logLik)
model_spec_lik$logLik

est <- stats::setNames(c(0,0,0,0,0,0), c("B_price", "B_time", "B_timeB", 
                                         "B_change", "ASC_B", "SIGMA_B"))

availabilities <- mixl::generate_default_availabilities(Train, model_spec_lik$num_utility_functions)

debugonce(mixl::estimate)
model <- mixl::estimate(model_spec_lik, est, Train, availabilities = availabilities, nDraws = 20)

