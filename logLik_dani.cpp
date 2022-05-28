#include <Rcpp.h>

#define _MIXED_MNL 1

#ifdef _OPENMP
#include <omp.h>
#endif
#include "inst/include/mixl/utility_function.h"

using Rcpp::RObject;
using Rcpp::DataFrame;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Rcpp::CharacterVector;
using Rcpp::_;
using Rcpp::stop;

// [[Rcpp::export]]
NumericVector logLik(NumericVector betas, DataFrame data,
                     int Nindividuals, NumericMatrix availabilities,
                     Nullable<NumericMatrix> nullableDraws, int nDraws,
                     NumericMatrix P, NumericVector weights, int num_threads=1, bool p_indices=false) {
  
  Rprintf("hello from dani\n");
  
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif  
  
  NumericMatrix draws;
  if (nullableDraws.isNotNull()) {
    draws = nullableDraws.get();
  }
  
  UF_args v(data, Nindividuals, availabilities, draws, nDraws, weights, P, p_indices);
  
  NumericVector LL(v.Nindividuals);
  std::fill(v.P.begin(), v.P.end(), 0);
  
  //double begin = omp_get_wtime();
  
  utilityFunction(betas, v);
  
  const double lognDraws = std::log(v.nDraws);
  
  //TODO: parallelise this as well
  //pragma omp parallel for
  for (int i=0; i<v.Nindividuals; i++) {
    double s = 0;
    for (int draw=0; draw<v.nDraws; draw++) {
      v.P(i,draw) = exp(v.P(i,draw));
      s += v.P(i,draw);
    }
    
    LL[i] = std::log(s) - lognDraws;
    
  }
  
  return LL;
}

//idea - preprocess the c++ code to declare all the variables at compilation time.
//or - through r, check the names in the utility function, that they are in the data, and return error if not. Then desugarise and compile
//need to distinquish between betas, random-coeefs and parameters

void utilityFunction(NumericVector betas, UF_args& v)
{
  
  if (!(v.data.containsElementNamed("ID") && v.data.containsElementNamed("CHOICE"))) {
    stop("Both ID and CHOICE columns need to be present in the data");
  }
  
  //delcare the variables that you will be using from the dataframe
  const NumericVector row_ids = v.data["ID"];
  const NumericVector choice = v.data["CHOICE"];
  
  NumericVector count;
  
  if (v.include_probability_indices){
    count = v.data["count"];
  }
  /////////////////////////////////////
  
  
  //betas
  double ASC_B = betas["ASC_B"];
  double SIGMA_B = betas["SIGMA_B"];
  double B_price = betas["B_price"];
  double B_time = betas["B_time"];
  double B_change = betas["B_change"];
  double B_timeB = betas["B_timeB"];
  
  
  //data
  const NumericVector data_price_A = v.data["price_A"];
  const NumericVector data_time_A = v.data["time_A"];
  const NumericVector data_change_A = v.data["change_A"];
  const NumericVector data_price_B = v.data["price_B"];
  const NumericVector data_time_B = v.data["time_B"];
  
  /////////////////////////////////////
  
#pragma omp parallel
{
  std::vector<double> utilities(2);  //specify here the number of alternatives
  
#pragma omp for
  for (int i=0; i < v.data.nrows(); i++) {
    
    int individual_index = row_ids[i]-1; //indexes should be for c, ie. start at 0
    //Rcpp::Rcout << "indv: " << individual_index << std::endl;
    for (int d=0; d<v.nDraws; d++) {
      
#ifdef _MIXED_MNL
      int draw_index = individual_index * v.nDraws + d; //drawsrep give the index of the draw, based on id, which we dont want to carry in here.
      NumericMatrix::ConstRow draw = v.draws(draw_index, _);
#endif
      
      std::fill(std::begin(utilities), std::end(utilities), 0.0);
      
      /////////////////////////
      
      
      double ASC_B_RND 	= ASC_B 	+ draw[0] * SIGMA_B;
      
      utilities[0] =             B_price * data_price_A[i] / 1000 + B_time * data_time_A[i] / 60 + B_change * data_change_A[i]; 
      utilities[1] = ASC_B_RND + B_price * data_price_B[i] / 1000 + B_timeB * data_time_B[i] / 60;
      
      
      /////////////////////////
      
      //dont edit beflow this line
      for (unsigned k=0; k < utilities.size(); ++k) {
        utilities[k] = std::min(700.0, std::max(-700.0, utilities[k])); //trip utilities to +- 700 for compuational reasons
        utilities[k] = exp(utilities[k]); //take the exponential of each utility
      }
      
      double chosen_utility = utilities[choice[i]-1]; //this -1 is needed if the choices start at 1 (as they should)
      
      double sum_utilities = 0.0;
      NumericMatrix::ConstRow  choices_avail = v.availabilities( i , _ );
      for (unsigned k=0; k < utilities.size(); ++k) {
        sum_utilities += utilities[k] * choices_avail[k];
      }
      
      double log_p_choice = std::log((chosen_utility / sum_utilities))  * v.weights[i];
      
      if (v.include_probability_indices){
        
        double p_indic_total = 0;
        //note: not a hybrid choice model
        log_p_choice += (1/count[i]) * std::log(p_indic_total) * v.weights[i];
      }
      
#pragma omp atomic 
      v.P(individual_index, d) += log_p_choice; //sum up the draws as we go along.
      
    }
  }
}
}