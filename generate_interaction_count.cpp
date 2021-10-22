#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame interaction_count(DataFrame interaction_table, NumericVector q_value_thresholds) {
  int length_out = q_value_thresholds.length();
  IntegerVector n_significant = IntegerVector(length_out);
  IntegerVector n_negative = IntegerVector(length_out);
  int n_significant_current = 0;
  int n_negative_current = 0;
  int pos = 0;
  int i = 0;
  NumericVector q_values = interaction_table["q.value"];
  NumericVector sim_scores = interaction_table["sim.score"];
  while(i < length_out){
    /*
    Rcout << pos << "\t" << i << std::endl;
    Rcout << "q-value test"<< "\t" << q_values[pos] << std::endl;
    Rcout << "q-value threshold"<< "\t" << q_value_thresholds[pos] << std::endl;
    **/
    if(pos >= interaction_table.nrow()){
      break;
    }
    if(q_values[pos] < q_value_thresholds[i]){
      n_significant_current++;
      if(sim_scores[pos] < 0){
        n_negative_current++;
      }
      pos++;
    }
    else{
      n_significant[i] = n_significant_current;
      n_negative[i] = n_negative_current;
      i++;
    }
    }
  // Fills up the results for the q-values not covered by the main loop
  for (;i < length_out;i++){
    n_significant[i] = n_significant_current;
    n_negative[i] = n_negative_current;
  }
  return DataFrame::create(Named("threshold") = q_value_thresholds,Named("significant") = n_significant, Named("negative") = n_negative);
}
