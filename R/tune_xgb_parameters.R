tune_xgb_parameters <- function(training_prop, parameters, motif_counts) { # a function which does
                                                                           # gridwise parameter tuning for you geared towards
                                                                           # xgboost for motif finding, will provide parameter
                                                                           # results along with good parameter options

  trainingsize <- floor(nrow(motif_counts)*training_prop) # number of samples for training
  indexes <- sample(1:nrow(motif_counts), size = trainingsize) # pull the index for the samples
  trainingset <- motif_counts[indexes,] # pull the samples set aside for training
  validationset <- motif_counts[-indexes,] # pull the samples for testing

  # xgboost does require that the data be in xgboost format, these two lines do that

  trainingset2 <- xgboost::xgb.DMatrix(data = as.matrix(trainingset[,-ncol(trainingset)]), label = as.matrix(trainingset[,ncol(trainingset)]))
  validationset2 <- xgboost::xgb.DMatrix(data = as.matrix(validationset[,-ncol(validationset)]), label = as.matrix(validationset[,ncol(validationset)]))

  # analyzing the list of parameters to see which parameters need to be tested/examined

  if(length(parameters) > 0){ # as long as there is a list of parameters

    number_col_for_output <- length(parameters) + 5 # this is to figure out how large the output dataframe needs to be

    xgboost_parameter_metrics_dataframe <- data.frame(matrix(nrow=0,ncol=number_col_for_output)) # to save each run of the parameter tuning results

    number_of_tuning_rounds <- 1 # need to figure out just how many things we are going to be testing

    for(i in 1:length(parameters)){
      number_of_tuning_rounds <- number_of_tuning_rounds*length(parameters[[i]]) # number of parameter comparisons/xgboost runs
    }

    all_parameters_tested <- list()

    for(i in 1:length(parameters)){

      n <- 1

      while(n <= number_of_tuning_rounds){

        current_parameter_vector <- parameters[[i]]
        length_current_parameter_vector <- length(current_parameter_vector)
        number_of_repeats <- number_of_tuning_rounds/length_current_parameter_vector

          if(i < 2) {

            for(j in 1:length(current_parameter_vector)){

              for(k in 1:number_of_repeats){

                new_addition <- list(current_parameter_vector[j])
                names(new_addition) <- names(parameters)[[i]]
                all_parameters_tested[[n]] <- new_addition
                n = n + 1

              }

            }

          } else{

            for(j in 1:length(current_parameter_vector)){

              for(k in 1:number_of_repeats){

                new_addition <- list(current_parameter_vector[j])
                names(new_addition) <- names(parameters)[[i]]
                all_parameters_tested[[n]] <- append(all_parameters_tested[[n]],new_addition)
                n = n + 1

              }

            }

          }



        }



      }

    return(all_parameters_tested)

  } else{ # if no parameters were provided then a warning is given

    print("Warning: No parameters provided")

  }

}
