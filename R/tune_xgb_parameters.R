#' Tune Parameters
#'
#' @param training_prop the proportion of data that should be sampled for the training set to do cross fold validation.
#' @param parameters a list of parameters and corresponding values you want to investigate. This should match what the parameters are called in XGBoost.
#' @param nrounds how many iterations of xgboost you want to perform
#' @param motif_counts your motif counts data frame
#'
#' @return a data frame with the parameter values tested and some evaluation metrics
#' @export

tune_xgb_parameters <- function(training_prop, parameters, nrounds, motif_counts) { # a function which does
                                                                           # gridwise parameter tuning for you geared towards
                                                                           # xgboost for motif finding, will provide parameter
                                                                           # results along with good parameter options

  #### i want to add in a parameter check where if somebody gives a parameter that doesn't
  #### exist then it won't run

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

    number_of_tuning_rounds <- 1 # need to figure out just how many things we are going to be testing

    for(i in 1:length(parameters)){
      number_of_tuning_rounds <- number_of_tuning_rounds*length(parameters[[i]]) # number of parameter comparisons/xgboost runs
    }

    xgboost_parameter_metrics_dataframe <- data.frame(matrix(nrow=number_of_tuning_rounds,ncol=number_col_for_output)) # to save each run of the parameter tuning results

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

                xgboost_parameter_metrics_dataframe[n,i] <- current_parameter_vector[j]
                names(xgboost_parameter_metrics_dataframe)[i] <- names(parameters)[[i]]

                n = n + 1

              }
            }

          } else{

            for(j in 1:length(current_parameter_vector)){

              for(k in 1:number_of_repeats){

                new_addition <- list(current_parameter_vector[j])
                names(new_addition) <- names(parameters)[[i]]
                all_parameters_tested[[n]] <- append(all_parameters_tested[[n]],new_addition)

                xgboost_parameter_metrics_dataframe[n,i] <- current_parameter_vector[j]
                names(xgboost_parameter_metrics_dataframe)[i] <- names(parameters)[[i]]

                n = n + 1

              }
            }
          }
        }
      }

    # return(all_parameters_tested)



    for(i in 1:length(all_parameters_tested)){

      toydataset_xgboostmodel <- xgboost::xgboost(data=trainingset2, early_stopping_rounds = 25, verbose = 0, nrounds = nrounds, objective = "binary:logistic", params = all_parameters_tested[[i]]) # , objective = "binary:logistic"

      pred <- predict(toydataset_xgboostmodel, validationset2)
      pred2 <- round(pred)

      genes_and_pred <- cbind(row.names(validationset),pred2)

      confusion.table <- table(observed=validationset[,ncol(validationset)],pred2)

      toydataset.accuracy <- sum(diag(confusion.table))/sum(confusion.table)

      xgboost_parameter_metrics_dataframe[i,(length(parameters) + 1)] <- toydataset.accuracy
      names(xgboost_parameter_metrics_dataframe)[(length(parameters) + 1)] <- "accuracy"

      if(ncol(confusion.table) > 1){
        toydataset.specificity <- confusion.table[1,1]/sum(confusion.table[1,])
        toydataset.sensitivity <- confusion.table[2,2]/sum(confusion.table[2,])
      } else{
        if(sum(pred2) == 0){
          toydataset.specificity <- confusion.table[1,1]/sum(confusion.table[1,])
          toydataset.sensitivity <- "NA"
        } else{
          toydataset.sensitivity <- confusion.table[1,1]/sum(confusion.table[1,])
          toydataset.specificity <- "NA"
        }
      }

      xgboost_parameter_metrics_dataframe[i,(length(parameters) + 2)] <- toydataset.sensitivity
      names(xgboost_parameter_metrics_dataframe)[(length(parameters) + 2)] <- "sensitivity"

      xgboost_parameter_metrics_dataframe[i,(length(parameters) + 3)] <- toydataset.specificity
      names(xgboost_parameter_metrics_dataframe)[(length(parameters) + 3)] <- "specificity"

      roc_full_resolution_xgb <- pROC::roc(validationset[,ncol(validationset)], pred, levels = c(0,1), direction = "<")

      auc_full_resolution_xgb <- pROC::auc(roc_full_resolution_xgb)

      xgboost_parameter_metrics_dataframe[i,(length(parameters) + 4)] <- auc_full_resolution_xgb
      names(xgboost_parameter_metrics_dataframe)[(length(parameters) + 4)] <- "auc"

      toydataset.variable.imp <- xgboost::xgb.importance(colnames(trainingset[,-ncol(trainingset)]), toydataset_xgboostmodel)

      number_of_features_identified <- nrow(toydataset.variable.imp)

      xgboost_parameter_metrics_dataframe[i,(length(parameters) + 5)] <- number_of_features_identified
      names(xgboost_parameter_metrics_dataframe)[(length(parameters) + 5)] <- "number_of_motifs_ided"

    }

    return(xgboost_parameter_metrics_dataframe)



  } else{ # if no parameters were provided then a warning is given

    print("Warning: No parameters provided")

  }

}
