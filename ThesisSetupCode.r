#### Initializiton ####
rm(list = ls())

library(data.table)
library(tidyverse)
library(rJava)
library(RNetLogo)

library(lhs)  # For maximin Latin hypercube sampling
library(ggplot2)
library(plotly)  # For beautiful plotting
library(caret)
library(randomForest)
library(factoextra)
library(e1071)
library(TSrepr)  # for evaluating predictive power

require(gridExtra)

options(warn = -1)

log_entry <- function(){
    write(paste0( "model =",nl.model,"\n"
             ,"nofrep =",nofrep,"\n"
             ,paste0(c("metarep =",metarep),collapse = " "),"\n"
             ,"ntree =",ntree,"\n"
             #,"mtry =",mtry,"\n"
             ,"mtry.multiplier =",mtry.multiplier,"\n"
             ,"nperm =",nperm,"\n"
             ,"iteration_budget =",iteration_budget,"\n"
             ,"unlabeled_ins =",unlabeled_ins,"\n"
             ,"unlabeled.type =",unlabeled.type,"\n"
             ,paste0(c("test_ins =",test_ins),collapse = " "),"\n"
             ,"train_ins_Ad =",train_ins_Ad,"\n"
             ,"selected_ins =",selected_ins,"\n"
             ,"p =",p,"\n"
             ,"oob_allowance =",oob_allowance,"\n"  
             ,paste0(c("seed.focus =",seed.focus),collapse = " "),"\n"             
             ,"error_type =",error_type,"\n" 
             ,"eliminaton.type =",elimination.type,"\n"
             ,"elimination_start_iter =",elimination_start_iter,"\n"    
             ,"selection_metric =",selection_metric,"\n" 
             ,"Date =", Sys.Date()
             )
      ,ReadMe, append=TRUE, sep = "\n" )   
}

log_entry.sampling <- function(){
    write(paste0( "model =",nl.model,"\n"
             ,"nofrep =",nofrep,"\n"
             ,paste0(c("metarep =",metarep),collapse = " "),"\n"
             ,"ntree =",ntree,"\n"
             ,"mtry =",mtry,"\n"
             ,"nperm =",nperm,"\n"
             ,"iteration_budget =",iteration_budget,"\n"
             ,"unlabeled_ins =",unlabeled_ins,"\n"
             ,"unlabeled.type =",unlabeled.type,"\n"
             ,paste0(c("test_ins =",test_ins),collapse = " "),"\n"
             ,"train_ins_Ad =",train_ins_Ad,"\n"
             ,"selected_ins =",selected_ins,"\n"
             ,"h =",h,"\n"
             ,paste0(c("seed.focus =",seed.focus),collapse = " "),"\n"             
             ,"error_type =",error_type,"\n" 
             ,"sample.type =",sample.type,"\n"
             ,"selection_metric =",selection_metric,"\n" 
             ,"Date =", Sys.Date()
             )
      ,ReadMe, append=TRUE, sep = "\n" )   
}

run_log_entry <- function(){
    fwrite( data.table(seed = i, rep = r, iter = iter, start = Sys.time())
           ,paste0(outputs.path,model.type,"_run_logs_",sample.type,".csv"), append = TRUE)
}

run_step_log_entry <- function(log){
    fwrite( data.table(seed = i, rep = r, iter = iter, start = Sys.time(), step = log)
           ,paste0(outputs.path,model.type,"_run_step_logs_",sample.type,".csv"), append = TRUE)
}

# run_model <- function(feature_names,feature_values){ # both should be in character list format both should be in character list format
run_model <- function(feature_values) {
    k = length(feature_names)
    for (i in 1:k) {
        NLCommand(paste0("set ", feature_names[i], " ", feature_values[i]), nl.obj = nl.model)
    }
    NLCommand("setup", nl.obj = nl.model)
    NLDoCommand(30, "go", nl.obj = nl.model)
    result <- NLReport(output_name, nl.obj = nl.model)
    return(result)
}

# run_replicas <- function(nofrep,feature_names,feature_values) {
run_replicas <- function(nofrep, feature_values) {
    replicas = matrix(NA, ncol = nofrep, nrow = 1)  # Save the result of each replication
    for (i in 1:nofrep) {
        # replicas[i]= run_model(feature_names,feature_values)
        replicas[i] = run_model(feature_values)
    }
    aggregated_result = mean(replicas)
    return(aggregated_result)
}

# run_ABM = function(nofrep,nofinstances,unlabeledset,featurenames = feature_names){
run_ABM = function(nofrep, nofinstances, unlabeledset) {
    # unlabeledset = setcolorder(unlabeledset,featurenames)
    unlabeledset = setcolorder(unlabeledset, feature_names)
    for (i in 1:nofinstances) {
        # unlabeledset[i, output := run_replicas(nofrep,featurenames,
        # as.matrix(unlabeledset[i,]))]
        unlabeledset[i, `:=`(output, run_replicas(nofrep, as.matrix(unlabeledset[i,])))]
        NLQuit(all = TRUE)
        NLStart(nl.path, gui = FALSE,nl.jarname='netlogo-6.0.4.jar', nl.obj=nl.model)
        NLLoadModel (model.path, nl.obj=nl.model)
    }
    return(unlabeledset)
}

# error functions on test data
rmse_func <- function(actual, predicted) {
    error = predicted - actual
    return(sqrt(mean(error^2)))
}

mape_func <- function(actual, predicted) {
    return((abs(actual - predicted)/actual) * 100)
}

bias_func <- function(actual, predicted) {
    return((actual - predicted)/actual)
}

# error functions on train data
obb_error_func <- function(model) {
    if (model$type == "regression") {
        oob_error = model$mse[model$ntree]
    } else if (model$type == "classification") {
        oob_error = model$err.rate
    }
    return(oob_error)
}

# prediction functions
get_test_predictions <- function(model, testset, errortype) {
    
    predictedLabels <- predict(model, testset)
    predictedLabels <- cbind(testset, predictedLabels)
    setnames(predictedLabels, "predictedLabels", "pred_output")
    
    predictedLabels[, `:=`(MAPE = mapply(function(x, y) mape_func(x, y), output, pred_output), 
                           RMSE = mapply(function(x, y) rmse_func(x, y), output, pred_output), 
                           BIAS = mapply(function(x, y) bias_func(x, y), output, pred_output))]
    
    performance_table = predictedLabels[, .(mae  = mae(output,pred_output),
                                            rmse = rmse(output,pred_output),
                                            mape = mape(output,pred_output))
                                        ,.(size)]
    
    
    return(list(predictedLabels, performance_table))    
}

# Adaptive sample selection function with an uncertainty measure depending on 'selection_metric'
sample_selection <- function(selected_ins, unlabeled_set, model,selection_metric) {
    ind_pred <- t(predict(model, unlabeled_set, predict.all = TRUE)$individual) %>% 
        data.table()  # predictions by each tree in the forest
    ind_pred_eval = data.table()  
 #   ind_pred_eval[, `:=`(idx, 1:.N)]
    
    if (selection_metric == "sd") {
        # standard deviation calculation
        s_dev = sapply(ind_pred, sd) %>% data.table()
        setnames(s_dev, ".", "sd")
        ind_pred_eval = cbind(ind_pred_eval, s_dev)
        
        ind_pred_eval[, `:=`(idx, 1:.N)]
        ind_pred_eval = ind_pred_eval[order(-sd)][1:selected_ins]
    } else if (selection_metric == "range") {
        # range calculation
        range = sapply(ind_pred, range) %>% t() %>% data.table()
        range = range[, .(range = abs(range[, 1] - range[, 2]))]
        setnames(range, "range.V1", "range")
        ind_pred_eval = cbind(ind_pred_eval, range)
        
        ind_pred_eval[, `:=`(idx, 1:.N)]
        ind_pred_eval = ind_pred_eval[order(-range)][1:selected_ins]
    } else if (selection_metric == "coefvar") {
        #coeff variance calculation
        s_dev = sapply(ind_pred, sd) %>% data.table()
        setnames(s_dev, ".", "sd")
        s_mean = sapply(ind_pred, mean) %>% data.table()
        setnames(s_mean, ".", "mean")
        coeff_var = cbind(s_dev,s_mean) 
        coeff_var = coeff_var[,.(c_var = (sd / mean)* 100)]
        ind_pred_eval = cbind(ind_pred_eval, coeff_var)
        
        ind_pred_eval[, `:=`(idx, 1:.N)]
        ind_pred_eval = ind_pred_eval[order(-c_var)][1:selected_ins]
    }
    
    unlabeled_set[, `:=`(idx, 1:.N)]
    train_candidates = unlabeled_set[ind_pred_eval$idx]
    
    return(train_candidates)
}

# Random sample selection
random_sample_selection <- function(selected_ins, unlabeled_set) {
    
    unlabeled_set[, `:=`(idx, 1:.N)]
    
    train_candidate_idx = sample(unlabeled_set$idx, selected_ins, replace = FALSE, prob = NULL)
    train_candidates = unlabeled_set[idx %in% train_candidate_idx]
    
    return(train_candidates)
}

dt <- data.table(a = c('a','b','c'),b = c('1','2','3'),c = c('x','y','z'))
v <- c('A','B','C')

random_sample_selection(1,dt)

get_variable_importance = function(model) {
    
    imp <- importance(model, type = 1, scale = FALSE)
    names = colnames(t(imp))
    ranked_features = names[order(imp, decreasing = TRUE)]
    return(ranked_features)
}

#get_variable_importance <- function(model) {
#    importances <- importance(model, type = 1, scale = FALSE)
#    selected.vars <- order(importances, decreasing = TRUE)
#    ranked_features = feature_names[selected.vars]
#    ordered.importances <- importances[selected.vars]
#    
#    return(ranked_features)
#}

##unique to each seed and rep, updated iter by iter
get_importance_data = function(model,type,scaled,n){
    #n shows how many iterations' values will be evaluated (moving average n)
    iteration_imp = data.table(t(importance(model, type = type, scaled = scaled))) %>% melt()
    setnames(iteration_imp,"variable","feature")
    
    imp_history = if(nrow(importance_data) > 0){
        copy(importance_data[iter_no >= 1 & iter_no > (iter - n)])
    }else{
        data.table()
    }
          
    imp_history = rbind(imp_history, data.table(iter_no = iter
                                                ,iteration_imp
                                                ,cm_mean = as.numeric(NA)))
    imp_history[, cm_mean := mean(value),.(feature)]
    
    return(imp_history[iter_no == iter])
    
}

#importance_data = rbind(importance_data, get_importance_data(model,1,TRUE,3))

# in each iteration the result of this function is checked
#importance data should hold historic data for an seed-rep combination through iterations
importance_partition = function(grid_partition, iteration_importance_data, column_name) {
    # name of the column which grid separation performed will be based on
    iteration_importance_data[, `:=`(grid_no, as.numeric(NA))]
    iteration_importance_data[get(column_name) < 0, `:=`(grid_no, 0)]
    
    grid = 0
    order = 0
    
    while (is.na(iteration_importance_data[which.max(iteration_importance_data[[column_name]])]$grid_no)) {
        order = order + 1
        
        grid_up = grid + grid_partition
        
        idx = which(iteration_importance_data[[column_name]] < grid_up & iteration_importance_data[[column_name]] > grid)       
        iteration_importance_data[idx, `:=`(grid_no, order)]        
        grid = copy(grid_up)       
    }
    iteration_importance_data[order(grid_no), `:=`(grid_order, .GRP), .(grid_no)]
    
    return(iteration_importance_data)
}

start_elimination = function(importance_data_partioned){
    #partitioned_importance_data
    importance_data_partioned = importance_data_partioned[order(iter_no)]
    importance_data_partioned[,prev_grid_order := shift(grid_order,1, type = "lag"),.(feature)]
    importance_data_partioned[,order_diff :=  grid_order - prev_grid_order ]

    elimination_start_iter = ifelse(nrow(importance_data_partioned[order_diff == 0,.N,.(iter_no)][N == length(feature_names)]) > 0,
                                   importance_data_partioned[order_diff == 0,.N,.(iter_no)][N == length(feature_names)]$iter_no,
                                   0)
    return(elimination_start_iter)
}
# returns the iteration number where the elimination starts

feature_elimination = function(h, feature_list = columns_left, ranked_features) {
    
    numof_columns_left = length(feature_list) - (h)
    columns_left = ranked_features[1:numof_columns_left]
    
    eliminated_columns = c(eliminated_columns, setdiff(feature_list, columns_left))
    
    # update total_numof_eliminated_vars
    total_numof_eliminated_vars = length(feature_names) - length(columns_left)
    
    return(list(columns_left, total_numof_eliminated_vars, h, eliminated_columns))
}

refresh_sample_pool <- function(selected.seed, columns_left = feature_names) {
    set.seed(selected.seed)
    
    k = length(columns_left)
    unlabeled_pool = as.data.table(maximinLHS(n = unlabeled_ins, k = k, dup = 5))
    setnames(unlabeled_pool, c(paste0("V", 1:k)), columns_left)
    
    for (c in 1:k) {
        unlabeled_pool[[c]] = qunif(  unlabeled_pool[[c]]
                                    , feature_ranges[feature == colnames(unlabeled_pool)[c]]$min_range
                                    , feature_ranges[feature == colnames(unlabeled_pool)[c]]$max_range)
    }
    
    random_pool_all = data.table()
    eliminated_columns = setdiff(feature_names, columns_left)
    
    if (nofparams > k) {
        for (e in 1:length(eliminated_columns)) {
            random_pool = data.table(runif(  unlabeled_ins
                                           , feature_ranges[feature == eliminated_columns[e]]$min_range
                                           , feature_ranges[feature == eliminated_columns[e]]$max_range))
            setnames(random_pool, eliminated_columns[e])
            
            random_pool_all = cbind(random_pool_all, random_pool)
        }
        
        unlabeled_pool = cbind(unlabeled_pool, random_pool_all)
    }
    
    return(unlabeled_pool)
}

upload_training_set <- function(model.type,seed.list,data.size){
    training_set_all = data.table()
    for( i in seed.list){
        training_set.name= paste0(data.path,"training_set","_",model.type,"_",data.size,"_seed",i,".csv")
        training_set <- fread(training_set.name) 
    
        training_set_all = rbind(training_set_all,data.table(training_set, "seed" = i))
        rm(training_set,training_set.name)  
        }
   return(training_set_all)
}

write_importance.rf = function(seed,rep,iter,model,sample.type){
    importance_table = data.table()
    importance_table = rbind(importance_table, data.table(seed = seed, rep = rep, iter_no = iter
                                                      , scaled = "yes", type = "1",
                                                      t(importance(model,type=1,scaled = TRUE)))
                            , use.names = FALSE)
    importance_table = rbind(importance_table, data.table(seed = seed, rep = rep, iter_no = iter
                                                  , scaled = "yes", type = "2",
                                                  t(importance(model,type=2,scaled = TRUE)))
                            , use.names = FALSE)
    importance_table = rbind(importance_table, data.table(seed = seed, rep = rep, iter_no = iter
                                                  , scaled = "no", type = "1",
                                                  t(importance(model,type=1,scaled = FALSE)))
                            , use.names = FALSE)
    importance_table = rbind(importance_table, data.table(seed = seed, rep = rep, iter_no = iter
                                                  , scaled = "no", type = "2",
                                                  t(importance(model,type=2,scaled = FALSE)))
                            , use.names = FALSE)
    
    importance_table = melt(importance_table, id.vars = c("seed","rep","iter_no","scaled","type")
                           , measure.vars = colnames(importance_table[,.SD,.SDcols = -c("seed","rep","iter_no","scaled","type")]))
    setnames(importance_table,c("variable"),c("feature"))
    
    fwrite(importance_table,paste0(outputs.path,model.type,"_importance_table_",sample.type,".csv"), append = TRUE)

}

write_test_accuracy = function(seed, rep, iter, model, test_set, error_type) {
        
    test_predictions_Sub = get_test_predictions(model, test_set, error_type)
    
    predictedLabels_Sub = data.table(seed = i,rep = r,iter,test_predictions_Sub[[1]])
    fwrite(predictedLabels_Sub, paste0(outputs.path, PL.folder, model.type, "_","predictedLabels.", sample.type, "_seed_", i, "_iter_", iter, "_rep_", r, ".csv"))
    
    performance_table_Sub = data.table(iter, seed = i,rep = r, test_predictions_Sub[[2]])
    fwrite(performance_table_Sub, paste0(outputs.path, sample.folder, model.type, "_","performance_table_", sample.type, ".csv"), append = TRUE)
          
}

mtry_default = function(features){
     max(floor(length(features)/3), 1)
}
