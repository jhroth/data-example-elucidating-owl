# Jeremy Roth
rm(list=ls())
library(DevTreatRules)
library(glmnet)

load("Results/datasets_for_data_example_with_missingness_weights.RData")
## set parameters
name.outcome <- "no_BREAST_after_10_yr" #"no_CHD_after_10_yr"
bootstrap.CI.replications <- 100
vec.building.propensity.method <- c("ridge", "lasso", "logistic.regression")
vec.building.rule.method <- c("ridge", "lasso", "glm.regression")

## RUN SPLIT REGRESSION
total.combinations.split.regression <- length(vec.building.propensity.method) * length(vec.building.rule.method)
run.split.regression <- TRUE
if (run.split.regression == TRUE) {
    # initialize objects to store results
    list.BuildRules.split.regression <- vector("list", total.combinations.split.regression)
    list.EvaluateRules.split.regression <- vector("list", total.combinations.split.regression)
    mat.split.regression <- matrix(NA, nrow=total.combinations.split.regression, ncol=5)
    colnames(mat.split.regression) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
    rownames(mat.split.regression) <- rep(NA, total.combinations.split.regression)
    
    before.building.split.regression <- proc.time()
    row.number <- 1
    for (p in 1:length(vec.building.propensity.method)) {
        for (r in 1:length(vec.building.rule.method)) {
            print(row.number)
            one.propensity.method <- vec.building.propensity.method[p]
            one.rule.method <- vec.building.rule.method[r]
            
           # building rule on training set
            set.seed(row.number)
            one.split.regression.rule <- BuildRule(data=training.data,
                                                    study.design="observational",
                                                    prediction.approach="split.regression",
                                                    name.outcome=name.outcome,
                                                    type.outcome="binary",
                                                    desirable.outcome=TRUE,
                                                    name.treatment=name.treatment,
                                                    names.influencing.treatment=names.influencing.treatment,
                                                    names.influencing.rule=names.influencing.rule,
                                                    propensity.method=one.propensity.method,
                                                    rule.method=one.rule.method,
                                                    additional.weights=training.data$IPW.CC)
            
            ## evaluate on validation
            evaluate.split.regression.rule.on.validation <- EvaluateRule(data=validation.data,
                                                                          BuildRule.object=one.split.regression.rule,
                                                                          study.design="observational",
                                                                          name.outcome=name.outcome,
                                                                          type.outcome="binary",
                                                                          desirable.outcome=TRUE,
                                                                          clinical.threshold=0,
                                                                          name.treatment=name.treatment,
                                                                          names.influencing.treatment=names.influencing.treatment,
                                                                          names.influencing.rule=names.influencing.rule,
                                                                          propensity.method=one.propensity.method,
                                                                          additional.weights=validation.data$IPW.CC,
                                                                          bootstrap.CI=FALSE,
                                                                          bootstrap.CI.replications=bootstrap.CI.replications)
            
            one.rowname <- paste0("propensity_", one.propensity.method, "_rule_", one.rule.method)
            rownames(mat.split.regression)[row.number] <- one.rowname
            mat.split.regression[one.rowname, 
                                 c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
                as.numeric(evaluate.split.regression.rule.on.validation[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
            
            list.BuildRules.split.regression[[row.number]] <- one.split.regression.rule
            names(list.BuildRules.split.regression)[[row.number]] <- one.rowname
            list.EvaluateRules.split.regression[[row.number]] <- evaluate.split.regression.rule.on.validation
            names(list.EvaluateRules.split.regression)[[row.number]] <- one.rowname
            row.number <- row.number + 1
        }
    }
    total.time.building.split.regression <- proc.time() - before.building.split.regression
    save(list.BuildRules.split.regression, list.EvaluateRules.split.regression,
         file=paste0("Results/", "list_split.regression.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
    save(total.time.building.split.regression,
         mat.split.regression, 
         file=paste0("Results/", "mat_split.regression.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
}


## RUN DIRECT.INTERACTIONS
total.combinations.direct.interactions <- length(vec.building.propensity.method) * length(vec.building.rule.method)
run.direct.interactions <- TRUE
if (run.direct.interactions == TRUE) {
    # initialize objects to store results
    list.BuildRules.direct.interactions <- vector("list", total.combinations.direct.interactions)
    list.EvaluateRules.direct.interactions <- vector("list", total.combinations.direct.interactions)
    mat.direct.interactions <- matrix(NA, nrow=total.combinations.direct.interactions, ncol=5)
    colnames(mat.direct.interactions) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
    rownames(mat.direct.interactions) <- rep(NA, total.combinations.direct.interactions)
    
    before.building.direct.interactions <- proc.time()
    row.number <- 1
    for (p in 1:length(vec.building.propensity.method)) {
        for (r in 1:length(vec.building.rule.method)) {
            print(row.number)
            one.propensity.method <- vec.building.propensity.method[p]
            one.rule.method <- vec.building.rule.method[r]
            
            # building rule on training set
            set.seed(row.number)
            one.direct.interactions.rule <- BuildRule(data=training.data,
                                                       study.design="observational",
                                                       prediction.approach="direct.interactions",
                                                       name.outcome=name.outcome,
                                                       type.outcome="binary",
                                                       desirable.outcome=TRUE,
                                                       name.treatment=name.treatment,
                                                       names.influencing.treatment=names.influencing.treatment,
                                                       names.influencing.rule=names.influencing.rule,
                                                       propensity.method=one.propensity.method,
                                                       additional.weights=training.data$IPW.CC,
                                                       rule.method=one.rule.method)
            
            ## evaluate on validation
            evaluate.direct.interactions.rule.on.validation <- EvaluateRule(data=validation.data,
                                                                             BuildRule.object=one.direct.interactions.rule,
                                                                             study.design="observational",
                                                                             name.outcome=name.outcome,
                                                                             type.outcome="binary",
                                                                             desirable.outcome=TRUE,
                                                                             clinical.threshold=0,
                                                                             name.treatment=name.treatment,
                                                                             names.influencing.treatment=names.influencing.treatment,
                                                                             names.influencing.rule=names.influencing.rule,
                                                                             propensity.method=one.propensity.method,
                                                                             additional.weights=validation.data$IPW.CC,
                                                                             bootstrap.CI=FALSE,
                                                                             bootstrap.CI.replications=bootstrap.CI.replications)
            
            one.rowname <- paste0("propensity_", one.propensity.method, "_rule_", one.rule.method)
            rownames(mat.direct.interactions)[row.number] <- one.rowname
            mat.direct.interactions[one.rowname, 
                                                         c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
                as.numeric(evaluate.direct.interactions.rule.on.validation[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
            
            list.BuildRules.direct.interactions[[row.number]] <- one.direct.interactions.rule
            names(list.BuildRules.direct.interactions)[[row.number]] <- one.rowname
            list.EvaluateRules.direct.interactions[[row.number]] <- evaluate.direct.interactions.rule.on.validation
            names(list.EvaluateRules.direct.interactions)[[row.number]] <- one.rowname
            row.number <- row.number + 1
        }
    }
    total.time.building.direct.interactions <- proc.time() - before.building.direct.interactions
    save(list.BuildRules.direct.interactions, list.EvaluateRules.direct.interactions,
         file=paste0("Results/", "list_direct.interactions.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
    save(total.time.building.direct.interactions,
         mat.direct.interactions, 
         file=paste0("Results/", "mat_direct.interactions.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
}

## RUN OWL.FRAMEWORK
total.combinations.OWL.framework <- length(vec.building.propensity.method) * length(vec.building.rule.method)
run.OWL.framework <- TRUE
if (run.OWL.framework == TRUE) {
    # initialize objects to store results
    list.BuildRules.OWL.framework <- vector("list", total.combinations.OWL.framework)
    list.EvaluateRules.OWL.framework <- vector("list", total.combinations.OWL.framework)
    mat.OWL.framework <- matrix(NA, nrow=total.combinations.OWL.framework, ncol=5)
    colnames(mat.OWL.framework) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
    rownames(mat.OWL.framework) <- rep(NA, total.combinations.OWL.framework)
    
    before.building.OWL.framework <- proc.time()
    row.number <- 1
    for (p in 1:length(vec.building.propensity.method)) {
        for (r in 1:length(vec.building.rule.method)) {
            print(row.number)
            one.propensity.method <- vec.building.propensity.method[p]
            one.rule.method <- vec.building.rule.method[r]
            
            # building rule on training set
            set.seed(row.number)
            one.OWL.framework.rule <- BuildRule(data=training.data,
                                                 study.design="observational",
                                                 prediction.approach="OWL.framework",
                                                 name.outcome=name.outcome,
                                                 type.outcome="binary",
                                                 desirable.outcome=TRUE,
                                                 name.treatment=name.treatment,
                                                 names.influencing.treatment=names.influencing.treatment,
                                                 names.influencing.rule=names.influencing.rule,
                                                 propensity.method=one.propensity.method,
                                                 additional.weights=training.data$IPW.CC,
                                                 rule.method=one.rule.method)
            
            ## evaluate on validation
            evaluate.OWL.framework.rule.on.validation <- EvaluateRule(data=validation.data,
                                                                       BuildRule.object=one.OWL.framework.rule,
                                                                       study.design="observational",
                                                                       name.outcome=name.outcome,
                                                                       type.outcome="binary",
                                                                       desirable.outcome=TRUE,
                                                                       clinical.threshold=0,
                                                                       name.treatment=name.treatment,
                                                                       names.influencing.treatment=names.influencing.treatment,
                                                                       names.influencing.rule=names.influencing.rule,
                                                                       propensity.method=one.propensity.method,
                                                                       additional.weights=validation.data$IPW.CC,
                                                                       bootstrap.CI=FALSE,
                                                                       bootstrap.CI.replications=bootstrap.CI.replications)
            
            one.rowname <- paste0("propensity_", one.propensity.method, "_rule_", one.rule.method)
            rownames(mat.OWL.framework)[row.number] <- one.rowname
            mat.OWL.framework[one.rowname, 
                              c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
                as.numeric(evaluate.OWL.framework.rule.on.validation[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
            
            list.BuildRules.OWL.framework[[row.number]] <- one.OWL.framework.rule
            names(list.BuildRules.OWL.framework)[[row.number]] <- one.rowname
            list.EvaluateRules.OWL.framework[[row.number]] <- evaluate.OWL.framework.rule.on.validation
            names(list.EvaluateRules.OWL.framework)[[row.number]] <- one.rowname
            row.number <- row.number + 1
        }
    }
    total.time.building.OWL.framework <- proc.time() - before.building.OWL.framework
    save(list.BuildRules.OWL.framework, list.EvaluateRules.OWL.framework,
         file=paste0("Results/", "list_OWL.framework.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
    save(total.time.building.OWL.framework,
         mat.OWL.framework, 
         file=paste0("Results/", "mat_OWL.framework.", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
}


## EVALUATE COMPETING RULE: TREAT EVERYONE
run.treat.everyone <- TRUE
total.combinations.treating.everyone <- length(vec.building.propensity.method)
if (run.treat.everyone == TRUE) {
    B.treat.everyone <- rep(1, nrow(validation.data))
    mat.treat.everyone <- matrix(NA, nrow=total.combinations.treating.everyone, ncol=5)
    colnames(mat.treat.everyone) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
    rownames(mat.treat.everyone) <- rep(NA, total.combinations.treating.everyone)
    list.EvaluateRules.treat.everyone <- vector("list", total.combinations.treating.everyone)
    
    before.building.treat.everyone <- proc.time()
    set.seed(123)
    row.number <- 1
    for (p in 1:length(vec.building.propensity.method)) {
        print(row.number)
        set.seed(row.number)
        one.propensity.method <- vec.building.propensity.method[p]
        evaluate.treat.everyone.rule.on.validation <- EvaluateRule(data=validation.data,
                                                                    BuildRule.object=NULL,
                                                                    B=B.treat.everyone,
                                                                    study.design="observational",
                                                                    name.outcome=name.outcome,
                                                                    type.outcome="binary",
                                                                    desirable.outcome=TRUE,
                                                                    separate.propensity.estimation=FALSE,
                                                                    clinical.threshold=0,
                                                                    name.treatment=name.treatment,
                                                                    names.influencing.treatment=names.influencing.treatment,
                                                                    names.influencing.rule=names.influencing.rule,
                                                                    propensity.method=one.propensity.method,
                                                                    additional.weights=validation.data$IPW.CC,
                                                                    bootstrap.CI=TRUE,
                                                                    bootstrap.CI.replications=bootstrap.CI.replications)
        one.rowname <- paste0("propensity_", one.propensity.method)
        rownames(mat.treat.everyone)[row.number] <- one.rowname
        mat.treat.everyone[one.rowname, 
                           c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
            as.numeric(evaluate.treat.everyone.rule.on.validation[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
        list.EvaluateRules.treat.everyone[[row.number]] <- evaluate.treat.everyone.rule.on.validation
        names(list.EvaluateRules.treat.everyone)[[row.number]] <- one.rowname
        row.number <- row.number + 1
    }
    total.time.building.treat.everyone <- proc.time() - before.building.treat.everyone
    save(list.EvaluateRules.treat.everyone,
         file=paste0("Results/", "list_treat_everyone_", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
    save(total.time.building.treat.everyone, mat.treat.everyone,
         file=paste0("Results/", "mat_treat_everyone_", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
} 

## EVALUATE COMPETING RULE: TREAT NO ONE
run.treat.noone <- TRUE
total.combinations.treating.noone <- length(vec.building.propensity.method)
if (run.treat.noone == TRUE) {
    B.treat.noone <- rep(0, nrow(validation.data))
    mat.treat.noone <- matrix(NA, nrow=total.combinations.treating.noone, ncol=5)
    colnames(mat.treat.noone) <- c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")
    rownames(mat.treat.noone) <- rep(NA, total.combinations.treating.noone)
    list.EvaluateRules.treat.noone <- vector("list", total.combinations.treating.noone)
    
    before.building.treat.noone <- proc.time()
    set.seed(123)
    row.number <- 1
    for (p in 1:length(vec.building.propensity.method)) {
        print(row.number)
        set.seed(row.number)
        one.propensity.method <- vec.building.propensity.method[p]
        evaluate.treat.noone.rule.on.validation <- EvaluateRule(data=validation.data,
                                                                 BuildRule.object=NULL,
                                                                 B=B.treat.noone,
                                                                 study.design="observational",
                                                                 name.outcome=name.outcome,
                                                                 type.outcome="binary",
                                                                 desirable.outcome=TRUE,
                                                                 separate.propensity.estimation=FALSE,
                                                                 clinical.threshold=0,
                                                                 name.treatment=name.treatment,
                                                                 names.influencing.treatment=names.influencing.treatment,
                                                                 names.influencing.rule=names.influencing.rule,
                                                                 propensity.method=one.propensity.method,
                                                                 additional.weights=validation.data$IPW.CC,
                                                                 bootstrap.CI=FALSE,
                                                                 bootstrap.CI.replications=bootstrap.CI.replications)
        one.rowname <- paste0("propensity_", one.propensity.method)
        rownames(mat.treat.noone)[row.number] <- one.rowname
        mat.treat.noone[one.rowname, 
                        c("Positives", "Negatives", "ATE in Positives", "ATE in Negatives", "ABR")] <- 
            as.numeric(evaluate.treat.noone.rule.on.validation[c("n.test.positives", "n.test.negatives", "ATE.test.positives", "ATE.test.negatives", "ABR")])
        list.EvaluateRules.treat.noone[[row.number]] <- evaluate.treat.noone.rule.on.validation
        names(list.EvaluateRules.treat.noone)[[row.number]] <- one.rowname
        row.number <- row.number + 1
    }
    total.time.building.treat.noone <- proc.time() - before.building.treat.noone
    save(list.EvaluateRules.treat.noone,
         file=paste0("Results/", "list_treat_noone_", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
    save(total.time.building.treat.noone, mat.treat.noone,
         file=paste0("Results/", "mat_treat_noone_", name.outcome, "_", name.treatment, "_", Sys.Date(), ".RData"))
}


