# code chunks come from: https://github.com/galenos-project/LSR3_taar1_A/blob/main/util_u1/util.R 
### full code at: https://github.com/galenos-project/LSR3_taar1_A/tree/main 


forest_subgroup <- function(modelsumm, moderator, outcome, moderator_text) {
  # this uses GGplot2 to draw a forest plot for the subgroup analyses, and returns the plot 
  
  title <- paste0("Effect of TAAR1 Agonists on ",outcome, " by ", moderator_text)           
  
  model <- modelsumm
  colnames(model) <- c('moderator','k','SMD','se','p','ci_l','ci_u','symbol','size','summary','fontfaace','fontsize','d1','d2')
  model$order <- as.numeric(rownames(model))
  model$estimate_lab = paste0(sprintf('%.2f',model$SMD)," (", sprintf('%.2f',model$ci_l,2),",", sprintf('%.2f',model$ci_u,2),")")
  model <- model %>%
    arrange(order) %>%
    mutate(moderator = factor(model[["moderator"]], levels = unique(model[["moderator"]])))
  lnth <- nrow(model)+1
  
  axis_min <- min(floor(min(model$ci_l, model$ci_u)),-2)
  axis_max <- max(ceiling(max(model$ci_l, model$ci_u)),1)
  span2 <- 1 + (axis_max - axis_min)
  span1 <- span2 * 0.8
  span3 <- span2 * 0.5
  r1 <- span1
  l2 <- span1 + 1
  r2 <- span1 + span2 + 1
  l3 <- span1 + span2 + 2
  r3 <- span1 + span2 + span3 + 2
  
  cf <- span2/lnth
  
  top_margin <- 2
  
  poly1 <- subset(model, model$moderator == "Overall estimate")
  upp <- 1 + ((poly1$SMD - poly1$ci_l)/(cf *2))
  lop<- 1 - ((poly1$SMD - poly1$ci_l)/(cf *2))
  dfp <- data.frame(x = c(poly1$SMD, poly1$ci_u, poly1$SMD, poly1$ci_l), y = c(lop, 1, upp, 1))
  
  model <- model %>%
    arrange(order) %>%
    mutate(moderator = factor(moderator, levels = unique(moderator)))
  p_mid <- model %>%
    ggplot(aes(y = fct_rev(moderator))) +
    theme_classic() +
    geom_point(aes(x = SMD), shape = model$symbol, size = model$size) +
    geom_linerange(aes(xmin = ci_l, xmax = ci_u)) +
    labs(x = "SMD Effect size") +
    coord_cartesian(ylim = c(0, lnth + top_margin), xlim = c(axis_min-1, axis_max+1)) +
    geom_vline(xintercept = 0, linetype = "solid") +
    geom_vline(xintercept = poly1$SMD, linetype = "dashed") +
    annotate("text", x = axis_min-1, y = lnth + top_margin - 0.5, label = "TAAR1 Agonist\nworse", hjust = 0) +
    annotate("text", x = axis_max+1, y = lnth + top_margin - 0.5, label = "TAAR1 Agonist\nbetter", hjust = 1) +
    geom_polygon(data = dfp, aes(x = x, y = y), fill = "grey") +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank())
  
  p_left <-
    model %>%
    ggplot(aes(y = fct_rev(moderator))) +
    geom_text(aes(x = 0, label = moderator), hjust = 0, size = model$fontsize) +
    geom_text(aes(x = r1, label = k), hjust = 1, size = model$fontsize) +
    annotate("text", x = r1, y = lnth  + top_margin - 0.5, label = "Number of\nexperimental contrasts", hjust=1) +
    theme_void() +
    coord_cartesian(ylim = c(0, lnth + top_margin), xlim = c(0, span1))
  
  p_right <-
    model %>%
    ggplot() +
    geom_text(aes(x = span3, y = fct_rev(moderator), label = estimate_lab),size = model$fontsize, hjust = 1) +
    coord_cartesian(ylim = c(0, lnth + top_margin), xlim = c(0, span3)) +
    theme_void()
  
  layout <- c(
    area(t = 0, l = 0, b = 30, r = r1),
    area(t = 1, l = l2, b = 30, r = r2),
    area(t = 0, l = l3, b = 30, r = r3)
  )
  
  p_left + p_mid + p_right + plot_layout(design = layout) + 
    plot_annotation(
      title = title, 
      theme = theme(plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),  # Add margin below title
                    plot.margin = margin(t = 40, r = 10, b = 10, l = 10)
      )
    )
}


subgroup_analysis <- function(df, experiment_type, outcome, moderator, rho_value) {
  # this returns a table of effect sizes etc by moderator, for passing to 'forest_subgroup'
  # for plotting
  
  # Ensure the moderator is a character string for later conversion to symbol
  moderator <- as.character(moderator)
  
  df2 <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>%  
    filter(!is.na(SMDv)) %>%
    filter(!is.na(!!sym(moderator))) # Filter out NA values in moderator column
  
  # Convert character to factor if necessary
  if (is.character(df2[[moderator]])) {
    df2[[moderator]] <- factor(df2[[moderator]])}
  
  # List of factors to consider
  factors_to_consider <- c("Strain", "StudyId", "ExperimentID_I")
  
  # Create the random effects formula
  random_formula <- create_formula(factors_to_consider, df2)
  
  if (is.null(random_formula)) {
    cat("Insufficient levels for random effects grouping. Skipping meta-analysis.\n")
    return(NULL)
  }
  # Add a check for the number of levels in the moderator variable
  if (length(levels(df2[[moderator]])) >= 1) {
    #  message("In this iteration of the review, there was insufficient data to perform subgroup analysis for this variable (data for one subgroup only)")
    #  return(NULL)
    #}
    
    
    df2 <- df2 %>% mutate(effect_id = row_number()) # add effect_id column
    
    #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
    
    VCVM_SMD <- vcalc(vi = SMDv,
                      cluster = StudyId, 
                      subgroup= ExperimentID_I,
                      obs=effect_id,
                      data = df2, 
                      rho = rho_value) 
    
    # ML model on df2 with subgroup
    subgroup_analysis <- rma.mv(
      yi = SMD,
      V = VCVM_SMD,
      random = random_formula,
      data = df2,
      mods = as.formula(paste("~", moderator, "-1")),
      method = 'REML',
      test = "t",
      dfs = "contain"
    )
    
    #subgroup_analysis_predict <- predict(subgroup_analysis)
    
    ## ML model on df2 without subgroup
    overall_estimate_rma <- rma.mv(yi = SMD,
                                   V = VCVM_SMD,
                                   random = random_formula, # nested levels
                                   test = "t", # use t- and F-tests for making inferences
                                   data = df2,
                                   rho = rho_value,
                                   dfs="contain", # improve degree of freedom estimation for t- and F-distributions
                                   control=list(optimizer="nlminb"))
    
    #overall_estimate_rma_predict <- predict(overall_estimate_rma)
    
    
    k_subgroups <- df2 %>%
      group_by(df2[[moderator]]) %>%
      count() %>%
      pull(n)
    
    
    subgroup_analysis_plotdata <- data.frame(levels(df2[[moderator]]), k_subgroups, subgroup_analysis$beta, subgroup_analysis$se, subgroup_analysis$pval, subgroup_analysis$ci.lb, subgroup_analysis$ci.ub)
    colnames(subgroup_analysis_plotdata) <- c(moderator, "k", "SMD", "se","p", "ci_l", "ci_u") #, "pi.lb", "pi.ub")
    subgroup_analysis_plotdata$symbol <- 15
    subgroup_analysis_plotdata$size <- (1/subgroup_analysis_plotdata$se)
    subgroup_analysis_plotdata$summary <- FALSE
    subgroup_analysis_plotdata$fontfaace <- "plain"
    subgroup_analysis_plotdata$fontsize <- 3.88
    subgroup_analysis_plotdata=rbind(subgroup_analysis_plotdata, c("Overall estimate",  overall_estimate_rma$k, overall_estimate_rma$beta, overall_estimate_rma$se, overall_estimate_rma$pval, overall_estimate_rma$ci.lb,overall_estimate_rma$ci.ub, 18,1,TRUE,"bold",5)) #overall_estimate_rma_predict$pi.lb, overall_estimate_rma_predict$pi.ub))
    
    
    
    
    rownames(subgroup_analysis_plotdata) <- 1:nrow(subgroup_analysis_plotdata)
    subgroup_analysis_plotdata$k <- as.numeric(subgroup_analysis_plotdata$k)
    subgroup_analysis_plotdata$SMD <- as.numeric(subgroup_analysis_plotdata$SMD)
    subgroup_analysis_plotdata$se <- as.numeric(subgroup_analysis_plotdata$se)
    subgroup_analysis_plotdata$ci_l <- as.numeric(subgroup_analysis_plotdata$ci_l)
    subgroup_analysis_plotdata$ci_u <- as.numeric(subgroup_analysis_plotdata$ci_u)
    subgroup_analysis_plotdata$p <- as.numeric(subgroup_analysis_plotdata$p)
    subgroup_analysis_plotdata$symbol <- as.numeric(subgroup_analysis_plotdata$symbol)
    subgroup_analysis_plotdata$size <- as.numeric(subgroup_analysis_plotdata$size)
    subgroup_analysis_plotdata$fontsize <- as.numeric(subgroup_analysis_plotdata$fontsize)
    
    subgroup_analysis_plotdata$d1 <- (subgroup_analysis_plotdata$SMD - subgroup_analysis_plotdata$ci_l)/1.92
    subgroup_analysis_plotdata$d2 <- subgroup_analysis_plotdata$ci_u - subgroup_analysis_plotdata$SMD
    
    return(list(plotdata = subgroup_analysis_plotdata, 
                analysis = subgroup_analysis))
  }
}




create_formula <- function(factor_names, data) {
  distinct_levels <- sapply(factor_names, function(factor) n_distinct(data[[factor]]))
  
  if (all(distinct_levels < 5)) {
    # If all factors have fewer than 5 distinct levels, return NULL
    return(NULL)
  }
  
  selected_factors <- factor_names[distinct_levels >= 5]
  
  if (length(selected_factors) == 0) {
    # If none of the factors has enough levels, use a default grouping variable
    formula_str <- "~1"
  } else {
    formula_str <- paste("~ 1 |", paste(selected_factors, collapse = "/"))
  }
  
  formula_obj <- as.formula(formula_str)
  return(formula_obj)
}
