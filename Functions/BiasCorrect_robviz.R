# Original function is from R package robvis https://github.com/mcguinlu/robvis/tree/main
# Modified by Chin Yang Shapland on 19/09/24 to visualize estimates, RoB and corrected estimates

# TODO add argument that prevents subgroup summary estimates - in this case, each subgroup will need to be one row closer to each other
# TODO Finish description of function

#' Bias direction plots
#'
#' @description Used to
#'
#' @param dat Dataframe
#' @param vi Vector containing the sampling variances (normally defined as the column within the dataset, i.e. dat$vi). Note: either vi or sei must be set.
#' @param sei Vector containing the corresponding standard errors (normally defined as the column within the dataset, i.e. dat$sei). Note: either vi or sei must be set.
#' @param title Graph title
#' @param legend_cex Expansion factor for figure legend.
#' @param grouping Variable of the provided dataset by which the resulting plot will be stratified. Often will study design or overall risk-of-bias level.
#' @param grouping_levels Ordering of grouping variable. Note: the levels will be plotted
#'   in order, starting at the bottom of the graph (i.e. the last item in the
#'   vector will be placed at the top of the graph)
#' @param label_subgroup_summary Annotation text for subgroup label
#' @param atransf Transformation applied only when drawing estimates and
#'   confidence intervals. Use identity for SMDs; use exp for log-ratio measures.
#' @param external_subgroup_summaries Optional data frame with externally fitted
#'   original subgroup summaries. Expected columns: type, yi, ci.lb, ci.ub, label.
#' @param external_subgroup_summaries_adj Optional data frame with externally
#'   fitted bias-adjusted subgroup summaries. Same expected columns.
#' @param ... Other arguments to pass to metafor::forest
#'
#' @keywords internal

rob_direction <-
  function(dat,
           sei = NULL,
           title = NULL,
           legend_cex = 0.9,
           grouping = "type",
           grouping_levels = c("MR","NRSI","Obs","RCT"),
           label_subgroup_summary = "RE Model for Subgroup",
           atransf = identity,
           external_subgroup_summaries = NULL,
           external_subgroup_summaries_adj = NULL,
           ...) {
    
    ### calculate log risk ratios and corresponding sampling variances (and use
    ### the 'slab' argument to store study labels as part of the data frame)
    
    if (("study" %in% colnames(dat)) == FALSE) {
      dat$study <- paste("Study", 1:nrow(dat))
    }
    
    rob_levels = c("Low","Moderate","High","Critical")
    
    dat <- dat %>%
      dplyr::mutate(type = factor(type, levels = grouping_levels)) %>%
      dplyr::mutate(overall = factor(overall, levels = rob_levels)) %>%
      dplyr::arrange(type, overall, dplyr::desc(study))
    
    #dat[is.na(dat)] <- "None"
    
    
    # Use this to define the gaps between different groups
    # Will be important when adding argument to prevent subgroup analyses
    offset_n <- 4
    
    dat_rob_vec <- dat %>%
      dplyr::mutate(row_n = 1:dplyr::n()) %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(n=dplyr::n(),max = max(row_n), min = min(row_n)) %>%
      dplyr::mutate(offset = seq(2,length(unique(.$type))*offset_n,by=offset_n)) %>%
      dplyr::mutate(min = min+offset, max =max+offset, heading = max+1, stats = min-1.25) %>%
      dplyr::mutate(min = ifelse(n==1,min-1,min),
                    max = ifelse(n==1,max-1,max),
                    heading = ifelse(n==1,heading-1,heading))
    
    if (length(unique(dat$type))==1) {
      dat_rob_vec <- dat_rob_vec %>%
        dplyr::mutate(dplyr::across(c(min, max, heading),~.-1))
    }
    
    rows <- c()
    
    for (i in 1:nrow(dat_rob_vec)) {
      
      rows <-c(rows, dat_rob_vec$min[i]:dat_rob_vec$max[i])
      
    }
    
    arg <- list(...)
    
    # Use a symmetric SMD scale around the line of no effect. This keeps the
    # original and bias-corrected panels comparable and ensures 0 is printed on
    # the bottom axis. The caller can override this with x_min/x_max/at.
    data_ci_min <- min(
      dat$yi - 1.96 * sqrt(dat$vi),
      dat$yi_adj - 1.96 * sqrt(dat$vi_adj),
      na.rm = TRUE
    )
    data_ci_max <- max(
      dat$yi + 1.96 * sqrt(dat$vi),
      dat$yi_adj + 1.96 * sqrt(dat$vi_adj),
      na.rm = TRUE
    )
    if (!is.null(external_subgroup_summaries)) {
      data_ci_min <- min(data_ci_min, external_subgroup_summaries$ci.lb, na.rm = TRUE)
      data_ci_max <- max(data_ci_max, external_subgroup_summaries$ci.ub, na.rm = TRUE)
    }
    if (!is.null(external_subgroup_summaries_adj)) {
      data_ci_min <- min(data_ci_min, external_subgroup_summaries_adj$ci.lb, na.rm = TRUE)
      data_ci_max <- max(data_ci_max, external_subgroup_summaries_adj$ci.ub, na.rm = TRUE)
    }
    symmetric_limit <- ceiling(max(abs(c(data_ci_min, data_ci_max)), na.rm = TRUE))
    
    if (is.null(arg$x_min)) {
      x_min = -symmetric_limit
    } else {
      x_min <- arg$x_min
    }
    
    if (is.null(arg$x_max)) {
      x_max = symmetric_limit
    } else {
      x_max <- arg$x_max
    }
    
    forest_alim <- c(x_min, x_max)
    
    if (is.null(arg$at)) {
      axis_at <- seq(x_min, x_max, by = 1)
      if (!0 %in% axis_at) {
        axis_at <- sort(unique(c(axis_at, 0)))
      }
    } else {
      axis_at <- arg$at
    }
    
    # Separate the forest/CIs, printed estimate text, and risk-of-bias columns.
    # Without an explicit plotting limit (`alim` below), long CIs and subgroup
    # polygons can extend into the risk-of-bias grid because the x-axis also
    # includes the ROB columns.
    effect_text_x <- x_max + 1.2
    textpos <- c(x_min, effect_text_x)
    # Reserve vertical headroom above the "Studies" header. The plot title and
    # risk-of-bias header otherwise sit on the PDF device boundary and can be
    # clipped.
    y_max <- max(rows)+10
    title_y <- y_max - 0.8
    rob_header_y <- y_max - 3.2
    domain_header_y <- y_max - 4.2
    forest_header_y <- y_max - 6
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Deal with adding rob data
    dat_adj <- dat %>%
      dplyr::mutate(dplyr::across(-c(result_id,study,type, yi_adj, vi_adj), clean_data))
    
    dat <- dat %>%
      dplyr::mutate(dplyr::across(-c(result_id,study,type,yi,vi), clean_data))
    
    # Combine direction and type for each plotted domain. The package version
    # only loops over D1-D7; TAAR1 animal studies contain D1-D10.
    for (j in paste0("d",1:10)) {
      for (i in 1:nrow(dat)) {
        dat[i,paste0(j,"d")] <- paste0(dat[i,paste0(j,"d")],dat[i,paste0(j,"t")])
      }
    }
    
    # Put the ten domain columns to the right of the printed estimate text.
    # The PDF output is widened in taar1-agonists.R, so this can use real
    # horizontal space rather than squeezing the risk-of-bias grid into the CIs.
    x_pos <- seq(effect_text_x + 2.2, by = 0.38, length.out = 10)
    
    x_overall_pos <- max(x_pos) + .5
    
    # Convenience vector, specifying x-axis positions for all risk of bias columns
    header_row <- c(x_pos, x_overall_pos)
    
    legend_pos <- x_max+(max(header_row)-min(header_row))/2
    
    # New right-hand x-axis limit
    new_x_lim <- x_overall_pos + .8
    
    # Setting colours (changed)
    rob_colours <- c()
    rob_colours$na_colour <- "#cccccc"
    rob_colours$low_colour <- "#02C100"
    rob_colours$concerns_colour <- "#E2DF07"
    rob_colours$high_colour <- "#BF0000"
    rob_colours$critical_colour <- "#820000"
    rob_colours$ni_colour <- "#4EA1F7"
    
    judgements<-c("Very high risk",  #changed
                  "High risk",
                  "Moderate risk",
                  "Low risk")
    cols <- c(
      c = rob_colours$critical_colour, #changed
      h = rob_colours$high_colour,
      m = rob_colours$concerns_colour,
      l = rob_colours$low_colour,
      n = rob_colours$ni_colour,
      x = "transparent"
    )
    
    # Direction/type icons drawn inside each ROB square.
    # clean_data() reduces absolute direction and bias type to one letter each.
    # In the current triangulate output the direction keys are usually:
    # - a = away from null
    # - t = towards null
    # - u = unpredictable
    # Older package internals also used l/r for left/right absolute direction,
    # so those keys are retained as fallbacks.
    #
    # Proportional bias changes the estimate multiplicatively, so it is shown
    # with horizontal icons. Additive bias shifts the estimate, so it is shown
    # with vertical icons.
    syms <- c(ua = "?",
              up = "?",
              aa = "^",
              ta = "v",
              ap = ">",
              tp = "<",
              lp = "<",
              rp = ">",
              la = "v",
              ra = "^",
              l = "v",
              r = "^",
              xx = "",
              x = "")
    
    shapes <- c(c = 15,
                v = 15,
                h = 15,
                m = 15,
                l = 15,
                n = 15,
                x = 15)
    
    
    rob_psize = 3
    tsize <- rob_psize * 0.3
    
    # Use plain ASCII model labels so non-interactive PDF devices do not need
    # CID fonts for plotmath symbols such as tau.
    mlab_text <- function(text, res) {
      paste0(
        text,
        " (p ",
        .pval(res$pval, digits = 2, showeq = TRUE, sep = " "),
        "; I2 = ",
        formatC(res$I2, digits = 1, format = "f"),
        "%, tau2 = ",
        formatC(res$tau2, digits = 2, format = "f"),
        ")"
      )
    }
    
    external_label <- function(row, fallback) {
      if ("label" %in% names(row) && !is.na(row$label) && nzchar(row$label)) {
        return(row$label)
      }
      fallback
    }
    
    par(cex=0.9, mai=c(1.5,0.1,0.8,0.1))
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Make forest plot
    
    res_all <- rma(yi, vi, data=dat)
    
    par(fig=c(0.01,0.76,0,1))
    ### set up forest plot (with 2x2 table counts added; the 'rows' argument is
    ### used to specify in which rows the outcomes will be plotted)
    metafor::forest(res_all, 
                    #x = dat$yi,
                    #vi = dat$vi,
                    #sei = sei,
                    xlim=c(x_min, new_x_lim),
                    alim=forest_alim,
                    at=axis_at,
                    atransf=atransf,
                    slab = paste0(" ", dat$study),
                    cex=1.2,
                    ylim=c(-1.5, y_max),
                    rows=rows,
                    textpos = textpos,
                    mlab = "RE Model for All Studies",
                    header=FALSE,
                    ...
    )
    
    ### set font expansion factor (as in forest() above) and use a bold font
    op <- graphics::par(font=2)
    
    ### switch to bold italic font
    graphics::par(font=2)
    
    ### add text for the subgroups
    graphics::text(x_min, -2, pos=4, mlab_text(" ", res_all), cex = 1.2)
    graphics::text(x_min, forest_header_y, pos=4, "Studies", cex = 1.2)
    
    for (i in 1:nrow(dat_rob_vec)) {
      
      graphics::text(x_min, dat_rob_vec$heading[i], pos=4, dat_rob_vec$type[i], cex = 1.2)
    }
    
    ### set par back to the original settings
    graphics::par(op)
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Add risk of bias data
    
    headers <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "O")
    
    graphics::par(font = 2)
    # Need to add handling of top here
    graphics::text(mean(header_row), rob_header_y, labels = "Risk of Bias", cex=1.2)
    graphics::text(header_row, domain_header_y, labels = headers, cex=1.2)
    graphics::par(op)
    
    # Plot domain points
    for (j in 1:length(x_pos)) {
      graphics::points(
        x = rep(x_pos[j], length(rows)),
        y = rows,
        pch = shapes[dat[[paste0("d", j,"j")]]],
        col = scales::alpha(cols[dat[[paste0("d", j,"j")]]],0.6),
        cex = rob_psize
      )
      graphics::text(x_pos[j], rows, syms[dat[[paste0("d", j,"d")]]], cex = tsize)
    }
    
    graphics::points(
      rep(x_overall_pos, length(rows)),
      rows,
      pch = 15,
      col = scales::alpha(cols[dat[["overall"]]],0.6),
      cex = rob_psize
    )
    # graphics::text(x_overall_pos, rows, syms[dat[["overall"]]], cex = tsize)
    graphics::par(op)
    
    # #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    #
    
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Add sub-group, summary polygons & text
    
    rma_flexi <- function(x) {
      metafor::rma(
        yi,
        vi,
        subset = (type == x),
        data = dat,
        #method = "DL"  ### CHANGE to have "REML" the default rma
      )
    }
    
    
    res <- purrr::map(dat_rob_vec$type, rma_flexi)
    
    if (length(unique(dat$type))>1) {
      
      ### add summary polygons for the three subgroups
      for (i in 1:nrow(dat_rob_vec)) {
        type_i <- as.character(dat_rob_vec$type[i])
        ext_i <- NULL
        if (!is.null(external_subgroup_summaries)) {
          ext_i <- external_subgroup_summaries[external_subgroup_summaries$type == type_i, , drop = FALSE]
          if (nrow(ext_i) == 0) {
            ext_i <- NULL
          } else {
            ext_i <- ext_i[1, , drop = FALSE]
          }
        }
        
        if (!is.null(ext_i)) {
          metafor::addpoly(
            x = ext_i$yi,
            ci.lb = ext_i$ci.lb,
            ci.ub = ext_i$ci.ub,
            rows = dat_rob_vec$stats[i],
            cex = 1.2,
            textpos = textpos,
            atransf = atransf,
            annotate = FALSE,
            # The external model label is drawn explicitly below the polygon.
            # Keep addpoly()'s mlab empty so the same label is not printed twice.
            mlab = ""
          )
          
          graphics::par(font = 2)
          graphics::text(x_min, dat_rob_vec$stats[i]-1, pos=4, external_label(ext_i, "External subgroup model"), cex = 1.2)
          graphics::par(op)
          
          annotate_poly(ext_i$yi,
                        ext_i$ci.lb,
                        ext_i$ci.ub,
                        atransf = atransf,
                        textpos = textpos,
                        rows = dat_rob_vec$stats[[i]])
          next
        }
        
        if (length(res[[i]]$slab) == 1) {
          next
        }
        
        metafor::addpoly(
          res[[i]],
          #fonts = 1,
          row = dat_rob_vec$stats[i],
          cex = 1.2,
          textpos=textpos,
          atransf = atransf,
          annotate = F,
          mlab = label_subgroup_summary
        )
        
        graphics::par(font = 2)
        # Need to add handling of top here
        graphics::text(x_min, dat_rob_vec$stats[i]-1, pos=4, mlab_text(" ", res[[i]]), cex = 1.2)
        graphics::par(op)
        
        annotate_poly(res[[i]]$b,
                      res[[i]]$ci.lb,
                      res[[i]]$ci.ub,
                      atransf = atransf,
                      textpos = textpos,
                      rows = dat_rob_vec$stats[[i]])
        
      }
    }
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    
    if(!is.null(title)){
      graphics::par(font = 2)
      graphics::text(x_min, title_y, pos=4, bquote(bold(underline(.(title)))), cex = 1.2)
      graphics::par(op)
    }
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    
    graphics::legend(
      legend_pos-3.5,                         #changed
      -4,                                   #changed
      c(judgements),
      pch = c(15,15,15,15,16,50),           #changed
      xjust = 0.5,
      col = c(cols[1:4],"white","white"),   #changed
      xpd = TRUE,
      title = parse(text = "bold(\"Extent of bias\")"),
      title.adj = 0.05,
      cex = legend_cex,
      pt.cex = legend_cex-.1,
      y.intersp = 0.9,                     #changed
      box.col = "white",
    )
    
    graphics::legend(
      legend_pos-1.5,                    #changed
      -4,                                #changed    
      c("v   ^   Additive","<   >   Proportional", "?       Unpredictable"),
      xjust = 0.5,
      xpd = TRUE,
      adj = 0.15,
      title = parse(text = "bold(\"Type of bias\")"),
      title.adj = 0.05,
      cex = legend_cex,
      pt.cex = legend_cex,
      y.intersp = 0.9,
      box.col = "white"
    )
    #
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Make forest plot with adjusted estimates
    
    par(fig=c(0.72,1,0,1), new=TRUE)
    
    res_adj <- rma(yi_adj, vi_adj, data=dat_adj)
    
    ### set up forest plot (with 2x2 table counts added; the 'rows' argument is
    ### used to specify in which rows the outcomes will be plotted)
    metafor::forest(res_adj, 
                    #x = dat$yi,
                    #vi = dat$vi,
                    #sei = sei,
                    xlim=c(x_min, x_max + 2),
                    alim=forest_alim,
                    at=axis_at,
                    atransf=atransf,
                    slab = NA,
                    cex=1.2,
                    ylim=c(-1.5, y_max),
                    rows=rows,
                    textpos = c(x_min, x_max + 1.2),
                    mlab = NA,
                    header=c(" ", "Bias-corrected Est. [95% CI]"),
                    ...
    )
    
    ### set font expansion factor (as in forest() above) and use a bold font
    op <- graphics::par(font=2)
    
    graphics::text(x_min-.4, -2, pos=4, mlab_text(" ", res_adj), cex = 1.2)
    
    ### set par back to the original settings
    graphics::par(op)
    
    
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
    # Add sub-group, summary polygons & text
    
    rma_flexi_adj <- function(x) {
      metafor::rma(
        yi_adj,
        vi_adj,
        subset = (type == x),
        data = dat_adj,
        #method = "DL"  ### CHANGE to have "REML" the default rma
      )
    }
    
    
    res_adj <- purrr::map(dat_rob_vec$type,  rma_flexi_adj)
    
    if (length(unique(dat$type))>1) {
      
      ### add summary polygons for the three subgroups
      for (i in 1:nrow(dat_rob_vec)) {
        type_i <- as.character(dat_rob_vec$type[i])
        ext_i <- NULL
        if (!is.null(external_subgroup_summaries_adj)) {
          ext_i <- external_subgroup_summaries_adj[external_subgroup_summaries_adj$type == type_i, , drop = FALSE]
          if (nrow(ext_i) == 0) {
            ext_i <- NULL
          } else {
            ext_i <- ext_i[1, , drop = FALSE]
          }
        }
        
        if (!is.null(ext_i)) {
          metafor::addpoly(
            x = ext_i$yi,
            ci.lb = ext_i$ci.lb,
            ci.ub = ext_i$ci.ub,
            rows = dat_rob_vec$stats[i],
            cex = 1.2,
            textpos = textpos,
            atransf = atransf,
            annotate = FALSE,
            mlab = ""
          )
          
          graphics::par(font = 2)
          # Center the adjusted model label on the 0 SMD reference line so it
          # stays inside the bias-corrected forest panel.
          graphics::text(0, dat_rob_vec$stats[i]-1, external_label(ext_i, "External adjusted subgroup model"), cex = 1.2)
          graphics::par(op)
          
          annotate_poly(ext_i$yi,
                        ext_i$ci.lb,
                        ext_i$ci.ub,
                        atransf = atransf,
                        textpos = c(x_min, x_max + 1.2),
                        rows = dat_rob_vec$stats[[i]])
          next
        }
        
        if (length(res_adj[[i]]$slab) == 1) {
          next
        }
        
        metafor::addpoly(
          res_adj[[i]],
          row = dat_rob_vec$stats[i],
          cex = 1.2,
          textpos=textpos,
          atransf = atransf,
          annotate = F,
          mlab = ""
        )
        
        graphics::par(font = 2)
        # Need to add handling of top here
        graphics::text(x_min-.4, dat_rob_vec$stats[i]-1, pos=4, mlab_text(" ", res_adj[[i]]), cex = 1.2)
        graphics::par(op)
        
        annotate_poly(res_adj[[i]]$b,
                      res_adj[[i]]$ci.lb,
                      res_adj[[i]]$ci.ub,
                      atransf = atransf,
                      textpos = c(x_min, x_max + 1.2),
                      rows = dat_rob_vec$stats[[i]])
        
      }
    }   
}
