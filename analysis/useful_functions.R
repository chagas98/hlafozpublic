#These sets of functions was added to the analysis from different resources.
#Each one has its resource mentioned.

# addfiglab & line2user --------------------------------------------------------
#workaround to add figure labels with corrplot:
#https://stackoverflow.com/questions/38439211/figure-labels-add-text-on-graphs-in-the-same-location-despite-figure-size

## ----line2user ---------------------------------------------------------------
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
  y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
  switch(side,
         `1` = grconvertY(-line * y_off, 'npc', 'user'),
         `2` = grconvertX(-line * x_off, 'npc', 'user'),
         `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
         `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
         stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

## ----addfiglab ---------------------------------------------------------------
addfiglab <- function(lab, xl = par()$mar[2], yl = par()$mar[3]) {

  text(x = line2user(xl, 2), y = line2user(yl, 3),
       lab, xpd = NA, font = 1, cex = 1.5, adj = c(0, 1))
}


# ggcheck_the_qq, ggcheck_thespreadlevel & ggcheck_the_model -------------------

#Set of functions extract from Applied Statistics for Experimental Biology
#Copyright 2018 Jeffrey A. Walker
#Draft: 2021-09-25
#Link: https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/

## ----ggcheck_the_qq, warning = FALSE------------------------------------------
ggcheck_the_qq = function(m1,
                          line = "robust",
                          n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  #sigma_m1_res <- sigma(m1)

  normal_qq <- ppoints(n) %>%
    qnorm()
  sample_qq <- m1_res[order(m1_res)]

  # mean + sd
  parametric_slope <- sd(sample_qq)
  parametric_intercept <- mean(sample_qq)

  # quartiles
  m1_quartiles <- quantile(m1_res, c(.25, .75))
  qnorm_quartiles <- qnorm( c(.25, .75))
  m1_diff <- m1_quartiles[2] - m1_quartiles[1]
  qnorm_diff <- qnorm_quartiles[2] - qnorm_quartiles[1] # = 1.349
  quartile_slope <- m1_diff/qnorm_diff
  quartile_intercept <- median(m1_quartiles) # median of quartiles not quantiles

  # robust uses MASS:rlm (default arguments?)
  qq_observed <- data.table(normal_qq = normal_qq,
                            sample_qq = sample_qq)
  m2 <- rlm(sample_qq ~ normal_qq, data = qq_observed)
  robust_intercept <- coef(m2)[1]
  robust_slope <- coef(m2)[2]

  # re-sample ribbon
  set.seed(1)
  resample_qq_model <- numeric(n_boot*n)
  Y <- simulate(m1, n_boot)
  fd <- model.frame(m1) %>%
    data.table
  inc <- 1:n
  for(sim_i in 1:n_boot){
    # parametric bound
    fd[, (1) := Y[,sim_i]]
    m1_class <- class(m1)[1]
    if(m1_class == "lm"){
      ff <- lm(formula(m1), data = fd)
    }
    if(m1_class == "lmerModLmerTest" | m1_class == "lmerMod"){
      ff <- lmer(formula(m1), data = fd)
    }
    y_res <- residuals(ff)
    resample_qq <- y_res[order(y_res)]
    resample_qq_model[inc] <- resample_qq
    inc <- inc + n

    # robust bound
    qq_resampled <- data.table(normal_qq = normal_qq,
                               resample_qq = resample_qq)
    m2_resample <- rlm(resample_qq ~ normal_qq, data = qq_resampled)

  }

  qq_sim <- data.table(normal_qq = normal_qq,
                       resample_qq_model = resample_qq_model)

  qq_ci_model <- qq_sim[, .(median = median(resample_qq_model),
                            lower = quantile(resample_qq_model, 0.025),
                            upper = quantile(resample_qq_model, 0.975)),
                        by = normal_qq]
  m2_boot <- rlm(median ~ normal_qq, data = qq_ci_model)
  robust_intercept_boot <- coef(m2_boot)[1]
  robust_slope_boot <- coef(m2_boot)[2]

  ggplot(data = qq_observed,
         aes(x = normal_qq, y = sample_qq)) +

    # ribbon
    geom_ribbon(data = qq_ci_model,
                aes(ymin = lower,
                    ymax = upper,
                    y = median,
                    fill = "band"),
                fill = "gray",
                alpha = 0.6) +
    # draw points
    geom_point() +

    # robust
    geom_abline(aes(intercept = robust_intercept,
                    slope = robust_slope,
                    color = "robust"),
                show.legend = FALSE,
                size = 0.75) +
    # robust_boot
    # geom_abline(aes(intercept = robust_intercept_boot,
    #                 slope = robust_slope_boot,
    #                 color = "robust boot"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    xlab("Normal Quantiles") +
    ylab("Sample Quantiles") +

    scale_color_manual(values = gradient_grey[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL

}


## ----ggcheck_the_spreadlevel----------------------------------------------------
ggcheck_the_spreadlevel <- function(m1,
                                    n_boot = 200){
  n <- nobs(m1)
  m1_res <- residuals(m1)
  m1_scaled <- m1_res/sd(m1_res)
  m1_root <- sqrt(abs(m1_scaled))
  m1_fitted <- fitted(m1)

  m2 <- lm(m1_root ~ m1_fitted)
  m2_intercept <- coef(m2)[1]
  m2_slope <- coef(m2)[2]

  plot_data <- data.table(
    m1_res = sqrt(abs(m1_scaled)),
    fitted = m1_fitted
  )

  ggplot(data = plot_data,
         aes(x = fitted, y = m1_res)) +

    # ribbon
    # geom_ribbon(data = qq_ci_model,
    #             aes(ymin = lower,
    #                 ymax = upper,
    #                 y = median,
    #                 fill = "band"),
    #             fill = "gray",
    #             alpha = 0.6) +
    # # draw points
    geom_point() +

    geom_smooth(method = lm) +
    # robust
    # geom_abline(aes(intercept = robust_intercept,
    #                 slope = robust_slope,
    #                 color = "robust"),
    #             show.legend = TRUE,
    #             size = 0.75) +
    # robust_boot
    # geom_abline(aes(intercept = robust_intercept_boot,
    #                 slope = robust_slope_boot,
    #                 color = "robust boot"),
    #             show.legend = TRUE,
  #             size = 0.75) +
  xlab("Fitted") +
    ylab("root abs-scaled-residual") +

    scale_color_manual(values = gradient_grey[c(1:2,5:6)]) +
    theme_minimal_grid() +
    NULL
}


## ----ggcheck_the_model--------------------------------------------------------
ggcheck_the_model <- function(m1){
  gg1 <- ggcheck_the_qq(m1)
  gg2 <- ggcheck_the_spreadlevel(m1)
  cowplot::plot_grid(gg1, gg2, nrow = 1)
}


