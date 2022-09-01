################################################################################
## Robust inferences from the HLA project (Foz do Iguaçu - Brazil):           ##
## Correlation (between and within 2020-2021)                                 ##
## Logistic Regression (outcome ~ alleles + covariates)                       ##
## Negative Binomial Model (continuous admission variables ~ alleles + covar) ##
##                                                                            ##
## sc.assis.2017@aluno.unila.edu.br - Aug 2022                                ##
################################################################################

# Data preparation -------------------------------------------------------------

  ##Selecting important variables for statistical modelling and setting up the new
  ##names for plot

transform_data <- function(datafinal, sel_vars, type = NULL) {
  #Mutate variables with binary code - required for downstream analysis
  new_dataset <- datafinal %>%
    dplyr::select(all_of(names(sel_vars)), Num_Sample) %>%
    dplyr::mutate(
      Age = as.numeric(Age),
      Outcome_icu = ifelse(Outcome_icu == "Alta", 0, 1),
      smok_bv = case_when(
        smok_bv == "Não" ~ 0,
        smok_bv == "Ex" ~ 1,
        smok_bv == "Sim" ~ 2,
        smok_bv == "Ignorado" ~ 0),
      sex = case_when(sex == "Masculino" ~ 0,
                      sex == "Feminino" ~ 1),
      BacCoinfec = case_when(BacCoinfec == "Não" ~ 0,
                             BacCoinfec == "Sim" ~ 1)) %>%
    dplyr::mutate(
      allele1 = ifelse(allele1 %in% alleles_MAF$allele,
                              "binned", allele1),
      allele2 = ifelse(allele2 %in% alleles_MAF$allele,
                              "binned", allele2))

  new_dataset <- new_dataset %>%
    dplyr::mutate_at(vars(smok_bv, sum_comorb), as.factor)

  if (type == "corr") {
    new_dataset <- new_dataset %>%
      subset(select = -c(allele1, allele2, Num_Sample)) %>%
      mutate(SamplesGroupYear = ifelse(SamplesGroupYear == "2020", 0, 1))
  }

  if (type == "lreg") {
    new_dataset <- new_dataset %>%
      mutate(SamplesGroupYear = as.factor(SamplesGroupYear))
  }


  if(type == "glm") {
    new_dataset <- new_dataset %>%
      tidyr::gather(., "Name", "allele", allele1, allele2)
  }

  return(new_dataset)

}

  ## Check distribution of residuals and Residuals-Fitted values

modelchecking <- function(cm, modeltype){

  if (modeltype == "lm"){

    qqres <-  ggplot(cm, aes(sample = rstudent(cm))) +
      geom_qq() +
      stat_qq_line() +
      xlab("Theoretical Quantiles") +
      ylab("Studentized Residuals") +
      ggtitle("Q-Q Plot") +
      theme_bw()

    denres <- ggplot(cm) +
      geom_density(aes(x = residuals(cm))) +
      xlab("Residuals") +
      ylab("Density") +
      ggtitle("Density") +
      theme_bw()
    return(qqres+denres)
  }

  if (modeltype == "glm"){
    print("ok")
    cm_simres <- simulateResiduals(fittedModel = cm, n = 250)

    plot1 <- ggplotify::as.ggplot(~plot(cm_simres, asFactor = FALSE), envir=environment())

    cm_simres_refit <- simulateResiduals(fittedModel = cm,
                                         n = 250,
                                         refit = TRUE)
    plot2 <- ggplotify::as.ggplot(~testDispersion(cm_simres_refit), envir=environment())

    plot3 <- (plot1) / plot2 +  patchwork::plot_layout(nrow = 2, heights = c(2/5, 3/5))

    return(plot3)
  }
}


# Correlation ----------------------------------------------------------------

  # Function: corr_values()
  # Split the dataset if necessary (i.e for each year);
  # Apply correlation analysis with correlation(), setting method as "auto" to
  # selecting the most relevant method for different type of data(see ?correlation).

corr_values <- function (response = NULL, exclude = NULL, designmatrix){

  if (!is.null(response)) {
    designmatrix <- designmatrix %>%
      dplyr::select(-all_of(exclude))

    var_dm <- split(designmatrix[, -match(response, names(designmatrix))],
                    designmatrix[[response]])

    corr_results <- lapply(var_dm, function(x) correlation::correlation(x,
                                                                       method = "auto",
                                                                       robust = TRUE,
                                                                       partial = TRUE,
                                                                       redundant = TRUE,
                                                                       include_factors = TRUE,
                                                                       p_adjust = "bonferroni"))
  } else {

    var_dm <- designmatrix %>%
      dplyr::select(-all_of(exclude))

    corr_results <- correlation::correlation(var_dm,
                                             method = "auto",
                                             robust = TRUE,
                                             partial = TRUE,
                                             redundant = TRUE,
                                             include_factors = TRUE,
                                             p_adjust = "bonferroni")

  }

  return(corr_results)
}

  # Prepare data to plotting correlations from corr_values results
corr_plots <- function(corrvalues,sel_vars, colors){

  cor_matrixvalues <- summary(corrvalues, redundant = TRUE) %>%
    column_to_rownames("Parameter") %>%
    data.matrix(rownames.force = NA)

  newcol <- unname(unlist(sel_vars[rownames(cor_matrixvalues)]))
  newrow <- unname(unlist(sel_vars[colnames(cor_matrixvalues)]))

  names_dim <- list(newcol, newrow)

  dimnames(cor_matrixvalues) <- names_dim

  cor.p <-matrix(corrvalues$p, nrow = nrow(cor_matrixvalues), ncol = ncol(cor_matrixvalues),
                 dimnames = names_dim)

  #plot <- ggcorrplot::ggcorrplot(cor_matrixvalues, hc.order = TRUE, outline.color = "white", p.mat = cor.p, lab = TRUE)

  cor.p[lower.tri(cor.p, diag = T)] <- 1

  plot <- corrplot.mixed(cor_matrixvalues, order= "original",
                         mar=c(0,0,2,0), tl.col = 'black',
                         p.mat = cor.p, insig = "label_sig",
                         sig.level = c(.001, .01, .05), pch.cex=1,
                         tl.cex = .5, number.font=1.8, number.cex=.6,
                         upper.col = colors, tl.pos = "lt",
                         lower.col = colors)

  return(plot)
}

# Logistic Regression --------------------------------------------------------#

  # Logistic Regression Function
logiting <- function(matrix, DepVar, year = NULL, CoVars = NULL, model = NULL) {


  if(!is.null(year))  matrix <- matrix %>% filter(SamplesGroupYear == year)


  hla_class <- HIBAG::hlaAllele(matrix$Num_Sample,
                         matrix$allele1,
                         matrix$allele2,
                         locus = "B", na.rm = FALSE)

  formula <- as.formula(paste0(DepVar, " ~ h", ifelse(is.null(CoVars), "", " + "), CoVars))

  print(DepVar)

  logit_model <- HIBAG::hlaAssocTest(hla_class, formula, model = "additive",
                              data = matrix, showOR = TRUE)

  logit_model <- rownames_to_column(logit_model,"allele")

  return(logit_model)

  }

# Negative Binomial Regression ------------------------------------------------#
# calcular p valor de lm(), mais significantes realizar analise de contraste


glm_emm <- function(modelglm, dataglm, specs){

  model_emm <- emmeans::emmeans(modelglm,
                                specs = specs,
                                type="response")

  return(model_emm)
}

glm_effect_rat <- function(model_emm){

  fin_res <- data.frame()

  for (meth in c("eff", "pairwise")) {

  contrast_res <- emmeans::contrast(
    model_emm,
    method = meth,
    adjust = "none") %>%
    summary(infer = TRUE)


  fin_res <- rbind(fin_res, contrast_res)

  }

  return(fin_res)

}

plot_effects <- function(effect_rat, xlab, plim) {

  effect_rat$contrast <- gsub("[()]|effect|\\s", "", as.character(effect_rat$contrast))

  effect_rat <- effect_rat %>%
    dplyr::filter(p.value < plim) %>%
    dplyr::mutate(p_round = p_round(p.value,digits = 3),
                  p_pretty = p_format(p_round,
                                      digits = 3,
                                      accuracy = 1e-04,
                                      add.p = TRUE)
                  ) %>%
    data.frame()

  c <- effect_rat$contrast
  levelsc <- rev(c[order(nchar(c),c)])
  effect_rat$contrast <- factor(effect_rat$contrast, levels = levelsc)

  plot_eff <- ggplot(data = effect_rat, aes(x = ratio, y = contrast)) +

    geom_point() +
    xlab(xlab) +
    ylab("") +

  # add error bars
    geom_errorbar(aes(x = ratio,
                    xmin = asymp.LCL,
                    xmax = asymp.UCL),
                width = 0.02) +
  # add zero effect line
    geom_vline(xintercept = 1,
               linetype = "dashed") +
  # add p-values
    geom_text(aes(x = ratio,
                  y = contrast,
                  label = p_pretty),
              vjust = -0.5) +
    ggpubr::theme_pubr() +
    ggpubr::rremove("ylab")+
    scale_x_continuous(position="top") # move x-axis to "top"

  return(plot_eff)
  }

plot_response <- function(dataraw, data_emm, y_var, palette, ylabel, legendtitle) {

  data_emm <- data.frame(data_emm)

  plot_res <- ggplot(data = dataraw,
                     aes(x = allele,
                         y = get(y_var)
    )) +
  # add jittered points
  geom_jitter(aes(colour = Age),
              width = 0.3) +

  # add layer containing means
  geom_point(data = data_emm,
             aes(y = response),
             size = 2.5) +

  # add layer containing error bars
  geom_errorbar(data = data_emm,
                aes(y = response,
                    ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.05) +

  # change the title of the y-axis
  labs(y = ylabel,
       colour = legendtitle) +

  # change the theme
  theme_pubr() +
  # these theme modifications need to be added after re-setting
  # the theme
  theme(legend.position= "right",
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
  ) +

  scale_colour_distiller(palette = palette)
  # change the colors palette for the points
 # scale_color_brewer(brewer.pal(n = 12, name = "Purples"))

  return(plot_res)
}


plot_response_int <- function(dataraw, data_emm, y_var, ylabel, legendtitle) {

  data_emm <- data.frame(data_emm)

  plot_int <- ggstripchart(
    data = dataraw,
    x = "allele",
    y = "ICU_days",
    color = "SamplesGroupYear",
    fill = "SamplesGroupYear",
    #  position = position_dodge(width = dodge_width)
    position = position_jitterdodge(dodge.width = .4,
                                    jitter.width = .2)
    ) +
    rremove("xlab") +
    # add layer containing means
    geom_point(data = data_emm,
               aes(y = response,
                   fill = factor(SamplesGroupYear)),
               size = 2.5,
               position = position_dodge(width = .4)
               ) +
    # add layer containing error bars
    geom_errorbar(data = data_emm,
                  aes(y = response,
                      ymin = asymp.LCL,
                      ymax = asymp.UCL,
                      fill = factor(SamplesGroupYear)),
                  width = 0.05,
                  position = position_dodge(width = .4)
                  ) +
    theme_pubr() +
    labs(y = ylabel,
         color = legendtitle) +
    #newtheme after theme pubr
    theme(legend.position= "right",
          axis.title.x = element_blank(),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
    ) +
    scale_colour_manual(values = c(gradient_orange[6], gradient_purple[6])
                        ) +
    scale_fill_discrete(guide="none")

  return(plot_int)
}


#---------------------------USEFUL REFERENCES----------------------------------#
# https://towardsdatascience.com/get-a-grip-when-to-add-covariates-in-a-linear-regression-f6a5a47930e5\
# https://towardsdatascience.com/the-difference-between-correlation-and-regression-134a5b367f7c
# https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/
# http://www.sthda.com/english/articles/40-regression-analysis/167-simple-linear-regression-in-r/
# http://r-statistics.co/Linear-Regression.html
################################################################################

