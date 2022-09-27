################################################################################
###########################  DATASETS FREQUENCIES  #############################
################################################################################

#Collect HLA frequencies from Paraná State with immunotation package
#from DataCleaningMutate
Redome_freq <- control_freq_redome %>%
  dplyr::rename(Freq = "allele_frequency"
  ) %>%
  dplyr::mutate(n_Redome = round(
    as.numeric(Freq) * (2*(as.numeric(gsub(",","",sample_size)))),0)
  ) %>%
  dplyr::select("allele", "Freq", "n_Redome")

#Get two alleles (genotypes) from each individual
allele_freq <- data_final %>%
  dplyr::select(allele1, allele2,
                Outcome_icu,
                SamplesGroupYear,
                SamplesGroup
  ) %>%
  dplyr::mutate(Outcome_icu = ifelse(
    Outcome_icu == "Óbito", 'Obito',  "Alta")
  ) %>%
  tidyr::gather(., "Name", "allele",1:2)

#Allele freq. between years
allele_freqyears <- frequencer(data = allele_freq,
                               var = "SamplesGroupYear",
                               redome = Redome_freq)

#Allele freq. between outcomes in both years
allele_freqgrouped <- frequencer(data = allele_freq,
                                 var = "SamplesGroup",
                                 redome = Redome_freq)

#Allele freq. between outcomes
allele_freqoutcome<- frequencer(data = allele_freq,
                                var = "Outcome_icu",
                                redome = Redome_freq)

# Relative Frequency between years
stat_alleleyears <- stats(data = data_final,
                          var = "SamplesGroupYear",
                          redome = Redome_freq)

# Relative Frequency between outcomes
stat_alleleoutcome <- stats(data = data_final,
                            var = "Outcome_icu",
                            redome = Redome_freq)

# Relative Frequency between outcomes in both years
stat_allelegrouped<- stats(data = data_final,
                           var = "SamplesGroup",
                           redome = Redome_freq)

#Minor Alleles Frequency
alleles_MAF <- stat_alleleyears %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  filter(Freq < 0.05 | n_2021 < 5 | n_2020 < 5)

##################  Summary Table 1 -  Clinical Variable  ######################

##Setup
gtsummary::set_gtsummary_theme(
  theme_gtsummary_compact(set_theme = TRUE))
#https://cran.r-project.org/web/packages/gtsummary/vignettes/themes.html

  ##Template - Design Summary Table
template_table <- function(data, div, translation) {
  tmp <- data %>%
    gtsummary::tbl_summary(by = all_of(div), # stratification variable.
                           include = c(Age, age_bin, sex,
                                       selfrep_race, categorical_comorb,
                                       smok_bv, alcoh_bv, BacCoinfec, HospitalPeriod_days,
                                       ICU_days, VentilatorySupport_days,
                                       FirstSymptons_days),
                           percent = "column",
                           statistic = list(
                             all_continuous() ~ "{mean} ({sd})",
                             all_categorical() ~ "{n} ({p}%)",
                             Age ~ "{median} ({p25}, {p75})",
                             VentilatorySupport_days ~ "{median} ({p25}, {p75})",
                             HospitalPeriod_days ~ "{median} ({p25}, {p75})",
                             ICU_days ~ "{median} ({p25}, {p75})"
                           ),
                           digits = list(all_continuous() ~ c(1, 1),
                                         all_categorical() ~ 1),
                           label = translation) %>%
    bold_labels() %>%
    italicize_levels()

  return(tmp)

}

  ##Translation
lista_labels = list(
  Age ~ "Idade (em anos)",
  age_bin ~ "Faixa-etária (em anos)",
  sex ~ "Sex",
  selfrep_race ~ "Raça Autodeclarada",
  categorical_comorb ~ "Número de Comorbidades",
  smok_bv ~ "Tabagismo",
  alcoh_bv ~ "Alcoolismo",
  BacCoinfec ~ "Coinfec. Bacteriana",
  HospitalPeriod_days ~ "Dias de Internação",
  ICU_days ~ "Dias de UTI",
  VentilatorySupport_days ~ "Dias de Suporte Ventilatório",
  FirstSymptons_days ~ "Período Primeiros Sintomas")

  ##Collect between outcomes
tab_stats <- template_table(data = data_final,
                            div = "SamplesGroup",
                            translation = lista_labels)

  ##Collect between years
tab_statsYears <- template_table(data = data_final,
                                 div = "SamplesGroupYear",
                                 translation = lista_labels)

##Between outcomes in 2020
pval20 <- data_final %>%
  dplyr::filter(SamplesGroupYear %in% "2020") %>%
  template_table(., "Outcome_icu", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**Alta20 vs. Obito20**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Between outcomes in 2021
pval21 <- data_final %>%
  dplyr::filter(SamplesGroupYear %in% "2021") %>%
  template_table(., "Outcome_icu", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**Alta21 vs. Obito21**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Between incidence in years
pvalYears <- data_final %>%
  template_table(., "SamplesGroupYear", lista_labels) %>%
  gtsummary::add_p(Age ~ "wilcox.test") %>%
  gtsummary::modify_header(p.value ~ "**2020 vs. 2021**") %>%
  gtsummary::modify_column_hide(all_stat_cols())

  ##Gather all tables
summary_table_final <-
  gtsummary::tbl_merge(list(tab_stats, tab_statsYears,
                            pval20, pval21, pvalYears)) %>%
  gtsummary::modify_header(label = "") %>%
  gtsummary::modify_spanning_header(
    list(
      gtsummary::all_stat_cols() ~ "**Desfechos**",
      gtsummary::starts_with("p.value") ~ "**p-valores**"
    )
  ) %>%
  gtsummary::as_flex_table() %>%
  flextable::set_table_properties(width = 1, layout = "autofit")


############################ FIGURE 1 ##########################################

count <- data_final %>%
  dplyr::select(allele1, allele2, SamplesGroupYear, Outcome_icu) %>%
  tidyr::gather(., "Name", "allele",1:2) %>%
  dplyr::mutate(faceet_color =
                  case_when(
                    SamplesGroupYear == "2020" ~ gradient_orange[6],
                    SamplesGroupYear == "2021" ~ gradient_purple[8]
                  ))

p_count <- ggplot(count, aes(x = allele, fill = as.character(SamplesGroupYear))) +
  geom_bar(width = 0.8, position = position_dodge(0.8, preserve = "single"),
           colour = "black") +
  geom_hline(aes(yintercept = 5), linetype = "dashed", size = 0.6) +
  labs(x = "",
       y = "Frequência (n)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(gradient_orange[6], gradient_purple[8])) +
  theme_pubr() +
  ggtitle('A')

p_count2 <- ggplot(count, aes(x = allele, fill = Outcome_icu)) +
  geom_bar(width = 0.8, position = position_dodge(0.8, preserve = "single"),
           colour = "black") +
  labs(x = "",
       y = "Frequência (n)") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(gradient_orange[6], gradient_purple[8])) +
  facet_wrap(~SamplesGroupYear, scales = "free", ncol = 1) +
  theme_pubr() +
  theme(
    legend.position= "right",
    legend.direction="vertical",
    legend.title = element_blank(),
    #legend.box.background = element_rect(size = 0.2), # Box for legend
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7))+
  ggtitle('B')

sortbold <- sort(unique(p_count$data$allele))

bold.labels <- ifelse(sortbold %in% alleles_MAF$allele,
                      "plain", "bold")

theme_count <-theme(
  legend.position= "right",
  legend.direction="vertical",
  legend.title = element_blank(),
  #legend.box.background = element_rect(size = 0.2), # Box for legend
  legend.text = element_text(size = 8),
  axis.text.x = element_text(size = 8, face = bold.labels),
  axis.text.y = element_text(size = 7))

p_count <- p_count +
  theme_count

fig1 <- (p_count / p_count2)

#Save Plot
png("images/fig1_count.png", width = 10, height = 7, units = 'in', res = 300)
print(fig1)
while (!is.null(dev.list()))  dev.off()

#############################  FIGURE S1  #######################################

bold.labels1 <- ifelse(allele_freqoutcome$allele %in% alleles_MAF$allele,
                       "plain", "bold")

#Middle
g.mid<-ggplot(allele_freqoutcome, aes(x=1,y=allele)) +
  geom_text(aes(label=allele), size = 4, color = "black",
            fontface = bold.labels1) +
  geom_segment(aes(x=0.95,
                   xend=0.96,
                   yend=allele)) +
  geom_segment(aes(x=1.04,
                   xend=1.05,
                   yend=allele)) +
  ggtitle("") +
  ylab(NULL) +
  scale_x_continuous(expand=c(0,0),
                     limits=c(0.94,1.065)
                     )+
  theme(axis.title.x = element_text(size = 11),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(0,-1,1,-1), "mm"))+
  xlab("Frequência \n Relativa (%)")

#Pattern Theme
pattern_theme <-
  theme_bw()+
  theme(panel.grid.minor = element_blank(), # No minor grid lines
        panel.grid.major.x = element_blank(), # No major x-axis grid lines
        panel.grid.major.y = element_blank(),
        legend.position = c(0.55, 0.95),
        legend.direction="vertical",
        legend.title = element_blank(),
        #legend.box.background = element_rect(size = 0.2), # Box for legend
        legend.key.size = unit(4, unit = 'mm'),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

#Left1
allele_freq_graph1 <- allele_freqyears %>%
  ggplot(aes(y = Freq, x = allele,
             fill = factor(SamplesGroupYear))
         ) +
  geom_bar(width = 0.8,
           position = position_dodge(0.8, preserve = "single"),
           stat = "identity",colour = "black"
           ) +
  geom_hline(aes(yintercept = 0.05),
             linetype = "dashed",
             size = 0.6) + # Dashed line at y = 1
  scale_fill_manual(values=c(gradient_orange[6],
                             gradient_orange[8],
                             gradient_grey[4])) +
  pattern_theme +
  theme(plot.margin = unit(c(1,1.1,1,0), "mm")) +
  scale_y_reverse(expand = c(0,0),limits = c(0.2,0)) +
  coord_flip() +
  ggtitle('A')

#Left2
allele_freq_graph2 <- allele_freqoutcome %>%
  ggplot(aes(y = Freq, x = allele,
             fill = factor(Outcome_icu))) +
  geom_bar(width = 0.8, position = position_dodge(0.8, preserve = "single"),
           stat = "identity",colour = "black")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,0.2)) +
  scale_fill_manual(values=c(gradient_purple[6],
                             gradient_purple[8],
                             gradient_grey[4]))+
  pattern_theme +
  theme(plot.margin = unit(c(1,1,1,-1), "mm"))+
  coord_flip()

#Left
plotfreq1 <- (allele_freq_graph1 + g.mid + allele_freq_graph2) +
  patchwork::plot_layout(ncol = 3, widths = c(3/9,3/9,3/9))

#Right1
allele_freq_graph3 <- allele_freqgrouped %>%
  filter(SamplesGroup %in% c("Alta20", "Óbito20", "Redome")) %>%
  ggplot(aes(y = Freq, x = allele, fill = factor(SamplesGroup))) +
  geom_bar(width = 0.8, position = position_dodge(0.8, preserve = "single"),
           stat = "identity", colour = "black") +
  scale_fill_manual(values=c(gradient_orange[6],
                             gradient_orange[8],
                             gradient_grey[4]))+
  pattern_theme +
  theme(plot.margin = unit(c(1,1.1,1,0), "mm")) +
  scale_y_reverse(expand = c(0,0),limits = c(0.3,0)) +
  coord_flip() +
  ggtitle('B')

#Right2
allele_freq_graph4 <- allele_freqgrouped %>%
  filter(SamplesGroup %in% c("Alta21", "Óbito21", "Redome")) %>%
  ggplot(aes(y = Freq, x = allele, fill = factor(SamplesGroup))) +
  geom_bar(width = 0.8, position = position_dodge(0.8, preserve = "single"),
           stat = "identity",colour = "black")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,0.3)) +
  scale_fill_manual(values=c(gradient_purple[6],
                             gradient_purple[8],
                             gradient_grey[4]))+
  pattern_theme +
  theme(plot.margin = unit(c(1,1,1,-1), "mm"))+
  coord_flip()

#Right
plotfreq2 <- (allele_freq_graph3 + g.mid + allele_freq_graph4) +
  patchwork::plot_layout(ncol = 3, widths = c(3/9,3/9,3.2/9))


fig2 <- (plotfreq1 | plot_spacer() | plotfreq2)  +
  plot_layout(widths = c(4, 0.02 ,4.2))

#Save Plot
png("images/figs1_frequencies.png", width = 7, height = 10, units = 'in', res = 300)
print(fig2)
while (!is.null(dev.list()))  dev.off()

############################### TABLE 2 ########################################

#Bind frequencies (n) for each allele
stat_all <-  cbind(stat_alleleyears["allele"],
                   stat_alleleoutcome["n_Redome"],
                   stat_alleleyears["n_2020"],
                   stat_alleleyears["n_2021"],
                   stat_allelegrouped["n_Alta20"],
                   stat_allelegrouped["n_Alta21"],
                   stat_allelegrouped["n_Óbito20"],
                   stat_allelegrouped["n_Óbito21"],
                   stat_alleleoutcome["n_Alta"],
                   stat_alleleoutcome["n_Óbito"])

# Get fisher row wise comparison

fisher_res <- rbind(fisher_tester(table = stat_all, groups = c("Alta20", "Óbito20"), plim = 1),
                    fisher_tester(table = stat_all, groups = c("Alta21", "Óbito21"), plim = 1),
                    fisher_tester(table = stat_all, groups = c("2021", "2020"), plim = 1),
                    fisher_tester(table = stat_all, groups = c("Óbito20", "Óbito21"), plim = 1)) %>%
  as.data.frame()

#Write Results

write.xlsx(fisher_res,
           file="results/fishertest.xlsx",
           sheetName="FisherResults", row.names=FALSE) #Fisher Results

write.xlsx(stat_all,
           file="results/allele_freq.xlsx",
           sheetName="Frequencies", row.names=FALSE) #Fisher Results

#############################  FIGURE 2  #######################################

 ##Select variables to compare
selvar <-
  list(Age = "Idade",
       sex = "Sexo",
       smok_bv = "Tabag." ,
       BacCoinfec = "Coinf. \n Bac.",
       HospitalPeriod_days = "Dias \n de \n Internação",
       ICU_days = "Dias \n de UTI",
       VentilatorySupport_days = "Suporte. \n Vent.",
       FirstSymptons_days = "Primeiros \n Sintomas",
       SamplesGroupYear = "Ano \n de Coleta",
       Outcome_icu = "Desfecho",
       sum_comorb = "Núm. de \n Comorb.",
       allele1 = "allele1",
       allele2 = "allele2")

selvar_fac <-
  list(Age = "Idade",
       sex = "Fem.",
       smok_bv.0 = "Não \n Tabag." ,
       smok_bv.1 = "Ex \n Tabag." ,
       smok_bv.2 = "Tabag.",
       BacCoinfec = "Coinf. \n Bac.",
       HospitalPeriod_days = "Dias \n de \n Internação",
       ICU_days = "Dias \n de UTI",
       VentilatorySupport_days = "Suporte. \n Vent.",
       FirstSymptons_days = "Primeiros \n Sintomas",
       SamplesGroupYear = "2021",
       Outcome_icu = "Óbito",
       sum_comorb = "Comorb.",
       allele1 = "allele1",
       allele2 = "allele2")

  ##Datasets for correlation
cor_dataset <- transform_data(datafinal = data_final,
                              sel_vars = selvar,
                              type = "corr")


fig2_corr_ab <- corr_values("SamplesGroupYear",
                       designmatrix = cor_dataset)

fig2_corr_c <- corr_values(designmatrix = cor_dataset)

 ##Save Plot
png("images/fig2_correlation.png", width = 10, height = 3, units = 'in', res = 300)
par(mfrow = c(1,3), xpd = NA)
corr_plots(fig2_corr_ab$`0`, selvar_fac, gradient_orange)
addfiglab("A")
corr_plots(fig2_corr_ab$`1`, selvar_fac, gradient_purple)
addfiglab("B")
corr_plots(fig2_corr_c, selvar_fac, gradiente_diverg)
addfiglab("C")
while (!is.null(dev.list()))  dev.off()

################################################################################
###########################  DATASETS MODELS ###################################
################################################################################


#Data for Logistic Regression
logit_dataset <- transform_data(data_final, selvar, type = "lreg")


#Datasets for Negative Binomial Regression
glm_dataset <- transform_data(data_final, selvar,type = "glm")


#############################  TABLE 4  ########################################

#outcome/years
logit20 <- logiting(logit_dataset, year = 2020, "Outcome_icu", model = "additive")
logit21 <- logiting(logit_dataset, year = 2021, "Outcome_icu", model = "additive")
logit20adj <- logiting(logit_dataset, year = 2020, "Outcome_icu", CoVars = "Age", model = "additive")
logit21adj <- logiting(logit_dataset, year = 2021, "Outcome_icu", CoVars = "Age", model = "additive")

#between years
logityear <- logiting(logit_dataset, "SamplesGroupYear", model = "additive")
logityearadj <- logiting(logit_dataset, "SamplesGroupYear", CoVars = "Age", model = "additive")


write.xlsx(logit20, file="results/logisticregression.xlsx",
           sheetName="outcome2020", row.names=FALSE)
write.xlsx(logit21, file="results/logisticregression.xlsx",
           sheetName="outcome2021", append=TRUE, row.names=FALSE)
write.xlsx(logit20adj, file="results/logisticregression.xlsx",
           sheetName="outcome2020adj", append=TRUE, row.names=FALSE)
write.xlsx(logit21adj, file="results/logisticregression.xlsx",
           sheetName="outcome2021adj", append=TRUE, row.names=FALSE)
write.xlsx(logityear, file="results/logisticregression.xlsx",
           sheetName="2020-2021", append=TRUE, row.names=FALSE)
write.xlsx(logityearadj, file="results/logisticregression.xlsx",
           sheetName="2020-2021adj", append=TRUE, row.names=FALSE)

#############################  FIGURE 4-5  #####################################
#interpret: file:///home/schagas/Downloads/NegativeBinomialRegressionGuide.pdf
## Split 2020 & 2021

glm_dataset0 <- glm_dataset[glm_dataset$SamplesGroupYear == "2020",]

glm_dataset1 <- glm_dataset[glm_dataset$SamplesGroupYear == "2021",]


#MODELS
glmmod <-  MASS:::glm.nb(ICU_days ~ allele*SamplesGroupYear, data = glm_dataset)

#2020
glmmod0 <-  MASS:::glm.nb(ICU_days ~ allele + Age, data = glm_dataset0)

#2021
glmmod1 <-  MASS:::glm.nb(ICU_days ~ allele + Age, data = glm_dataset1)


glmmodint0 <-  MASS:::glm.nb(ICU_days ~ allele*Outcome_icu, data = glm_dataset0)
glmmodint1 <-  MASS:::glm.nb(ICU_days ~ allele*Outcome_icu, data = glm_dataset1)

h <- MASS:::glm.nb(ICU_days ~ allele*SamplesGroupYear, data = glm_dataset)

#Responses
response0 <- glm_emm(glmmod0,
                     glm_dataset0,
                     "allele")
response1 <- glm_emm(glmmod1,
                     glm_dataset1,
                     "allele")

responseint0 <- glm_emm(glmmodint0,
                       glm_dataset0,
                       c("allele", "Outcome_icu"))
responseint1 <- glm_emm(glmmodint1,
                        glm_dataset1,
                        c("allele", "Outcome_icu"))



#Effects Ratio
effect_rat0 <- glm_effect_rat(response0)
effect_rat1 <- glm_effect_rat(response1)
effect_ratint0 <- glm_effect_rat(responseint0)
effect_ratint1 <- glm_effect_rat(responseint1)


#Plot Effects
plot_eff0 <- plot_effects(effect_rat0, "Razão de Efeito", 0.05) + #IRR is exp(coef(B))
  ggtitle("A")

plot_eff1 <- plot_effects(effect_rat1, "Razão de Efeito", 0.05) +
  ggtitle("B")

plot_eff3 <- plot_effects(effect_ratint0, "Razão de Efeito", 0.05) +
  ggtitle("A")

plot_eff4 <- plot_effects(effect_ratint1, "Razão de Efeito", 0.05) +
  ggtitle("B")



#Plot Responses
plot_res0 <- plot_response(dataraw = glm_dataset0,
                           data_emm = response0,
                           y_var = "ICU_days",
                           palette = "Oranges",
                           ylabel = "Dias de UTI",
                           strat = "Age",
                           legend = "Idade")

plot_res1 <- plot_response(dataraw = glm_dataset1,
                           data_emm = response1,
                           y_var = "ICU_days",
                           palette = "Purples",
                           ylabel = "Dias de UTI",
                           strat = "Age",
                           legend = "Idade")

plot_resint0 <- plot_response_int(dataraw = glm_dataset0,
                                 data_emm = responseint0,
                                 y_var = "ICU_days",
                                 ylabel = "Dia de UTI",
                                 strat = "Outcome_icu",
                                 legend = "Desfecho",
                                 palette = gradient_orange[c(6,8)])

plot_resint1 <- plot_response_int(dataraw = glm_dataset1,
                                  data_emm = responseint1,
                                  y_var = "ICU_days",
                                  ylabel = "Dia de UTI",
                                  strat = "Outcome_icu",
                                  legend = "Desfecho",
                                  palette = gradient_purple[c(6,8)])

#Save Plot
fig3a <- (plot_eff0 / plot_res0)
fig3b <- (plot_eff1 / plot_res1)
fig4a <- (plot_eff3 / plot_resint0)
fig4b <- (plot_eff4 / plot_resint1)


fig3 <- (fig3a | fig3b)  +
  patchwork::plot_layout(ncol = 2, widths = c(2/5,3/5))

fig4 <- (fig4a | fig4b)  +
  patchwork::plot_layout(ncol = 2, widths = c(2/5,3/5))

fig4
png("images/fig3_nbmodel.png", width = 14, height = 6, units = 'in', res = 300)
print(fig3)
while (!is.null(dev.list()))  dev.off()

png("images/fig4_nbmodel.png", width = 14, height = 6, units = 'in', res = 300)
print(fig4)
while (!is.null(dev.list()))  dev.off()


################################################################################
#############################  SUPPLEMENTARY  ##################################
################################################################################

############################# OUTLIERS ########################################
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

outliersplot <- function(var, data) {

   outliers <- data %>%
     dplyr::select(var, SamplesGroupYear, Num_Sample, allele1, allele2) %>%
     group_by(SamplesGroupYear) %>%
     mutate(outlier = ifelse(is_outlier(get(var)),
                             paste("Sample", Num_Sample, "\n", "=", get(var)), NA_character_),
            SamplesGroupYear = as.character(SamplesGroupYear)) %>%
     ggplot(aes(x = SamplesGroupYear, y = get(var),
                color = SamplesGroupYear)) +
     geom_boxplot() +
     geom_text(aes(label = outlier), hjust = -.1,  size = 2) +
     scale_color_manual(values = c(gradient_orange[6], gradient_purple[6])) +
     theme_pubr() +
     xlab("") +
     ylab(var) +
     theme(legend.position = "right",
           legend.direction="vertical",
           legend.title = element_blank())

   return(outliers)
}

fig_out1 <- outliersplot("ICU_days", data0) + outliersplot("HospitalPeriod_days", data0)

fig_out2 <- outliersplot("ICU_days", data_final) + outliersplot("HospitalPeriod_days", data_final)

fig_out <- (fig_out1 / fig_out2) + patchwork::plot_annotation(tag_levels = "A")

fig_out
png("images/fig_outsupp.png", width = 10, height = 10, units = 'in', res = 300)
print(fig_out)
while (!is.null(dev.list()))  dev.off()

############################# EPIDEMIOLOGICAL WEEKS ############################
#(ie, comparing the period of weeks 8 to 43 in 2020 vs the period from week 44 in 2020 to week 21 in 2021)
#https://www.hpsc.ie/notifiablediseases/resources/epidemiologicalweeks/

epi_week <- function(date) {paste(lubridate::epiweek(date), lubridate::isoyear(date), sep = "/")}

#https://www.thelancet.com/journals/lanres/article/PIIS2213-2600%2821%2900287-3/fulltext
wave1 <- as.Date(c("2020/8/1","2020/43/1"), "%Y/%U/%u")
wave1_mid <- wave1[1] + floor((wave1[2]-wave1[1])/2)

wave2 <-  as.Date(c("2020/43/1","2021/21/1"), "%Y/%U/%u")
wave2_mid <- wave2[1] + floor((wave2[2]-wave2[1])/2)


out_diff <-  setdiff(data0$Num_Sample, data_final$Num_Sample)
data0$seq <- seq_len(nrow(data0))


time_epiw <- data0 %>%
  dplyr::select(
    Num_Sample, seq, allele1, allele2, HospitalPeriod_days, SamplesGroupYear,
    Hospital_admission, ICU_admission, ICU_discharge) %>%
  dplyr::mutate(
    SamplesGroupYear = as.character(SamplesGroupYear),
    `Alta Hospitalar` =as.Date(as.Date(Hospital_admission) + HospitalPeriod_days),
    `Admissão Hospitalar`= as.Date(Hospital_admission),
    `Entrada UTI` = as.Date(ICU_admission),
    `Saída UTI` = as.Date(ICU_discharge)) %>%
  gather("event", "date", "Alta Hospitalar":"Saída UTI")


time_epiwout <- time_epiw %>%
  mutate(outlier = ifelse(HospitalPeriod_days > 65 & event == "Alta Hospitalar",
                          Num_Sample, NA_character_))
time_epiwout
Samplescollection_plot <- time_epiwout %>%
  ggplot(aes(x = date, y = seq, group = Num_Sample, col = event)) +
  annotate("rect", fill = "lightgrey", alpha = 0.5,
           xmin = wave1[1], xmax = wave1[2],
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "grey", alpha = 0.5,
           xmin = wave2[1], xmax = wave2[2],
           ymin = -Inf, ymax = Inf) +
  geom_text(
    aes(x = wave1_mid, y = 70, label = "Primeira Onda"),
    size = 3, vjust = 0, hjust = 0, color = "black"
  )+
  geom_text(
    aes(x = wave2_mid, y = 70, label = "Segunda Onda"),
    size = 3, vjust = 0, hjust = 0, color = "black"
  )+
  #geom_rect(data = rect1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha = 0.4)
  geom_line()  +
  geom_point(shape = 15, size = 1.2) +
  scale_x_date(date_breaks = "10 week", date_labels = "%U/%Y") +
  scale_color_manual(values = c(gradient_orange[c(5,8)], gradient_purple[c(5,8)])) +
  geom_text(aes(label = outlier), hjust = -.3,  size = 3, color = "black", show.legend = FALSE) +
  labs(y = "Amostras") +
  theme_pubr() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()#,
   # legend.position = "right"
  )

Samplescollection_plot

png("images/fig_collectionsupp.png", width = 15, height = 10, units = 'in', res = 300)
print(Samplescollection_plot)
while (!is.null(dev.list()))  dev.off()

oi <- data0 %>%
  filter(Num_Sample %in% out_diff) %>%
  dplyr::select(Num_Sample, allele1, allele2, SamplesGroupYear, HospitalPeriod_days, Outcome_icu) %>%
  gather("name", "allele", 2:3) %>%
  group_by(SamplesGroupYear, allele, Outcome_icu) %>%
  summarise(n = n())

###############################  MODEL CHECK  ##################################


glmmod0_check <- modelchecking(glmmod0, "glm")
glmmod1_check <- modelchecking(glmmod1, "glm")
glmmodint0_check <- modelchecking(glmmodint0, "glm")
glmmodint1_check <- modelchecking(glmmodint1, "glm")


#More Robust Testing
#plot(performance::check_model(glmmod0))
#plot(performance::check_distribution(glmmod0))
#plot(performance::binned_residuals(glmmod0))
