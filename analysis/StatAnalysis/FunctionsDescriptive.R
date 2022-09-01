
#Options
options(scipen = 999) #no scientific notation


## Alleles Frequencies and Fisher Tests ----------------------------------------

  # Function: frequencer()
  # Calculate the frequencies between "var" categorical variables
  # Add Redome frequencies to plot allele frequencies

frequencer <- function(data, var, redome){

    freqtab <- allele_freq %>%
      janitor::tabyl(.data[[var]], allele) %>%
      janitor::adorn_percentages() %>%
      tidyr::gather(., "allele", "Freq", starts_with("B*")) %>%
      dplyr::bind_rows(dplyr::select(Redome_freq, -n_Redome)) #%>%

    freqtab[[var]] <- ifelse(is.na(freqtab[[var]]),"Redome", as.character(freqtab[[var]]))


  return(freqtab)

}

  # Function: stats()
  # Give the count data and frequencies for each allele in function of "var".

stats <- function(data, var, redome) {

  stats <- data_final %>%
    dplyr::select(allele1, allele2, var) %>%
    tidyr::gather(., "Name", "allele",1:2) %>%
    dplyr::group_by(.data[[var]], allele) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(allele_frequency = prop.table(n)) %>%
    tidyr::pivot_wider(id_cols = allele,
                       names_from = var,
                       values_from = c("allele_frequency", "n")) %>%
    dplyr::full_join(redome, by = "allele")%>%
    dplyr::select(sort(names(.))) %>%
    replace(is.na(.),0)

  return(stats)

  }

  # Function: fisher_tester()
  # Calculate the fisher's exact tests for any combinations of "groups"'s elements.
  # fisher_tester execute row wise tests through a 2x2 contigency table where
  # Group 1 (allele) x Group 2 (allele) | Group 1 (no-allele) x Group 2 (no-allele)
  # Remove the zeros from the alleles inserted by Redome data
  # Execute multiple comparisons corrections with Bonferroni approach
  # Save the output only the results p < "plim"

fisher_tester <- function(table, groups, plim){

  combinations <- data.frame(combn(groups, 2))

  final <- data.frame()

  if("n_Redome" %in% groups){
    table <- filter_at(
      table, vars(-c(n_Redome,allele)),
      any_vars(. > 0))
  } else {
    table <- filter_at(
      table, vars(-c(n_Redome,allele)),
      any_vars(. > 0))
  }

  for (x in seq(combinations)){

    vars <- as.vector(combinations[x])[[1]]
    var1 <- vars[1]
    var2 <- vars[2]

    names(table) <- gsub(pattern = "n_", replacement = "", x = names(table))

    tab_fis <- table %>%
      column_to_rownames(var="allele") %>%
      dplyr::select(var1, var2) %>%
      as.matrix(labels = TRUE)

    tab_frame <- table %>%
      dplyr::select("allele", var1, var2) %>%
      mutate(`Outros 1 (n)` =  sum(table[[var1]]) - table[[var1]],
             `Outros 2 (n)` =  sum(table[[var2]]) - table[[var2]]) %>%
      rename(`Grupo 1 (n)` = var1,
             `Grupo 2 (n)` = var2)

    colnames(tab_fis) <- c("Grupo 1 (n)", "Grupo 2 (n)")

    fis_result <-  rstatix::row_wise_fisher_test(tab_fis, alternative = "two.sided", p.adjust.method = "bonferroni") %>%
      dplyr::select(-n) %>%
      dplyr::rename(allele = "group") %>%
      dplyr::mutate(`Grupo 1` = vars[1],
             `Grupo 2`= vars[2])

    fis_count_result <- fis_result %>%
      dplyr::full_join(tab_frame, by = "allele") %>%
      dplyr::filter(p <= plim)

    if (nrow(fis_count_result) > 0) {final <- rbind(final,fis_count_result)}
  }

   return(final)

}




