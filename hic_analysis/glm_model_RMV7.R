install.packages("reshape2")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("lme4")
install.packages("sjPlot")
install.packages("sjmisc")
install.packages("sjlabelled")
library(lme4)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(reshape2)
library(ggplot2)
library(dplyr)

load("/Users/tylim/Downloads/chromosome_interactions_matrix.RDS")

rmv7_hic_cm <- 
  dplyr::bind_rows(
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv7_hic_data[[4]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_hic_data[[5]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_hic_data[[6]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV7_rare"),
    
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv7_hic_data[[1]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_hic_data[[2]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_hic_data[[3]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV7_domi")
  ) %>%
  dplyr::rename(
    "X" = Var1,
    "Y" = Var2
  ) %>%
  dplyr::mutate(
    Replicate = paste0(Variant,"_",Replicate)
  )

rmv7_porec_cm <- 
  dplyr::bind_rows(
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv7_porec_data[[4]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_porec_data[[5]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_porec_data[[6]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV7_rare"),
    
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv7_porec_data[[1]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_porec_data[[2]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv7_porec_data[[3]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV7_domi")
  ) %>%
  dplyr::rename(
    "X" = Var1,
    "Y" = Var2
  ) %>%
  dplyr::mutate(
    Replicate = paste0(Variant,"_",Replicate)
  )

# RMV7 characteristics
rmv7_characteristics_df <-
  read.csv("Desktop/ICL/PJ2- Nic/S_pneumoniae_RMV7_dominant_characteristics.csv")

#code to measure prop_reads
porec_v7_contact_density <-
  rmv7_porec_cm %>%
  dplyr::mutate(distance_from_X = dplyr::if_else((X == min(X) & Y == max(Y)) | (X == max(X) & Y == min(Y)),
                                                 1,
                                                 abs(X-Y))
  ) %>%
  dplyr::mutate(distance_classification = dplyr::case_when(
    distance_from_X == 0 ~ "<1kb",
    distance_from_X == 1 ~ "1-2kb",
    TRUE ~ ">2kb"
  )
  ) %>%
  dplyr::group_by(Variant,Replicate,X) %>%
  dplyr::mutate(total_reads = sum(ContactFrequency)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(Variant,X,Replicate,distance_classification) %>%
  dplyr::mutate(prop_reads = sum(ContactFrequency)/total_reads) %>%
  dplyr::ungroup() %>%
  dplyr::select(Variant,X,Y,distance_classification,prop_reads,total_reads,Replicate,ContactFrequency) %>%
  dplyr::distinct() %>%
  dplyr::filter(ContactFrequency > 0) %>%
  dplyr::mutate(distance_classification = factor(distance_classification,
                                                 levels = c("<1kb","1-2kb",">2kb"))) %>%
  dplyr::mutate(Position = (X-1)*1000+500) %>%
  dplyr::left_join(rmv7_characteristics_df,by = c("Position")) %>% 
  dplyr::mutate(Position_Y = (Y-1)*1000+500) %>%
  dplyr::left_join(rmv7_characteristics_df %>%
                     rename_with(~paste0(.x, "_Y"), everything()),
                   by = c("Position_Y"))


# RMV7 annotation
rmv7_annotation.df <- data.frame(
  "start" = c(1792101,1601016,1268197,1218320,1347217,1409274,1631773),
  "end" = c(1804961,1619555,1296621,1223200,1352097,1414154,1636653),
  "Classification" = c("PRCI_uvrA","PRCI_dnaN","Tn916","rRNA","rRNA","rRNA","rRNA")
) %>% 
  dplyr::mutate(
    adjusted_start = start + 500,
    adjusted_end = end - 500
  )

ggplot(hic_v7_contact_density %>%
         dplyr::filter(distance_classification=="<1kb"),
       aes(x = MluCI,
           y=ContactFrequency)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10() +
  facet_wrap(~Replicate) +
  theme_bw()


# Annotated RMV7
genome_length <- 2112955
position_of_origin <- 1597872
annotated_rmv7_porec_df <-
  hic_v7_contact_density %>%
  dplyr::left_join(rmv7_annotation.df,
                   join_by(Position>=start,Position<=end)) %>%
  dplyr::select(-c(start,end)) %>%
  # Rescale everything between 0 and 1 for model fitting
  dplyr::mutate(GC = (GC+GC_Y)/200,
                MluCI = dplyr::case_when(MluCI == 0 | MluCI_Y == 0 ~ "None",
                                         MluCI == 1 | MluCI_Y == 1 ~ "Single",
                                         MluCI >= 5 & MluCI < 10 | MluCI_Y >= 5 & MluCI_Y < 10 ~ "Medium",
                                         MluCI < 5 | MluCI_Y < 5 ~ "Low",
                                         TRUE ~ "High"),
                NlaIII = dplyr::case_when(NlaIII == 0 | NlaIII_Y == 0 ~ "None",
                                          NlaIII == 1 | NlaIII_Y == 1 ~ "Single",
                                          NlaIII < 5 | NlaIII_Y < 5 ~ "Low",
                                          TRUE ~ "High")
  ) %>%
  dplyr::mutate(Classification = dplyr::if_else(is.na(Classification),"Unclassified",Classification)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(Distance_from_origin = min(abs(position_of_origin - Position),
                                           abs(genome_length - Position + position_of_origin),
                                           abs(genome_length - position_of_origin + Position))/genome_length
  ) %>%
  dplyr::ungroup()


write.csv(
  annotated_rmv7_porec_df,
  file = "Desktop/ICL/PJ2- Nic/RMV7_porec_characteristics_new.csv",
  row.names = FALSE,
  quote = FALSE
) 

ggplot(annotated_rmv7_hic_df,
       aes(x = ContactFrequency)) +
  geom_density() +
  scale_x_continuous(trans = "log10") +
  theme_bw()

annotated_RMV7_hic_df <-
  read.csv(file = "Desktop/ICL/PJ2- Nic/RMV7_hic_characteristics_new.csv") %>%
  dplyr::mutate(distance_classification = factor(distance_classification,
                                                 levels = c("<1kb","1-2kb",">2kb")),
                MluCI = factor(MluCI,
                               levels = c("None","Single","Low","Medium","High")),
                NlaIII = factor(NlaIII,
                                levels = c("None","Single","Low","High")),
                Classification = factor(Classification,
                                        levels = c("Unclassified","PRCI_uvrA","PRCI_dnaN","Tn916","rRNA"))
  )


# Just fit model to RMV7 
rmv7_porec_contact_density_model <-
  lme4::glmer(
    ContactFrequency ~ distance_classification + GC + MluCI + NlaIII + Classification*Variant + log(Distance_from_origin)*Variant + (1|Replicate) + (1|Position),
    data = annotated_RMV7_porec_df %>%
      dplyr::mutate(Position = factor(Position)),
    family = poisson(link = "log"),
    nAGQ = 0,
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  )

summary(rmv7_porec_contact_density_model)
sjPlot::plot_model(rmv7_porec_contact_density_model, 
                   show.values = TRUE, 
                   value.offset = .35,
                   value.size = 3,
                   title = "Model Coefficients of PoreC RMV7")
sjPlot::plot_model(rmv7_hic_contact_density_model, 
                   show.values = TRUE, 
                   value.offset = .35,
                   value.size = 3,
                   title = "Model Coefficients of HiC RMV7")

#prediction model 
annotated_rmv7_porec_with_prediction_df <-
  annotated_RMV7_porec_df %>%
  dplyr::mutate(log_poisson_prediction = fitted(rmv7_porec_contact_density_model)) %>%
  dplyr::mutate(residual = ContactFrequency - log_poisson_prediction)

ggplot(annotated_rmv7_porec_with_prediction_df,
       aes(x = log(Distance_from),
           y = ContactFrequency,
           colour = distance_classification)) +
  geom_point(alpha = 0.25) +
  facet_wrap(~Replicate, scales = 'free') +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_bw() +
  theme(strip.background = element_rect(NA))

ggplot(annotated_RMV8_porec_df,
       aes(x = log(Distance_from_origin),
           y = ContactFrequency,
           colour = distance_classification)) +
  geom_point() +
  facet_wrap(~Replicate) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

p <- ggplot(annotated_rmv7_porec_with_prediction_df,
       aes(x = Position,
           y = residual,
           colour = distance_classification)) +
  geom_point() +
  facet_grid(~Replicate) + 
  theme_bw()

p + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

#Obtain fixed effects from model 
fe_hic_7 <- tidy(rmv7_hic_contact_density_model, effects = "fixed")
fe_porec_7 <- tidy(rmv7_porec_contact_density_model, effects = "fixed")


# Add a column to distinguish between models
fe_hic_7$model <- "Hic model"
fe_porec_7$model <- "Porec model"

# Combine into a single data frame
fixed_effects_df_7 <- rbind(fe_hic_7, fe_porec_7)

write.csv(  
  fixed_effects_df_7,
  file = "Desktop/ICL/PJ2- Nic/fixed_effects_df_rmv7.csv",
  row.names = FALSE,
  quote = FALSE)

ggplot(fixed_effects_df_7, aes(x = term, y = estimate, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                position = position_dodge(width = 0.7), width = 0.25) +
  theme_bw() +
  scale_x_discrete(limits = c("(Intercept)", "distance_classification1-2kb", "distance_classification>2kb", "GC", "MluCISingle", "MluCISingle", "MluCIMedium", "MluCIHigh", "NlaIIISingle", "NlaIIILow", "NlaIIIHigh", "ClassificationPRCI_uvrA", "ClassificationPRCI_dnaN", "ClassificationTn916", "ClassificationrRNA", "VariantRMV7_rare", "log(Distance_from_origin)", "ClassificationPRCI_uvrA:VariantRMV7_rare", "ClassificationPRCI_uvrA:VariantRMV7_rare", "ClassificationTn916:VariantRMV7_rare", "ClassificationrRNA:VariantRMV7_rare")) + 
  labs(title = "Comparison of Fixed Effects Between Models",
       x = "Fixed Effect Term",
       y = "Estimate",
       fill = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##anything below are old codes##
annotated_rmv7_porec_with_prediction_df$model <- "Porec model"
annotated_rmv7_hic_with_prediction_df$model <- "Hic model"
obs_pred_combined_model_rmv7 <- rbind(annotated_rmv7_porec_with_prediction_df, annotated_rmv7_hic_with_prediction_df)

write.csv(  
  obs_pred_combined_model_rmv7,
  file = "Desktop/ICL/PJ2- Nic/obs_vs_pred_df_rmv7.csv",
  row.names = FALSE,
  quote = FALSE)

ggplot(obs_pred_combined_model,
       aes(x = Classification,
           y = ContactFrequency,
           colour = Variant)) +
  geom_point() +
  facet_grid(~model) + 
  theme_bw()

#filter overpredicted values 
filtered_annotated_rmv8_hic_with_prediction_df <- 
  annotated_rmv8_hic_with_prediction_df %>%
  dplyr::filter(log_poisson_prediction > 7.5 & log(ContactFrequency) < 4)

lowmap_annotated_rmv8_hic_with_prediction_df  <- 
  annotated_rmv8_hic_with_prediction_df %>%
  dplyr::filter(log(ContactFrequency) < log_poisson_prediction) %>%
  dplyr::mutate(log_ContactFrequncy = log(ContactFrequency))

write.csv(
  lowmap_annotated_rmv8_hic_with_prediction_df,
  file = "Desktop/ICL/PJ2- Nic/lowmap_annotated_rmv8_hic_with_prediction.csv",
  row.names = FALSE,
  quote = FALSE
) 



# (old) RMV7 model 
rmv7_porec_contact_density_model <-
  lme4::glmer(
    ContactFrequency ~ GC + MluCI + NlaIII + distance_classification + Classification + Distance_from_origin + (1|Replicate),
    data = annotated_rmv7_porec_df, 
    family = poisson(link = "log"),
    nAGQ = 0,
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  )

# (old) simple RMV7 model 
rmv7_hic_contact_density_model_1 <-
  glm(
    ContactFrequency ~ distance_classification + GC + MluCI + NlaIII + Classification + Distance_from_origin,
    data = annotated_rmv7_hic_df_1 %>% dplyr::filter(Variant == "RMV7_domi" & ContactFrequency > 10),
    family = poisson(link = "log")
  )
