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
rmv8_hic_cm <- 
  dplyr::bind_rows(
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv8_hic_data[[4]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_hic_data[[5]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_hic_data[[6]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV8_rare"),
    
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv8_hic_data[[1]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_hic_data[[2]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_hic_data[[3]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV8_domi")
  ) %>%
  dplyr::rename(
    "X" = Var1,
    "Y" = Var2
  ) %>%
  dplyr::mutate(
    Replicate = paste0(Variant,"_",Replicate)
  )

rmv8_porec_cm <- 
  dplyr::bind_rows(
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv8_porec_data[[4]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_porec_data[[5]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_porec_data[[6]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV8_rare"),
    
    dplyr::bind_rows(
      tibble(
        reshape2::melt(rmv8_porec_data[[1]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_porec_data[[2]]@matrix,
                       value.name = "ContactFrequency")
      ),
      tibble(
        reshape2::melt(rmv8_porec_data[[3]]@matrix,
                       value.name = "ContactFrequency")
      ),
      .id = "Replicate"
    ) %>%
      dplyr::mutate(Variant = "RMV8_domi")
  ) %>%
  dplyr::rename(
    "X" = Var1,
    "Y" = Var2
  ) %>%
  dplyr::mutate(
    Replicate = paste0(Variant,"_",Replicate)
  )

# RMV8 characteristics
rmv8_characteristics_df <-
  read.csv("Desktop/ICL/PJ2- Nic/S_pneumoniae_RMV8_dominant_characteristics.csv")

#code to measure prop_reads
porec_v8_contact_density <-
  rmv8_porec_cm %>%
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
  dplyr::left_join(rmv8_characteristics_df,by = c("Position")) %>% 
  dplyr::mutate(Position_Y = (Y-1)*1000+500) %>%
  dplyr::left_join(rmv8_characteristics_df %>%
                     rename_with(~paste0(.x, "_Y"), everything()),
                   by = c("Position_Y"))

#RMV8n annotation 
rmv8_annotation.df<-data.frame(
  "start" = c(25288,1998066,16401,1777850,1879619,1943684),
  "end" = c(58975,2010058,21281,1782730,1884499,1948564),
  "Classification" = c("prophage","PRCI_mal","rRNA","rRNA","rRNA","rRNA")
) %>% 
  dplyr::mutate(
    adjusted_start = start + 500,
    adjusted_end = end - 500
  )

ggplot(hic_v8_contact_density %>%
         dplyr::filter(distance_classification=="<1kb"),
       aes(x = MluCI,
           y=ContactFrequency)) +
  geom_point() +
  geom_smooth() +
  scale_y_log10() +
  facet_wrap(~Replicate) +
  theme_bw()


# Annotated RMV8
genome_length <- 2147252
position_of_origin <- 64
annotated_RMV8_hic_df <-
  porec_v8_contact_density %>%
  dplyr::left_join(rmv8_annotation.df,
                   join_by(Position>=start,Position<end)) %>%
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
  annotated_RMV8_porec_df,
  file = "Desktop/ICL/PJ2- Nic/RMV8_porec_characteristics_new.csv",
  row.names = FALSE,
  quote = FALSE
) 

ggplot(annotated_RMV8_porec_df,
       aes(x = ContactFrequency)) +
  geom_density() +
  scale_x_continuous(trans = "log10") +
  theme_bw()


annotated_rmv8_porec_df <-
  read.csv(file = "Desktop/ICL/PJ2- Nic/RMV8_porec_characteristics_new.csv") %>%
  dplyr::mutate(distance_classification = factor(distance_classification,
                                                 levels = c("<1kb","1-2kb",">2kb")),
                MluCI = factor(MluCI,
                               levels = c("None","Single","Low","Medium","High")),
                NlaIII = factor(NlaIII,
                                levels = c("None","Single","Low","High")),
                Classification = factor(Classification,
                                        levels = c("Unclassified","prophage","PRCI_mal","rRNA"))
  )


# Fit model to RMV8 
rmv8_porec_contact_density_model <-
  lme4::glmer(
    ContactFrequency ~ distance_classification + GC + MluCI + NlaIII + Classification*Variant + log(Distance_from_origin)*Variant + (1|Replicate) + (1|Position),
    data = annotated_rmv8_porec_df %>%
      dplyr::mutate(Position = factor(Position)),
    family = poisson(link = "log"),
    nAGQ = 0,
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  )


summary(rmv8_porec_contact_density_model)
sjPlot::plot_model(rmv8_porec_contact_density_model, 
                   show.values = TRUE, 
                   value.offset = .35,
                   value.size = 3,
                   title = "Model Coefficients of PoreC RMV8")

sjPlot::plot_model(rmv8_hic_contact_density_model, 
                   show.values = TRUE, 
                   value.offset = .35,
                   value.size = 3,
                   title = "Model Coefficients of HiC RMV8")


#prediction model
rmv8_hic_contact_density_model_data <- model.frame(rmv8_hic_contact_density_model)

annotated_rmv8_hic_with_prediction_df <-
  rmv8_hic_contact_density_model_data %>%
  dplyr::mutate(log_poisson_prediction = fitted(rmv8_hic_contact_density_model)) %>%
  dplyr::mutate(residual = ContactFrequency - log_poisson_prediction) %>% 
  dplyr::mutate(Position = as.numeric(as.character(Position)))

#plot CF vs fitted/predicted data
ggplot(annotated_rmv8_porec_with_prediction_df,
       aes(x = ContactFrequency,
           y = log_poisson_prediction,
           colour = distance_classification)) +
  geom_point(alpha = 0.25) +
  facet_wrap(~Replicate, scales = 'free') +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_bw() +
  theme(strip.background = element_rect(NA))

#plot residual vs Position
ggplot(annotated_rmv8_porec_with_prediction_df,
            aes(x = Position,
                y = residual,
                colour = distance_classification)) +
  geom_point() +
  theme_bw() +
  facet_grid(~Replicate) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust= 1))

# Extract fixed effects
fe_hic_8 <- tidy(rmv8_hic_contact_density_model, effects = "fixed")
fe_porec_8 <- tidy(rmv8_porec_contact_density_model, effects = "fixed")

# Add a column to distinguish between models
fe_hic_8$model <- "Hic model"
fe_porec_8$model <- "Porec model"

# Combine into a single data frame
fixed_effects_df_8 <- rbind(fe_hic_8, fe_porec_8)

write.csv(  
  fixed_effects_df_8,
file = "Desktop/ICL/PJ2- Nic/fixed_effects_df_rmv8.csv",
            row.names = FALSE,
            quote = FALSE)

# Plot the fixed effects 
ggplot(fixed_effects_df, aes(x = term, y = estimate, fill = model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error),
                position = position_dodge(width = 0.7), width = 0.25) +
  theme_bw() +
  scale_x_discrete(limits = c("(Intercept)", "distance_classification1-2kb", "distance_classification>2kb", "GC", "MluCISingle", "MluCISingle", "MluCIMedium", "MluCIHigh", "NlaIIISingle", "NlaIIILow", "NlaIIIHigh", "Classificationprophage", "ClassificationPRCI_mal", "ClassificationrRNA", "VariantRMV8_rare", "log(Distance_from_origin)", "Classificationprophage:VariantRMV8_rare", "ClassificationPRCI_mal:VariantRMV8_rare", "ClassificationrRNA:VariantRMV8_rare")) + 
  labs(title = "Comparison of Fixed Effects Between Models",
       x = "Fixed Effect Term",
       y = "Estimate",
       fill = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##anything below is old codes##
annotated_rmv8_porec_with_prediction_df$model <- "Porec model"
annotated_rmv8_hic_with_prediction_df$model <- "Hic model"
obs_pred_combined_model_v8 <- rbind(annotated_rmv8_porec_with_prediction_df, annotated_rmv8_hic_with_prediction_df)

write.csv(  
  obs_pred_combined_model_v8,
  file = "Desktop/ICL/PJ2- Nic/obs_vs_pred_df_rmv8.csv",
  row.names = FALSE,
  quote = FALSE)

ggplot(obs_pred_combined_model_v8,
       aes(x = NlaIII,
           y = ContactFrequency,
           colour = Variant,Classification )) +
  geom_point() +
  facet_grid(~model) + 
  theme_bw()

#old models for single variant 
rmv8_hic_contact_density_model <-
  lme4::glmer(
    ContactFrequency ~ GC + MluCI + NlaIII + Classification + Distance_from_origin + (1|Replicate),
    data = annotated_rmv8_hic_df_01kb, 
    family = poisson(link = "log"),
    nAGQ = 0,
    control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  )

rmv8_hic_contact_density_model_1 <-
  glm(
    ContactFrequency ~ distance_classification + GC + MluCI + NlaIII + Classification + Distance_from_origin,
    data = annotated_rmv8_hic_df_1 %>% dplyr::filter(Variant == "RMV8_domi" & ContactFrequency > 10),
    family = poisson(link = "log")
  )

annotated_rmv8_hic_with_prediction_df_1 <-
  annotated_rmv8_hic_df_1 %>%
  dplyr::filter(Variant == "RMV8_rare" & ContactFrequency > 10) %>%
  dplyr::mutate(log_poisson_prediction = predict(rmv8_hic_contact_density_model_1))

# Create a new dataframe for prediction, ensuring it matches the model data exactly
model_data_1kb <- model.frame(rmv8_hic_contact_density_model_1kb)

annotated_rmv8_hic_with_prediction_df_1kb <-
  model_data_1kb %>%
  dplyr::mutate(log_poisson_prediction = predict(rmv8_hic_contact_density_model_1kb))

ggplot(annotated_rmv8_hic_with_prediction_df_1kb,
       aes(x = log(ContactFrequency),
           y = log_poisson_prediction,
           colour = Replicate)) +
  geom_point(position = position_jitter(), alpha = 0.25) +
  theme_bw()

ggplot(annotated_rmv8_hic_with_prediction_df,
       aes(x = distance_classification,
           y = log(ContactFrequency),
           colour = Replicate)) +
  geom_point(position = position_jitter(), alpha = 0.25) +
  theme_bw()

ggplot(annotated_rmv8_hic_with_prediction_df_1kb,
       aes(x = log(ContactFrequency),
           y = log_poisson_prediction,
           colour = Replicate)) +
  geom_point() +
  theme_bw()


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

