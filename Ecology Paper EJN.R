

# Esteban J Nogueira - Ecology paper

# "Impacts of thermal disturbances on the spatial structure of reef fish communities are influenced by biogeography factors"


#### Loop ####
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(betapart)
library(openxlsx)
library(vegan)

process_site_beta <- function(arquivo_excel, sheet_name, N = 1000,
                              out_dir = NULL,
                              hellinger_transform = TRUE) {
  
  df <- read_excel(arquivo_excel, sheet = sheet_name)
  
  if (!"Local_census" %in% names(df)) {
    stop("Não encontrei coluna 'Local_census' na aba: ", sheet_name)
  }
  
  df <- df %>%
    mutate(
      Year = as.factor(Year),
      Local_census = paste(Dive_Site, Year, Local_census, sep = "_")
    )
  
  sub1 <- df %>%
    select(Dive_Site, Year, Area, Local_census) %>%
    distinct()
  
  # MSA
  msa_by_year <- tapply(sub1$Area, sub1$Year, sum)
  msa <- min(msa_by_year, na.rm = TRUE)
  
  years <- sort(unique(sub1$Year))
  

  species_cols <- setdiff(names(df), c("Dive_Site","Year","Area","Local_census"))
  

  select_transects_one_iter <- function(iter_id) {
    sel <- lapply(years, function(yy) {
      df_y <- sub1 %>% filter(Year == yy) %>% slice_sample(prop = 1)
      
      cum_area <- cumsum(df_y$Area)
      keep_n <- which(cum_area >= msa)[1]
      if (is.na(keep_n)) keep_n <- nrow(df_y)
      
      df_y %>%
        slice(1:keep_n) %>%
        transmute(iter = iter_id, Year = yy, Local_census)
    })
    bind_rows(sel)
  }
  
  selected_transects <- bind_rows(lapply(seq_len(N), select_transects_one_iter))
  
  results_beta <- vector("list", N)
  
  for (i in seq_len(N)) {
    chosen <- selected_transects %>% filter(iter == i)
    
    database_red <- df %>% semi_join(chosen, by = c("Year","Local_census"))
    
    comm_long <- database_red %>%
      select(Year, all_of(species_cols)) %>%
      pivot_longer(cols = all_of(species_cols),
                   names_to = "species", values_to = "abundance")
    
    comm_year <- comm_long %>%
      group_by(Year, species) %>%
      summarise(abundance = sum(as.numeric(abundance)), .groups = "drop") %>%
      pivot_wider(names_from = species, values_from = abundance, values_fill = 0) %>%
      arrange(Year)
    
    yy <- comm_year$Year
    

    # 1) Sørensen (PA): sor/sim/nestedness
   
    CommPA <- decostand(as.matrix(comm_year %>% select(-Year)), method = "pa")
    beta_pa <- beta.pair(betapart.core(CommPA))
    
    sor_mat <- as.matrix(beta_pa$beta.sor)
    sim_mat <- as.matrix(beta_pa$beta.sim)
    nes_mat <- as.matrix(beta_pa$beta.sne)
    

    # 2) Bray abundance: total + balance + gradient 

    CommAB <- as.matrix(comm_year %>% select(-Year))
    
    CommAB <- decostand(CommAB, method = "hellinger") 
    
    
    beta_ab <- beta.pair.abund(CommAB, index.family = "bray")
    
    bray_tot_mat <- as.matrix(beta_ab$beta.bray)
    bray_bal_mat <- as.matrix(beta_ab$beta.bray.bal)
    bray_gra_mat <- as.matrix(beta_ab$beta.bray.gra)
    
    total_df <- data.frame(
      iter = i,
      site1 = yy[row(sor_mat)],
      site2 = yy[col(sor_mat)],
      sor_fam = as.vector(sor_mat),
      sim_fam = as.vector(sim_mat),
      nes_fam = as.vector(nes_mat),
      bray_total = as.vector(bray_tot_mat),
      bray_balance = as.vector(bray_bal_mat),
      bray_gradient = as.vector(bray_gra_mat)
    ) %>%
      filter(site1 != site2) %>%
      mutate(pair = map2_chr(site1, site2, ~ paste(sort(c(.x, .y)), collapse = "_"))) %>%
      distinct(iter, pair, .keep_all = TRUE) %>%
      select(-pair)
    
    results_beta[[i]] <- total_df
  }
  
  matrix_total_values <- bind_rows(results_beta) %>%
    mutate(Locality = sheet_name)
  
  total_comparison <- matrix_total_values %>%
    mutate(combination = paste(site2, site1, sep = "_")) %>%
    group_by(Locality, combination) %>%
    summarise(
      sor_fam        = mean(sor_fam, na.rm=TRUE),
      sim_fam        = mean(sim_fam, na.rm=TRUE),
      nes_fam        = mean(nes_fam, na.rm=TRUE),
      bray_total     = mean(bray_total, na.rm=TRUE),
      bray_balance   = mean(bray_balance, na.rm=TRUE),
      bray_gradient  = mean(bray_gradient, na.rm=TRUE),
      .groups="drop"
    )
  
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.xlsx(matrix_total_values, file.path(out_dir, paste0(sheet_name, "_Total_iter.xlsx")))
    write.xlsx(total_comparison,   file.path(out_dir, paste0(sheet_name, "_Media.xlsx")))
    write.xlsx(selected_transects, file.path(out_dir, paste0(sheet_name, "_SelectedTransects.xlsx")))
  }
  
  list(
    selected_transects = selected_transects,
    matrix_total_values = matrix_total_values,
    total_comparison = total_comparison,
    msa = msa
  )
}


arquivo_excel <- "C:/" #Insert user directory
abas <- c("Rocas", "Noronha", "Trindade", "SpSp", "SantaCatarina", "Arraial")

res_list <- lapply(abas, function(s) process_site_beta(arquivo_excel, s, N = 1000, hellinger_transform = TRUE))

SorBra_big <- bind_rows(lapply(res_list, `[[`, "matrix_total_values"))
SorBra_meanPairs <- bind_rows(lapply(res_list, `[[`, "total_comparison"))

write.xlsx(SorBra_big, 
           "C:/") #Insert user directory

write.xlsx(SorBra_meanPairs, "C:/") #Insert user directory


#### Spatial Beta ####
library(readxl)
library(betapart)
library(vegan)
library(openxlsx)
library(dplyr)


processar_ilha <- function(arquivo, nome_ilha) {
  
  
  dados <- read_excel(arquivo)
  
  
  anos <- sort(unique(dados$Year))
  
  for (ano in anos) {
    # ----------- Presence/abscence (Sorensen) -----------
    IlhaPa <- decostand(dados[dados$Year == ano, c(-1, -2, -3, -4)], method = "pa")
    Ilha.core <- betapart.core(IlhaPa)
    IlhaBeta <- beta.pair(IlhaPa, index.family = "sorensen")
    
    data.frameIlha <- data.frame(
      Sorensen = round(as.numeric(IlhaBeta$beta.sor), 2),
      Simpson = round(as.numeric(IlhaBeta$beta.sim), 2),
      Aninhamento = round(as.numeric(IlhaBeta$beta.sne), 2)
    )
    
    # Results
    write.xlsx(
      data.frameIlha,
      file = paste0(
        "C:/", #Insert user directory
        nome_ilha, "Beta", ano, ".xlsx"
      )
    )
    
    # ----------- Abundance (Bray-Curtis) -----------
    resultado_AB <- beta.pair.abund(
      dados[dados$Year == ano, c(-1, -2, -3, -4)],
      index.family = "bray"
    )
    
    data.frameIlhaAB <- data.frame(
      Bray = round(as.numeric(resultado_AB$beta.bray), 2),
      Balanceada = round(as.numeric(resultado_AB$beta.bray.bal), 2),
      Gradiente = round(as.numeric(resultado_AB$beta.bray.gra), 2)
    )
    
    # Results
    write.xlsx(
      data.frameIlhaAB,
      file = paste0(
        "C:/", #Insert user directory
        nome_ilha, "Bray", ano, ".xlsx"
      )
    )
    
    message("Processado: ", nome_ilha, " - Ano: ", ano)
  }
}

# Localities
arquivos <- list(
  Arraial = "C:/",
  SpSp = "C:/",
  Rocas = "C:/",
  Noronha = "C:/",
  SantaCatarina = "C:/",
  Trindade = "C:/" #Insert user directory
)


for (ilha in names(arquivos)) {
  processar_ilha(arquivos[[ilha]], ilha)
}


#### nMDS ####
library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
library(readxl)
library(scales)
# Presence absence
arquivo_excel <- "C:/" 
abas <- c("Rocas", "Noronha", "Trindade", "SpSp", "SantaCatarina", "Arraial", "Juntos Media", "Media espacial", "SorBra espacial")


walk(abas, function(sheet_name) {
  df <- read_excel(arquivo_excel, sheet = sheet_name)
  assign(paste0(sheet_name, ""), df, envir = .GlobalEnv)
})



#JACCARD#
novo <- decostand(`Juntos Media`[, 3:ncol(`Juntos Media`)], method = "pa")
novo[is.na(novo)] <- 0


dados_dist <- vegdist(novo, method = "jaccard")

#nMDS
set.seed(123)
nmds_result <- metaMDS(dados_dist)


nmds_plot <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_plot$Ilhas <- `Juntos Media`$Dive_Site
nmds_plot$Anos <- `Juntos Media`$Year

stress_val <- round(nmds_result$stress, 3)
nmds_plot <- nmds_plot %>%
  group_by(Ilhas) %>%
  arrange(Anos) %>%
  mutate(point_type = case_when(
    Anos == min(Anos) ~ "first",
    Anos == max(Anos) ~ "last",
    TRUE ~ "middle"
  )) %>%
  ungroup()
site_cols <- setNames(hue_pal()(length(unique(nmds_plot$Ilhas))), unique(nmds_plot$Ilhas))

# Plot
Jac<- ggplot(nmds_plot, aes(x = NMDS1, y = NMDS2, color = Ilhas, group = interaction(Ilhas), shape =point_type )) +
  geom_point( size = 3.4) +
  geom_path(aes(group = Ilhas), alpha = 0.7, size = 0.9) +
  geom_text(data = filter(nmds_plot, point_type !="middle"),
            aes(label = Anos), size = 4, vjust = 1.5, show.legend = FALSE, colour = "black") +
  scale_shape_manual(values = c(first = 17, middle = 1, last = 15))+
  theme_minimal(base_size = 14) +
  labs(title = paste0("Jaccard Assemblages nMDS (Stress = ", stress_val, ")",linetype = ""),
       x = "nMDS1",
       y = "nMDS2",
       color = "Locallity",
       shape = "Years") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(0, "mm"),
        legend.key.height = unit(0, "mm"),
        legend.background = element_rect(fill = "transparent", colour = NA))

Jac

# Save the plot
ggsave("C:/", # User directory
       plot = Jac, width = 30, height = 20, units = "cm")



#BRAY-CURTIS#



dados_dist <- vegdist(`Juntos Media`[, 3:ncol(`Juntos Media`)], method = "bray")

#nMDS
set.seed(123)
nmds_result <- metaMDS(dados_dist)


nmds_plot <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_plot$Ilhas <- `Juntos Media`$Dive_Site
nmds_plot$Anos <- `Juntos Media`$Year

stress_val <- round(nmds_result$stress, 3)
nmds_plot <- nmds_plot %>%
  group_by(Ilhas) %>%
  arrange(Anos) %>%
  mutate(point_type = case_when(
    Anos == min(Anos) ~ "first",
    Anos == max(Anos) ~ "last",
    TRUE ~ "middle"
  )) %>%
  ungroup()

site_cols <- setNames(hue_pal()(length(unique(nmds_plot$Ilhas))), unique(nmds_plot$Ilhas))
# Plot
bray<-ggplot(nmds_plot, aes(x = NMDS1, y = NMDS2, color = Ilhas, group = interaction(Ilhas), shape =point_type )) +
  geom_point( size = 3.4) +
  geom_path(aes(group = Ilhas), alpha = 0.7, size = 0.9) +
  geom_text(data = filter(nmds_plot, point_type !="middle"),
            aes(label = Anos), size = 4, vjust = 1.5, show.legend = FALSE, colour = "black") +
  scale_shape_manual(values = c(first = 17, middle = 1, last = 15))+
  theme_minimal(base_size = 14) +
  labs(title = paste0("Bray-Curtis Assemblages nMDS (Stress = ", stress_val, ")",linetype = ""),
       x = "nMDS1",
       y = "nMDS2",
       color = "Locallity",
       shape = "Years") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(0, "mm"),
        legend.key.height = unit(0, "mm"),
        legend.background = element_rect(fill = "transparent", colour = NA))

bray
# Save the plot
ggsave("C:/", #User directory
       plot = bray, width = 30, height = 30, units = "cm")

#### Linear Models ####

library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(lmtest)
library(car)
library(purrr)


#  SØRENSEN – MODEL


modelos_Sor <- `Media espacial` %>%
  group_by(Locality, Type) %>%
  nest() %>%
  mutate(
    modelo = map(data, ~ lm(Sor ~ Year, data = .x)),
    
   
    tidy = map(modelo, ~ tidy(.x, conf.int = TRUE)),
    
   
    glance = map(modelo, glance),
    
    
    shapiro = map(modelo, ~ shapiro.test(residuals(.x))),
    breusch_pagan = map(modelo, bptest),
    durbin_watson = map(modelo, dwtest),
    
    
    cooks = map(modelo, cooks.distance),
    
    max_cooks = map_dbl(
      cooks,
      ~ max(.x, na.rm = TRUE)
    )
  )

for(i in seq_along(modelos_Sor$modelo)) {
  
  local <- modelos_Sor$Locality[i]
  tipo <- modelos_Sor$Type[i]
  
  modelo <- modelos_Sor$modelo[[i]]
  
  png(
    filename = paste0(
      "Diagnostic_Sorensen_",
      local,
      "_",
      tipo,
      ".png"
    ),
    width = 1600,
    height = 1600,
    res = 200
  )
  
  par(mfrow = c(2, 2))
  plot(modelo)
  
  par(mfrow = c(1, 1))
  
  dev.off()
}

#   SØRENSEN TABLE


resultados_Sor <- modelos_Sor %>%
  mutate( 
    coef = map(
      tidy,
      ~ .x %>%
        filter(term == "Year")
    ),
    
    stats = map(
      glance,
      ~ .x
    )
  ) %>%
  select(
    Locality, Type, coef, stats, shapiro,
    breusch_pagan, durbin_watson, max_cooks
  ) %>%
  mutate(
    estimate = map_dbl(coef, ~ .x$estimate),
    std.error = map_dbl(coef, ~ .x$std.error),
    statistic = map_dbl(coef, ~ .x$statistic),
    p.value = map_dbl(coef, ~ .x$p.value),
    CI_low = map_dbl(coef, ~ .x$conf.low),
    CI_high = map_dbl(coef, ~ .x$conf.high),
    
    R2 = map_dbl(stats, ~ .x$r.squared),
    R2_adj = map_dbl(stats, ~ .x$adj.r.squared),
    
    Shapiro_W = map_dbl(
      shapiro,
      ~ unname(.x$statistic)
    ),
    
    Shapiro_p = map_dbl(
      shapiro,
      ~ .x$p.value
    ),
    
    BP_statistic = map_dbl(
      breusch_pagan,
      ~ unname(.x$statistic)
    ),
    
    BP_p = map_dbl(
      breusch_pagan,
      ~ .x$p.value
    ),
    
    DW_statistic = map_dbl(
      durbin_watson,
      ~ unname(.x$statistic)
    ),
    
    DW_p = map_dbl(
      durbin_watson,
      ~ .x$p.value
    )
  ) %>%
  select(
    Locality, Type, estimate, std.error, statistic, p.value, CI_low, CI_high, 
    R2, R2_adj, Shapiro_W, Shapiro_p, BP_statistic, BP_p, DW_statistic, DW_p, max_cooks
  )


resultados_Sor
write.xlsx(resultados_Sor, 
           "C:/") #User directory

# BRAY-CURTIS – MODEL

modelos_Bray <- `Media espacial`%>%
  group_by(Locality, Type) %>%
  nest() %>%
  mutate(
    modelo = map(data, ~ lm(Bray ~ Year, data = .x)),
    
    tidy = map(
      modelo,
      ~ tidy(.x, conf.int = TRUE)
    ),
    
    glance = map(
      modelo,
      glance
    ),
    
    shapiro = map(
      modelo,
      ~ shapiro.test(residuals(.x))
    ),
    
    breusch_pagan = map(
      modelo,
      bptest
    ),
    
    durbin_watson = map(
      modelo,
      dwtest
    ),
    
    cooks = map(
      modelo,
      cooks.distance
    ),
    
    max_cooks = map_dbl(
      cooks,
      ~ max(.x, na.rm = TRUE)
    )
  )


# DIAGNOSTIC BRAY-CURTIS


for(i in seq_along(modelos_Bray$modelo)) {
  
  local <- modelos_Bray$Locality[i]
  tipo <- modelos_Bray$Type[i]
  
  modelo <- modelos_Bray$modelo[[i]]
  
  png(
    filename = paste0(
      "Diagnostic_Bray_",
      local,
      "_",
      tipo,
      ".png"
    ),
    width = 1600,
    height = 1600,
    res = 200
  )
  
  par(mfrow = c(2, 2))
  plot(modelo)
  par(mfrow = c(1, 1))
  
  dev.off()
}

# BRAY-CURTIS TABLE


resultados_Bray <- modelos_Bray %>%
  mutate(
    coef = map(
      tidy,
      ~ .x %>%
        filter(term == "Year")
    ),
    
    stats = map(
      glance,
      ~ .x
    )
  ) %>%
  select(
    Locality, Type, coef, stats, shapiro,
    breusch_pagan, durbin_watson, max_cooks
  ) %>%
  mutate(
    estimate = map_dbl(coef, ~ .x$estimate),
    std.error = map_dbl(coef, ~ .x$std.error),
    statistic = map_dbl(coef, ~ .x$statistic),
    p.value = map_dbl(coef, ~ .x$p.value),
    CI_low = map_dbl(coef, ~ .x$conf.low),
    CI_high = map_dbl(coef, ~ .x$conf.high),
    
    R2 = map_dbl(stats, ~ .x$r.squared),
    R2_adj = map_dbl(stats, ~ .x$adj.r.squared),
    
    Shapiro_W = map_dbl(
      shapiro,
      ~ unname(.x$statistic)
    ),
    
    Shapiro_p = map_dbl(
      shapiro,
      ~ .x$p.value
    ),
    
    BP_statistic = map_dbl(
      breusch_pagan,
      ~ unname(.x$statistic)
    ),
    
    BP_p = map_dbl(
      breusch_pagan,
      ~ .x$p.value
    ),
    
    DW_statistic = map_dbl(
      durbin_watson,
      ~ unname(.x$statistic)
    ),
    
    DW_p = map_dbl(
      durbin_watson,
      ~ .x$p.value
    )
  ) %>%
  select(
    Locality, Type, estimate, std.error, statistic, p.value, CI_low, CI_high,
    R2, R2_adj, Shapiro_W, Shapiro_p, BP_statistic, BP_p, DW_statistic,
    DW_p, max_cooks
  )


resultados_Bray
write.xlsx(resultados_Bray, 
           "C:/") #User directory

# Graphic
dados_mean <- `Media espacial` %>%
  group_by(Locality, Type, Year) %>%
  summarise(
    Sor = mean(Sor, na.rm = TRUE),
    Bray = mean(Bray, na.rm = TRUE),
    .groups = "drop"
  )

# Sorensen
regressoesSor <- `Media espacial` %>%
  group_by(Locality) %>%
  do(tidy(lm(Sor ~ Year, data = .))) %>%
  ungroup()

regressoesSor
inclinações <- regressoesSor %>%
  filter(term == "Year") %>%
  select(Locality, estimate, std.error, p.value)
regressoes_glance <- `Media espacial` %>%
  group_by(Locality) %>%
  do(glance(lm(Sor ~ Year, data = .))) %>%
  ungroup()

regressoes_glance %>%
  select(Locality, r.squared, adj.r.squared, p.value)

df_mean <- `Media espacial` %>%
  group_by(Locality, Type, Year) %>%
  summarise(Sor_mean = mean(Sor, na.rm = TRUE), .groups = "drop")
SorLinear <-ggplot(`Media espacial`, aes(x = Year, y = Sor, color = Locality,fill= Locality, linetype = Type)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(method = "lm", se = T, size = 1.2) +
  geom_point(data = df_mean, 
             aes(x = Year, y = Sor_mean, color = Locality,), 
             size = 2.5) +
  labs(x = "Year", 
       y = "Sørensen Dissimilarity", 
       color = "Locality",
       linetype = "",
       shape = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90")
  )
SorLinear
# Bray-Curtis

regressoesBray <- `Media espacial` %>%
  group_by(Locality) %>%
  do(tidy(lm(Bray ~ Year, data = .))) %>%
  ungroup()
regressoesBray
inclinações <- regressoesBray %>%
  filter(term == "Year") %>%
  select(Locality, estimate, std.error, p.value)
regressoes_glance <- `Media espacial` %>%
  group_by(Locality) %>%
  do(glance(lm(Bray ~ Year, data = .))) %>%
  ungroup()

regressoes_glance %>%
  select(Locality, r.squared, adj.r.squared, p.value)



df_mean <- `Media espacial` %>%
  group_by(Locality, Type, Year) %>%
  summarise(Sor_mean = mean(Bray, na.rm = TRUE), .groups = "drop")
BrayLinear<-ggplot(`Media espacial`, aes(x = Year, y = Bray, color = Locality,fill= Locality, linetype = Type)) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_smooth(method = "lm", se = T, size = 1.2) +
  geom_point(data = df_mean, 
             aes(x = Year, y = Sor_mean, color = Locality,), 
             size = 2.5) +
  labs(x = "Year", 
       y = "Bray-Curtis Dissimilarity", 
       color = "Locality",
       linetype = "",
       shape = "") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.major.y = element_line(color = "grey90")
  )

BrayLinear

ggsave("C:/", #User directory
       plot = BrayLinear, width = 30, height = 20, units = "cm")
ggsave("C:/", #User directory
       plot = SorLinear, width = 30, height = 20, units = "cm")

Linear<-(BrayLinear / SorLinear) + plot_layout(guides = "collect") +plot_annotation(tag_levels = "A") & theme(legend.position = "right")

# Save the plot
ggsave("C:/", #User directory
       plot = Linear, width = 30, height = 30, units = "cm")

#### GAM ####
# Modelo GAM para variação temporal

library(dplyr)
library(mgcv)
library(DHARMa)
library(gratia)
library(lmtest)
library(ggplot2)
library(patchwork)

df_scaled_spa <- `Media espacial` %>%
  mutate(
    Locality = factor(trimws(Locality)),
    Type = factor(Type), Year_sc = as.numeric(scale(Year)),
    DHW_sc = as.numeric(scale(mean_DHW)),
    MHW_sc = as.numeric(scale(MHW_mean_intensity_cumulative)),
    SST_sc = as.numeric(scale(mean_sst))
  ) %>%
  droplevels()


gam_bray <- gam(
  Bray ~ s(Year, Locality, bs = "fs", k = 3) +
    s(DHW_sc, k = 3) +
    s(MHW_sc, k = 3) +
    s(SST_sc, k = 3),
  data = df_scaled_spa,
  method = "REML"
)

summary(gam_bray)

gam.check(gam_bray)

gam_sor <- gam(
  Sor ~ s(Year, Locality, bs = "fs", k = 3) +
    s(DHW_sc, k = 3) +
    s(MHW_sc, k = 3) +
    s(SST_sc, k = 3),
  data = df_scaled_spa,
  method = "REML"
)

summary(gam_sor)

gam.check(gam_sor)


set.seed(123)

sim_bray <- simulateResiduals(
  fittedModel = gam_bray,
  n = 1000
)

plot(sim_bray)

testUniformity(sim_bray)
testDispersion(sim_bray)
testOutliers(sim_bray)


set.seed(123)

sim_sor <- simulateResiduals(
  fittedModel = gam_sor,
  n = 1000
)

plot(sim_sor)

testUniformity(sim_sor)
testDispersion(sim_sor)
testOutliers(sim_sor)

# Predictions

pred_grid <- df_scaled_spa %>%
  group_by(Locality, Type) %>%
  summarise(min_year = min(Year),
            max_year = max(Year),
            DHW_sc = mean(DHW_sc, na.rm = TRUE),
            MHW_sc = mean(MHW_sc, na.rm = TRUE),
            SST_sc = mean(SST_sc, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Locality, Type) %>%
  do(
    expand.grid(
      Year = seq(.$min_year, .$max_year, length.out = 100),
      Locality = .$Locality,
      DHW_sc = .$DHW_sc,
      MHW_sc = .$MHW_sc,
      SST_sc = .$SST_sc
    )
  ) %>%
  ungroup()
pred <- predict(gam_bray, newdata = pred_grid, se.fit = TRUE, type = "response")
pred_grid$fit <- pred$fit
pred_grid$upper <- pred$fit + 1.96 * pred$se.fit
pred_grid$lower <- pred$fit - 1.96 * pred$se.fit
plot(pred)
df_year_mean_Sor <- df_scaled_spa %>%
  group_by(Year, Locality, Type) %>%
  summarise(
    Sor_mean = mean(Sor, na.rm = TRUE),
    Sor_sd = sd(Sor, na.rm = TRUE),     
    .groups = "drop"
  )
df_year_mean_Bray <- df_scaled_spa %>%
  group_by(Year, Locality, Type) %>%
  summarise(
    Bray_mean = mean(Bray, na.rm = TRUE),
    Bray_sd = sd(Bray, na.rm = TRUE),     
    .groups = "drop"
  )

# Spatial Graphics
# Temporal trends (GAM)

Sorjunto<-ggplot(pred_grid, aes(x = Year, y = fit, color = Locality, fill = Locality, linetype = Type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) + 
  geom_line(size = 1.2) +
  geom_point(
    data = df_year_mean_Sor,
    aes(x = Year, y = Sor_mean, color = Locality),
    size = 3,
    alpha = 0.9
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Year", y = "Bray (predicted)", 
       title = "Temporal Sørensen Variation per Locality (GAM)", linetype ="") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 12))

Bray<-ggplot(pred_grid, aes(x = Year, y = fit, color = Locality, fill = Locality, linetype = Type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) + 
  geom_line(size = 1.2) +
  geom_point(
    data = df_year_mean_Bray,
    aes(x = Year, y = Bray_mean, color = Locality),
    size = 3,
    alpha = 0.9
  ) +
  theme_minimal(base_size = 14) +
  labs(x = "Year", y = "Bray (predicted)", 
       title = "Temporal Bray-Curtis Variation per Locality (GAM)", linetype ="") +
  theme(legend.position = "right",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 12))

Brayy<- bray + theme(legend.position = "none")
GamBray<-(Brayy / Bray) + plot_layout(guides = "collect") +plot_annotation(tag_levels = "A") & theme(legend.position = "right")
GamSor<-(Jac / Sorjunto) + plot_layout(guides = "collect") +plot_annotation(tag_levels = "A") & theme(legend.position = "right")

# Save the plot
ggsave("C:/", #User directory
       plot = GamBray, width = 30, height = 30, units = "cm")
plot(combined_nmds)

ggsave("C:/", #User directory
       plot = GamSor, width = 30, height = 30, units = "cm")
#### glmm ####

library(ggplot2)
library(dplyr)
library(broom)
library(stringr)
library(lme4)
library(ggeffects)
library(car)
library(performance)
library(lmerTest)
library(mgcv)
library(MuMIn)
library(readxl)
library(tidyverse)
library(dplyr)
library(DHARMa)
library(glmmTMB)
library(writexl)

arquivo_excel <- "C:/" #User directory
abas <- c("Media SorBra Espacial", "SorBra espacial", "Temp", "Media espacial")

# Import data frame
walk(abas, function(sheet_name) {
  df <- read_excel(arquivo_excel, sheet = sheet_name)
  assign(paste0(sheet_name, ""), df, envir = .GlobalEnv)
})

df_spa <- `Media espacial`

df_spa <- df_spa %>%
  mutate( Locality = factor(trimws(Locality)),
          Type = factor(Type),
          Year = as.numeric(Year)
  ) %>%
  droplevels()


str(df_spa)

table(df_spa$Locality)
table(df_spa$Type)
table(df_spa$Year)


adjust_beta <- function(y) { n <- sum(!is.na(y))
y <- ifelse(
  y <= 0 | y >= 1,
  (y * (n - 1) + 0.5) / n,
  y
)
return(y)
}

df_spa <- df_spa %>%
  mutate(
    Bray_beta = adjust_beta(Bray),
    Sor_beta  = adjust_beta(Sor)
  )



range(df_spa$Bray_beta, na.rm = TRUE)
range(df_spa$Sor_beta, na.rm = TRUE)


df_spa <- df_spa %>%
  mutate(
    DHW_sc = as.numeric(scale(mean_DHW)),
    MHW_sc = as.numeric(scale(MHW_mean_intensity_cumulative)),
    SST_sc = as.numeric(scale(mean_sst))
  )

# BRAY-CURTIS

bray_baseline <- glmmTMB(
  Bray_beta ~ Type + (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

bray_MHW <- glmmTMB(
  Bray_beta ~ MHW_sc +
    (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

bray_DHW <- glmmTMB(
  Bray_beta ~ DHW_sc +
    (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

bray_SST <- glmmTMB(
  Bray_beta ~ SST_sc +
    (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_full <- glmmTMB(Bray_beta ~ SST_sc + DHW_sc+ MHW_sc
                     + (1|Locality), data = df_spa,
                     family = beta_family(link = "logit"),
                     REML = FALSE)
# SØRENSEN

sor_baseline <- glmmTMB(
  Sor_beta ~ Type + (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

sor_MHW <- glmmTMB(
  Sor_beta ~ MHW_sc + (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

sor_SST <- glmmTMB(
  Sor_beta ~ SST_sc + (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

sor_DHW <- glmmTMB(
  Sor_beta ~ DHW_sc + (1 | Locality),
  data = df_spa,
  family = beta_family(link = "logit"),
  REML = FALSE
)

sor_full <- glmmTMB(
  Sor_beta ~ SST_sc  +DHW_sc + MHW_sc
  + (1|Locality), data = df_spa, family = beta_family(link = "logit"),
  REML = FALSE
)

options(na.action = "na.fail")

AICc_bray <- model.sel(
  bray_baseline, bray_MHW, bray_full, bray_DHW,bray_SST
)

AICc_sor <- model.sel(
  sor_baseline, sor_MHW, sor_DHW, sor_SST, sor_full
)

AICc_bray
AICc_sor

summary(bray_baseline)
summary(sor_full)

coef_bray_MHW <- broom.mixed::tidy(
  bray_full, effects = "fixed",
  conf.int = TRUE
)

coef_sor_MHW <- broom.mixed::tidy(
  sor_full, effects = "fixed",
  conf.int = TRUE
)

coef_bray_MHW
coef_sor_MHW

resultado_modelos <- bind_rows(
  
  broom.mixed::tidy(
    bray_baseline, effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(Response = "Bray-Curtis",
           Model = "Baseline"),
  
  broom.mixed::tidy(
    bray_full, effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(Response = "Bray-Curtis",
           Model = "Thermal"),
  
  broom.mixed::tidy(
    sor_baseline, effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(Response = "Sørensen",
           Model = "Baseline"),
  
  broom.mixed::tidy(
    sor_full, effects = "fixed",
    conf.int = TRUE
  ) %>%
    mutate(Response = "Sørensen",
           Model = "Thermal")
)

resultado_modelos %>%
  filter(term %in% c("TypeIsland", "Thermal"))

tabela_AICc <- bind_rows(
  
  data.frame(
    Response = "Bray-Curtis", Model = "Baseline",
    AICc = AICc_bray$models[[1]]$AICc
  ),
  
  data.frame(
    Response = "Bray-Curtis", Model = "Baseline + Thermal",
    AICc = AICc_bray$models[[2]]$AICc
  ),
  
  data.frame(
    Response = "Sørensen", Model = "Baseline",
    AICc = AICc_sor$models[[1]]$AICc
  ),
  
  data.frame(
    Response = "Sørensen", Model = "Baseline + Thermal",
    AICc = AICc_sor$models[[2]]$AICc
  )
)

tabela_bray <- as.data.frame(AICc_bray) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(Response = "Bray-Curtis")

tabela_sor <- as.data.frame(AICc_sor) %>%
  tibble::rownames_to_column("Model") %>%
  mutate(Response = "Sørensen")

tabela_AICc <- bind_rows(
  tabela_bray,tabela_sor
) %>% select(Response, Model, df, logLik, AICc, delta, weight)

tabela_AICc

pred_MHW_bray <- ggpredict(
  bray_MHW,terms = "Thermal")

pred_MHW_sor <- ggpredict(
  sor_MHW,terms = "Thermal")

plot(pred_MHW_bray) +
  labs(x = "Thermal intensity (standardized)",
       y = "Predicted Bray-Curtis dissimilarity"
  ) +theme_classic()

plot(pred_MHW_sor) +
  labs(x = "Thermal intensity (standardized)",
       y = "Predicted Sørensen dissimilarity"
  ) +theme_classic()


set.seed(123)

res_bray_thermal <- simulateResiduals(
  fittedModel = bray_full,
  n = 1000
)

set.seed(123)

res_sor_thermal <- simulateResiduals(
  fittedModel = sor_full,
  n = 1000
)

# Bray-Curtis

plot(res_bray_thermal)

testUniformity(res_bray_thermal)
testDispersion(res_bray_thermal)
testOutliers(res_bray_thermal)

# Sørensen

plot(res_sor_thermal)

testUniformity(res_sor_thermal)
testDispersion(res_sor_thermal)
testOutliers(res_sor_thermal)

res_bray_base <- simulateResiduals(
  fittedModel = bray_baseline,
  n = 1000
)

set.seed(123)

res_sor_base <- simulateResiduals(
  fittedModel = sor_baseline,
  n = 1000
)

# Bray-Curtis

plot(res_bray_base)

testUniformity(res_bray_base)
testDispersion(res_bray_base)
testOutliers(res_bray_base)

# Sørensen

plot(res_sor_base)

testUniformity(res_sor_base)
testDispersion(res_sor_base)
testOutliers(res_sor_base)


plotResiduals(res_bray_thermal,
              df_spa$DHW_sc,
              main = "Bray-Curtis: DHW"
)

plotResiduals(res_bray_thermal,
              df_spa$MHW_sc,main = "Bray-Curtis: MHW"
)

plotResiduals(res_bray_thermal,
              df_spa$SST_sc,main = "Bray-Curtis: SST"
)

plotResiduals(res_bray_thermal,
              df_spa$Locality,main = "Bray-Curtis: Locality"
)

plotResiduals(res_sor_thermal,
              df_spa$DHW_sc,main = "Sørensen: DHW"
)

plotResiduals(res_sor_thermal,
              df_spa$MHW_sc,main = "Sørensen: MHW"
)

plotResiduals(res_sor_thermal,
              df_spa$SST_sc,main = "Sørensen: SST"
)

plotResiduals(res_sor_thermal,
              df_spa$Locality,main = "Sørensen: Locality"
)

AICc_bray_df <- as.data.frame(AICc_bray)
AICc_sor_df  <- as.data.frame(AICc_sor)

coef_bray <- broom.mixed::tidy(
  bray_full,
  effects = "fixed",
  conf.int = TRUE
)

coef_sor <- broom.mixed::tidy(
  sor_full,
  effects = "fixed",
  conf.int = TRUE
)

pred_bray <- as.data.frame(
  ggeffects::ggpredict(bray_full, terms = "DHW_sc")
)

pred_sor <- as.data.frame(
  ggeffects::ggpredict(sor_full, terms = "Type")
)
write_xlsx(
  list(
    AICc_Bray = AICc_bray_df,
    AICc_Sorensen = AICc_sor_df,
    Coef_Bray = coef_bray,
    Coef_Sorensen = coef_sor
  ),
  "C:/" # User directory
)

png("residuos_painel_dharma.png", width = 1200, height = 1600, res = 150)

par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))


#DHARMA Panel

plotQQunif(gam_bray, main = "Bray-Curtis GAM")
plotResiduals(gam_bray)
plotQQunif(gam_sor, main = "Sørensen GAM")
plotResiduals(gam_sor)


plotQQunif(bray_full, main = "Envrionmental Bray-Curtis GLMM")
plotResiduals(bray_full, main = "")
plotQQunif(sor_full, main = "Sørensen GLMM")
plotResiduals(sor_full, main = "")


plotQQunif(bray_baseline, main = "Bray-Curtis GLMM baseline")
plotResiduals(bray_baseline, main = "")
plotQQunif(res_bray_bio, main = "Biogeographical GLMM")
plotResiduals(res_bray_bio, main = "")

dev.off()

#### glmm biogeographic ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(DHARMa)
library(performance)
library(broom.mixed)
library(ggeffects)
library(MuMIn)
library(purrr)
library(openxlsx)

df_bio <- df_spa %>%
  mutate(
    Locality = factor(Locality),
    Type = factor(Type)
  )


bio_vars <- `Media espacial` %>%
  group_by(Locality, Type) %>%
  summarise(
    Richness  = mean(Richness, na.rm = TRUE),
    Isolation = mean(Isolation, na.rm = TRUE),
    Area      = mean(Area, na.rm = TRUE),
    Distance  = mean(Distance, na.rm = TRUE),
    Reef_30   = mean(Reef_30, na.rm = TRUE),
    Latitude  = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    .groups = "drop"
  )



bio_vars <- bio_vars %>%
  mutate(
    Isolation_sc = as.numeric(scale(Isolation)),
    Area_sc      = as.numeric(scale(Area)),
    Reef_30_sc   = as.numeric(scale(Reef_30)),
    Richness_sc  = as.numeric(scale(Richness)),
    Distance_sc = as.numeric(scale(Distance))
  )


df_bio <- df_bio %>%
  select(-any_of(c(
    "Isolation_sc",
    "Area_sc",
    "Reef_30_sc",
    "Richness_sc",
    "Distance_sc"
  ))) %>%
  left_join(
    bio_vars %>%
      select(
        Locality,
        Isolation_sc,
        Area_sc,
        Reef_30_sc,
        Richness_sc,
        Distance_sc
      ),
    by = "Locality"
  )


bio_cor <- bio_vars %>%
  select(
    Isolation_sc,
    Area_sc,
    Reef_30_sc,
    Richness_sc,
    Distance_sc
  ) %>%
  cor(
    method = "spearman",
    use = "complete.obs"
  )

round(bio_cor, 2)

library(corrplot)

corrplot(
  bio_cor,
  method = "number",
  type = "upper",
  tl.col = "black",
  tl.srt = 45
)



bray_bio_null <- glmmTMB(
  Bray_beta ~
    Type +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)
summary(bray_bio_null)


bray_bio_iso_mhw <- glmmTMB(
  Bray_beta ~
    Type +
    MHW_sc +
    Isolation_sc +
    MHW_sc:Isolation_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)

summary(bray_bio_iso_mhw)


bray_bio_dis_mhw <- glmmTMB(
  Bray_beta ~
    Type +
    MHW_sc +
    Distance_sc +
    MHW_sc:Distance_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_bio_reef_mhw <- glmmTMB(
  Bray_beta ~
    Type +
    MHW_sc +
    Reef_30_sc +
    MHW_sc:Reef_30_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_bio_iso_sst <- glmmTMB(
  Bray_beta ~
    Type +
    SST_sc +
    Isolation_sc +
    SST_sc:Isolation_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)

summary(bray_bio_iso)


bray_bio_dis_sst <- glmmTMB(
  Bray_beta ~
    Type +
    SST_sc +
    Distance_sc +
    SST_sc:Distance_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_bio_reef_sst <- glmmTMB(
  Bray_beta ~
    Type +
    SST_sc +
    Reef_30_sc +
    SST_sc:Reef_30_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_bio_iso_dhw <- glmmTMB(
  Bray_beta ~
    Type +
    DHW_sc +
    Isolation_sc +
    DHW_sc:Isolation_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)

summary(bray_bio_iso)


bray_bio_dis_dhw <- glmmTMB(
  Bray_beta ~
    Type +
    DHW_sc +
    Distance_sc +
    DHW_sc:Distance_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


bray_bio_reef_dhw <- glmmTMB(
  Bray_beta ~
    Type +
    DHW_sc +
    Reef_30_sc +
    DHW_sc:Reef_30_sc +
    (1 | Locality) +
    (1 | Year),
  data = df_bio,
  family = beta_family(link = "logit"),
  REML = FALSE
)


options(na.action = "na.fail")

AICc_bray_bio <- MuMIn::model.sel(
  bray_bio_null,
  bray_bio_iso_dhw,bray_bio_iso_mhw,bray_bio_iso_sst,
  bray_bio_dis_mhw,bray_bio_dis_dhw,bray_bio_dis_sst,
  bray_bio_reef_mhw,bray_bio_reef_dhw,bray_bio_reef_sst
)

AICc_sor_bio <- MuMIn::model.sel(
  sor_bio_null,
  sor_bio_iso,
  sor_bio_area,
  sor_bio_reef
)

AICc_bray_bio

best_bray_bio <- get.models(
  AICc_bray_bio,
  subset = 1
)[[1]]

summary(best_bray_bio)

check_convergence(best_bray_bio)
check_singularity(best_bray_bio)
check_collinearity(best_bray_bio)


set.seed(123)

res_bray_bio <- simulateResiduals(
  best_bray_bio,
  n = 1000
)

plot(res_bray_bio)

testUniformity(res_bray_bio)
testDispersion(res_bray_bio)
testOutliers(res_bray_bio)