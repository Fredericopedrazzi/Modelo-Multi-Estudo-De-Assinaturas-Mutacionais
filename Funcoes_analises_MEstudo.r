###########################################################
#    Manipulação dos dados para graficar as Exposições       #<---- código essencial para a classe que vai executar
###########################################################

#S <- length(dados$A_list)
#E_list_medio <- list()

#for (s in 1:S) {
    
#  E_list_medio[[s]] <- dados$E_list_medio[[s]]
      
#}  

#############################################
#     Histograma Assinaturas Mutacionais
############################################


Assinaturas_Mutacionais <- function(assinatura,
                                    cores,
                                    largura = 12,
                                    altura = 4,
                                    pausa = 0.5) {

  # Ajustar tamanho do gráfico no Jupyter
  options(repr.plot.width = largura,
          repr.plot.height = altura)

  # Converter para formato longo
  P_long <- melt(assinatura)
  colnames(P_long) <- c("Amostra", "Assinatura", "Valor")

  unique_assinaturas <- unique(P_long$Assinatura)

  # Salvar parâmetros gráficos
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(1, 1))

  for (ass in unique_assinaturas) {

    dados_ass <- P_long[P_long$Assinatura == ass, ]
    m <- max(dados_ass$Valor, na.rm = TRUE)

    barplot(
      dados_ass$Valor,
      col  = rep(cores, length.out = nrow(dados_ass)),
      ylim = c(0, m * 1.1),
      main = paste("Assinatura:", ass),
      ylab = "Valor",
      xlab = "Amostras",
      border = NA
    )

    abline(v = c(19.5, 38.5, 57.5, 77, 96), lty = 2)
    box()
    
    if (pausa > 0) Sys.sleep(pausa)
  }

  invisible(P_long)
}


extrair_diagonal_A <- function(A_list) {

  # Extrair a diagonal de cada matriz
  diags <- lapply(A_list, diag)

  # Empilhar em data.frame
  dfs <- as.data.frame(do.call(rbind, diags))

  # Nomear linhas e colunas
  rownames(dfs) <- paste0("Estudo ", seq_along(A_list))
  colnames(dfs) <- paste0("N", seq_len(ncol(dfs)))

  return(dfs)
}


#########################################
#            Boxplot Exposição
#########################################

Boxplot_exposicoes <- function(E_list,
                                        idx = 1,
                                        ylim = c(0, 0.075),
                                        titulo = "Distribuição das exposições por assinatura",
                                        base_size = 14) {

  suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(dplyr)
    library(RColorBrewer)
  })

  # Todas as assinaturas possíveis
  todas_assinaturas <- paste0(
    "Assinatura_",
    1:max(sapply(E_list, nrow))
  )

  # Matriz de exposições
  E_mat <- E_list[[idx]]

  rownames(E_mat) <- paste0("Assinatura_", 1:nrow(E_mat))
  colnames(E_mat) <- paste0("Paciente_", 1:ncol(E_mat))

  # Long format
  E_df <- melt(E_mat)
  colnames(E_df) <- c("Assinatura", "Paciente", "Exposicao")

  # Fatores com níveis fixos
  E_df$Assinatura <- factor(E_df$Assinatura, levels = todas_assinaturas)
  E_df$Paciente   <- factor(E_df$Paciente)

  # Remover zeros
  E_df <- dplyr::filter(E_df, Exposicao != 0)

  # Paleta dinâmica
  cores <- colorRampPalette(
    brewer.pal(8, "Set2")
  )(length(todas_assinaturas))

  # Plot
  p <- ggplot(E_df, aes(x = Assinatura, y = Exposicao, fill = Assinatura)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
    theme_minimal(base_size = base_size) +
    labs(
      title = titulo,
      x = "Assinaturas",
      y = "Valor de exposição"
    ) +
    scale_fill_manual(values = cores, drop = FALSE) +
    coord_cartesian(ylim = ylim) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  return(p)
}
 

# Ajusta tamanho do gráfico
options(repr.plot.width = 12, repr.plot.height = 8)


###########################################################
#              Boxplot Medida Math Score 
###########################################################
plot_box_MATHSCORE  <- function(Maior,
                                Menor,
                                coluna = 3,
                                titulo = "Math Score",
                                base_size = 14) {

  # Verificar se os inputs são dataframes
  if (!is.data.frame(Maior)) {
    stop("O primeiro argumento deve ser um data.frame (maior medida Math Score)")
  }
  if (!is.data.frame(Menor)) {
    stop("O segundo argumento deve ser um data.frame (menor medida Math Score)")
  }

  # Criar dataframes padronizados
  Maior_data <- data.frame(
    MATHSCORE = Maior[[coluna]],
    Group = "Estudo 2: Heterogênios"
  )

  Menor_data <- data.frame(
    MATHSCORE = Menor[[coluna]],
    Group = "Estudo 1: Homogênios"
  )

  # Combinar
  data <- rbind(Menor_data, Maior_data)

  # Plotagem
  library(ggplot2)
  m <- ggplot(data, aes(x = Group, y = MATHSCORE, fill = Group)) +
    geom_boxplot() +
    labs(
      title = titulo,
      x = "Grupos Math Score",
      y = "MATH SCORE Valor"
    ) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "none")

  return(m)
}


##############################################
#            Boxplot Pontuação TMB 
##############################################
plot_box_TMB  <- function(Maior,
                                Menor,
                                coluna = 4,
                                titulo = "TMB",
                                base_size = 14) {

  # Verificar se os inputs são dataframes
  if (!is.data.frame(Maior)) {
    stop("(maior medida TMB)")
  }
  if (!is.data.frame(Menor)) {
    stop("(menor medida Math Score)")
  }

  # Criar dataframes padronizados
  Maior_data <- data.frame(
    TMB = Maior[[coluna]],
    Group = "Estudo 2: Maiores medidas TMB"
  )

  Menor_data <- data.frame(
    TMB = Menor[[coluna]],
    Group = "Estudo 1: Menores medidas TMB"
  )

  # Combinar
  data <- rbind(Menor_data, Maior_data)

  # Plotagem
  library(ggplot2)
  t <- ggplot(data, aes(x = Group, y = TMB, fill = Group)) +
    geom_boxplot() +
    labs(
      title = titulo,
      x = "Grupos TMB",
      y = "Valor TMB"
    ) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "none")

  return(t)
}

##############################################
#            Boxplot Pontuação TMB 
##############################################


plot_box_FUMO  <- function(Maior,
                           Menor,
                           coluna = 2,
                           titulo = "Pack_Smoked_years",
                           base_size = 14) {
  
  # Verificar se os inputs são dataframes
  if (!is.data.frame(Maior)) {
    stop("(maior medida FUMO)")
  }
  if (!is.data.frame(Menor)) {
    stop("(menor medida Math FUMO)")
  }
  
  # Criar dataframes padronizados
  Maior_data <- data.frame(
    FUMO = Maior[[coluna]],
    Group = "Estudo 2: Maiores medidas FUMO"
  )
  
  Menor_data <- data.frame(
    FUMO = Menor[[coluna]],
    Group = "Estudo 1: Menores medidas FUMO"
  )
  
  # Combinar
  data <- rbind(Menor_data, Maior_data)
  
  # Plotagem
  library(ggplot2)
  f <- ggplot(data, aes(x = Group, y = FUMO, fill = Group)) +
    geom_boxplot() +
    labs(
      title = titulo,
      x = "Grupos FUMO",
      y = "Valor FUMO"
    ) +
    theme_minimal(base_size = base_size) +
    theme(legend.position = "none")
  
  return(f)
}

#######################################################
#               Função Barplot Exposição 
#######################################################

Barplot_exposição <- function(
  E_list, idx = 1,
  titulo = "Porcentagem de cada assinatura por paciente",
  base_size = 14
) {

  suppressPackageStartupMessages({
    library(reshape2)
    library(ggplot2)
    library(RColorBrewer)
  })

  # ====== Matriz de exposições ======
  E_mat <- E_list[[idx]]

  # ====== Reduzir barcode TCGA (PACIENTE) ======
  colnames(E_mat) <- sub(
    "(TCGA-[A-Za-z0-9]+-[A-Za-z0-9]+).*",
    "\\1",
    colnames(E_mat)
  )

  # ====== Nomes das assinaturas ======
  if (is.null(rownames(E_mat))) {
    rownames(E_mat) <- paste0("Assinatura_", 1:nrow(E_mat))
  }

  # ====== Long format ======
  E_df <- melt(E_mat)
  colnames(E_df) <- c("Assinatura", "Paciente", "Exposicao")

  # ====== Remover zeros ======
  E_df <- subset(E_df, Exposicao != 0)

  # ====== Percentual por paciente ======
  total_por_paciente <- aggregate(Exposicao ~ Paciente, data = E_df, sum)
  colnames(total_por_paciente)[2] <- "Total"

  E_df <- merge(E_df, total_por_paciente, by = "Paciente")
  E_df$Percentual <- E_df$Exposicao / E_df$Total * 100

  # ====== Paleta ======
  todas_assinaturas <- paste0(
    "Assinatura_",
    1:max(sapply(E_list, nrow))
  )

  pal <- colorRampPalette(brewer.pal(8, "Set2"))
  cores <- pal(length(todas_assinaturas))
  names(cores) <- todas_assinaturas

  cores_presentes <- cores[unique(E_df$Assinatura)]

  # ====== Plot ======
  ggplot(E_df, aes(x = Paciente, y = Percentual, fill = Assinatura)) +
    geom_bar(stat = "identity", color = "black", width = 1) +
    theme_minimal(base_size = base_size) +
    labs(
      title = titulo,
      x = "Pacientes (TCGA reduzido)",
      y = "Porcentagem da exposição",
      fill = "Assinatura"
    ) +
    scale_fill_manual(values = cores_presentes) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "top"
    )
}



##################################################################################
# Função que compara as assnaturas com assinaturas Cosmic, Distancia do cosseno 
##################################################################################
#comparar_assinaturas_COSMIC <- function(P_mat,
#                                        arquivo_COSMIC,
#                                        limiar_dist = 0.8) {
#
#  suppressPackageStartupMessages({
#    library(proxy)
#  })

  # -----------------------------
  # Ler COSMIC
  # -----------------------------
#  dt <- read.table(arquivo_COSMIC, header = TRUE)
#  dt_mat <- as.matrix(dt[, ])  # remove coluna de nomes

  # -----------------------------
  # Assinaturas inferidas
  # -----------------------------
#  P1 <- P_mat

#  if (is.null(colnames(P1))) {
#    colnames(P1) <- paste0("Sig", seq_len(ncol(P1)))
#  }

  # -----------------------------
  # Função distância do cosseno
  # -----------------------------
#  cosine_dist_proxy <- function(a, b) {
#    as.numeric(proxy::dist(t(a), t(b), method = "cosine"))
#  }

  # -----------------------------
  # Tabela de resultados
  # -----------------------------
#  resultado_dist <- data.frame(
#    Assinatura_P1 = colnames(P1),
#    Melhor_COSMIC = NA,
#    Distancia     = NA,
#    Status        = NA,
 #   stringsAsFactors = FALSE
 # )

  # -----------------------------
  # Loop principal
  # -----------------------------
#  for (j in seq_len(ncol(P1))) {

#    assinatura_j <- P1[, j]

#    distancias <- apply(
#      dt_mat,
#      2,
#      function(x) cosine_dist_proxy(assinatura_j, x)
#    )

#    melhor_idx  <- which.min(distancias)
#    melhor_dist <- distancias[melhor_idx]

#    resultado_dist$Melhor_COSMIC[j] <- colnames(dt_mat)[melhor_idx]
#    resultado_dist$Distancia[j]     <- melhor_dist

#    if (melhor_dist <= limiar_dist) {
#      resultado_dist$Status[j] <- "Detectada no COSMIC"
#    } else {
#      resultado_dist$Status[j] <- "Não encontrada"
#    }
#  }

#  return(resultado_dist)
#}


##################################################################################
# Função que compara as assnaturas com assinaturas Cosmic, Similaridade do cosseno 
##################################################################################

comparar_assinaturas_COSMIC <- function(P_mat,
                                        arquivo_COSMIC,
                                        limiar_sim = 0.8) {

  suppressPackageStartupMessages({
    library(proxy)
  })

  # -----------------------------
  # Ler COSMIC
  # -----------------------------
  dt <- read.table(arquivo_COSMIC, header = TRUE)
  dt_mat <- as.matrix(dt)

  # -----------------------------
  # Assinaturas inferidas
  # -----------------------------
  P1 <- P_mat

  if (is.null(colnames(P1))) {
    colnames(P1) <- paste0("Sig", seq_len(ncol(P1)))
  }

  # -----------------------------
  # Similaridade cosseno
  # -----------------------------
  cosine_sim_proxy <- function(a, b) {
    1 - as.numeric(proxy::dist(t(a), t(b), method = "cosine"))
  }

  # -----------------------------
  # Tabela de resultados
  # -----------------------------
  resultado_sim <- data.frame(
    Assinatura_P1 = colnames(P1),
    Melhor_COSMIC = NA,
    Similaridade  = NA,
    Status        = NA,
    stringsAsFactors = FALSE
  )

  # -----------------------------
  # Loop principal
  # -----------------------------
  for (j in seq_len(ncol(P1))) {

    assinatura_j <- P1[, j]

    similares <- apply(
      dt_mat,
      2,
      function(x) cosine_sim_proxy(assinatura_j, x)
    )

    melhor_idx  <- which.max(similares)
    melhor_sim  <- similares[melhor_idx]

    resultado_sim$Melhor_COSMIC[j] <- colnames(dt_mat)[melhor_idx]
    resultado_sim$Similaridade[j]   <- melhor_sim

    if (melhor_sim >= limiar_sim) {
      resultado_sim$Status[j] <- "Aceita no COSMIC"
    } else {
      resultado_sim$Status[j] <- "Não encontrada"
    }
  }

  return(resultado_sim)
}


###################################
#         Função Sobrevida
##################################

# ============================================
# 0. Pacotes necessários
# ============================================
library(survival)
library(ggplot2)

# ============================================
# 1. Função Kaplan-Meier
# ============================================
kaplan_meier = function (
    data,
    vars = list("feature"="feature", "event"="event", "time"="time"),
    title = "Kaplan-Meier (Survival)",
    xlab = "Time (days)",
    ylab = "Survival Probability",
    min_samples_by_feature = 0,
    min_features = 1,
    confidence.interval = FALSE
) {

    if (!(vars$feature %in% colnames(data))) {
        data$feature = "all"
    }

    data = data[,c(vars$feature, vars$event, vars$time)]
    colnames(data) = c("feature", "event", "time")

    data = na.omit(data)
    data = data[data$time >= 0, ]

    n_features = length(unique(data$feature))

    out = list(
        pval = 1,
        fit = NA,
        plot = ggplot() + theme_void()
    )

    if (n_features < min_features) return(out)
    if (!all(table(data$feature) >= min_samples_by_feature)) return(out)

    fit = survfit(Surv(time = time, event = event) ~ feature, data = data)

    strata = ""
    if (n_features == 1) {
        strata = paste0("all [", fit$n, "]")
    } else {
        names(fit$strata) = gsub("feature=", "", names(fit$strata))
        strata = rep(paste0(names(fit$strata), " [", fit$strata, "]"), fit$strata)
    }

    fit = data.frame(
        time = fit$time,
        surv = fit$surv,
        lower = fit$lower,
        upper = fit$upper,
        n.censor = fit$n.censor,
        strata = strata
    )

    fit = rbind(
        fit,
        data.frame(
            time = 0,
            surv = 1,
            lower = 1,
            upper = 1,
            n.censor = 0,
            strata = unique(strata)
        )
    )

    out$fit = fit

    censor_df = data.frame(
        time = fit$time[fit$n.censor > 0],
        surv = fit$surv[fit$n.censor > 0],
        strata = fit$strata[fit$n.censor > 0]
    )

    pval_label = ""
    if (n_features > 1) {
        sdiff = survdiff(Surv(time, event) ~ feature, data = data)
        pval = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
        
        out$pval = pval
        
        pval_label = ifelse(
            pval < 0.001,
            "Log-rank p < 0.001",
            paste0("Log-rank p = ", round(pval, 3))
        )
    }

    out$plot = ggplot(fit, aes(x = time, y = surv, color = strata, fill = strata)) +
        geom_step(linewidth = 1)

    if (confidence.interval) {
        out$plot = out$plot +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = 0)
    }

    out$plot = out$plot +
        geom_point(
            data = censor_df,
            aes(x = time, y = surv, color = strata),
            shape = 3,
            size = 2
        ) +
        scale_x_continuous(expand = c(0,0.01), name = xlab, n.breaks = 10) +
        scale_y_continuous(expand = c(0,0), limits = c(0,1), name = ylab,
                           breaks = seq(0,1,0.25)) +
        labs(title = title, color = "", fill = "") +
        annotate("text", x = max(fit$time)*0.1, y = 0.1, label = pval_label) +
        theme_minimal(base_size = 14) +
        theme(
            legend.position="top",
            panel.grid.minor = element_blank(),
            axis.line.x.bottom = element_line(linewidth = 0.5, color = "black"),
            axis.line.y.left = element_line(linewidth = 0.5, color = "black"),
            axis.text = element_text(colour = "black")
        )

    return(out)
}

#scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)), name = xlab, n.breaks = 10) +
#scale_y_continuous(expand = expansion(mult = c(0.02, 0.02)), limits = c(0,1), name = ylab,
#                   breaks = seq(0,1,0.25))



#======================
# Forest Plot
#======================


# USAGE:
# df = lung
# df$sex = as.character(df$sex)
# forest_plot(
#     df,
#     vars = list("feature"=c("age", "sex"), "event"="status", "time"="time"),
#     vars_must_have = c(),
#     title = "Lung",
#     plot = TRUE,
#     multivariate = FALSE,
#     sort_by_HR = TRUE,
#     standardize_continuous = TRUE
# )

# 2025-10-09: LI -> round to ceiling; LS -> round to floor

forest_plot = function (
    data,
    vars = list("feature"="feature", "event"="event", "time"="time"),
    vars_must_have = vars$feature, # used in filters
    title = "",
    plot = TRUE,
    filter_NA = TRUE,
    multivariate = TRUE,
    sort_by_HR = FALSE,
    show_label = TRUE,
    pval_max = 1.0,
    standardize_continuous = FALSE,
    collapse_factors = FALSE,
    min_factor_count = 5,
    interactions = FALSE # TRUE -> surv ~ A * B * C (variables interact with each other); FALSE -> surv ~ A + B + C
) {
    out = list(
        plot = ggplot() + theme_void(),
        fit = NA,
        pval = 1,
        c_index = NA,
        aic = NA,
        logLik = NA,
        removed_vars = 0,
        hazard = NA,
        diagnostics = c()
    )

    ############################################################
    
    data = data[,c(vars$event, vars$time, vars$feature)]
    colnames(data)[c(1,2)] = c("event", "time")

    ############################################################

    # force event column to 0/1 numeric
    data$event = as.numeric(as.character(data$event))
    if (all(is.na(data$event))) {
        out$diagnostics = c(out$diagnostics, "All event values are NA")
        return(out)
    }
    if (any(data$time < 0, na.rm = TRUE)) {
        out$diagnostics = c(out$diagnostics, "Column time has negative values")
        data$time[data$time < 0] = NA
    }

    ############################################################

    # force factor columns to character
    for (i in vars$feature) {
        if (is.factor(data[[i]])) {
            data[[i]] = as.character(data[[i]])
        }
    }

    ############################################################
    ############################################################
    ############################################################

    # NA's

    if (filter_NA) {
        for (i in vars$feature) {
            all_na = all(is.na(data[[i]]))
            if (all_na) {
                if (i %in% vars_must_have) {
                    out$diagnostics = c(out$diagnostics, paste0("there is a mandatory column that only has NA's (", i, ")"))
                    return(out)
                } else {
                    # remove column
                    out$diagnostics = c(out$diagnostics, paste0("Column ", i, " was removed (All values are NA)"))
                    data = data[,!(colnames(data) %in% i)]
                    vars$feature = vars$feature[!(vars$feature==i)]
                    out$removed_vars = out$removed_vars + 1
                }
            }
        }
        
        data = na.omit(data)
    }

    ############################################################
    ############################################################
    ############################################################

    # check for infinite values

    for (i in vars$feature) {
        if (any(is.infinite(data[,i]))) stop(paste0("Infinite values in column: ", i))
    }

    ############################################################
    ############################################################
    ############################################################

    # standardize continuous

    if (standardize_continuous) {
        for (i in vars$feature) {
            if (is.numeric(data[[i]]) & length(unique(data[[i]])) > 2) {
                data[[i]] = as.numeric(scale(data[[i]]))
                out$diagnostics = c(out$diagnostics, paste0("Standardized: ", i))
            }
        }
    }

    ############################################################
    ############################################################
    ############################################################

    # collapse rare levels
        
    if (collapse_factors) {
        for (i in vars$feature) {
            if (is.character(data[[i]])) {
                tab = table(data[[i]])
                rare_levels = names(tab)[tab < min_factor_count]
                if (length(rare_levels) > 0) {
                    data[[i]] = forcats::fct_lump_min(factor(data[[i]]), min = min_factor_count)
                    data[[i]] = as.character(data[[i]])
                    out$diagnostics = c(out$diagnostics, paste("Rare levels collapsed in: ", i))
                }
            }
        }
    }
    
    ############################################################
    ############################################################
    ############################################################

    # check if data has at least 2 categories
    for (i in vars$feature) {
        n_events = length(unique(data[[i]]))
        
        if (n_events < 2) {
            if (i %in% vars_must_have) {
                out$diagnostics = c(out$diagnostics, paste("there is a mandatory column that has only one event (", i, ")"))
                return(out)
            }
            
            # remove column
            out$diagnostics = c(out$diagnostics, paste("Column ", i, " was removed (Less than 2 categories)"))
            data = data[,!(colnames(data) %in% i)]
            vars$feature = vars$feature[vars$feature != i]
            out$removed_vars = out$removed_vars + 1
        }
    }

    if (length(vars$feature) == 0) {
        out$diagnostics = c(out$diagnostics, "All columns was removed")
        return(out)
    }

    ############################################################
    ############################################################
    ############################################################

    # check if data has at least 1 event

    # all data
    has_event = sum(data$event, na.rm = TRUE) > 0
    if (!has_event) {
        out$diagnostics = c(out$diagnostics, "There are no events remaining")
        return(out)
    }

    # by feature (only needed by coxph, it needs to have at least 1 event by category)
    for (i in vars$feature) {
        if (is.character(data[[i]])) {
            tab = table(data[[i]], data$event)
            if (any(tab[, "1"] == 0)) {
                if (i %in% vars_must_have) {
                    out$diagnostics = c(out$diagnostics, paste("there is a category of a mandatory column that has no events (", i, ")"))
                    return(out)
                }

                data = data[,!(colnames(data) %in% i)]
                vars$feature = vars$feature[vars$feature != i]
                out$removed_vars = out$removed_vars + 1
                out$diagnostics = c(out$diagnostics, paste("Column ", i, " was removed (No events)"))
            }
        }
    }
    
    if (length(vars$feature) == 0) {
        out$diagnostics = c(out$diagnostics, "All columns was removed")
        return(out)
    }
    
    ############################################################
    ############################################################
    ############################################################

    # correct categorical variable names
    # transform character to factor
    for (i in vars$feature) {
        if (is.character(data[[i]])) {
            data[[i]] = paste0("\n[", data[[i]], "]")
            data[[i]] = factor(data[[i]])
        }
    }
    
    ############################################################
    ############################################################
    ############################################################

    if (multivariate) { # MULTIVARIATE
        # cox model
        if (interactions) {
            form = as.formula(paste0("Surv(time, event) ~ ", paste0(vars$feature, collapse = " * ")))
        } else {
            form = as.formula(paste0("Surv(time, event) ~ ", paste0(vars$feature, collapse = " + ")))
        }
        suppressWarnings({
            fit = coxph(
                form,
                data = data
            )
        })
        
        s = summary(fit)

        ########################################################

        # # check cox_model result (merge ou exclude categories to solve some of the issues)
        confint_upper = s$conf.int[,"upper .95"]
        confint_lower = s$conf.int[,"lower .95"]
    
        # has inf in conf.int (low qty in some or multiple groups)
        has_inf_CI = any(is.infinite(confint_upper)) || any(is.infinite(confint_lower))
        
        has_na_CI = any(is.na(confint_upper)) || any(is.na(confint_lower)) # has NA's in conf.int
        has_na_coef = any(is.na(s$coefficients[,"coef"])) # has NA in coef
        
        if (has_inf_CI | has_na_CI | has_na_coef) return(out)
        # if (has_inf_CI) return(out)
        
        ########################################################

        if (plot) {
            df_plot = data.frame(
                feature = rownames(s$coef),
                HR = s$coef[, "exp(coef)"],
                LI = s$conf.int[, "lower .95"],
                LS = s$conf.int[, "upper .95"],
                p = s$coef[, "Pr(>|z|)"]
            )

            # append baseline
            for (i in vars$feature) {
                if (is.factor(data[[i]])) {
                    pos = grep(i, gsub("(.*)\n.*", "\\1", df_plot$feature))[1]
                    if (pos > 1) {
                        df_plot = rbind(
                            df_plot[1:(pos-1),],
                            data.frame(
                                feature = paste0(i, levels(data[[i]])[1]),
                                HR = 1,
                                LI = 1,
                                LS = 1,
                                p = NA
                            ),
                            df_plot[pos:nrow(df_plot),]
                        )
                    }
                }
            }

            df_plot$feature = factor(df_plot$feature, levels = rev(df_plot$feature))

            # pval max
            df_plot = df_plot[is.na(df_plot$p) | df_plot$p<=pval_max,]
            if (nrow(df_plot)==0) return(out)

            # append label
            if (show_label) {
                df_plot$label = paste0(
                    ifelse(
                        df_plot$p<0.001,
                        "p<0.001",
                        paste0("p=",round(df_plot$p, 3))
                    ), "\n",
                        "[", ceiling(df_plot$LI*100)/100, "~", floor(df_plot$LS*100)/100, "]"
                )
                df_plot$label_pos = max(df_plot$LS, na.rm = TRUE)
                df_plot$label = ifelse(is.na(df_plot$p), "", df_plot$label)
            }

            # sort by HR
            if (sort_by_HR) {
                df_plot$feature = with(df_plot, reorder(feature, HR))
            }

            # PLOT
            p = ggplot(df_plot, aes(x = feature, y = HR, ymin = LI, ymax = LS))

            # add pval
            if (show_label) {
                p = p + geom_text(mapping=aes(x=feature, y=label_pos, label=label), nudge_y = 0.2)
            }
            
            p = p +
                geom_pointrange(size = 0.8) +
                geom_hline(yintercept = 1, linetype = "dashed") +
                coord_flip() +
                scale_y_log10() +
                theme_minimal(base_size = 20) +
                labs(
                    x = "",
                    y = "Hazard ratio",
                    title = title,
                    caption = paste0(
                        "# Events: ", s$nevent,
                        "; Global p-value (Log-Rank): ", ifelse(s$logtest["pvalue"]<=0.001, "<= 0.001", round(s$logtest["pvalue"], 3)),
                        "\nAIC: ", round(AIC(fit), 3),
                        "; Concordance Index: ", round(s$concordance[1], 3)
                    )
                )
            
            out$plot = p
        }
        
        ########################################################
        
        out$fit = fit
        out$pval = s$logtest["pvalue"]
        out$c_index = unname(s$concordance[1])
        out$aic = extractAIC(fit)[2]
        out$logLik = as.numeric(logLik(fit))
        out$hazard = exp(coef(fit))
        
    ############################################################
    ############################################################
    ############################################################
        
    } else { # UNIVARIATE

        df_plot = lapply(vars$feature, function(var) {
            f = as.formula(paste("Surv(time, event) ~", var))
            m = coxph(f, data = data)
            s = summary(m)
            df_aux = data.frame(
                feature = rownames(s$coef),
                HR = s$coef[, "exp(coef)"],
                LI = s$conf.int[, "lower .95"],
                LS = s$conf.int[, "upper .95"],
                p = s$coef[, "Pr(>|z|)"],
                nevent = s$nevent
            )

            if (is.factor(data[[var]])) {
                rbind(
                    data.frame(
                        feature = paste0(var, levels(data[[var]])[1]),
                        HR = 1,
                        LI = 1,
                        LS = 1,
                        p = NA,
                        nevent = s$nevent
                    ),
                    df_aux
                )
            } else {
                df_aux
            }
        })

        df_plot = bind_rows(df_plot)

        # filter analysis
        df_plot = df_plot[!is.infinite(df_plot$LS),]
        if (nrow(df_plot)==0) return(out)

        # pval max
        df_plot = df_plot[is.na(df_plot$p) | df_plot$p<=pval_max,]
        if (nrow(df_plot)==0) return(out)

        # append label
        if (show_label) {
            df_plot$label = paste0(
                ifelse(
                    df_plot$p<0.001,
                    "p<0.001",
                    paste0("p=",round(df_plot$p, 3))
                ), "\n",
                    "[", ceiling(df_plot$LI*100)/100, "~", floor(df_plot$LS*100)/100, "]"
            )
            df_plot$label_pos = max(df_plot$LS, na.rm = TRUE)
            df_plot$label = ifelse(is.na(df_plot$p), "", df_plot$label)
        }

        df_plot$feature = factor(df_plot$feature, levels = rev(df_plot$feature))

        # sort by HR
        if (sort_by_HR) {
            df_plot$feature = with(df_plot, reorder(feature, HR))
        }

        out$fit = df_plot
        rownames(out$fit) = NULL
        
        ########################################################

        if (plot) {
            p = ggplot(df_plot, aes(x = feature, y = HR, ymin = LI, ymax = LS))

            # add pval
            if (show_label) {
                p = p + geom_text(mapping=aes(x=feature, y=label_pos, label=label), nudge_y = 0.2)
            }
            
            p = p +
                geom_pointrange(size = 0.8) +
                geom_hline(yintercept = 1, linetype = "dashed") +
                coord_flip() +
                scale_y_log10() +
                theme_minimal(base_size = 20) +
                labs(
                    x = "",
                    y = "Hazard ratio",
                    title = title,
                    caption = paste0(
                        "# Events: ", df_plot$nevent[1]
                    )
                )

            out$plot = p
        }
    }
    
    ############################################################
    ############################################################
    ############################################################

    if (plot) {
        out$plot = out$plot +
            theme(
                plot.title = element_text(hjust = 0.5),
                panel.grid.minor = element_blank(),
                # text = element_text(face="bold", colour = "black"),

                panel.grid.major = element_blank(),
                # panel.grid.minor = element_blank(),

                axis.line.x.bottom = element_line(color = "grey"),

                axis.ticks = element_line(color = "grey"),
                axis.ticks.length = unit(5, "pt"),
            )            
    }    

    ############################################################
    ############################################################
    ############################################################
    
    return(out)
}




#===================================================================
# Função para coluna de exposição como variavel continua na tablela
# E também grafico Klapan dividido em quartis de exposição
#===================================================================

km_quartis_exposure <- function(
  E_list_medio,
  which_list,
  Bcr,
  dados_survival
) {

  # ===== Exposição =====
  nomes_raw <- E_list_medio[[which_list]][Bcr, ]

  final_df <- data.frame(
    TCGA = sub(
      ".*(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}).*",
      "\\1",
      names(nomes_raw)
    ),
    exposure = as.numeric(nomes_raw),
    stringsAsFactors = FALSE
  )

  final_df <- final_df[!is.na(final_df$TCGA), ]

  # ===== Clínico =====
  dados <- dados_survival[
    !is.na(dados_survival$feature) &
    !is.na(dados_survival$time),
  ]

  extract_patient_barcode <- function(x) {
    pat <- "(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4})"
    regmatches(x, regexpr(pat, x, perl = TRUE))
  }

  dados$TCGA <- extract_patient_barcode(dados$bcr_patient_barcode)

  dados <- merge(dados, final_df, by = "TCGA", all = FALSE)

# ===== Checagem de variabilidade =====
if(length(unique(dados$exposure)) < 4){
  message("Menos de 4 valores distintos. Quartis não confiáveis.")
  return(NULL)
}

    
 # ===== Quartis =====
  dados$quartil <- factor(
    ntile(dados$exposure, 4),
    levels = 1:4,
    labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")
  )

  # ===== Evento =====
  dados$event <- ifelse(
    tolower(dados$vital_status) == "dead",
    1, 0
  )

  # ===== KM =====
  surv_obj <- Surv(dados$time, dados$event)
  fit_km <- survfit(surv_obj ~ quartil, data = dados)

  plot(
    fit_km,
    col  = c("#2166ac", "#67a9cf", "#fdae61", "#b2182b"),
    lwd  = 2,
    main = paste("KM – Lista", which_list, "| Assinatura", Bcr)
  )

  legend(
    "bottomleft",
    legend = levels(dados$quartil),
    col    = c("#2166ac", "#67a9cf", "#fdae61", "#b2182b"),
    lwd    = 2,
    bty    = "n"
  )

  # ===== Log-rank =====
  lr <- survdiff(surv_obj ~ quartil, data = dados)

  p_val <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

  logrank_table <- data.frame(
    Quartil  = names(lr$n),
    N        = lr$n,
    Observed = lr$obs,
    Expected = lr$exp
  )

  # ===== Retorno =====
  return(list(
    fit_km        = fit_km,
    logrank_obj   = lr,
    logrank_table = logrank_table,
    chisq         = lr$chisq,
    p_value       = p_val,
    dados         = dados
  ))
}


#=====================================================================================
# Função para coluna de exposição como variavel continua na tabela para forest plot
#=====================================================================================


run_forest_signature <- function(
  E_list_medio,
  which_list,
  Bcr,
  dados_survival,
  title_prefix = "Forest Plot"
) {

  # ==============================
  # 1. Exposição
  # ==============================
  nomes_raw <- E_list_medio[[which_list]][Bcr, ]

  final_df <- data.frame(
    TCGA = sub(
      ".*(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}).*",
      "\\1",
      names(nomes_raw)
    ),
    exposure = as.numeric(nomes_raw),
    stringsAsFactors = FALSE
  )

  final_df <- final_df[!is.na(final_df$TCGA), ]

  # ==============================
  # 2. Clínico
  # ==============================
  dados <- dados_survival[
    !is.na(dados_survival$feature) &
    !is.na(dados_survival$time),
  ]

  extract_patient_barcode <- function(x) {
    pat <- "(TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4})"
    regmatches(as.character(x),
               regexpr(pat, as.character(x), perl = TRUE))
  }

  dados$TCGA <- extract_patient_barcode(dados$bcr_patient_barcode)

  # ==============================
  # 3. Merge
  # ==============================
  dados <- merge(dados, final_df, by = "TCGA", all = FALSE)

  # ===== Checagem de variabilidade =====
if(length(unique(dados$exposure)) < 4){
  message("Assinatura não existe para essa Lista (Estudo).")
  return(NULL)
}
  

  # ==============================
  # 4. Evento
  # ==============================
  dados$event <- ifelse(
    tolower(dados$vital_status) == "dead",
    1, 0
  )

  # ==============================
  # 5. Forest – exposição contínua
  # ==============================
  res_cont <- forest_plot(
    data = dados,
    vars = list(
      feature = "exposure",
      event   = "event",
      time    = "time"
    ),
    title = paste(
      title_prefix,
      "| Lista", which_list,
      "| Assinatura", Bcr,
      "– Exposição como variável contínua"
    ),
    multivariate = FALSE,
    standardize_continuous = TRUE,
    plot = TRUE
  )

  # ==============================
  # 6. Forest – quartis
  # ==============================
  dados$quartil <- factor(
    ntile(dados$exposure, 4),
    levels = 1:4,
    labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")
  )

  dados$quartil <- as.character(dados$quartil)

  res_quartil <- forest_plot(
    data = dados,
    vars = list(
      feature = "quartil",
      event   = "event",
      time    = "time"
    ),
    title = paste(
      title_prefix,
      "| Lista", which_list,
      "| Assinatura", Bcr,
      "– Quartis de exposição"
    ),
    multivariate = FALSE,
    plot = TRUE
  )

  # ==============================
  # 7. Retorno organizado
  # ==============================
  return(list(
    continuous = list(
      plot        = res_cont$plot,
      fit         = res_cont$fit,
      diagnostics = res_cont$diagnostics
    ),
    quartile = list(
      plot        = res_quartil$plot,
      fit         = res_quartil$fit,
      diagnostics = res_quartil$diagnostics
    ),
    dados = dados
  ))
}
