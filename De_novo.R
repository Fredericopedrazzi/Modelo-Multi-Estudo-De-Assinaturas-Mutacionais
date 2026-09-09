library(stringi)  
library(Rcpp)
library(RcppArmadillo)


counts_list <- lapply(data, function(df) {
  numeric_cols <- df[, sapply(df, is.numeric), drop = FALSE]
  as.matrix(numeric_cols)
})

num_jotas <- sapply(counts_list, nrow)

# ---------------------------
# Hiper e iniciais
# ---------------------------
Ap <- matrix(rexp(i * n_signatures, rate = lp), i, n_signatures)
Bp <- matrix(rgamma(i * n_signatures, shape = ap, rate = bp), i, n_signatures)

P  <- matrix(0, nrow = i, ncol = n_signatures)
new_Ap <-  matrix(0, nrow = i, ncol = n_signatures)
alphaP <- matrix(0, nrow = i, ncol = n_signatures)
rhoP   <- matrix(0, nrow = i, ncol = n_signatures)

# Se existir P fixo (somente se há colunas fixas)
if (exists("dados") && !is.null(dados$P) && length(fix_colsP) > 0) {
  P[, fix_colsP] <- dados$P[, fix_colsP, drop = FALSE]
}

# Inicialização de P ~ Gamma(Ap+1, Bp)
for (m in 1:n_signatures) {
  if (m %in% fix_colsP) next
  for (k in 1:i) {
    P[k, m] <- rgamma(1, shape = Ap[k, m] + 1, rate = Bp[k, m])
  }
}


# ---------------------------
# Estruturas por estudo
# ---------------------------
M_list <- list()
E_list <- list()
Z_list <- list()
Fi_list <- list()
Ae_list <- list()
Be_list <- list()
PE_list <- list()
j_list <- list()
new_Ae_list <- list()
alphaE_list <- list()
rhoE_list <- list()
A_list <- list()

for (s in 1:S) {
  j_list[[s]] <- num_jotas[s]
  Z_list[[s]] <- array(0, dim = c(i, j_list[[s]], n_signatures))
  E_list[[s]] <- matrix(0, nrow = n_signatures, ncol = j_list[[s]])
  Fi_list[[s]] <- array(0, dim = c(i, j_list[[s]], n_signatures))
  new_Ae_list[[s]]  <- matrix(0, nrow = n_signatures, ncol = j_list[[s]])
  alphaE_list[[s]]  <- matrix(0, nrow = n_signatures, ncol = j_list[[s]])
  rhoE_list[[s]]    <- matrix(0, nrow = n_signatures, ncol = j_list[[s]])
  A_list[[s]] <- matrix(0, nrow = n_signatures, ncol = n_signatures)
  
  # Se tiver E/A fixos
  if (exists("dados") && !is.null(dados$E_list) && length(dados$E_list) >= s &&
      !is.null(dados$E_list[[s]]) && length(fix_colsE) > 0) {
    E_list[[s]][fix_colsE, ] <- dados$E_list[[s]][fix_colsE, , drop = FALSE]
  }
  if (exists("dados") && !is.null(dados$A_list) && length(dados$A_list) >= s &&
      !is.null(dados$A_list[[s]]) && length(fix_colsA) > 0) {
    A_list[[s]][fix_colsA, fix_colsA] <- dados$A_list[[s]][fix_colsA, fix_colsA, drop = FALSE]
  }
}

# ---------------------------
# Inicialização por estudo
# ---------------------------
for (s in 1:S) {
  
  # Hiper de E
  Ae_list[[s]] <- matrix(rexp(n_signatures * j_list[[s]], rate = le), n_signatures, j_list[[s]])
  Be_list[[s]] <- matrix(rgamma(n_signatures * j_list[[s]], shape = ae, rate = be), n_signatures, j_list[[s]])
  
  # E ~ Gamma
  for (m in 1:n_signatures) {
    if (m %in% fix_colsE) next
    for (g in 1:j_list[[s]]) {
      E_list[[s]][m, g] <- rgamma(1, shape = Ae_list[[s]][m, g] + 1, rate = Be_list[[s]][m, g])
    }
  }
  
  # M = t(counts) (k x g)
  M_list[[s]] <- t(counts_list[[s]])
  
  # PE e Fi (proteger com eps)
  PE_list[[s]] <- P %*% E_list[[s]]
  for (m in 1:n_signatures) {
    num <- (P[, m, drop = FALSE] %*% E_list[[s]][m, , drop = FALSE])
    Fi_list[[s]][, , m] <- num / PE_list[[s]]
  }
  
  # Z via multinomial
  for (k in 1:i) {
    Z_list[[s]][k, , ] <- t(sapply(1:j_list[[s]], function(g) {
      prob <- as.vector(Fi_list[[s]][k, g, ])
      if (anyNA(prob)) prob[is.na(prob)] <- 0
      s_prob <- sum(prob)
      if (s_prob <= 0) return(rep(0, length(prob)))
      prob <- prob / s_prob
      rmultinom(1, size = M_list[[s]][k, g], prob = prob)
    }))
  }
  
  # Chute inicial de A (usa probf de fato)
  X <- pmax(P %*% E_list[[s]], eps)
  log_X <- log(X)
  lgamma_term <- lgamma(M_list[[s]] + 1)
  
  prob1 <-   sum(-X + M_list[[s]] * log_X - lgamma_term)
  prob2 <-   sum(-X + M_list[[s]] * log_X - lgamma_term)
  
  denom <- prob1 + prob2
  probf <- if (is.finite(denom) && denom != 0) prob1 / denom else 0.5
  probf <- max(min(probf, 1), 0)
  
  for (m in 1:n_signatures) {
    if (m %in% fix_colsA) next
    A_list[[s]][m, m] <- rbinom(1, 1, prob = 0.5)  #probf
  }
}

# ---------------------------
# Esquema de temperamento
# ---------------------------
incremento_y <- 1 / 34
y <- -2
psi <- 10^y

decressimo_l_1 <- 0.05
decressimo_l_2 <- 0.0048
decressimo_l_3 <- 0.00045
l <- 1

# ---------------------------

W_list <- vector("list", length(M_list))

for (h in seq_along(M_list)) {
  
  x <- M_list[[h]]
  
  # Soma das colunas
  soma_col <- colSums(x)
  
  # Cria matriz replicando as somas
  mat_temp <- matrix(soma_col,
                     nrow = i,
                     ncol = ncol(x),
                     byrow = TRUE)
  
  # Calcula a mediana da matriz
  med_val <- median(mat_temp)
  
  # Substitui valores positivos pela mediana
  mat_temp[mat_temp > 0] <- med_val
  
  W_list[[h]] <- mat_temp
}

#Rcpp::sourceCpp("teste.cpp")  
Rcpp::sourceCpp("Gibbs_Sampler_D.cpp")  

# j_list como vetor inteiro
j_list <- as.integer(num_jotas)

# Chamada do burn_in
resultado_1 <- burn_in_8(
  num_iterations, P, Ap, Bp, E_list, Ae_list, Be_list, A_list, Z_list, Fi_list,
  M_list, W_list, S, i, n_signatures, j_list,
  fix_colsP, fix_colsA, fix_colsE, le, 
  incremento_y = 1/34,
  decressimo_l_1 = 0.05,
  decressimo_l_2 = 0.0048,
  decressimo_l_3 = 0.00045
)   


find_top5_frequent_combinations <- function(A_store, n_signatures) {
  combination_count <- list()
  combination_iters <- list()
  
  num_iterations <- length(A_store[[1]])  # número de iterações armazenadas
  S <- length(A_store)                    # número de estudos
  
  for (iter in seq_len(num_iterations)) {
    combo_key_parts <- c()
    
    for (s in seq_len(S)) {
      current_matrix <- A_store[[s]][[iter]]
      matrix_key <- paste(as.vector(current_matrix), collapse = "_")
      combo_key_parts <- c(combo_key_parts, matrix_key)
    }
    
    combo_key <- paste(combo_key_parts, collapse = "|")  # Delimitador entre estudos
    
    if (!combo_key %in% names(combination_count)) {
      combination_count[[combo_key]] <- 1
      combination_iters[[combo_key]] <- list(iter)
    } else {
      combination_count[[combo_key]] <- combination_count[[combo_key]] + 1
      combination_iters[[combo_key]] <- c(combination_iters[[combo_key]], iter)
    }
  }
  
  # Top 5 combinações mais frequentes
  sorted_keys <- names(sort(unlist(combination_count), decreasing = TRUE))
  top_keys <- head(sorted_keys, 5)
  
  top_combinations <- list()
  
  for (k in top_keys) {
    combo_parts <- strsplit(k, "\\|")[[1]]
    matrices_list <- list()
    
    for (s in seq_len(length(combo_parts))) {
      mat_vals <- as.numeric(unlist(strsplit(combo_parts[s], "_")))
      matrices_list[[s]] <- matrix(mat_vals, nrow = n_signatures, ncol = n_signatures)
    }
    
    top_combinations[[k]] <- list(
      matrices = matrices_list,
      frequency = combination_count[[k]],
      iterations = combination_iters[[k]]
    )
  }
  
  return(top_combinations)
}

# === Exemplo de uso com sua saída ===
A_store <- resultado_1$A_store
top5_combos <- find_top5_frequent_combinations(A_store, n_signatures)

idx <- 1
for (combo_key in names(top5_combos)) {
  #  cat("Rank", idx, "- Frequência:", top5_combos[[combo_key]]$frequency, "\n")
  
  for (s in seq_along(top5_combos[[combo_key]]$matrices)) {
    #    cat("Estudo s =", s, "\n")
    #  print(top5_combos[[combo_key]]$matrices[[s]])
    #    cat("\n")
  }
  
  idx <- idx + 1
  #  cat("------------------------------------------------------------\n\n")
}  


#resultado_1$P_store         # cubo com todos os P armazenados
#resultado_1$A_store         # lista de listas de A armazenados
#resultado_1$E_store         # lista de listas de E armazenados
#resultado_1$store_range     # intervalo [start_store, end_store] 


# Função para pegar P e E da PRIMEIRA aparição de uma combinação de A mais frequente
analyze_top_combination_first <- function(resultado, top5_combos, rank_choice = 1) {
  
  P_store <- resultado$P_store
  A_store <- resultado$A_store
  E_store <- resultado$E_store
  store_range <- resultado$store_range
  
  S <- length(A_store)            # número de estudos
  
  # Ordenar pela frequência
  ordered_keys <- names(sort(sapply(top5_combos, function(x) x$frequency), decreasing = TRUE))
  
  if (rank_choice > length(ordered_keys)) {
    stop("Rank escolhido fora do intervalo das top combinações disponíveis.")
  }
  
  combo_key  <- ordered_keys[rank_choice]
  combo_info <- top5_combos[[combo_key]]
  
  # Iterações relativas onde essa combinação apareceu
  iteracoes_rel <- unlist(combo_info$iterations)
  
  # A primeira aparição é simplesmente:
  iter_primeira_rel <- iteracoes_rel[1]
  
  # Converter para índice absoluto
  iter_primeira_abs <- iter_primeira_rel + store_range[1] - 1
  
  # Selecionar o P e E dessa iteração
  P_primeiro <- P_store[,,iter_primeira_rel]
  
  E_primeiro <- lapply(seq_len(S), function(s) {
    E_store[[s]][[iter_primeira_rel]]
  })
  
  return(list(
    rank = rank_choice,
    combo_key = combo_key,
    frequencia_total = combo_info$frequency,
    primeira_iter_rel = iter_primeira_rel,
    primeira_iter_abs = iter_primeira_abs,
    
    A_primeiro = combo_info$matrices,   # A do par
    
    P_primeiro = P_primeiro,
    E_primeiro = E_primeiro
  ))
}    

#------------------------------Segunda Parte ------------------------------------------------------

# -----------------------------
# 0) Descobrir ranks válidos
# -----------------------------
n_top <- length(top5_combos)  # número de combinações disponíveis
rank_choices <- 1:n_top
cat("Rank choices disponíveis:", rank_choices, "\n")

# -----------------------------
# 1) Preparar armazenamento
# -----------------------------
bic_values <- numeric(length(rank_choices))
resultados_lista <- vector("list", length(rank_choices))

# -----------------------------
# 2) Loop automático sobre ranks válidos
# -----------------------------
for (r in seq_along(rank_choices)) {
  
  rank_choice <- rank_choices[r]
  
  # -----------------------------
  # 2a) Obter resultado para este rank
  # -----------------------------
  resultado_par <- analyze_top_combination_first(
    resultado = resultado_1,
    top5_combos = top5_combos,
    rank_choice = rank_choice
  )
  
  # -----------------------------
  # 2b) Determinar número de assinaturas real para este rank
  # -----------------------------
  n_signatures_rank <- ncol(resultado_par$P_primeiro)  # antes de remover zeros
  
  # Construir listas E, A e P
  E_list <- lapply(1:S, function(s) {
    vec <- resultado_par$E_primeiro[[s]]
    if (is.null(vec)) stop(paste("E_primeiro[[", s, "]] é NULL"))
    matrix(vec, nrow = n_signatures_rank, ncol = j_list[s], byrow = FALSE)
  })
  
  A_list <- lapply(resultado_par$A_primeiro, function(v) {
    v <- as.numeric(unlist(strsplit(as.character(v), "")))
    matrix(v, nrow = n_signatures_rank, ncol = n_signatures_rank, byrow = FALSE)
  })
  
  P <- resultado_par$P_primeiro
  
  # -----------------------------
  # 2c) Detectar assinaturas vazias e filtrar
  # -----------------------------
  assinatura_vazia <- rep(TRUE, n_signatures_rank)
  for (s in 1:S) {
    assinatura_vazia <- assinatura_vazia & (rowSums(A_list[[s]]) == 0)
  }
  keep <- which(!assinatura_vazia)
  
  # Filtrar listas e atualizar n_signatures
  A_list <- lapply(A_list, function(A) A[keep, keep, drop = FALSE])
  E_list <- lapply(E_list, function(E) E[keep, , drop = FALSE])
  P <- P[, keep, drop = FALSE]
  n_signatures <- length(keep)  # número real de assinaturas deste rank
  
  # -----------------------------
  # 2d) Inicialização de Ap, Bp, Ae_list, Be_list, Fi_list, Z_list
  # -----------------------------
  Ap <- matrix(rexp(i * n_signatures, rate = lp), i, n_signatures)
  Bp <- matrix(rgamma(i * n_signatures, shape = ap, rate = bp), i, n_signatures)
  
  Ae_list <- vector("list", S)
  Be_list <- vector("list", S)
  PE_list <- vector("list", S)
  Fi_list <- vector("list", S)
  Z_list <- vector("list", S)
  
  for (s in 1:S) {
    Ae_list[[s]] <- matrix(rexp(n_signatures * j_list[[s]], rate = le), n_signatures, j_list[[s]])
    Be_list[[s]] <- matrix(rgamma(n_signatures * j_list[[s]], shape = ae, rate = be), n_signatures, j_list[[s]])
    
    PE_list[[s]] <- P %*% E_list[[s]]
    Fi_list[[s]] <- array(0, dim = c(i, j_list[[s]], n_signatures))
    Z_list[[s]] <- array(0, dim = c(i, j_list[[s]], n_signatures))
    
    for (m in 1:n_signatures) {
      num <- (P[, m, drop = FALSE] %*% E_list[[s]][m, , drop = FALSE])
      Fi_list[[s]][, , m] <- num / PE_list[[s]]
    }
    
    for (k in 1:i) {
      Z_list[[s]][k, , ] <- t(sapply(1:j_list[[s]], function(g) {
        prob <- as.vector(Fi_list[[s]][k, g, ])
        prob[is.na(prob)] <- 0
        if(sum(prob) <= 0) return(rep(0, length(prob)))
        prob <- prob / sum(prob)
        rmultinom(1, size = M_list[[s]][k, g], prob = prob)
      }))
    }
  }
  
  fix_colsA <- as.integer(1:n_signatures)
  
  # -----------------------------
  # 2e) Rodar burn_in
  # -----------------------------
  resultado_2 <- burn_in_9(
    num_iterations, P, Ap, Bp, E_list, Ae_list, Be_list, A_list, Z_list, Fi_list,
    M_list, W_list, S, i, n_signatures, j_list,
    fix_colsP, fix_colsA, fix_colsE, le, 
    incremento_y = 1/34,
    decressimo_l_1 = 0.05,
    decressimo_l_2 = 0.0048,
    decressimo_l_3 = 0.00045
  )
  
  # -----------------------------
  # 2f) Calcular BIC
  # -----------------------------
  n_iter <- dim(resultado_2$P_store)[3]
  liks <- numeric(n_iter)
  
  for (iter in 1:n_iter) {
    P_iter <- resultado_2$P_store[,,iter]
    E_iter_list <- lapply(1:S, function(s) {
      matrix(resultado_2$E_store[[s]][[iter]], nrow = n_signatures, ncol = j_list[s], byrow = FALSE)
    })
    E_iter_AE_list <- lapply(1:S, function(s) A_list[[s]] %*% E_iter_list[[s]]) #<----- ok
    
    lik_total_iter <- 0
    for (s in 1:S) {
      P_A1 <- P_iter %*% A_list[[s]] * W_list[[s]][1,1]
      termo_log <- log(P_A1 %*% E_iter_AE_list[[s]])
      termo_log[is.infinite(termo_log) | is.na(termo_log)] <- 0
      lgamma_term <- lgamma(M_list[[s]] + 1)
      lik_total_iter <- lik_total_iter + sum((-P_A1 %*% E_iter_AE_list[[s]]) + M_list[[s]] * termo_log - lgamma_term)
    }
    liks[iter] <- lik_total_iter
  }
  
  max_lik <- max(liks)
  Vp1 <- nrow(P_iter) * ncol(P_iter) + sum(sapply(1:S, function(s) sum(diag(A_list[[s]])*j_list[[s]])))
  Vp2 <- log(nrow(P_iter)*sum(j_list))
  BIC <- Vp1*Vp2  - 2 * max_lik
  
  bic_values[r] <- BIC
  resultados_lista[[r]] <- list(
    resultado_2 = resultado_2,
    P_medio = Reduce("+", lapply(1:n_iter, function(iter) resultado_2$P_store[,,iter])) / n_iter,
    E_list_medio = lapply(1:S, function(s) {
      Reduce("+", lapply(1:n_iter, function(iter) {
        matrix(resultado_2$E_store[[s]][[iter]], nrow = n_signatures, ncol = j_list[s], byrow = FALSE)
      })) / n_iter
    }),
    A_list = A_list,
    BIC = BIC
  )
}

# -----------------------------
# 3) Selecionar resultado com menor BIC
# -----------------------------
best_index <- which.min(bic_values)
best_result <- resultados_lista[[best_index]]

cat("Melhor rank_choice:", rank_choices[best_index], "com BIC =", bic_values[best_index], "\n")   

# -----------------------------
# 4) Coletar resultados principais
# -----------------------------

dados <- list(
  P_medio      = best_result$P_medio,
  E_list_medio = lapply(seq_along(best_result$E_list_medio), function(s) {
                    best_result$A_list[[s]] %*% best_result$E_list_medio[[s]]
                  }),
  A_list       = best_result$A_list,
  BIC          = best_result$BIC,
  rank_choice  = rank_choices[best_index]
)


#==========================================================
#    Manipulação dos dados para graficar as Exposições 
#==========================================================

S <- length(dados$A_list)
E_list_medio <- list()

for (s in 1:S) {
    
  E_list_medio[[s]] <- dados$E_list_medio[[s]]

}  
                     