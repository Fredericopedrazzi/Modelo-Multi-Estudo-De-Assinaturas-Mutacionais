library(labeling)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(proxy)
library(IRdisplay)  

# ---------------------------
# Parâmetros
# ---------------------------
i <- 96
n_signatures <- 50
S <- 2
num_iterations <- 13000

#--------------------
# Hiper Multi-Estudo
#--------------------
ap = 5 
bp = 0.05
ae = 5
be = 0.01 
lp = 0.5
le = 0.1  

# Fix columns como vetores inteiros (vazios se não houver colunas fixas)
fix_colsP <- integer(0)   
fix_colsA <- integer(0)    
fix_colsE <- integer(0)


eps <- .Machine$double.xmin  

df1 <- read.table("countsMATH_Homogeneos_COAD.txt", header = FALSE, sep = "\t")

df2 <- read.table("countsMATH_Heterogeneos_COAD.txt", header = FALSE, sep = "\t")

data <- list(df1, df2)

source("De_novo.R")  


#============================
# Acesso a resultados
# ===========================

#dados$P_medio     
#dados$E_list_medio                   
#dados$A_list       
#dados$BIC          
#dados$rank_choice  
