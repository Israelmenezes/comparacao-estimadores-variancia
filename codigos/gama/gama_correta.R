#Simulação da distribuição gama otimizada:

#EMV:
calculo_var <- function(par){
  par[1]/(par[2])^2 
}

#log função de verossimilhança: 
logLik_gama <-  function(p ,n ,sum_amostra ,sum_logAmostra ){
  alfa <- p[1]
  beta <- p[2]
  
  - ((n*alfa*log(beta))-(n*lgamma(alfa))+((alfa-1)*sum_logAmostra)-(beta*sum_amostra))
}

set.seed(484)
tamanhos_amostrais <- c(10,20,40,60,100)
#install.packages("openxlsx")
library(openxlsx)
resultados_tabela <- list() 

for (n in tamanhos_amostrais) {
#var_amostral:
amostra_varAmostral_n <- numeric(5000)
#var_Amostral_n_menos1:
amostra_varAmostral_n_menos1 <- numeric(5000)
#var_EMV:
amostra_varEMV_n <- numeric(5000)
#var_EMV_Bootstrap
amostra_varEMV_Bootstrap_n <- numeric(5000)

#EQM_diferenC'a_quadratica = dq:
amostra_dq_varEMV_n <- numeric(5000)
amostra_dq_varAmostral_n <- numeric(5000)
amostra_dq_varAmostral_n_menos1 <- numeric(5000)
amostra_dq_varEMV_Bootstrap_n <- numeric(5000)

for ( i in  1:5000) {
  #cv = 0,3
  y <- rgamma(n, shape = 11.111, rate = 1.111)
  #Variance_n-1:
  amostra_varAmostral_n_menos1[i] <- var(y)
  #var_n:
  amostra_varAmostral_n[i] <- var(y) * (n-1)/n
  
  #chute inicial com base em cada amostra:
  
  alfa_0 <-  (mean(y))^2/var(y)
  beta_0  <- mean(y)/var(y)
  p_chapeu <- c(alfa_0,beta_0)
  
  
  #somatorio do log de y p/ optim:
  sly <- sum(log(y))
  
  #somatorio de y p/ optim:
  sy <- sum(y)
  
  xx <- optim(p_chapeu,fn = logLik_gama, n = n, sum_amostra = sy, sum_logAmostra = sly,
  gr = NULL, 
  method = "L-BFGS-B",
  lower = 0.00001 ,
  hessian = FALSE
  )
  
  #EMVr
  amostra_varEMV_n[i] <- calculo_var(xx$par)
  
  #EMV_Bootstrap:
  amostra_varb_n <- numeric(1000)
  

  for (b in 1:1000){
    w <- rgamma(n,xx$par[1],xx$par[2])
    
    #chute inicial com base em cada amostra:
    
    alfa <-  (mean(w))^2/var(w)
    beta <- mean(w)/var(w)
    p_chapeu <- c(alfa,beta)
    
    #somatorio de w p/ o segundo optim:
    sw <- sum(w)
    
    #somatorio do log de w p/ o segundo optim:
    slw <- sum(log(w))
    
    
    yy <- optim(p_chapeu,fn = logLik_gama, n = n,  sum_amostra = sw, sum_logAmostra = slw,
    gr = NULL, 
    method = "L-BFGS-B",
    lower = 0.00001 ,
    hessian = FALSE
    )
    
    #var_b:
    amostra_varb_n[b] <- calculo_var(yy$par)
    
  }
  
  #var_EMV_Bootstrap:
  var_EMV_Bootstrap <- 2*calculo_var(xx$par)-mean(amostra_varb_n)
  
  amostra_varEMV_Bootstrap_n[i] <- var_EMV_Bootstrap
  
  
  #EQM_DiferenC'as^2: 
  dq_varEMV_n <-(calculo_var(xx$par) - 9)^2
  amostra_dq_varEMV_n[i] <- dq_varEMV_n

  dq_varAmostral_n <- ((var(y) * (n-1)/n) - 9)^2 
  amostra_dq_varAmostral_n[i] <- dq_varAmostral_n

  dq_varAmostral_n_menos1 <- (var(y) - 9)^2 
  amostra_dq_varAmostral_n_menos1[i] <- dq_varAmostral_n_menos1
    
  dq_varEMV_Bootstrap_n <- (var_EMV_Bootstrap - 9)^2
  amostra_dq_varEMV_Bootstrap_n[i] <- dq_varEMV_Bootstrap_n
  
  
  if (i%%250 == 0)
  {print(i)
    print(format(Sys.time(), "%H:%M:%S"))
  }
}



#Vies das variancias:
#####
#amostral:
vies_varAmostral_n <- mean(amostra_varAmostral_n) - 9
print(vies_varAmostral_n)
#Amostral_n_menos1:
vies_varAmostral_n_menos1 <- mean(amostra_varAmostral_n_menos1) - 9
print(vies_varAmostral_n_menos1)
#EMV
vies_varEMV_n <- mean(amostra_varEMV_n) - 9
print(vies_varEMV_n)
#EMV_Bootstrap:
vies_varEMV_Bootstrap_n <- mean(amostra_varEMV_Bootstrap_n) - 9
print(vies_varEMV_Bootstrap_n)
#####

#EQMs das variancias:
#####
#amostral_n
EQM_var_amostral <- mean(amostra_dq_varAmostral_n)
print(EQM_var_amostral)
#Amostral_n_menos1:
EQM_var_Amostral_n_menos1 <- mean(amostra_dq_varAmostral_n_menos1)
print(EQM_var_Amostral_n_menos1)
#EMV
EQM_var_EMV <- mean(amostra_dq_varEMV_n)
print(EQM_var_EMV)
#EMV_Bootstrap:
EQM_var_EMV_Bootstrap <- mean(amostra_dq_varEMV_Bootstrap_n)
print(EQM_var_EMV_Bootstrap)

resultados_tabela[[as.character(n)]] <- list(
  vies_varAmostral_n        = vies_varAmostral_n,
  vies_varAmostral_n_menos1 = vies_varAmostral_n_menos1,
  vies_varEMV_n             = vies_varEMV_n,
  vies_varEMV_Bootstrap_n   = vies_varEMV_Bootstrap_n,
  EQM_var_amostral          = EQM_var_amostral,
  EQM_var_Amostral_n_menos1 = EQM_var_Amostral_n_menos1,
  EQM_var_EMV               = EQM_var_EMV,
  EQM_var_EMV_Bootstrap     = EQM_var_EMV_Bootstrap
)

}
#####

#Criação da tabela:
#####
wb <- createWorkbook()
addWorksheet(wb, "Resultados")
ws <- "Resultados"
 
#Estilos
st_titulo <- createStyle(
  fontSize = 12, fontName = "Arial", fontColour = "white",
  fgFill = "#2E4057", textDecoration = "bold",
  halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#FFFFFF"
)
st_header <- createStyle(
  fontSize = 10, fontName = "Arial", fontColour = "white",
  fgFill = "#048A81", textDecoration = "bold",
  halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#FFFFFF"
)
st_n_label <- createStyle(
  fontSize = 10, fontName = "Arial", fontColour = "white",
  fgFill = "#54798C", textDecoration = "bold",
  halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#FFFFFF"
)
st_metrica <- createStyle(
  fontSize = 10, fontName = "Arial", fontColour = "#2E4057",
  textDecoration = "bold", halign = "center", valign = "center",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC"
)
st_vies <- createStyle(
  fontSize = 10, fontName = "Arial", fgFill = "#EAF4FB",
  halign = "center", valign = "center", numFmt = "0.0000",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC"
)
st_eqm <- createStyle(
  fontSize = 10, fontName = "Arial", fgFill = "#FEF9EF",
  halign = "center", valign = "center", numFmt = "0.0000",
  border = "TopBottomLeftRight", borderColour = "#CCCCCC"
)
 
#Título e cabeçalhos:
writeData(wb, ws, "Estimadores da Variância — Viés e EQM", startRow = 1, startCol = 1)
mergeCells(wb, ws, rows = 1, cols = 1:6)
addStyle(wb, ws, st_titulo, rows = 1, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, ws, rows = 1, heights = 28)
 
headers <- c("n", "Métrica", "Var Amostral (n)", "Var Amostral (n-1)", "EMV", "EMV Bootstrap")
for (k in seq_along(headers)) writeData(wb, ws, headers[k], startRow = 2, startCol = k)
addStyle(wb, ws, st_header, rows = 2, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb, ws, rows = 2, heights = 26)
 
#Preencher linhas:
linha <- 3
for (nm in names(resultados_tabela)) {
  r <- resultados_tabela[[nm]]
 
  # Coluna n mesclada nas 2 linhas (Viés + EQM)
  writeData(wb, ws, paste0("n = ", nm), startRow = linha, startCol = 1)
  mergeCells(wb, ws, rows = linha:(linha + 1), cols = 1)
  addStyle(wb, ws, st_n_label, rows = linha:(linha + 1), cols = 1, gridExpand = TRUE)
 
  # Linha Viés
  writeData(wb, ws, "Viés", startRow = linha, startCol = 2)
  addStyle(wb, ws, st_metrica, rows = linha, cols = 2)
  vals_vies <- c(r$vies_varAmostral_n, r$vies_varAmostral_n_menos1,
                 r$vies_varEMV_n, r$vies_varEMV_Bootstrap_n)
  for (k in seq_along(vals_vies)) {
    writeData(wb, ws, vals_vies[k], startRow = linha, startCol = 2 + k)
    addStyle(wb, ws, st_vies, rows = linha, cols = 2 + k)
  }
  setRowHeights(wb, ws, rows = linha, heights = 20)
 
  # Linha EQM
  writeData(wb, ws, "EQM", startRow = linha + 1, startCol = 2)
  addStyle(wb, ws, st_metrica, rows = linha + 1, cols = 2)
  vals_eqm <- c(r$EQM_var_amostral, r$EQM_var_Amostral_n_menos1,
                r$EQM_var_EMV, r$EQM_var_EMV_Bootstrap)
  for (k in seq_along(vals_eqm)) {
    writeData(wb, ws, vals_eqm[k], startRow = linha + 1, startCol = 2 + k)
    addStyle(wb, ws, st_eqm, rows = linha + 1, cols = 2 + k)
  }
  setRowHeights(wb, ws, rows = linha + 1, heights = 20)
 
  linha <- linha + 2
}
 
#Largura das colunas e salvar:
setColWidths(wb, ws, cols = 1,   widths = 10)
setColWidths(wb, ws, cols = 2,   widths = 10)
setColWidths(wb, ws, cols = 3:6, widths = 18)
 
saveWorkbook(wb, "simulação_gama_CV_0,3.xlsx", overwrite = TRUE)
cat("Tabela salva em: simulação_gama_CV_0,3\n")

getwd()
#####