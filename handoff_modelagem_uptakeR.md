# Relatório de Handoff: Modelagem de Transporte Reativo e Calibração no R

Este documento consolida as metas, decisões metodológicas, problemas superados, implementações de código e referências para a calibração hidráulica (NaCl) e reativa ($NH_4-N$) em riachos, aplicando os conceitos de **Soetaert & Meysman (2012)** [1], **Runkel (2007)** [74], **Tromboni et al. (2017)** [258] e a abordagem de identificabilidade de **Bonanno et al. (2022)** [379].

---

## 1. Objetivos do Projeto
O objetivo principal é simular o transporte e o uptake biológico de nitrogênio amoniacal ($NH_4-N$) em um trecho de rio de **220 metros**, utilizando o R para prototipagem rápida através do acoplamento dos pacotes `ReacTran` [3] (equações de transporte advectivo-dispersivo com armazenamento transitório), `deSolve` [3] (resolução numérica de equações diferenciais rígidas via `ode.1D`) e `FME` [19] (calibração e estimativa de parâmetros).

Os objetivos específicos incluem:
1. **Calibração Hidráulica (Tratamento Conservativo):** Utilizar dados de uma injeção de sal (NaCl) para calibrar os parâmetros de dispersão ($D$), troca hidráulica ($\alpha$) e área da zona de armazenamento ($A_s$ ou `area2`) [95].
2. **Calibração Reativa (Tratamento Ativo):** Utilizar uma injeção instantânea (*slug*) de cloreto de amônio ($NH_4Cl$) para estimar os coeficientes de decaimento de primeira ordem no canal principal ($\lambda$) e na zona hiporreica ($\lambda_s$) [94].
3. **Partição Biológica de Processos:** Isolar a remoção do nutriente que ocorre na coluna d'água (via perifíton/produtores primários) daquela ocorrida na zona hiporreica (via comunidades microbianas nos sedimentos) [82, 259].
4. **Extração da Tríade Métrica de Uptake:** Calcular parâmetros independentes da hidrologia local: comprimento de uptake ($S_w$), velocidade de uptake ($v_f$) e taxa areal de uptake ($U$) [128, 177].

---

## 2. Decisões Metodológicas e Resultados

### A. Calibração Hidráulica (Curva de NaCl)
A calibração do modelo físico-hidráulico foi concluída com sucesso usando uma **estratégia Global-Local** [39, 399]:
* **Passo 1 (Busca Global):** Uso do método `"Pseudo"` (algoritmo de busca aleatória controlada de Price) para varrer o espaço de parâmetros e evitar mínimos locais ou falsas convergências [39, 384].
* **Passo 2 (Refinamento Local):** Aplicação do algoritmo Levenberg-Marquardt (`"Marq"`) a partir do ponto de ótimo global para estimar erros padrão e significância estatística exata [39, 423].

**Parâmetros Calibrados Obtidos:**
* Coeficiente de Dispersão ($D$): **$0,2136 \, m^2/s$** ($p < 2e-16$) [Histórico]
* Taxa de Troca ($\alpha$): **$0,001522 \, s^{-1}$** ($p < 2e-16$) [Histórico]
* Área da Zona de Armazenamento ($A_s$ / `area2`): **$0,5586 \, m^2$** ($p < 2e-16$) [Histórico]
* **Razão $A_s/A$:** **$33,25\%$** (calculada como $A_s / (vazao\_Q / veloc\_v)$) [Histórico]

*Interpretação Física:* O valor de $A_s/A = 33,25\%$ indica uma expressiva zona de armazenamento hidráulico no trecho estudado (como piscinas naturais, remansos ou trocas com o sedimento de fundo), o que difere de sistemas mais rápidos como o Barra Pequena ($A_s/A \approx 4\%$) [273], sugerindo alta capacidade de retenção de água e solutos [41].

### B. Segmentação Dinâmica da Curva (Rising, Peak, Tail)
Seguindo a abordagem de **Bonanno et al. (2022)** [390], a curva de breakthrough (BTC) do nutriente foi categorizada em seções para contornar a equifinalidade:
1. **Rising Limb & Peak (Aumento e Pico):** Domínio advectivo-dispersivo [398]. Controlado pela velocidade ($v$), área do canal ($A$), dispersão ($D$) e uptake no canal ($\lambda$) [415, 437].
2. **Tail (Cauda):** Domínio do armazenamento transitório [386]. É a cauda que contém quase $100\%$ do conteúdo de informação sobre a taxa de troca ($\alpha$), área de armazenamento ($A_s$) e decaimento hiporreico ($\lambda_s$) [416, 446].

---

## 3. Problemas Enfrentados e Soluções Técnicas

### Problema 1: Falta de Identificabilidade e Alta Correlação ($-0,99$)
Ao calibrar $\lambda$ e $\lambda_s$ com ajustes locais diretos, os parâmetros apresentaram correlação extrema ($-0,9989$) e p-values insignificantes ($p > 0,7$), ficando "presos" no chute inicial [Histórico].
* **Causa:** O otimizador local não conseguia distinguir se a remoção ocorria no canal principal ou na zona hiporreica, compensando um parâmetro com o outro [432].
* **Solução:** Implementação do método de **Bonanno et al. (2022)** [381]. Substituição do ajuste local direto pela **Amostragem de Hipercubo Latino (LHS)** com 10.000 simulações para mapear o espaço de parâmetros globalmente antes do refinamento [395, 448].

### Problema 2: Resultados `NA %` na Partição de Runkel (2007)
A partição de processos entre Canal e Zona Hiporreica retornava `NA %` para os percentuais de atribuição [Histórico].
* **Causa:** Como o trecho simulado era curto (220 m) e a velocidade baixa ($0,07 \, m/s$), a integral de massa do pulso reativo (`mass_B`) e conservador (`mass_A`) era truncada antes de toda a cauda do nutriente ter saído do trecho, resultando em diferença de massa zero (divisão por zero) [117, Histórico].
* **Solução:** Aumento do vetor de tempo de simulação (`times_sim`) para **7.500 segundos**, garantindo que a cauda retorne totalmente ao nível de background e o balanço de massa seja fechado sem erros numéricos [137, 446].

### Problema 3: Erro de Sintaxe `Latinhyper` do Pacote `FME`
O R retornava o erro `Error in 1:npar : argument of length 0` ao rodar a amostragem global [Histórico].
* **Causa:** A função `Latinhyper` exige receber os limites dos parâmetros na forma de uma **matriz nomeada** bidimensional, mas o script utilizava uma lista [403, Histórico].
* **Solução:** O intervalo de busca foi reconfigurado em uma matriz estruturada com row names correspondentes aos parâmetros [Histórico]:
  ```R
  param_range <- matrix(c(0, 0.005, 0, 0.01), nrow = 2, byrow = TRUE)
  rownames(param_range) <- c("lambda", "lambda_s")
  ```

---

## 4. Código R Consolidado (Versão Corrigida e Otimizada)

```R
################################################################################
# MODELAGEM DE TRANSPORTE REATIVO E PARTIÇÃO DE UPSTAKE (NH4-N)
# Baseado em Runkel (2007), Tromboni et al. (2017) e Bonanno et al. (2022)
################################################################################

library(ReacTran)
library(deSolve)
library(FME)

# --- 1. CONFIGURAÇÕES FÍSICAS E HIDRÁULICAS (CALIBRADAS VIA NaCl) ---
riverlen   <- 220       # m (Comprimento do trecho)
vazao_Q    <- 0.1431    # m3/s (Vazão)
veloc_v    <- 0.07248   # m/s (Velocidade)
largura_w  <- 5.42      # m (Largura molhada)
area       <- vazao_Q / veloc_v
prof_mc    <- area / largura_w

# Parâmetros físicos consolidados na calibração de NaCl
D_calibrado      <- 0.2136    # Coeficiente de dispersão (m2/s)
alpha_calibrado  <- 0.001522  # Coeficiente de troca (1/s)
area2_calibrado  <- 0.5586    # Área de armazenamento transitório (As, m2)

# Configuração da grade numérica (Finite Differences)
dx   <- 5
Nb   <- riverlen / dx
Vol  <- area * dx
Yini <- rep(0, 2 * Nb) # Concentrações iniciais (C e Cs) zeradas acima do background

# --- 2. DADOS EXPERIMENTAIS DO PULSO (BREAKTHROUGH CURVE) ---
# Dados categorizados em rising, peak e tail para DYNIA
slug_N <- data.frame(
  time = c(0, 0, 2281, 2321, 2350, 2376, 2418, 2452, 2487, 2520, 2543, 2575, 
           2610, 2647, 2712, 2856, 3102, 3288, 3558, 3710, 4338, 5027),
  NH4_N_ug_L = c(0.1657, 0.1406, 18.5401, 23.3402, 27.2334, 35.3367, 40.7488, 48.0632, 
                 53.5047, 58.0615, 64.9630, 71.3188, 78.1761, 82.7033, 100.7681, 
                 127.2238, 121.8855, 113.0374, 81.3614, 66.7473, 25.0287, 10.9677),
  limb = c(rep("bg", 2), rep("rising", 13), "peak", rep("tail", 6))
)

# Estequiometria do pulso: 116.9g de NH4Cl dissolvem em balde -> ~30,610 mg de N puro
massa_N_mg <- 116.9 * (14.01 / 53.49) * 1000  # ~30610 mg
t_pulse_N  <- 5                              # Duração da injeção em segundos
Camb       <- mean(slug_N$NH4_N_ug_L[1:2])    # Concentração de fundo (background)

# --- 3. FUNÇÃO DO MODELO MATEMÁTICO ---
Rivermodel_Nutrient <- function(time, state, pars) {
  lambda   <- pars["lambda"]   
  lambda_s <- pars["lambda_s"] 
  
  C  <- state[1:Nb]
  Cs <- state[(Nb+1):(2 * Nb)]
  
  # Concentração de entrada Upstream (Cup) em ug/L
  Cup <- if(time <= t_pulse_N) (massa_N_mg * 1000) / (vazao_Q * 1000 * t_pulse_N) else 0
  
  # Transporte de advecção-dispersão no canal principal
  tranC <- tran.volume.1D(C = C, C.up = Cup, flow = vazao_Q, Disp = D_calibrado, V = Vol)
  exchange <- alpha_calibrado * (Cs - C)
  
  # Equações diferenciais de balanço de massa (C no canal, Cs no sedimento)
  dC  <- tranC$dC + exchange - lambda * C
  dCs <- - exchange * (area / area2_calibrado) - lambda_s * Cs
  
  list(c(dC, dCs))
}

# --- 4. CUSTO DO MODELO ---
ModelCost_N <- function(p) {
  out <- ode.1D(y = Yini, func = Rivermodel_Nutrient, times = slug_N$time, parms = p, nspec = 2)
  sim_at_end <- data.frame(time = out[, 1], NH4_N_ug_L = out[, Nb + 1])
  return(modCost(model = sim_at_end, obs = slug_N[, 1:2]))
}

# --- 5. EXECUÇÃO DA AMOSTRAGEM GLOBAL (LHS) & CALIBRAÇÃO (BONANNO 2022) ---
# Criação correta da matriz de parâmetros para Latin Hypercube Sampling (LHS)
param_range <- matrix(c(0, 0.005,  # Intervalo de busca para lambda (canal)
                        0, 0.010), # Intervalo de busca para lambda_s (hiporreico)
                      nrow = 2, byrow = TRUE)
rownames(param_range) <- c("lambda", "lambda_s")
colnames(param_range) <- c("min", "max")

cat("Executando amostragem LHS de 10.000 simulações para identificabilidade...\n")
set.seed(123) # Reprodutibilidade
LHS_samples <- Latinhyper(param_range, n = 10000)
costs_LHS   <- apply(LHS_samples, 1, function(p) ModelCost_N(p)$model)
best_p_LHS  <- LHS_samples[which.min(costs_LHS), ]

# Refinamento Local Marquardt robusto pós-LHS
Fit_Final <- modFit(f = ModelCost_N, p = best_p_LHS, method = "Marq", lower = c(0, 0))
summary(Fit_Final)

# --- 6. PARTIÇÃO DE MASSA (RUNKEL 2007) ---
# vetor de tempo longo (7.500s) evita erros numéricos de NA % na cauda
times_sim <- seq(0, 7500, by = 10) 

run_and_mass <- function(p) {
  out <- ode.1D(y = Yini, func = Rivermodel_Nutrient, times = times_sim, parms = p, nspec = 2)
  conc_final <- out[, Nb + 1]
  sum(vazao_Q * conc_final * 10) # Integral sob a curva
}

pars_est <- Fit_Final$par
mass_A <- run_and_mass(c(lambda = 0, lambda_s = 0))        # Simulação A: Conservador
mass_B <- run_and_mass(pars_est)                          # Simulação B: Reativo Completo
mass_C <- run_and_mass(c(lambda = pars_est["lambda"], lambda_s = 0)) # Simulação C: Só Canal
mass_D <- run_and_mass(c(lambda = 0, lambda_s = pars_est["lambda_s"])) # Simulação D: Só SZ

total_uptake <- mass_A - mass_B
pct_mc_raw   <- (mass_A - mass_C) / total_uptake * 100
pct_sz_raw   <- (mass_A - mass_D) / total_uptake * 100

# Normalização proporcional para fechamento de balanço de massa (100%)
atribuicao_mc <- (pct_mc_raw / (pct_mc_raw + pct_sz_raw)) * 100
atribuicao_sz <- (pct_sz_raw / (pct_mc_raw + pct_sz_raw)) * 100

# --- 7. EXIBIÇÃO DA TRÍADE MÉTRICA DE UPTAKE ---
Sw_mc <- veloc_v / pars_est["lambda"]
vf_mc <- pars_est["lambda"] * prof_mc 
U_mc  <- (vazao_Q * Camb) / (Sw_mc * largura_w) 

# Outputs consolidados
cat("--- RESULTADOS FINAIS NH4-N ---\n")
cat("Atribuição Canal Principal (Perifíton):", round(atribuicao_mc, 1), "%\n")
cat("Atribuição Zona Hiporreica (Microbiano):", round(atribuicao_sz, 1), "%\n")
cat("Sw (Comprimento de Uptake):", round(Sw_mc, 2), "m\n")
cat("vf (Velocidade de Uptake):", round(vf_mc * 1000 * 60, 4), "mm/min\n")
cat("U (Taxa Areal):", round(U_mc * 3600, 4), "ug/m2/h\n")
```

---

## 5. Referências Bibliográficas

*   **[1, 12, 19] Soetaert, K. & Meysman, F. (2012).** *Reactive transport in aquatic ecosystems: Rapid model prototyping in the open source software R*. Environmental Modelling & Software, 32, 49-60.
*   **[74, 82, 137] Runkel, R. L. (2007).** *Toward a transport-based analysis of nutrient spiraling and uptake in streams*. Limnology and Oceanography: Methods, 5, 50-62.
*   **[172, 210] Covino, T. P., McGlynn, B. L. & McNamara, R. A. (2010).** *Tracer Additions for Spiraling Curve Characterization (TASCC): Quantifying stream nutrient uptake kinetics from ambient to saturation*. Limnology and Oceanography: Methods, 8, 484-498.
*   **[258, 259, 311] Tromboni, F., Dodds, W. K., Neres-Lima, V., Zandona, E. & Moulton, T. P. (2017).** *Heterogeneity and scaling of photosynthesis, respiration, and nitrogen uptake in three Atlantic Rainforest streams*. Ecosphere, 8(9), e01959.
*   **[379, 384, 400] Bonanno, E., Blöschl, G. & Klaus, J. (2022).** *Exploring tracer information in a small stream to improve parameter identifiability and enhance the process interpretation in transient storage models*. Hydrology and Earth System Sciences, 26, 6003–6028.
