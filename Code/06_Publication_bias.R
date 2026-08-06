#06 Publication bias: egger's tests 
## Any pathogens
funnelpathogens <- rma.glmm(
  measure = "PLO",
  data = filtered_anypathogens,
  xi = Use,
  ni = Case
)
egger_pathogens <- regtest(
  x = funnelpathogens$yi,
  sei = sqrt(funnelpathogens$vi),
  model = "lm",
  predictor = "sei"
)
print(egger_pathogens)


## Any viruses
funnelVirus <- rma.glmm(
  measure = "PLO",
  data = filtered_anyviruses,
  xi = Use,
  ni = Case
)

egger_virus <- regtest(
  x = funnelVirus$yi,
  sei = sqrt(funnelVirus$vi),
  model = "lm",
  predictor = "sei"
)
print(egger_virus)


## Any bacteria
funnelBacteria <- rma.glmm(
  measure = "PLO",
  data = filtered_anybacteria,
  xi = Use,
  ni = Case,
  method = "FE"
)

egger_bacteria <- regtest(
  x = funnelBacteria$yi,
  sei = sqrt(funnelBacteria$vi),
  model = "lm",
  predictor = "sei"
)
print(egger_bacteria)


## Mono-infection
funnelmono <- rma.glmm(
  measure = "PLO",
  data = filtered_mono,
  xi = Use,
  ni = Case,
  method = "FE"
)

egger_mono <- regtest(
  x = funnelmono$yi,
  sei = sqrt(funnelmono$vi),
  model = "lm",
  predictor = "sei"
)
print(egger_mono)
