# Constants
P_VALUE_THRESHOLD  = 0.05
BIAS_THRESHOLD     = 0.1
SHAPIRO_MIN_SAMPLE = 3
SHAPIRO_MAX_SAMPLE = 5000
LJUNG_BOX_LAG      = 5

# Helper function to get p-values from matrix or vector
getPValues<-function(pVals) {
  if (is.matrix(pVals)) {
    pVals[, 1]
  } else {
    pVals}}

# Helper function to check convergence test
checkConvergenceTest<-function(pVals, threshold = P_VALUE_THRESHOLD) {
  if (is.null(pVals)) {
    return(list(passed = NA, mean = NA))
  }
  vals = getPValues(pVals)
  list(
    passed = all(vals > threshold, na.rm = TRUE),
    mean = mean(vals, na.rm = TRUE))}

# get convergence diagnostics
getConvergenceDiagnostics<-function(jb) {
  result = list(
    Geweke_passed = NA,
    Geweke_mean   = NA,
    Heidel_passed = NA,
    Heidel_mean   = NA)
  
  if (!is.null(jb$pars)) {
    geweke = checkConvergenceTest(jb$pars$Geweke.p)
    result$Geweke_passed = geweke$passed
    result$Geweke_mean   = geweke$mean
    
    heidel = checkConvergenceTest(jb$pars$Heidel.p)
    result$Heidel_passed = heidel$passed
    result$Heidel_mean   = heidel$mean
  }
  
  return(result)}

# get residual diagnostics
getResidualDiagnostics<-function(jb) {
  result = list(
    n_residuals      = NA,
    mean_residual    = NA,
    sd_residual      = NA,
    Normality_passed = NA,
    Normality_p      = NA,
    Autocorr_passed  = NA,
    Autocorr_p       = NA,
    Runs_passed      = NA,
    Runs_p           = NA,
    Bias_passed      = NA)
  
  if (is.null(jb$diag) || !"residual" %in% names(jb$diag)) {
    return(result)}
  
  residuals  = jb$diag$residual
  nResiduals = length(residuals)
  
  result$n_residuals   = nResiduals
  result$mean_residual = mean(residuals, na.rm = TRUE)
  result$sd_residual   = sd(residuals,   na.rm = TRUE)
  
  # Normality test (Shapiro-Wilk)
  if (nResiduals >= SHAPIRO_MIN_SAMPLE && nResiduals <= SHAPIRO_MAX_SAMPLE) {
    shapiroTest = shapiro.test(residuals)
    result$Normality_passed = shapiroTest$p.value > P_VALUE_THRESHOLD
    result$Normality_p = shapiroTest$p.value}
  
  # Autocorrelation test (Ljung-Box)
  ljungBox = Box.test(residuals, lag = LJUNG_BOX_LAG, type = "Ljung-Box")
  result$Autocorr_passed = ljungBox$p.value > P_VALUE_THRESHOLD
  result$Autocorr_p = ljungBox$p.value
  
  # Runs test
  runsTest = jbrunstest(jb)
  result$Runs_passed = runsTest$runs.p > P_VALUE_THRESHOLD
  result$Runs_p = runsTest$runs.p
  
  # Bias test
  result$Bias_passed = abs(result$mean_residual) < BIAS_THRESHOLD
  
  return(result)}

# get terminal parameter estimates
getTerminalParameters<-function(jb) {
  result = list(
    BBmsy_terminal = NA,
    FFmsy_terminal = NA)
  
  if (is.null(jb$kobe) || !is.data.frame(jb$kobe) || nrow(jb$kobe) == 0) {
    return(result)}
  
  terminalKobe = jb$kobe[nrow(jb$kobe), ]
  
  if ("stock" %in% names(terminalKobe)) {
    result$BBmsy_terminal = terminalKobe$stock[1]}
  if ("harvest" %in% names(terminalKobe)) {
    result$FFmsy_terminal = terminalKobe$harvest[1]}
  
  return(result)}

# Create empty result data frame
createEmptyResult<-function(id, scenario) {
  data.frame(
    Stock = id,
    Scenario = scenario,
    Geweke_passed    = NA,
    Geweke_mean      = NA,
    Heidel_passed    = NA,
    Heidel_mean      = NA,
    n_residuals      = NA,
    mean_residual    = NA,
    sd_residual      = NA,
    Normality_passed = NA,
    Normality_p      = NA,
    Autocorr_passed  = NA,
    Autocorr_p       = NA,
    Runs_passed      = NA,
    Runs_p           = NA,
    Bias_passed      = NA,
    BBmsy_terminal   = NA,
    FFmsy_terminal   = NA,
    stringsAsFactors = FALSE)}

# Function to run diagnostics for a single stock/scenario
runDiags<-function(jb) {
    
  # get all diagnostics
  convergence = getConvergenceDiagnostics(jb)
  residuals   = getResidualDiagnostics(jb)
  terminal    = getTerminalParameters(jb)
    
  # Combine all results
  result = cbind(as.data.frame(convergence), 
                 as.data.frame(residuals), 
                 as.data.frame(terminal))
    
  return(result)}


# Ensure all data frames have the same columns
chkCols<-function(diagsList) {
  allCols = unique(unlist(lapply(diagsList, names)))
  
  lapply(diagsList, function(x) {
    missingCols = setdiff(allCols, names(x))
    if (length(missingCols) > 0) {
      for (col in missingCols) {
        x[[col]] = NA}}
    x[, allCols, drop = FALSE]})}

# Calculate summary statistics for a scenario
calcSmry<-function(data) {
  data.frame(
    nStocks = nrow(data),
    nGewekePassed   = sum(data$Geweke_passed, na.rm = TRUE),
    nHeidelPassed   = sum(data$Heidel_passed, na.rm = TRUE),
    nNormalPassed   = sum(data$Normality_passed, na.rm = TRUE),
    nNormalTotal    = sum(!is.na(data$Normality_passed)),
    nAutocorrPassed = sum(data$Autocorr_passed, na.rm = TRUE),
    nAutocorrTotal  = sum(!is.na(data$Autocorr_passed)),
    nRunsPassed     = sum(data$Runs_passed, na.rm = TRUE),
    nRunsTotal      = sum(!is.na(data$Runs_passed)),
    nBiasPassed     = sum(data$Bias_passed, na.rm = TRUE),
    nBiasTotal      = sum(!is.na(data$Bias_passed)))}

runAll<-function(jb){
  rtn = list()
  for (id in names(jb)) {
    rslt = runDiags(jb[[id]])
    if (!is.null(rslt))
      rtn[[id]] = rslt}
  
  # Standardize columns after each stock
  do.call(rbind, chkCols(rtn))}

jbICES=llply(histICES, function(x) x[["ICES"]]$fit)
jb1903=llply(histICES, function(x) x[["1903"]]$fit)

dgICES=runAll(jbICES)
dg1903=runAll(jb1903)
