relaimpo_shap <- function(x, y, weights, intercept) {
  
  data <- as.data.frame(cbind(y, x))
  
  if (intercept) {
    fit <- lm(reformulate(termlabels = ".", response = as.name(names(data)[1])),
              weights = weights, data = data)
  } else {
    fit <- lm(reformulate(termlabels = ".-1", response = as.name(names(data)[1])),
              weights = weights, data = data)
  }
  
  # relaimpo: relative importance equals r-squared for univariate case
  if (ncol(x) == 1) {
    result <- summary(fit)$r.squared
  } else {
    
    # relaimpo: "model must contain intercept"
    result <- relaimpo::calc.relimp(fit, type = "lmg", rela = FALSE)@lmg
    
  }
  
  return(result)
  
}
