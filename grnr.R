rm(list = ls())
if (!require("SixSigma")) install.packages("SixSigma")
library(SixSigma)
library(readxl)
library(ggplot2)
# function
varcomp_func <- function(a, b, n, appr, part, digits){
  varComp <- matrix(ncol = 2, nrow = 7)
  rownames(varComp) <- c("Total Gage R&R", "  Repeatability", 
                         "  Reproducibility", 
                         paste0("    ", appr),
                         paste0(part, ":", appr), 
                         "Part-To-Part", "Total Variation")
  varComp[2, 1] <- modelm[[1]][4, 3] # repeatability
  varComp[4, 1] <- max(c((modelm[[1]][2, 3] - modelm[[1]][3, 3])/(a * n), 0))  # Appraiser Reproducibility
  varComp[5, 1] <- max(c((modelm[[1]][3, 3] - modelm[[1]][4, 3])/n, 0))  # part:Appr Reproducibility
  varComp[3, 1] <- varComp[4, 1] + varComp[5, 1]  # total reproducibility
  varComp[6, 1] <- max(c((modelm[[1]][1, 3] - modelm[[1]][3, 3])/(b * n),0))  # part to part
  varComp[1, 1] <- varComp[2, 1] + varComp[3, 1]  # Total GRR
  varComp[7, 1] <- varComp[1, 1] + varComp[6, 1] 
  varComp[, 2] <- sqrt(varComp[, 1])
  colnames(varComp) <- c('VarComp','StdDev')
  varComp <- data.frame(round(varComp, digits = digits))
  
  return(varComp)
}
# -------------------

df <- read_xlsx('./gage.xlsx')
df$operator <- factor(df$operator)
df$subject <- factor(df$subject)
right <- df[df$side == 'r',]
left <- df[df$side == 'l',]

# fit anova with 2 factor interaction
fit <- aov(temp ~ subject*operator, data = df)
(gar <- round(summary(fit)[[1]], 3))
# gar <- data.frame(cbind(Source = row.names(gar), gar))
# writexl::write_xlsx(gar, '../original_anova.xlsx')

# treat interaction term as residuals
# algorithm for F-statistic won't be the same
# Original: "MS_residual" as denominator
# Revise: "MS_interaction" as denominator

# residual diagnostics
par(mfrow = c(2, 2))
plot(fit)
shapiro.test(df$temp)

# GR & R analysis
result <- ss.rr(temp, subject, operator, 
                35, 40, 5.15, data = df,
                signifstars = T)  # default: errorTerm = 'interaction'
(gar1 <- round(result$anovaTable[[1]], 3))
# gar1 <- data.frame(cbind(Source = row.names(gar1), gar1))
# writexl::write_xlsx(gar1, '../anova_correct.xlsx')
(varcomp <- round(cbind(result$varComp, StdDev = result$studyVar[,1]), 3))
(varcomp <- varcomp[, c(1, 3, 2)])
# (std_original <- data.frame(cbind(Source = rownames(result$varComp), varcomp)))
# writexl::write_xlsx(std_original, '../varcomp.xlsx')


#-------------------------
# right
fit_r <- aov(temp ~ subject*operator, data = right)
modelm <- summary(fit_r)
appr <- 'operator'
part <- 'subject'
a <- 10; b <- 3; n <- 2
std_right <- varcomp_func(a, b, n, appr, part, 3)

# left
fit_l <- aov(temp ~ subject*operator, data = left)
modelm <- summary(fit_l)
std_left <- varcomp_func(a, b, n, appr, part, 3)

# calculate mean & max
seperate <- cbind(right, left$temp)
seperate <- seperate[,-3]
colnames(seperate)[4:5] <- c('rt', 'lt')
seperate$means <- rowMeans(seperate[,c('rt','lt')])
seperate$max_ <- apply(seperate[,c('rt','lt')], 1, max)
head(seperate)

# mean
fit_mean <- aov(means ~ subject*operator, data = seperate)
modelm <- summary(fit_mean)
a <- 10; b <- 3; n <- 2
std_mean <- varcomp_func(a, b, n, appr, part, 3)

# max
fit_max <- aov(max_ ~ subject*operator, data = seperate)
modelm <- summary(fit_max)
std_max <- varcomp_func(a, b, n, appr, part, 3)

# combined stddev components (all, right, left, mean, max) 
combined_std <- cbind(std_original[, 3], 
                      Right = std_right[, 2],
                      Left = std_left[, 2], 
                      Mean = std_mean[, 2],
                      Max = std_max[, 2])
colnames(combined_std)[1] <- 'All'
pt_ratio <- round(((5.15*as.numeric(combined_std[1, ])) / 5) * 100, 2)
combined_std <- data.frame(rbind(combined_std, pt_ratio))
combined_std <- cbind(Source = c(std_original[, 1], 'P/T Ratio'), combined_std)
combined_std
# writexl::write_xlsx(combined_std, '../all_std.xlsx')


# calculate contribution of each variance components & pt-ratio
part2part <- round(as.numeric(combined_std[6,-1])**2/as.numeric(combined_std[7,-1])**2 * 100, 2)
gage <- 100 - part2part
ratios <- rbind(part2part, gage)
Source <- c('Part-to-Part','Gage')
ratios <- data.frame(cbind(Source, ratios))
colnames(ratios)[-1] <- colnames(combined_std)[-1]
ratios
# writexl::write_xlsx(ratios, '../ratios.xlsx')

# scatterplot of measurements
plot(left$temp, right$temp, xlab = 'Left', ylab = 'Right',
     main = 'Scatterplot of measurements right vs. left',
     pch = 19)

# ------------------------------------
# range UCL/CL/LCL
r <- aggregate(x = df$temp,
          # Specify group indicator
          by = list(df$subject, df$operator),      
          # Specify function (i.e. mean)
          FUN = range)
r$x[,1] <- r$x[,1] * (-1)
ar <- mean(rowSums(r$x))
d4 <- 2.28  # by table
d3 <- 0   # by table
uclr <- ar*d4; lclr <- ar*d3
# data.xrange <- aggregate(as.formula(paste('temp', "~", 'operator', "+", 'subject')),
#                          data = df,
#                          function(x) {
#                            max(x) - min(x)
#                          })
# ar <- mean(data.xrange[['temp']])

# xbar UCL/CL/LCL
m <- aggregate(x = df$temp,
               by = list(df$subject, df$operator),      
               FUN = mean)
xbar <- mean(df$temp)
a2 <- 0.73  # by table
uclm <- xbar + a2 * ar; lclm <- xbar - a2 * ar
