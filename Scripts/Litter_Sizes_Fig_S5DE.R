################################
##### Litter size analysis #####
################################

# Last updated: 13/10/22 by LVB

# Description: Compare litter sizes and pre-weaning deaths across mice carrying
# the m.5024C>T and m.5019A>G mutations.

#----- Packages -----#
require("tidyverse")
require("trending")

#----- Data -----#
df <- read_csv("litter_sizes/Litter_Sizes.csv")

#----- Plot setup -----#
mut_cols <- c("m.5019" = "dodgerblue4", "m.5024" = "firebrick")

plot_save <- function(p, filename, size = 1, ar = 1, dev = "jpeg"){
  allowed_devs <- c("eps", "ps", "tex", "pdf", "jpeg", 
                    "tiff", "png", "bmp", "svg")
  if (!(dev %in% allowed_devs))
    stop("Invalid device.")
  if (dev != "jpeg" & !str_detect(filename, paste0("\\.", dev)))
    filename <- paste0(filename, ".", dev)
  if (dev == "jpeg" & !str_detect(filename, "\\.jpg|\\.jpeg"))
    filename <- paste0(filename, ".jpg")
  w <- round(180 * size)
  h <- w/ar
  ggsave(filename = filename,
         plot = p,
         width = w,
         height = h,
         units = "mm",
         device = dev)
}

#------ Figure S5D -----#
# Q: Does heteroplasmy have an effect on weaning litter size?
df$lh <- qlogis (.01 * df$mother_h) # heteroplasmy log-odds

# Poisson model
glm_pois <- glm_model(wean_size ~ lh : mutation,
                      family = "quasipoisson")
glm_pois_fitted <- fit(glm_pois, filter(df, wean_size>0))
df_fits <- predict(glm_pois_fitted, simulate_pi = F, uncertain = T, alpha = 0.1)
df_fits <- cbind(filter(df, wean_size>0), df_fits[, 4:6])

summary(glm_pois_fitted[[1]])
# Call:
#   glm(formula = wean_size ~ lh:mutation, family = "quasipoisson", 
#       data = data)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.8818  -0.8315   0.0665   0.8737   2.6034  
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.70295    0.04042  42.135  < 2e-16 ***
#   lh:mutationm.5019 -0.20871    0.04744  -4.400 1.45e-05 ***
#   lh:mutationm.5024  0.10694    0.06200   1.725   0.0854 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasipoisson family taken to be 1.67127)
# 
# Null deviance: 828.40  on 349  degrees of freedom
# Residual deviance: 779.44  on 347  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 5

# Figure
p <- ggplot(df_fits, aes(x = mother_h)) +
  theme_classic() + theme(legend.position = "top") +
  geom_point(aes(y = wean_size, colour = mutation)) +
  geom_line(aes(y = estimate, colour = mutation)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = mutation), alpha = .2) +
  scale_x_continuous("Mother heteroplasmy") +
  scale_y_continuous("Litter size at weaning", breaks = seq(0, 12, 2)) +
  scale_fill_manual("", values = mut_cols, labels = c("m.5019A>G", "m.5024C>T")) +
  scale_colour_manual("", values = mut_cols, labels = c("m.5019A>G", "m.5024C>T"))
p
plot_save(p, "litter_sizes/Figure_S5D.jpg", ar = 1.5, size = .5)

#----- Figure S5E -----#
# Q: Are m.5019 and m.5024 pre-weaning losses different?
losses <- df$litter_size - df$wean_size
df$loss_group <- case_when(losses == 0 ~ "0",
                            losses < 4 ~ "1-3",
                            TRUE ~ "4+")

# Chi-squared test
testthis <- split(df$loss_group, df$mutation)
testthis <- lapply(testthis, function(x) (table(x) %>% as.data.frame))
testthis <- lapply(1:2, function(i) cbind(testthis[[i]], group = names(testthis)[i])) %>% 
  bind_rows
colnames(testthis)[1] <- "loss"
testthis <- pivot_wider(testthis, id_cols = "loss", names_from = "group", values_from = "Freq")
testthis
# # A tibble: 3 × 3
# loss  m.5019 m.5024
# <fct>  <int>  <int>
#   1 0         95    149
# 2 1-3       52     33
# 3 4+        12      9

chisqtest <- chisq.test(testthis[, -1], simulate.p.value = T, B = 1e5) 
chisqtest
# Pearson's Chi-squared test with simulated p-value (based on 1e+05 replicates)
# 
# data:  testthis[, -1]
# X-squared = 13.816, df = NA, p-value = 0.00077
chisqtest$p.value %>% round(3) # 0.001

# Figure
p <- ggplot(df, aes(x = loss_group, fill = mutation)) +
  theme_classic() + theme(legend.position = "top") +
  geom_bar(aes(y =  after_stat(100 * count / sum(count))),
           position = position_dodge2(width = 0.75, preserve = "single")) +
  scale_x_discrete("Deaths per litter") +
  scale_y_continuous("Fraction of litters (%)", limits = c(0, 50)) +
  scale_fill_manual("", values = mut_cols, labels = c("m.5019A>G", "m.5024C>T")) +
  annotate("text", x = 2.7, y = 42, label = "Chi-squared test\np-value = 0.001")
p
plot_save(p, "litter_sizes/Figure_S5E.jpg", size = 0.5, ar = 3/2)
