library(dplyr)
library(ggplot2)

expression <- read.table("proteins.tsv", header = TRUE, stringsAsFactors = FALSE)

ldl_predictions <- readRDS("TRIM25-C7orf50_predictions.rds")
ad_predictions <- readRDS("APOE-LDLR_predictions.rds")
pd_predictions <- readRDS("DARS1_SNCA_predictions.rds")

ldl_predictions <- as.data.frame(ldl_predictions)
ad_predictions <- as.data.frame(ad_predictions)
pd_predictions <- as.data.frame(pd_predictions)

names(ldl_predictions) <- c("a", "b", "c")
names(ad_predictions) <- c("a", "b", "c")
names(pd_predictions) <- c("a", "b", "c")

ldl_predictions$IID <- expression$IID
ad_predictions$IID <- expression$IID
pd_predictions$IID <- expression$IID

ldl_predictions <- ldl_predictions %>% mutate(quantile = ntile(c, 10))
ad_predictions <- ad_predictions %>% mutate(quantile = ntile(c, 10))
pd_predictions <- pd_predictions %>% mutate(quantile = ntile(c, 10))

ldl_predictions <- merge(ldl_predictions, subset(expression, select = c("IID", "TRIM25", "C7orf50")), by = "IID")
ad_predictions <- merge(ad_predictions, subset(expression, select = c("IID", "APOE", "LDLR")), by = "IID")
pd_predictions <- merge(pd_predictions, subset(expression, select = c("IID", "DARS1", "SNCA")), by = "IID")

ldl_correlations <- ldl_predictions %>% group_by(quantile) %>% 
  summarize(correlation = cor(TRIM25, C7orf50, use = "complete.obs", method = "pearson"))

ad_correlations <- ad_predictions %>% group_by(quantile) %>% 
  summarize(correlation = cor(APOE, LDLR, use = "complete.obs", method = "pearson"))

pd_correlations <- pd_predictions %>% group_by(quantile) %>% 
  summarize(correlation = cor(DARS1, SNCA, use = "complete.obs", method = "pearson"))

ldl_correlations %>% ggplot(aes(x = quantile, y = correlation)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Quantile of imputed co-expression", y = "Observed correlation", title = "TRIM25 and C7orf50") +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
  theme_light() +
  theme(text = element_text(family = "sans"))
ggsave(filename = "ldl_corplot.png", device = "png", width = 4, height = 3, dpi = 600)
  
ad_correlations %>% ggplot(aes(x = quantile, y = correlation)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Quantile of imputed co-expression", y = "Observed correlation", title = "APOE and LDLR") +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
  theme_light() +
  theme(text = element_text(family = "sans"))
ggsave(filename = "ad_corplot.png", device = "png", width = 4, height = 3, dpi = 600)

pd_correlations %>% ggplot(aes(x = quantile, y = correlation)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(x = "Quantile of imputed co-expression", y = "Observed correlation", title = "DARS1 and SNCA") +
  scale_x_continuous(breaks = 1:10, minor_breaks = NULL) +
  theme_light() +
  theme(text = element_text(family = "sans"))
ggsave(filename = "pd_corplot.png", device = "png", width = 4, height = 3, dpi = 600)

