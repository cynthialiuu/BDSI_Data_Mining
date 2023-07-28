# 1.1.1 Survival
library(tidyverse)

survivor = as.data.frame(Y)

survivor %>% 
  filter(Y > 60)
# 8 patients had survival times greater than 60

#1.1.1.2
survivor %>% 
  summarise(med_survival = median(Y))
# The median was 25.89

#1.1.1.3
survivor %>% 
  ggplot() +
  aes(y = Y) + geom_boxplot()
# 5 outliers

#1.1.1.4
is.na(survivor$Y)
# no NA values


#1.1.2 Pathways
#1.1.2.1
ncol(pathway.scores)

#1.1.2.2
hi = apply(pathway.scores, 2, var)
summary(hi)
#max = GNF2_MKI67 0.3072158 min = 0.0024670380.002467038

#1.1.2.3
paths = as.data.frame(pathway.scores)
paths %>% filter()

data_cor <- cor(pathway.scores[ , colnames(pathway.scores) != "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"], pathway.scores[ , "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"])
sort(data_cor, decreasing = T)[1] 

kegg_drug_cors <- apply(
  pathway.scores, 2, function(x) {
    cor(x, pathway.scores[, "KEGG_DRUG_METABOLISM_CYTOCHROME_P450"])
  })
sort(kegg_drug_cors, decreasing = TRUE)[2]

idx <- Y > 10
y_cors <- apply(
  pathway.scores, 2, function(x) {
    cor(Y[idx], x[idx])
  })
which(y_cors == max(y_cors))

pathways <- sample(colnames(pathway.scores), 10)
idx_lower <- Y < quantile(Y, 0.25)
idx_upper <- Y >= quantile(Y, 0.75)
avg_expressions_lower <- apply(pathway.scores[idx_lower, pathways], 2, mean)
avg_expressions_upper <- apply(pathway.scores[idx_upper, pathways], 2, mean)
plot_df <- as.data.frame(
  rbind(
    cbind(type = "lower",
          pathway = names(avg_expressions_lower),
          avg_expression = avg_expressions_lower),
    cbind(type = "upper",
          pathway = names(avg_expressions_upper),
          avg_expression = avg_expressions_upper)
  )
)
plot_df$avg_expression <- as.numeric(plot_df$avg_expression)
ggplot(plot_df, aes(x = pathway, y = avg_expression,
                    fill = type)) +
  geom_bar(stat = "identity", position = "dodge2") +
  theme(axis.text.x = element_text(angle = 30))


# 1.1.3 Imaging
scores = as.data.frame(pc_scores)
scores1 = scores %>% pivot_longer(cols=colnames(scores),
                        names_to='combo',
                        values_to='num') 
scores1$combo <- gsub("\\.[0-9]+","", scores1$combo)
scores1
scores1 %>% 
  group_by(combo) %>% 
  summarise(count = n() / 61)


two = as.data.frame(pc_scores)
two1 = two %>% select(starts_with("T2_ED"))

two1 %>% 
  summarise_all(var)


two2 = two1 %>% head(1) %>% t()
two2 = as.data.frame(two2)

paths %>% View()

two2 %>% 
  rename("TCGACS4942" = "TCGA-CS-4942") %>% 
  ggplot() +
  aes(x = TCGACS4942) + geom_boxplot()







