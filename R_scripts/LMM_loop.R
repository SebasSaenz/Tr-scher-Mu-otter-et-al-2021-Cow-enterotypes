#Analysis of Linear and Mixed linear models

#install packages
install.packages("ade4")
install.packages("MuMIn")
install.packages("lmerTest")
install.packages("Hmisc")

#load libraries
library(tidyverse)
library(ade4)
library(MuMIn)
library(lmerTest)
library(Hmisc)

#read file
df  <- read.csv("Biocrates_normalized_for_R_Metabotype.CSV", header = TRUE)


df_2 <- df[15:191] %>% select_if(colSums(.) != 0)
df_2 <- df[, colSums(df != 0, na.rm = TRUE) > 0]


#Transfor some variables to factors
df$Group <- factor(df$Group, levels = c("CON", "CAR"))
df$Age <- factor(df$Age, levels = c(3, 4, 5, 7))
df$Time.point <- factor(df$Time.point, levels = c("72hL", "p126", "-42",  "12hL", "p100",  "12hC", "72hC", "m42"))
df$animal <- as.character(df$animal)

#remove NA, it is neccesary for LinearMixModels
nona_df <- na.omit(df)

#check structure of data
str(df)

#define the columns to do the regression, change the dataframe if you are going
#to use linear moderls

y <- df[202]
x <- df_2[1:177]
animals <- df[3]


######################################################################
################### LINEAR MIX MODEL   #########################

#X <- summary(lme4::lmer(df$H..loge. ~ df$Ala + (1|df$Animal), data = df, REML = FALSE))
#p.z <- 2 * (1 - pnorm(abs(X$coefficients[2, 3])))


lmr_out <- data.frame(NULL)              # create object to keep results

#loop to apply the function to all the columns in the dataset
#this loop would give us the coefficients and p-vals
#Here the animlas are included as a random effect
for (i in 1:length(x)) {
  summary_lmr <- summary(lme4::lmer(y[,1] ~ x[, i] + df$Age + df$Group + df$Time.point + (1|df$animal)))    # run model
  lmr_out[i, 1] <- names(x)[i]           # print variable name
  lmr_out[i, 2] <- summary_lmr$coefficients[1,1]   # intercept
  lmr_out[i, 3] <- summary_lmr$coefficients[2,1]
  lmr_out[i, 4] <- 2 * (1 - pnorm(abs(summary_lmr$coefficients[2, 3])))
}
#set the names
names(lmr_out) <- c("y.variable", "intercept", "Coefficient", "p-val")


rsq_out <- data.frame(NULL)              # create object to keep results

#loop to apply the function to all the columns in the dataset
#this loop would give us Rsquares values
#Here the animlas are included as a random effect
for (i in 1:length(x)) {
  rsq_lmr <- r.squaredGLMM(lme4::lmer(y[,1] ~ x[, i] + df$Age + df$Group + df$Time.point + (1|df$animal)))    # run model
  rsq_out[i, 1] <- names(x)[i]           # print variable name
  rsq_out[i, 2] <- rsq_lmr[1,1]   # intercept
  rsq_out[i, 3] <- rsq_lmr[1,2]
}


#set the names
names(rsq_out) <- c("y.variable", "Rmarginal", "Rconditional")

#join both df
full_out <- full_join(lmr_out, rsq_out, by ="y.variable")

#adjust the p.vals 
full_out$p.adju <- p.adjust(full_out$`p-val`, method = "fdr", length(full_out$`p-val`))

#filter the p-vals <0.05
full_out_filter <- full_out %>%
  filter(p.adju < 0.05)

#save full table
write.csv(full_out, file = "Newentero/LMM_output_olse.csv")

#Plot

df_plot <- read.csv("LMM_output_clos_bifi_plot.csv", header = TRUE)


sort_samples <- df_plot  %>%
  filter(Taxa == "Bifidobacterium") %>%
  arrange((Value)) %>%
  pull(y.variable) %>%
  unique()

sort_samples_meta <- df_plot %>%
  arrange((Metabolite)) %>%
  pull(y.variable) %>%
  unique()

df_plot %>%
  mutate(y.variable = factor(y.variable, levels = sort_samples_meta)) %>%
  ggplot() +
  geom_bar(aes(x = y.variable,
               y = Value,
               fill = Taxa),
           stat = "identity") +
  coord_flip() +
  ylab("Percentage of variance explained in the Relative abundance") +
  xlab("") +
  scale_y_continuous(breaks = seq(-30, 30, by =5),
                     limits = c(-30, 30)) +
  scale_fill_manual(values = c("darkgreen", "blue")) +
  theme_classic()

########################################################################################
######################## CORRELATION ###################################################


#read file
cor_df <- read.csv("Biocrates_normalized_for_R_Metabotype.csv",
                   header = TRUE)
#filter the df by columns
cor_df_filter <- cbind(cor_df[5:7], cor_df[16:192])

#run correlation
correaltion_mat <- rcorr(as.matrix(cor_df_filter), type = "pearson")


##########Function for creating a nice and tidy dataframe  from cor value and p value
#http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#extract values using the created fucntion
cor_results <- flattenCorrMatrix(correaltion_mat$r, correaltion_mat$P)

#adjust the p.vals 
cor_results$p.adju <- p.adjust(cor_results$p, method = "fdr", length(cor_results$p))


#filter result by enterotype and p adjusted value
cor_results %>%
  filter(row == "E.Clos") %>%
  filter(p.adju < 0.05)
