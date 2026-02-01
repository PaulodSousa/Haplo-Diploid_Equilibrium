## Basic stuf
## Assuming equal sex ratio and mendelian inhertence 



geno.data <- c("0/0", "0/0", "1/1", "1/1", "0/0", "1", "0")

allele.freq <- function(geno.data) {
  
  # 1st step: count each genotype
  AA.f <- sum(geno.data %in% "0/0", na.rm = TRUE)
  aa.f <- sum(geno.data %in% "1/1", na.rm = TRUE)
  Aa.f <- sum(geno.data %in% c("1/0", "0/1"), na.rm = TRUE)
  A.m <- sum(geno.data %in% "0", na.rm = TRUE)
  a.m <- sum(geno.data %in% "1", na.rm = TRUE)
  
  # 2nd step: Calculate number of females and males and total sample number
  N.F <- AA.f + aa.f + Aa.f
  N.M <- A.m + a.m
  Total_Samples <- N.F + N.M
  
  # 3rd step: Calculate allele frequencies
  F.A <- (2*AA.f + Aa.f + A.m) / (N.F * 2 + N.M)
  F.a <- (2*aa.f + Aa.f + a.m) / (N.F * 2 + N.M)
  
  # 4th step: Calculate expected genotype frequencies (assuming an equal sex ratio)
  Exp.AA <- (F.A * F.A) * 0.5
  Exp.aa <- (F.a * F.a) * 0.5
  Exp.Aa <- 2 * (F.A * F.a) * 0.5
  Exp.A <- F.A  * 0.5
  Exp.a <- F.a  * 0.5
 
  results <- data.frame(
     N_samples = Total_Samples,
     Prop_Males = (N.M) / Total_Samples,
     Freq.Ref = F.A,
     Freq.Alt = F.a,
     Obs.AA = AA.f  / Total_Samples,
     Obs.aa = aa.f / Total_Samples,
     Obs.Het = Aa.f / Total_Samples,
     Obs.M.Ref = A.m / Total_Samples,
     Obs.M.Alt = a.m / Total_Samples,
     Exp.AA = Exp.AA,
     Exp.aa = Exp.aa,
     Exp.Het = Exp.Aa,
     Exp.M.Ref = Exp.A,
     Exp.M.Alt = Exp.a
     )
}

HWE.freq <- function(geno.data) {
  
  # 1st step: count each genotype
  AA.f <- sum(geno.data %in% "0/0", na.rm = TRUE)
  aa.f <- sum(geno.data %in% "1/1", na.rm = TRUE)
  Aa.f <- sum(geno.data %in% c("1/0", "0/1"), na.rm = TRUE)
  
  # 2nd step: Calculate number samples 
  Total_Samples <- AA.f + aa.f + Aa.f
  
  # 3rd step: Calculate allele frequencies
  F.A <- (2*AA.f + Aa.f) / (Total_Samples * 2)
  F.a <- (2*aa.f + Aa.f) / (Total_Samples * 2)
  
  # 4th step: Calculate expected genotype frequencies (assuming an equal sex ratio)
  Exp.AA <- (F.A * F.A)
  Exp.aa <- (F.a * F.a)
  Exp.Aa <- 2 * (F.A * F.a)
  
  results <- data.frame(
    N_samples = Total_Samples,
    Freq.Ref = F.A,
    Freq.Alt = F.a,
    Obs.AA = AA.f  / Total_Samples,
    Obs.aa = aa.f / Total_Samples,
    Obs.Het = Aa.f / Total_Samples,
    Obs.Hom = Aa.f + aa.f,
    Exp.AA = Exp.AA,
    Exp.aa = Exp.aa,
    Exp.Het = Exp.Aa,
    Exp.Hom = Exp.AA + Exp.aa
  )
}

dt <- as.data.frame(allele.freq(geno.data = geno.data))
dt



dt_list <- list()

set.seed(1234)

for (i in 1:10000) {
  # Randomly choose size of vector
  N <- sample(500:2000, 1)
  
  # Generate normalized probabilities
  p <- round(runif(2, min = 0, max = 1), 5)
  p <- p / sum(p)
  
  probs <- c(p[1]*p[1]*0.5, p[1]*p[2], p[2]*p[2]*0.5, p[1]*0.5, p[2]*0.5)
  
  # Sample genotypes of size N
  geno <- sample(c("0/0", "0/1", "1/1", "0", "1"), N, replace = TRUE, prob = probs)
  
  # Pad with NAs to length 30
  geno_padded <- c(geno, rep(NA, 2000 - N))
  
  # Save to list
  dt_list[[i]] <- geno_padded
}

# Combine all into a data frame
dt <- as.data.frame(do.call(rbind, dt_list))

remove(dt_list, geno, geno_padded, i, N, p)

results_list <- list()

results_list <- lapply(1:nrow(dt), function(i) {
  allele.freq(dt[i, ])
})
results_df <- as.data.frame(do.call(rbind, results_list))

library(dplyr)
results_df %>% summarise(Mean = round(mean(Exp.M.Alt), 3), SD= round(sd(Exp.M.Alt), 3),
                         Min = round(min(Exp.M.Alt), 3), Max = round(max(Exp.M.Alt), 3))

results_df$Exp.Hom <- results_df$Exp.AA + results_df$Exp.aa 

library(ggplot2)
Allele.freq.plot <- 
  ggplot(results_df) +
    geom_histogram(aes(x= Freq.Ref),color="darkblue", fill="#56B1F7")  +
    geom_histogram(aes(x= Freq.Alt), color="darkred", fill="#F8766D",  alpha= 0.5)  +
    #annotate(geom="text", x=0.15, y=760, label="10k simulations", size=3) +
    #annotate(geom="text", x=0.15, y=710, label="N [500 - 2000]", size=3) +
    ggtitle("") + xlab("Frequency") + ylab("Number of simulations") +
    theme(plot.title = element_text(),
          axis.title.y = element_text(size= 14),
          axis.title.x = element_text(size= 14),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white", color="grey7")
          )
Female.Geno.plot <-
  ggplot(results_df) +
    geom_histogram(aes(x= Exp.Hom),color="darkblue", fill="#56B1F7")  +
    geom_histogram(aes(x= Exp.Het), color="darkred", fill="#F8766D",  alpha= 0.5)  +
    #annotate(geom="text", x=0.05, y=2200, label="10k simulations", size=3) +
    #annotate(geom="text", x=0.05, y=2000, label="N [500 - 2000]", size=3) +
    ggtitle("") + xlab("Frequency") + ylab("Number of simulations") +
    theme(plot.title = element_text(),
          axis.title.y = element_text(size= 14),
          axis.title.x = element_text(size= 14),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white", color="grey7")
    )
Male.Geno.plot <- 
  ggplot(results_df) +
    geom_histogram(aes(x= Exp.M.Ref),color="darkblue", fill="#56B1F7")  +
    geom_histogram(aes(x= Exp.M.Alt), color="darkred", fill="#F8766D",  alpha= 0.5)  +
    #annotate(geom="text", x=0.05, y=690, label="10k simulations", size=3) +
    #annotate(geom="text", x=0.05, y=650, label="N [500 - 2000]", size=3) +
    ggtitle("") + xlab("Frequency") + ylab("Number of simulations") +
    theme(plot.title = element_text(),
          axis.title.y = element_text(size= 14),
          axis.title.x = element_text(size= 14),
          axis.text.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=10),
          plot.background = element_rect(fill = "white"), 
          panel.background = element_rect(fill = "white", color="grey7")
    )
library(gridExtra)
grid.arrange(Allele.freq.plot, Female.Geno.plot, Male.Geno.plot, nrow=3)




ggplot(results_df, aes(x= Exp.Het))+
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept= mean(Exp.Het)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.05, y=4000, label="10k simulations", size=3) +
  annotate(geom="text", x=0.05, y=3600, label="N [500 - 2000]", size=3) +
  ggtitle("") + xlab("Expected heterozygosity") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_df, aes(x= Exp.AA))+
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept= mean(Exp.AA)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.35, y=3000, label="10k simulations", size=3) +
  annotate(geom="text", x=0.35, y=2600, label="N [500 - 2000]", size=3) +
  ggtitle("") + xlab("Expected homozygosity") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_df, aes(x= Exp.M.Ref))+
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept= mean(Exp.M.Ref)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.05, y=900, label="10k simulations", size=3) +
  annotate(geom="text", x=0.05, y=850, label="N [500 - 2000]", size=3) +
  ggtitle("") + xlab("Expected males reference allele") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_df, aes(x= N_samples, y= Exp.Het, color= Prop_Males)) +
    geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
    annotate(geom="text", x=600, y=0.28, label="10k simulations", size=3) +
    annotate(geom="text", x=600, y=0.26, label="N [500 - 2000]", size=3) +
    ggtitle("") + xlab("Number of samples") + ylab("Expected Heterozygosity") +
    theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )


cor.test(results_df$N_samples, results_df$Exp.Het)
cor.test(results_df$N_samples, results_df$Prop_Males)
cor.test(results_df$Prop_Males, results_df$Exp.Het)

ggplot(results_df, aes(x= N_samples, y= Exp.AA, color= Prop_Males)) +
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  annotate(geom="text", x=600, y=0.53, label="10k simulations", size=3) +
  annotate(geom="text", x=600, y=0.51, label="N [500 - 2000]", size=3) +
  ggtitle("") + xlab("Number of samples") + ylab("Expected frequency of AA") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )
cor.test(results_df$N_samples, results_df$Exp.AA)
cor.test(results_df$Prop_Males, results_df$Exp.AA)

ggplot(results_df, aes(x= N_samples, y= Exp.M.Ref, color= Prop_Males)) +
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  annotate(geom="text", x=600, y=0.53, label="10k simulations", size=3) +
  annotate(geom="text", x=600, y=0.51, label="N [500 - 2000]", size=3) +
  ggtitle("") + xlab("Number of samples") + ylab("Expected males reference allele") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )
cor.test(results_df$N_samples, results_df$Exp.M.Ref)
cor.test(results_df$Prop_Males, results_df$Exp.M.Ref)

################################################################################
# Tests for sampling sizes and biases
set.seed(1234)

geno.original <- sample(c("0/0", "0/1", "1/1", "0", "1"), 100000, replace = TRUE, 
                        prob = c(1/8, 1/4, 1/8, 1/4, 1/4 ))

complete <- allele.freq(geno.original)
complete$Exp.Hom <- complete$Exp.AA + complete$Exp.aa

# Define probabilities for each character
char_probs <- c("0/0" = 1/4, "0/1" = 1/2, "1/1" = 1/4, "0" = 0, "1" = 0)

# Map these probabilities to the large vector
element_probs <- char_probs[geno.original]

geno.subs <- list()
for (i in 1:100000) {
  # size of sampling
  N <- sample(2:20, 1)
  
  #Sample 10 characters from the large vector using those mapped probabilities
  sub.females <- sample(geno.original, size = N, replace = TRUE, prob = element_probs)
  
  # probability of sampling genotype X is its frequency on the population
  #sub <- sample(geno.original, N, replace = TRUE)
 
  # Pad with NAs to length 20
  sub_padded <- c(sub.females, rep(NA, 20 - N))
  geno.subs[[i]] <- sub_padded
}
sub.dt <- as.data.frame(do.call(rbind, geno.subs))

remove(sub_padded, sub, i, N, geno.subs, sub.females, char_probs, element_probs)

results_list <- list()

results_list <- lapply(1:nrow(sub.dt), function(i) {
  #allele.freq(sub.dt[i, ])
  HWE.freq(sub.dt[i,]) 
})
results_sub_df <- as.data.frame(do.call(rbind, results_list))
remove(results_list, sub.dt)

results_sub_df$Exp.Hom <- results_sub_df$Exp.AA + results_sub_df$Exp.aa

library(ggplot2)
ggplot(results_sub_df, aes(x = Freq.Ref)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  geom_vline(data = complete, aes(xintercept = Freq.Ref), color = "red", linetype = "solid") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Frequency reference allele") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_sub_df, aes(x = Exp.Het)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  geom_vline(data = complete, aes(xintercept = Exp.Het), color = "red", linetype = "solid") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Expected Heterozygosity") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_sub_df, aes(x = Exp.Hom)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  geom_vline(data = complete, aes(xintercept = Exp.Hom), color = "red", linetype = "solid") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Expected Homozygosity") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )


ggplot(results_sub_df, aes(x = Exp.M.Ref)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  geom_vline(data = complete, aes(xintercept = Exp.M.Ref), color = "red", linetype = "solid") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Expected Male reference allele") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_sub_df, aes(x = Prop_Males)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  geom_vline(data = complete, aes(xintercept = Prop_Males), color = "red", linetype = "solid") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Proportion of Males") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_sub_df, aes(x = N_samples)) +
  geom_histogram(color="darkblue", fill="lightblue") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Number of samples") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

results_sub_df_diffs <- data.frame(
  N_Samples = results_sub_df$N_samples,
  #Prop_Males = results_sub_df$Prop_Males,
  Diff_Ref.Allele = complete$Freq.Ref - results_sub_df$Freq.Ref,
  Diff_Exp.Het = complete$Exp.Het - results_sub_df$Exp.Het,
  Diff_Exp.Hom = complete$Exp.Hom - results_sub_df$Exp.Hom#,
  #Diff_M.Ref.Allele = complete$Exp.M.Ref - results_sub_df$Exp.M.Ref
)

ggplot(results_sub_df_diffs, aes(x = N_Samples, y= Diff_Ref.Allele, color= Prop_Males)) + 
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Number of samples") + ylab("Delta of reference allele frequency") +
  theme_bw() +          # white background with grid lines
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        panel.grid.major = element_line(color = "grey80"),  # major grid lines
        panel.grid.minor = element_line(color = "grey90")   # minor grid lines
  )


ggplot(results_sub_df_diffs, aes(x = N_Samples, y= Diff_Exp.Het, color= Prop_Males)) + 
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Number of samples") + ylab("Delta of expected heterozygosity") +
  theme_bw() +          # white background with grid lines
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        panel.grid.major = element_line(color = "grey80"),  # major grid lines
        panel.grid.minor = element_line(color = "grey90")   # minor grid lines
  )

ggplot(results_sub_df_diffs, aes(x = N_Samples, y= Diff_Exp.Hom, color= Prop_Males)) + 
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Number of samples") + ylab("Delta of expected homozygosity") +
  theme_bw() +          # white background with grid lines
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        panel.grid.major = element_line(color = "grey80"),  # major grid lines
        panel.grid.minor = element_line(color = "grey90")   # minor grid lines
  )

ggplot(results_sub_df_diffs, aes(x = N_Samples, y= Diff_M.Ref.Allele, color= Prop_Males)) + 
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  ggtitle("Comparison between population and sample of 10") + 
  xlab("Number of samples") + ylab("Delta of Male reference allele frequency") +
  theme_bw() +          # white background with grid lines
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        panel.grid.major = element_line(color = "grey80"),  # major grid lines
        panel.grid.minor = element_line(color = "grey90")   # minor grid lines
  )

# Mean square error
mean((results_sub_df$Exp.M.Ref - complete$Exp.M.Ref)^2)
