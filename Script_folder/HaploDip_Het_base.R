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
     Obs.Hom = (AA.f + aa.f) / Total_Samples,
     Obs.Het = Aa.f / Total_Samples,
     Obs.M.Ref = A.m / Total_Samples,
     Exp.Hom = (Exp.AA + Exp.aa),
     Exp.Het = Exp.Aa,
     Exp.M.Ref = Exp.A
     )
}


dt <- as.data.frame(allele.freq(geno.data = geno.data))
dt



dt_list <- list()

for (i in 1:10000) {
  # Randomly choose size of vector
  N <- sample(500:2000, 1)
  
  # Generate normalized probabilities
  p <- round(runif(5, min = 0, max = 1), 3)
  p <- p / sum(p)
  
  # Sample genotypes of size N
  geno <- sample(c("0/0", "0/1", "1/1", "0", "1"), N, replace = TRUE, prob = p)
  
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

library(ggplot2)
ggplot(results_df, aes(x= Freq.Alt))+
    geom_histogram(color="darkblue", fill="lightblue")  +
    geom_vline(aes(xintercept= mean(Freq.Ref)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.15, y=760, label="10k simulations", size=3) +
  annotate(geom="text", x=0.15, y=710, label="N [500 - 2000]", size=3) +
  annotate(geom="text", x=0.15, y=660, label="Sex ratio [0 - 1]", size=3) +
  ggtitle("") + xlab("Frequency reference allele") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
        )

ggplot(results_df, aes(x= Exp.Het))+
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept= mean(Exp.Het)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.05, y=4000, label="10k simulations", size=3) +
  annotate(geom="text", x=0.05, y=3600, label="N [500 - 2000]", size=3) +
  annotate(geom="text", x=0.05, y=3200, label="Sex ratio [0 - 1]", size=3) +
  ggtitle("") + xlab("Expected heterozygosity") + ylab("Count") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )

ggplot(results_df, aes(x= Exp.Hom))+
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_vline(aes(xintercept= mean(Exp.Hom)),
             color="black", linetype="dashed", size=1) +
  annotate(geom="text", x=0.35, y=3000, label="10k simulations", size=3) +
  annotate(geom="text", x=0.35, y=2600, label="N [500 - 2000]", size=3) +
  annotate(geom="text", x=0.35, y=2200, label="Sex ratio [0 - 1]", size=3) +
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
  annotate(geom="text", x=0.05, y=810, label="Sex ratio [0 - 1]", size=3) +
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
    annotate(geom="text", x=600, y=0.09, label="10k simulations", size=3) +
    annotate(geom="text", x=600, y=0.08, label="N [500 - 2000]", size=3) +
    annotate(geom="text", x=600, y=0.07, label="Sex ratio [0 - 1]", size=3) +
    ggtitle("") + xlab("Number of samples") + ylab("Expected heterozygosity") +
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

ggplot(results_df, aes(x= N_samples, y= Exp.Hom, color= Prop_Males)) +
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  annotate(geom="text", x=600, y=0.47, label="10k simulations", size=3) +
  annotate(geom="text", x=600, y=0.45, label="N [500 - 2000]", size=3) +
  annotate(geom="text", x=600, y=0.43, label="Sex ratio [0 - 1]", size=3) +
  ggtitle("") + xlab("Number of samples") + ylab("Expected homozygosity") +
  theme(plot.title = element_text(),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_text(size= 14),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", color="grey7")
  )
cor.test(results_df$N_samples, results_df$Exp.Hom)
cor.test(results_df$Prop_Males, results_df$Exp.Hom)

ggplot(results_df, aes(x= N_samples, y= Exp.M.Ref, color= Prop_Males)) +
  geom_point() + scale_color_gradient(low = "#56B1F7", high = "#F8766D") +
  annotate(geom="text", x=600, y=0.6, label="10k simulations", size=3) +
  annotate(geom="text", x=600, y=0.56, label="N [500 - 2000]", size=3) +
  annotate(geom="text", x=600, y=0.52, label="Sex ratio [0 - 1]", size=3) +
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
