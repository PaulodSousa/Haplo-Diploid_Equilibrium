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


dt <- as.data.frame(allele.freq(geno.data = geno.data))
dt



dt_list <- list()

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
