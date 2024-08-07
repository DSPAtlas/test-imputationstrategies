setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# load data and library ---------------------------------------------------

load("201911_imputation data files.RData")
source("functions_3 proteome.R")

library(dplyr)
library(missForest)
library(purrr)
library(ggplot2)

# generate summary tables -------------------------------------------------

sum_ion <- SumTable(ion)

sum_ludo <- lapply(ion_ludo,SumTable)
names(sum_ludo) <- names(ion_ludo)

sum_lod <- lapply(ion_lod,SumTable)
names(sum_lod) <- names(ion_lod)

sum_nd <- lapply(ion_nd,SumTable)
names(sum_nd) <- names(ion_nd)

sum_rf <- lapply(ion_rf,SumTable)
names(sum_rf) <- names(ion_rf)



# generate result tables --------------------------------------------------

## lod

res_lod <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_lod) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_lod) <- rownames(time.table)
res_lod$time <- time.table$LOD

for(i in 1:90){
  
  res_lod$NRMSE[i] <- missForest::nrmse(ion_lod[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_lod$MAE.B.A[i] <- mean(abs(sum_lod[[i]]$FC.B.A - sum_lod[[i]]$ori.FC.B.A))
  res_lod$MAE.C.A[i] <- mean(abs(sum_lod[[i]]$FC.C.A - sum_lod[[i]]$ori.FC.C.A))
  res_lod$MAE.D.A[i] <- mean(abs(sum_lod[[i]]$FC.D.A - sum_lod[[i]]$ori.FC.D.A))
  res_lod$NRMSE.B.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.B.A,sum_lod[[i]]$FC.B.A)
  res_lod$NRMSE.C.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.C.A,sum_lod[[i]]$FC.C.A)
  res_lod$NRMSE.D.A[i] <- cal.nrmse(sum_lod[[i]]$ori.FC.D.A,sum_lod[[i]]$FC.D.A)
  res_lod$TP.B.A[i] <- sum(sum_lod[[i]]$adj.p.B.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$TP.C.A[i] <- sum(sum_lod[[i]]$adj.p.C.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$TP.D.A[i] <- sum(sum_lod[[i]]$adj.p.D.A < 0.05 & sum_lod[[i]]$species == "ECOLI")
  res_lod$FP.B.A[i] <- sum(sum_lod[[i]]$adj.p.B.A < 0.05 & sum_lod[[i]]$species == "HUMAN")
  res_lod$FP.C.A[i] <- sum(sum_lod[[i]]$adj.p.C.A < 0.05 & sum_lod[[i]]$species == "HUMAN")
  res_lod$FP.D.A[i] <- sum(sum_lod[[i]]$adj.p.D.A < 0.05 & sum_lod[[i]]$species == "HUMAN")

}

res_lod$FADR.B.A <- res_lod$FP.B.A/(res_lod$TP.B.A+res_lod$FP.B.A)
res_lod$FADR.C.A <- res_lod$FP.C.A/(res_lod$TP.C.A+res_lod$FP.C.A)
res_lod$FADR.D.A <- res_lod$FP.D.A/(res_lod$TP.D.A+res_lod$FP.D.A)

## nd

res_nd <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_nd) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_nd) <- rownames(time.table)
res_nd$time <- time.table$ND

for(i in 1:90){
  
  res_nd$NRMSE[i] <- missForest::nrmse(ion_nd[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_nd$MAE.B.A[i] <- mean(abs(sum_nd[[i]]$FC.B.A - sum_nd[[i]]$ori.FC.B.A))
  res_nd$MAE.C.A[i] <- mean(abs(sum_nd[[i]]$FC.C.A - sum_nd[[i]]$ori.FC.C.A))
  res_nd$MAE.D.A[i] <- mean(abs(sum_nd[[i]]$FC.D.A - sum_nd[[i]]$ori.FC.D.A))
  res_nd$NRMSE.B.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.B.A,sum_nd[[i]]$FC.B.A)
  res_nd$NRMSE.C.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.C.A,sum_nd[[i]]$FC.C.A)
  res_nd$NRMSE.D.A[i] <- cal.nrmse(sum_nd[[i]]$ori.FC.D.A,sum_nd[[i]]$FC.D.A)
  res_nd$TP.B.A[i] <- sum(sum_nd[[i]]$adj.p.B.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$TP.C.A[i] <- sum(sum_nd[[i]]$adj.p.C.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$TP.D.A[i] <- sum(sum_nd[[i]]$adj.p.D.A < 0.05 & sum_nd[[i]]$species == "ECOLI")
  res_nd$FP.B.A[i] <- sum(sum_nd[[i]]$adj.p.B.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  res_nd$FP.C.A[i] <- sum(sum_nd[[i]]$adj.p.C.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  res_nd$FP.D.A[i] <- sum(sum_nd[[i]]$adj.p.D.A < 0.05 & sum_nd[[i]]$species == "HUMAN")
  
}

res_nd$FADR.B.A <- res_nd$FP.B.A/(res_nd$TP.B.A+res_nd$FP.B.A)
res_nd$FADR.C.A <- res_nd$FP.C.A/(res_nd$TP.C.A+res_nd$FP.C.A)
res_nd$FADR.D.A <- res_nd$FP.D.A/(res_nd$TP.D.A+res_nd$FP.D.A)

## ludo

res_ludo <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_ludo) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_ludo) <- rownames(time.table)
res_ludo$time <- time.table$LUDO

for(i in 1:90){
  
  res_ludo$NRMSE[i] <- missForest::nrmse(ion_ludo[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_ludo$MAE.B.A[i] <- mean(abs(sum_ludo[[i]]$FC.B.A - sum_ludo[[i]]$ori.FC.B.A))
  res_ludo$MAE.C.A[i] <- mean(abs(sum_ludo[[i]]$FC.C.A - sum_ludo[[i]]$ori.FC.C.A))
  res_ludo$MAE.D.A[i] <- mean(abs(sum_ludo[[i]]$FC.D.A - sum_ludo[[i]]$ori.FC.D.A))
  res_ludo$NRMSE.B.A[i] <- cal.nrmse(sum_ludo[[i]]$ori.FC.B.A,sum_ludo[[i]]$FC.B.A)
  res_ludo$NRMSE.C.A[i] <- cal.nrmse(sum_ludo[[i]]$ori.FC.C.A,sum_ludo[[i]]$FC.C.A)
  res_ludo$NRMSE.D.A[i] <- cal.nrmse(sum_ludo[[i]]$ori.FC.D.A,sum_ludo[[i]]$FC.D.A)
  res_ludo$TP.B.A[i] <- sum(sum_ludo[[i]]$adj.p.B.A < 0.05 & sum_ludo[[i]]$species == "ECOLI")
  res_ludo$TP.C.A[i] <- sum(sum_ludo[[i]]$adj.p.C.A < 0.05 & sum_ludo[[i]]$species == "ECOLI")
  res_ludo$TP.D.A[i] <- sum(sum_ludo[[i]]$adj.p.D.A < 0.05 & sum_ludo[[i]]$species == "ECOLI")
  res_ludo$FP.B.A[i] <- sum(sum_ludo[[i]]$adj.p.B.A < 0.05 & sum_ludo[[i]]$species == "HUMAN")
  res_ludo$FP.C.A[i] <- sum(sum_ludo[[i]]$adj.p.C.A < 0.05 & sum_ludo[[i]]$species == "HUMAN")
  res_ludo$FP.D.A[i] <- sum(sum_ludo[[i]]$adj.p.D.A < 0.05 & sum_ludo[[i]]$species == "HUMAN")
  
}

res_ludo$FADR.B.A <- res_ludo$FP.B.A/(res_ludo$TP.B.A+res_ludo$FP.B.A)
res_ludo$FADR.C.A <- res_ludo$FP.C.A/(res_ludo$TP.C.A+res_ludo$FP.C.A)
res_ludo$FADR.D.A <- res_ludo$FP.D.A/(res_ludo$TP.D.A+res_ludo$FP.D.A)



## rf

res_rf <- data.frame(matrix(NA,nrow = 90, ncol = 17))
colnames(res_rf) <- c("time","NRMSE","MAE.B.A","MAE.C.A","MAE.D.A","NRMSE.B.A","NRMSE.C.A","NRMSE.D.A","TP.B.A","TP.C.A","TP.D.A","FP.B.A","FP.C.A","FP.D.A","FADR.B.A","FADR.C.A","FADR.D.A")
rownames(res_rf) <- rownames(time.table)
res_rf$time <- time.table$RF

for(i in 1:90){
  
  res_rf$NRMSE[i] <- missForest::nrmse(ion_rf[[i]],ion_miss[[i]],ion[rownames(ion_miss[[i]]),])
  res_rf$MAE.B.A[i] <- mean(abs(sum_rf[[i]]$FC.B.A - sum_rf[[i]]$ori.FC.B.A))
  res_rf$MAE.C.A[i] <- mean(abs(sum_rf[[i]]$FC.C.A - sum_rf[[i]]$ori.FC.C.A))
  res_rf$MAE.D.A[i] <- mean(abs(sum_rf[[i]]$FC.D.A - sum_rf[[i]]$ori.FC.D.A))
  res_rf$NRMSE.B.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.B.A,sum_rf[[i]]$FC.B.A)
  res_rf$NRMSE.C.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.C.A,sum_rf[[i]]$FC.C.A)
  res_rf$NRMSE.D.A[i] <- cal.nrmse(sum_rf[[i]]$ori.FC.D.A,sum_rf[[i]]$FC.D.A)
  res_rf$TP.B.A[i] <- sum(sum_rf[[i]]$adj.p.B.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$TP.C.A[i] <- sum(sum_rf[[i]]$adj.p.C.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$TP.D.A[i] <- sum(sum_rf[[i]]$adj.p.D.A < 0.05 & sum_rf[[i]]$species == "ECOLI")
  res_rf$FP.B.A[i] <- sum(sum_rf[[i]]$adj.p.B.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  res_rf$FP.C.A[i] <- sum(sum_rf[[i]]$adj.p.C.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  res_rf$FP.D.A[i] <- sum(sum_rf[[i]]$adj.p.D.A < 0.05 & sum_rf[[i]]$species == "HUMAN")
  
}

res_rf$FADR.B.A <- res_rf$FP.B.A/(res_rf$TP.B.A+res_rf$FP.B.A)
res_rf$FADR.C.A <- res_rf$FP.C.A/(res_rf$TP.C.A+res_rf$FP.C.A)
res_rf$FADR.D.A <- res_rf$FP.D.A/(res_rf$TP.D.A+res_rf$FP.D.A)


## export and save results

sum_result <- list("LOD" = res_lod,
                   "ND" = res_nd,
                   "RF" = res_rf,
                   "LUDO" = res_ludo)
library(openxlsx)
write.xlsx(sum_result,file = "20200329_3 proteome summary results.xlsx")

load("202003_summary and reults.RData")

# boxplot -----------------------------------------------------------------

## boxplot of running time

plot_color <- c("purple","blue","hotpink","green")

df_time <- data.frame("Time" = c(sum_result$LOD[,1],sum_result$ND[,1],sum_result$kNN[,1],sum_result$LLS[,1],sum_result$RF[,1],sum_result$SVD[,1],sum_result$BPCA[,1]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("kNN",90),rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("LOD","ND","kNN","LLS","RF","SVD","BPCA")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_time,aes(x=Method,y=Time,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("Run time (s)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)


## boxplot of NRMSE

df_nrmse <- data.frame("NRMSE" = c(sum_result$LOD[,2],sum_result$ND[,2],sum_result$LUDO[,2],sum_result$RF[,2]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_nrmse,aes(x=Method,y=NRMSE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("NRMSE")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_nrmse2 <- df_nrmse[271:630,]
df_nrmse2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_nrmse2,aes(x=Method,y=NRMSE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("NRMSE")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.45)
plot(p)



## boxplot of MAE B/A

df_maeB <- data.frame(
  "MAE" = c(sum_result$LOD[,3], sum_result$ND[,3], sum_result$LUDO[,3], sum_result$RF[,3]),
  "Method" = factor(c(rep("LOD", 90), rep("ND", 90), rep("LUDO", 90), rep("RF", 90)), levels = c("LOD", "ND", "LUDO", "RF")),
  "MV.ratio" = factor(rep(c(rep("10%MV", 30), rep("20%MV", 30), rep("30%MV", 30)), 4), levels = c("10%MV", "20%MV", "30%MV")),
  "MNAR.ratio" = factor(rep(c(rep("20%MNAR", 10), rep("50%MNAR", 10), rep("80%MNAR", 10)), 12), levels = c("20%MNAR", "50%MNAR", "80%MNAR"))
)


p <- ggplot(df_maeB,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (B/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeB2 <- df_maeB[271:630,]
df_maeB2$Method <- factor(c(rep("LUDO",90),rep("RF",90)),levels = c("LUDO","RF"))

p <- ggplot(df_maeB2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (B/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.18)
plot(p)

## boxplot of MAE C/A

df_maeC <- data.frame("MAE" = c(sum_result$LOD[,4],sum_result$ND[,4],sum_result$LUDO[,4],sum_result$RF[,4]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeC,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (C/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeC2 <- df_maeC[271:630,]
df_maeC2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeC2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (C/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.23)
plot(p)

## boxplot of MAE D/A

df_maeD <- data.frame("MAE" = c(sum_result$LOD[,5],sum_result$ND[,5],sum_result$LUDO[,5],sum_result$RF[,5]),
                      "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                      "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                      "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_maeD,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (D/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)

df_maeD2 <- df_maeD[271:630,]
df_maeD2$Method <- factor(c(rep("LLS",90),rep("RF",90),rep("SVD",90),rep("BPCA",90)),levels = c("kNN","LLS","RF","SVD","BPCA"))

p <- ggplot(df_maeD2,aes(x=Method,y=MAE,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color[4:7])
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("MAE (D/A)")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
#p <- p + ylim(0,0.5)
plot(p)


# Boxplot of fold changes -------------------------------------------------

plot_color <- c("purple","blue","hotpink","green")

i = 1

tmp <- do.call(rbind.data.frame, sum_lod[i])
df_lod <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "LOD")

tmp <- do.call(rbind.data.frame, sum_nd[i])
df_nd <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "ND")

tmp <- do.call(rbind.data.frame, sum_ludo[i])
df_ludo <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                    "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                    "Species" = rep(tmp$species,3),
                    "Method" = "LUDO")

tmp <- do.call(rbind.data.frame, sum_rf[i])
df_rf <- data.frame("value" = c(tmp$FC.B.A,tmp$FC.C.A,tmp$FC.D.A),
                     "Ratio" = c(rep("B/A",nrow(tmp)),rep("C/A",nrow(tmp)),rep("D/A",nrow(tmp))),
                     "Species" = rep(tmp$species,3),
                     "Method" = "RF")


df_fc <- rbind(df_lod,df_nd,df_ludo,df_rf)

p <- ggplot(df_fc,aes(x=Method,y=value,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1)# + geom_point(position = position_jitter(width = 0.1))# 
p <- p +  scale_color_manual(values = plot_color)
p <- p + facet_grid(Species~Ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("Protein Ratio")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
p <- p + ylim(0,3)
plot(p)





# ROC plot ----------------------------------------------------------------

#load("201911_summary and reults.RData")
load("201911_imputation data files.RData")
load("202003_summary and reults.RData")


library(pROC)

get_roc <- function(a){
  
  tmp <- a[,15:19]
  tmp$DE <- grepl("ECOLI|YEAST", tmp$species)
  #tmp$DE <- grepl("ECOLI", tmp$species)
  tmp$BG <- grepl("HUMAN", tmp$species)
  
  for(i in 1:3){
    
    tmp2 <- tmp[,c(i,6,7)]
    tmp2 <- tmp2[order(tmp2[,1],decreasing = F),]
    tmp2$TPR <- cumsum(as.numeric(tmp2$DE)) / sum(tmp$DE)
    tmp2$FPR <- cumsum(as.numeric(tmp2$BG)) / sum(tmp$BG)
    tmp <- cbind(tmp,tmp2[rownames(tmp),4:5])
    
  }
  
  colnames(tmp)[8:13] <- paste0(colnames(tmp)[8:13],".",c("B.A","B.A","C.A","C.A","D.A","D.A"))
  
  return(tmp)
}

get_auc <- function(df,rank = "adj.p.B.A"){
  
  return(as.numeric(auc(roc(df[order(df[,grep(rank,colnames(df))],decreasing = F),"DE"],rev(seq(1,nrow(df))),direction = "<"))))
  
}

roc_ion <- get_roc(sum_ion)
roc_lod <- lapply(sum_lod,get_roc)
roc_nd <- lapply(sum_nd,get_roc)
roc_ludo <- lapply(sum_ludo,get_roc)
roc_rf <- lapply(sum_rf,get_roc)


get_auc(roc_ion,rank = "adj.p.B.A")  #0.9146457
get_auc(roc_ion,rank = "adj.p.C.A")  #0.9634159
get_auc(roc_ion,rank = "adj.p.D.A")  #0.9719284

auc_lod <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_lod,get_auc, rank="adj.p.D.A"))
)

auc_nd <- data.frame("dataset" = rownames(time.table),
                     "AUC.B.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.B.A")),
                     "AUC.C.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.C.A")),
                     "AUC.D.A" = unlist(lapply(roc_nd,get_auc, rank="adj.p.D.A"))
)

auc_ludo <- data.frame("dataset" = rownames(time.table),
                      "AUC.B.A" = unlist(lapply(roc_ludo,get_auc, rank="adj.p.B.A")),
                      "AUC.C.A" = unlist(lapply(roc_ludo,get_auc, rank="adj.p.C.A")),
                      "AUC.D.A" = unlist(lapply(roc_ludo,get_auc, rank="adj.p.D.A"))
)


auc_rf <- data.frame("dataset" = rownames(time.table),
                     "AUC.B.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.B.A")),
                     "AUC.C.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.C.A")),
                     "AUC.D.A" = unlist(lapply(roc_rf,get_auc, rank="adj.p.D.A"))
)


auc_result <- list("LOD" = auc_lod,
                   "ND" = auc_nd,
                   "LUDO" = auc_ludo,
                   "RF" = auc_rf)
library(openxlsx)
#write.xlsx(auc_result,file = "201911_3 proteome auc results.xlsx")
write.xlsx(auc_result,file = "20200329_3 proteome auc results.xlsx")

## boxplot of AUCs

plot_color <- c("purple","blue","hotpink","green")

## B vs A

df_auc_B <- data.frame("AUC" = c(auc_result$LOD[,2],auc_result$ND[,2],auc_result$LUDO[,2],auc_result$RF[,2]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_B,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)



## C vs A

df_auc_C <- data.frame("AUC" = c(auc_result$LOD[,3],auc_result$ND[,3],auc_result$LUDO[,3],auc_result$RF[,3]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_C,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)


## D vs A

df_auc_D <- data.frame("AUC" = c(auc_result$LOD[,4],auc_result$ND[,4],auc_result$LUDO[,4],auc_result$RF[,4]),
                       "Method" = factor(c(rep("LOD",90),rep("ND",90),rep("LUDO",90),rep("RF",90)),levels = c("LOD","ND","LUDO","RF")),
                       "MV.ratio" = factor(c(rep("10%MV",30),rep("20%MV",30),rep("30%MV",30)),levels = c("10%MV","20%MV","30%MV")),
                       "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",10),rep("50%MNAR",10),rep("80%MNAR",10)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR"))
)

p <- ggplot(df_auc_D,aes(x=Method,y=AUC,color=Method))
p <- p + geom_boxplot(width = 0.5,lwd = 1,outlier.shape = NA) + geom_point(position = position_jitter(width = 0.1)) +  scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + xlab("Imputation methods") + ylab("AUC")
#p <- p + guides(fill=guide_legend(title="Legend_Title"))
p <- p + theme_bw(base_size = 16)# + labs(x="Isoform", y = "Log2 Isoform Intensity")
p <- p + theme(axis.text.x=element_text(angle = 45, vjust = 1,hjust = 1))
plot(p)


## plot AUC curve

idx <- 1
sum.pro <- as.numeric(unlist(lapply(roc_lod[seq(idx,idx+89,10)],nrow)))

roc.fun <- function(a,b = "TPR",idx = idx){
  
  tmp <- do.call(rbind.data.frame, a[seq(idx,idx+89,10)])
  
  if(b == "TPR"){
    return(c(tmp$TPR.B.A,tmp$TPR.C.A,tmp$TPR.D.A))
  }
  
  if(b == "FPR"){
    return(c(tmp$FPR.B.A,tmp$FPR.C.A,tmp$FPR.D.A))
  }
  
}

df1 <- data.frame("TPR" = roc.fun(roc_lod,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_lod,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "LOD")


df2 <- data.frame("TPR" = roc.fun(roc_nd,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_nd,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "ND")

df3 <- data.frame("TPR" = roc.fun(roc_ludo,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_ludo,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "LUDO")


df5 <- data.frame("TPR" = roc.fun(roc_rf,b = "TPR", idx = 1),
                  "FPR" = roc.fun(roc_rf,b = "FPR", idx = 1),
                  "group" = c(rep("B.A",sum(sum.pro)),rep("C.A",sum(sum.pro)),rep("D.A",sum(sum.pro))),
                  "MV.ratio" = factor(c(rep("10%MV",sum(sum.pro[1:3])),rep("20%MV",sum(sum.pro[4:6])),rep("30%MV",sum(sum.pro[7:9]))),levels = c("10%MV","20%MV","30%MV")),
                  "MNAR.ratio" = factor(c(rep(c(rep("20%MNAR",6689),rep("50%MNAR",6689),rep("80%MNAR",6689)),3)), levels = c("20%MNAR","50%MNAR","80%MNAR")),
                  "method" = "RF")

df.all <- rbind(df1,df2,df3,df5)

library(ggplot2)

p <- ggplot(df.all[df.all$group == "D.A",], aes(FPR, TPR, col = method))
p <- p + geom_line() + scale_color_manual(values = plot_color)
p <- p + facet_grid(MNAR.ratio ~ MV.ratio, scales="free")
p <- p + theme_bw(base_size = 14) + ggtitle("ROC-curve") + xlim(0, 0.1)
plot(p)





