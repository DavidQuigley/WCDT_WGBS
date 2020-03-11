library(survival)
##################
# Survival curves
#################
df_clinical = read.delim(file=fn_surv, sep=',', header=TRUE, stringsAsFactors=FALSE)
df_clinical$Months.to.death.or.last.contact = as.numeric(df_clinical$OS.mCRPC) / 365 * 12
group = rep("Non-CMP", dim(df_clinical)[1])
patient_ids_cmp = paste("DTB-", get.split.col( samples_cmp, "-", col=2), sep='')
group[ match.idx( df_clinical$Patient.ID, patient_ids_cmp )$idx.A ] = "CMP"

df_clinical$Group <- group
df_clinical$Group[df_clinical$Group == "CMP"] <- 1
df_clinical$Group[df_clinical$Group == "Non-CMP"] <- 2

times = df_clinical$Months.to.death.or.last.contact
#had.events = rep(1, length(times) ) #df_clinical$Event
had.events = df_clinical$Event
conditions = as.numeric( df_clinical$Group )
surv.all = survfit(Surv(times, had.events)~conditions)
cox = coxph(formula = Surv(times,had.events==1)~conditions)
p.val = as.numeric(summary(cox)$logtest[3])

pdf( fn_figure_surv, height=5, width=8)
plot(surv.all, lwd=3, col=c("black", "blue"), 
     mark.time=TRUE,
     main="Survival", 
     las=1,
     xlab="months", cex.axis=1.5, bty="n")
legend( 80, 1, c("CMP", "Non-CMP"), fill=c("black", "blue"))
dev.off()

