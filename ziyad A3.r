library(stringr) #assist with text manipulation
library(dplyr) # data manipulation
library(readr) # data input
library(readr) #assist with text manipulation
library(dplyr)

cancer=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/77_cancer_proteomes_CPTAC_itraq.csv")
clinc=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/clinical_data_breast_cancer.csv")
pam=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/PAM50_proteins.csv")


#-------------------------------------------------Handling Data ------------------------------------------------------------
for (i in 4:86) {
  temp=paste("TCGA-",substr(colnames(cancer)[i],0,2),"-",substr(colnames(cancer)[i],4,7),sep="")
  colnames(cancer)[i] = temp
}
cancer <- cancer[,!duplicated(colnames(cancer))]
cancer <- cancer[4:80]
#Subset the rows of proteomes table by removing rows with missing values
no_clinical <- na.omit(clinc)
no_cancer <- na.omit(cancer)
cancer=cancer[complete.cases(cancer), ]
head(cancer)




trans <- as.data.frame(t(cancer))
write.csv(trans,"E:/university/fourth year first term/biostat/assignments/assignment3/trans.csv")
read <- read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/trans.csv",header = T)
join <- inner_join(read, no_clinical, by = c("X" = "Complete.TCGA.ID"))
write.csv(join,"E:/university/fourth year first term/biostat/assignments/assignment3/JOIN.csv")
read2 <- read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/JOIN.csv",header = T)


#--------------------------------------------------- part 2----------------------------------------------------------------
cancer_2 <- read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/77_cancer_proteomes_CPTAC_itraq.csv",header = T)
for (i in 4:83) {
  result = paste("TCGA-",substr(colnames(cancer_2)[i],0,2),sep = "")
  result = paste(result,"-",sep = "")
  result = paste(result,substr(colnames(cancer_2)[i],4,7),sep = "")
  #print(result)
  colnames(cancer_2)[i] = result
}
cancer_2 <- cancer_2[,!duplicated(colnames(cancer_2))]
joi <- inner_join(cancer_2,pam,by = c("RefSeq_accession_number"= "RefSeqProteinID"))
# print(na.omit(joi[1:23,c(9)]))
#joi[1:78]
#head(joi)
trans2=as.data.frame(t(joi))
write.csv(trans2,"E:/university/fourth year first term/biostat/assignments/assignment3/Trans2.csv")
#head(trans2)
read2=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/Trans2.csv", check.names = FALSE)
#head(read2)
#dim(read2)
read2=read2[-c(2,3,4,84,85,86),-c(1)]
read2=read2[ , colSums(is.na(read2)) == 0]
#dim(read2)
write.csv(read2,"E:/university/fourth year first term/biostat/assignments/assignment3/Proteins.csv")
protein=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/Proteins.csv", check.names = FALSE)
colnames(protein)<-NULL
names(protein) <- lapply(protein[1, ], as.character)
protein= protein[-1,]
#tail(read2)
#colnames(protein)
#head(protein)
#dim(protein)
HER2=as.numeric(clinc$HER2.Final.Status)
#HER2
mydata2=data.frame(protein)
mydata2$HER2.Final.Status<-clinc$HER2.Final.Status[1:79]
write.csv(mydata2,"E:/university/fourth year first term/biostat/assignments/assignment3/DataWithHER2.csv")
DataWithHER2=read.csv("E:/university/fourth year first term/biostat/assignments/assignment3/DataWithHER2.csv", check.names = FALSE)
#print(mydata2)
#na.omit(joi[1:104,11])


#2.1
correlations <- cor(protein[4] %>% type.convert(as.is=TRUE),HER2[1:79])
#print(correlations)
for (i in 2:length(protein)) {
  #z=colSums(!is.na(joi[4:83]))
  #print(protein[i])
  correlations[i] <- cor(protein[i] %>% type.convert(as.is=TRUE),HER2[1:79])
  Cor=c(correlations[i],colnames(protein[i]),HER2[i])
  print(Cor)
}


datafram=data.frame(correlations,colnames(protein))
#2.2
datafram[order(datafram$correlations,decreasing = TRUE),]



#2.3
threshold=0.013
for (k in 1:26) {
  if (datafram$correlations[k]<threshold) {
    print(c(as.character(datafram$colnames.protein.[k]),datafram$correlations[k]))
  }  
}

#-------------------------------------part3---------------------------------------------
#3.1
positive=subset(DataWithHER2,DataWithHER2$HER2.Final.Status=="Positive")
#print(positive)
write.csv(positive,"E:/university/fourth year first term/biostat/assignments/assignment3/Positive.csv")
negative=subset(DataWithHER2,DataWithHER2$HER2.Final.Status=="Negative")
#print(negative[3])
positive=positive[,-c(1,2,29)]
#print(positive)
#length(positive)
negative=negative[,-c(1,2,29)]
#print(negative)
#length(negative)

mat= matrix(ncol = 0, nrow = 26)
SortedValues=data.frame(mat)
SortedValues$T_Test<-0
SortedValues$ProteinName<-1
for (i in 3:length(positive)) {
  Test=t.test(positive[,i],negative[,i])
  #Data[[length(Data) + 1]] <- list(A=c(Test$statistic),B=c(colnames(positive[i])))
  #my_list[[length(my_list) + 1]] <- Test$statistic
  #my_list2[[length(my_list2) + 1]] <-colnames(positive[i])
  SortedValues$T_Test[i] <- Test$statistic
  SortedValues$ProteinName[i]<- colnames(positive[i])
  print(c(Test$statistic,colnames(positive[i])))
}


#3.2
SortedValues[order(SortedValues$T_Test,decreasing = TRUE),]
#print(SortedValues)

#3.3
threshold1=0.001
for (k in 1:26) {
  if (SortedValues$T_Test[k]<threshold1) {
    print(c(as.character(SortedValues$ProteinName[k]),SortedValues$T_Test[k]))
  }  
}


#3.4
#there is similarity between protein(NP_000415,NP_000413,NP_000517,NP_005219,NP_058519,NP_058518
#NP_001116539,NP_061155,NP_001035932,....)
#it had tested the null hypothesis: 
#HO : X = Y
#(X − Y ) = 0

