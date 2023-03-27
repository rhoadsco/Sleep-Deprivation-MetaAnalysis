## Author: Cosette Rhoads, code originally from Dr. Hagenauer (hagenaue on Github)
## Date: July 28, 2022
## Last Updated: March 27, 2023
## UPDATE THIS WITH sessionInfo()
## Title: Meta-Analysis Code for Sleep Deprivation Project


### Section 1: Libraries
library(plyr)
library(metafor)
library(reshape)
library(multtest)
library(dplyr)

### Section 2: Functions
### None of these functions have been changed to be specific to this meta-analysis.
### All of this code was created by Dr. Megan Hagenauer and annotated versions can be found on her github:
### https://github.com/hagenaue/BrainDataAlchemy
ReadingInGemmaDE<-function(ResultSetFileNames){
    TempAnalysisResults<-read.delim("analysis.results.txt", sep="\t", stringsAsFactors = FALSE, comment.char = "#")
    TempResultsToJoin<-list(TempAnalysisResults)
    for(i in c(1:length(ResultSetFileNames))){
        TempResultsToJoin[[i]]<-read.delim(ResultSetFileNames[i], sep="\t", stringsAsFactors = FALSE, comment.char = "#")
    }
    TempResultsJoined<<-join_all(TempResultsToJoin, by="Element_Name")
    write.csv(TempResultsJoined, "TempResultsJoined.csv")
    rm(TempAnalysisResults, TempResultsToJoin)
    print("Outputted object: TempResultsJoined")
}

FilteringDEResults_GoodAnnotation<-function(TempResultsJoined){
    print("# of rows with missing NCBI annotation:")
    print(sum(TempResultsJoined$NCBI_ID==""|TempResultsJoined$NCBI_ID=="null"))
    print("# of rows with missing Gene Symbol annotation:")
    print(sum(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null"))
    print("# of rows mapped to multiple NCBI_IDs:")
    print(length(grep('\\|', TempResultsJoined$NCBI_ID)))
    print("# of rows mapped to multiple Gene Symbols:")
    print(length(grep('\\|', TempResultsJoined$Gene_Symbol)))
    TempResultsJoined_NoNA<-TempResultsJoined[(TempResultsJoined$Gene_Symbol==""|TempResultsJoined$Gene_Symbol=="null")==FALSE,]
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$Gene_Symbol)),]
    print("# of rows with good annotation")
    print(nrow(TempResultsJoined_NoNA_NoMultimapped))
    write.csv(TempResultsJoined_NoNA, "TempResultsJoined_NoNA.csv")
    write.csv(TempResultsJoined_NoNA_NoMultimapped, "TempResultsJoined_NoNA_NoMultimapped.csv")
    rm(TempResultsJoined_NoNA, TempResultsJoined_NoNA_NoMultimapped)
    print("Outputted object: TempResultsJoined_NoNA_NoMultimapped")
}

CollapsingDEResults_OneResultPerGene<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest, NamesOfFoldChangeColumns, NamesOfTstatColumns){
    print("Double check that the vectors containing the two fold change and tstat column names contain the same order as the group comparisons of interest - otherwise this function won't work properly!  If the order matches, proceed:")
    print("# of rows with unique NCBI IDs:")
    print(length(unique(TempResultsJoined_NoNA_NoMultimapped$NCBI_ID)))
    print("# of rows with unique Gene Symbols:")
    print(length(unique(TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol)))
    TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average<-list()
    for(i in c(1:length(NamesOfFoldChangeColumns))){
        TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average[[i]]<-tapply(NamesOfFoldChangeColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    }
    TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average)
    print("Dimensions of Fold Change matrix, averaged by gene symbol:")
    print(dim(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol))
    colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol)<-ComparisonsOfInterest
    write.csv(TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol.csv")
    TempResultsJoined_NoNA_NoMultimapped_Tstat_Average<-list()
    for(i in c(1:length(NamesOfFoldChangeColumns))){
        TempResultsJoined_NoNA_NoMultimapped_Tstat_Average[[i]]<-tapply(NamesOfTstatColumns[i][[1]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    }
    TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
    colnames(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol)<-ComparisonsOfInterest
    write.csv(TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol.csv")
    TempResultsJoined_NoNA_NoMultimapped_SE<-list()
    for(i in c(1:length(NamesOfFoldChangeColumns))){
        TempResultsJoined_NoNA_NoMultimapped_SE[[i]]<-NamesOfFoldChangeColumns[i][[1]]/NamesOfTstatColumns[i][[1]]
    }
    TempResultsJoined_NoNA_NoMultimapped_SE_Average<-list()
    for(i in c(1:length(NamesOfFoldChangeColumns))){
        TempResultsJoined_NoNA_NoMultimapped_SE_Average[[i]]<-tapply(TempResultsJoined_NoNA_NoMultimapped_SE[[i]], TempResultsJoined_NoNA_NoMultimapped$Gene_Symbol, mean)
    }
    TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol<-do.call(cbind, TempResultsJoined_NoNA_NoMultimapped_SE_Average)
    colnames(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)<-ComparisonsOfInterest
    write.csv(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, "TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol.csv")
    TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol)^2
    write.csv(TempResultsJoined_NoNA_NoMultimapped_SV, "TempResultsJoined_NoNA_NoMultimapped_SV.csv")
    TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, SE=TempResultsJoined_NoNA_NoMultimapped_SE_AveragedByGeneSymbol, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
    assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
    print(paste("Output: Named DEResults", GSE_ID, sep="_"))
    rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, TempResultsJoined_NoNA_NoMultimapped_FoldChange_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_FoldChange_Average, TempResultsJoined_NoNA_NoMultimapped_Tstat_AveragedByGeneSymbol, TempResultsJoined_NoNA_NoMultimapped_Tstat_Average)
}

AligningDatasets<-function(ListOfDEResults){
    MetaAnalysis_FoldChange_Dfs<-list()
    for(i in c(1:length(ListOfDEResults))){
        MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[1]]),ListOfDEResults[[i]][[1]], stringsAsFactors=FALSE)
    }
    print("MetaAnalysis_FoldChange_Dfs:")
    print(str(MetaAnalysis_FoldChange_Dfs))
    MetaAnalysis_FoldChanges<<-join_all(MetaAnalysis_FoldChange_Dfs, by="x", type="full")
    print("MetaAnalysis_FoldChanges:")
    print(str(MetaAnalysis_FoldChanges))
    MetaAnalysis_SV_Dfs<-list()
    for(i in c(1:length(ListOfDEResults))){
        MetaAnalysis_SV_Dfs[[i]]<-data.frame(x=row.names(ListOfDEResults[[i]][[4]]),ListOfDEResults[[i]][[4]], stringsAsFactors=FALSE)
    }
    print("MetaAnalysis_SV_Dfs:")
    print(str(MetaAnalysis_SV_Dfs))
    MetaAnalysis_SV<<-join_all(MetaAnalysis_SV_Dfs, by="x", type="full")
    print("MetaAnalysis_SV:")
    print(str(MetaAnalysis_SV))
    rm(MetaAnalysis_SV_Dfs, MetaAnalysis_FoldChange_Dfs)
}

RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
    MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
    print("Table of # of NAs per Row (Gene):")
    print(table(MetaAnalysis_FoldChanges_NAsPerRow))
    MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
    MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
    print("MetaAnalysis_FoldChanges_ForMeta:")
    print(str(MetaAnalysis_FoldChanges_ForMeta))
    metaOutput<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 6)
    for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
        effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
        var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
        skip_to_next <- FALSE
        tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next){}else{
            TempMeta<-rma(effect, var)
            metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
            metaOutput[i, 2]<-TempMeta$se #gives standard error
            metaOutput[i, 3]<-TempMeta$pval #gives pval
            metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
            metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
            metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
            rm(TempMeta)
        }
        rm(effect, var)
    }
    colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
    row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta$x
    metaOutput<<-metaOutput
    return(metaOutput)
    print("metaOutput:")
    print(str(metaOutput))
    print("Top of metaOutput:")
    print(head(metaOutput))
    print("Bottom of metaOutput")
    print(tail(metaOutput))
}
RunMetaAnalysisWTC<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
    MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges, 1, function(y) sum(is.na(y)))
    print("Table of # of NAs per Row (Gene):")
    print(table(MetaAnalysis_FoldChanges_NAsPerRow))
    MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
    MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
    print("MetaAnalysis_FoldChanges_ForMeta:")
    print(str(MetaAnalysis_FoldChanges_ForMeta))
    metaOutputWTC<-matrix(NA, length(MetaAnalysis_FoldChanges_ForMeta$x), 16)
    for(i in c(1:length(MetaAnalysis_FoldChanges_ForMeta$x))){
        effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-1])
        var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-1])
        skip_to_next <- FALSE
        tryCatch(TempMetaWTC<-rma(yi=effect~Predictor1+Predictor2, vi=var), error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next){}else{
            TempMetaWTC<-rma(yi=effect~Predictor1+Predictor2, vi=var)
            metaOutputWTC[i, c(1:3)]<-TempMetaWTC$b #gives estimate Log2FC
            metaOutputWTC[i, c(4:6)]<-TempMetaWTC$se #gives standard error
            metaOutputWTC[i, c(7:9)]<-TempMetaWTC$pval #gives pval
            metaOutputWTC[i, c(10:12)]<-TempMetaWTC$ci.lb #gives confidence interval lower bound
            metaOutputWTC[i, c(13:15)]<-TempMetaWTC$ci.ub #gives confidence interval upper bound
            metaOutputWTC[i, 16]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
            rm(TempMetaWTC)
        }
        rm(effect, var)
    }
    colnames(metaOutputWTC)<-c("Main_Log2FC_estimate", "Predictor1_Log2FC_estimate", "Predictor2_Log2FC_estimate", "Main_SE", "Predictor1_SE", "Predictor2_SE",  "Main_pval", "Predictor1_pval","Predictor2_pval", "Main_CI_lb","Predictor1_CI_lb","Predictor2_CI_lb", "Main_CI_ub", "Predictor1_CI_ub","Predictor2_CI_ub","Number_Of_Comparisons")
    row.names(metaOutputWTC)<-MetaAnalysis_FoldChanges_ForMeta$x
    metaOutputWTC<<-metaOutputWTC
    return(metaOutputWTC)
    print("metaOutputWTC:")
    print(str(metaOutputWTC))
    print("Top of metaOutputWTC:")
    print(head(metaOutputWTC))
    print("Bottom of metaOutputWTC")
    print(tail(metaOutputWTC))
}
FalseDiscoveryCorrection<-function(metaOutput){
    tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
    metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
    metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
    colnames(metaOutputFDR)[7]<-"FDR"
    metaOutputFDR<<-metaOutputFDR
    print("metaOutputFDR:")
    print(str(metaOutputFDR))
    write.csv(metaOutputFDR, "metaOutputFDR.csv")
    #a version of the output in order by p-value:
    metaOutputFDR_OrderbyPval<<-metaOutputFDR[order(metaOutputFDR[,3]),]
    write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval_wHDRFData.csv")
    print("Do we have any genes that are statistically significant following false discovery rate correction?")
    print(sum(metaOutputFDR[,7]<0.10, na.rm=TRUE))
    print("What are the top results?")
    print(head(metaOutputFDR[order(metaOutputFDR[,3]),]))
    rm(tempPvalAdjMeta, metaPvalAdj)
}
FalseDiscoveryCorrectionWTC<-function(metaOutputWTC){
    tempPvalAdjMeta<-mt.rawp2adjp(metaOutputWTC[,7], proc=c("BH"))
    metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
    metaOutputWTCFDR<-cbind(metaOutputWTC, metaPvalAdj[,2])
    colnames(metaOutputWTCFDR)[17]<-"Main_FDR"
    tempPvalAdjMeta<-mt.rawp2adjp(metaOutputWTC[,8], proc=c("BH"))
    metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
    metaOutputWTCFDR<-cbind(metaOutputWTCFDR, metaPvalAdj[,2])
    colnames(metaOutputWTCFDR)[18]<-"Predictor1_FDR"
    tempPvalAdjMeta<-mt.rawp2adjp(metaOutputWTC[,9], proc=c("BH"))
    metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
    metaOutputWTCFDR<-cbind(metaOutputWTCFDR, metaPvalAdj[,2])
    colnames(metaOutputWTCFDR)[19]<-"Predictor2_FDR"
    metaOutputWTCFDR<<-metaOutputWTCFDR
    print("metaOutputWTCFDR:")
    print(str(metaOutputWTCFDR))
    write.csv(metaOutputWTCFDR, "metaOutputWTCFDR.csv")
    metaOutputWTCFDR_OrderbyPval<<-metaOutputWTCFDR[order(metaOutputWTCFDR[,7]),]
    write.csv(metaOutputWTCFDR_OrderbyPval, "metaOutputWTCFDR_orderedByPval_wHDRFData.csv")
    print("Do we have any genes that are statistically significant following false discovery rate correction?")
    print(sum(metaOutputWTCFDR[,17]<0.10, na.rm=TRUE))
    print("What are the top results?")
    print(head(metaOutputWTCFDR[order(metaOutputWTCFDR[,7]),]))
    rm(tempPvalAdjMeta, metaPvalAdj)
}

MakeForestPlots<-function(GeneSymbol){
    pdf(paste("ForestPlot_", GeneSymbol, ".pdf", sep=""), height=5, width=8)
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$x==GeneSymbol,-1])
    forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-1],  xlim=c(-3, 3))
    mtext(paste(GeneSymbol), line=-1.5, cex=2)
    dev.off()
}



### Section 3: Cosette Rhoads' Meta-Analysis on Sleep Deprivation
### These datasets were found on Gemma and are identified by their Gemma ID number.
### To reproduce this process, look at the  outline below.

# OUTLINE FOR REPRODUCING SECTION 3
# Create a comment noting the Gemma ID and Authors.
# run setwd() to the path to the folder that contains the downloaded differential expression analysis for this study
# run list.files() will produce the names of the files in that folder. 
# there will be a file named "analysis.results.txt" as well as any number of resultset files, e.g. "resultset_1" and "resultset_2".
# use the function ReadingInGemmaDE and set the inputs to the resultset files: 
# ReadingInGemmaDE(ResultSetFileNames = c("resultset_1", "resultset_2")
# run FilteringDEResults_GoodAnnotation(TempResultsJoined)  
# before running the next function, identify the treatment groups you want to include in your meta-analysis and find the column in the object "TempResultsJoined"
# run the function: CollapsingDEResults_OneResultPerGene
# (GSE_ID="input the Gemma ID here", TempResultsJoined_NoNA_NoMultimapped, 
# ComparisonsOfInterest=c("treatment group 1 title", "treatment group 2 title, etc."), 
# NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$treatment group 1 fold change column name, TempResultsJoined_NoNA_NoMultimapped$treatment group 2 fold change column name), 
# NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$treatment group 1 t statcolumn name, TempResultsJoined_NoNA_NoMultimapped$treatment group 2 t statcolumn name)


# Meta-Analysis on Sleep Deprivation

#GSE6514, Mackiewicz, Miroslaw et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE6514/11418_GSE6514_diffExpAnalysis_92880")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID477377.data.txt"                                             
# [3] "resultset_ID477378.data.txt" 
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID477377.data.txt", "resultset_ID477378.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE6514", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("6 h SD vs. GSE6514 ctrl", "9 h SD vs. GSE6514 ctrl", "12 h SD vs. GSE6514 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_6.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_9.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_12.h), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_6.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_9.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_12.h))

#GSE33491, Hinard, Valérie et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE33491/13786_GSE33491_diffExpAnalysis_99952")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID480876.data.txt"    
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID480876.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE33491", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("6 h SD vs. GSE33491 ctrl", "6 h SD, 18 h RS vs. GSE33491 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_sleep.deprivation, TempResultsJoined_NoNA_NoMultimapped$FoldChange_sleep.deprivation.for.6.h.and.recovery.for.18.h), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_sleep.deprivation, TempResultsJoined_NoNA_NoMultimapped$Tstat_sleep.deprivation.for.6.h.and.recovery.for.18.h))

#GSE78215, Gerstner, Jason R et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE78215/12969_GSE78215_diffExpAnalysis_99294")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID480552.data.txt"                                             
# [3] "resultset_ID480553.data.txt"                                             
# [4] "resultset_ID480554.data.txt" 
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID480552.data.txt", "resultset_ID480553.data.txt", "resultset_ID480554.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE78215", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("6 h SD, 1 h RS vs. GSE78215 ctrl", "6 h SD, 2 h RS vs. GSE78215 ctrl", "5 h SD, 3 h RS vs. GSE78215 ctrl", "5 h SD, 6 h RS vs. GSE78215 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_1.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_2.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_3.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_6.h), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_1.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_2.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_3.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_6.h))

#GSE93041, Orozco-Solis, Ricardo et al. 
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE93041/13143_GSE93041_diffExpAnalysis_100030")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID480918.data.txt" 
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID480918.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE93041", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("12 h SD vs. GSE93041 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_sleep.deprivation), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_sleep.deprivation))

#GSE113754, Ingiosi, Ashley M et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE113754/15159_GSE113754_diffExpAnalysis_164400")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID494159.data.txt"                                             
# [3] "resultset_ID494160.data.txt" 
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID494159.data.txt", "resultset_ID494160.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE113754", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("5 h SD vs. GSE113754 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_Sleep.deprivation.5h.), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_Sleep.deprivation.5h.))

#GSE128770, Guo, Xiaofeng et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE128770/15623_GSE128770_diffExpAnalysis_131464")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID490841.data.txt"                                             
# [3] "resultset_ID490842.data.txt"                                             
# [4] "resultset_ID490843.data.txt"   
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID490841.data.txt", "resultset_ID490842.data.txt", "resultset_ID490843.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE128770", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("3 h SD vs. GSE128770 ctrl", "6 h SD vs. GSE128770 ctrl", "9 h SD vs. GSE128770 ctrl", "12 h SD vs. GSE128770 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_3.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_6.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_9.h, TempResultsJoined_NoNA_NoMultimapped$FoldChange_12.h), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_3.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_6.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_9.h, TempResultsJoined_NoNA_NoMultimapped$Tstat_12.h))

#GSE132076, Muheim, Christine M et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE132076/17237_GSE132076_diffExpAnalysis_176993")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID499390.data.txt"                                             
# [3] "resultset_ID499391.data.txt"  
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID499390.data.txt", "resultset_ID499391.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE132076", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("4 h SD vs. GSE132076 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_Sleep.deprivation), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_Sleep.deprivation))

#GSE144957, Bjorness, Theresa E et al.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE144957/17613_GSE144957_diffExpAnalysis_168449")
list.files()
# [1] "analysis.results.txt"                                                    
# [2] "resultset_ID495477.data.txt"                                             
# [3] "resultset_ID495478.data.txt"
ReadingInGemmaDE(ResultSetFileNames = c("resultset_ID495477.data.txt", "resultset_ID495478.data.txt"))
FilteringDEResults_GoodAnnotation(TempResultsJoined)
CollapsingDEResults_OneResultPerGene(GSE_ID="GSE144957", TempResultsJoined_NoNA_NoMultimapped, ComparisonsOfInterest=c("6 h SD vs. GSE144957 ctrl", "4 h SD, 2 h RS vs. GSE144957 ctrl"), NamesOfFoldChangeColumns=list(TempResultsJoined_NoNA_NoMultimapped$FoldChange_sleep.deprived.6.hrs., TempResultsJoined_NoNA_NoMultimapped$FoldChange_recovery.sleep.2.hrs.), NamesOfTstatColumns = list(TempResultsJoined_NoNA_NoMultimapped$Tstat_sleep.deprived.6.hrs., TempResultsJoined_NoNA_NoMultimapped$Tstat_recovery.sleep.2.hrs.))


### Section 4: Running the Meta-Analysis
### This section is where, using the functions in section 2, we run the meta-analysis on the data from section 3. 

# setwd to a file where you will store results
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results")
# define ListOfDEResults as the DEResults_GemmaID files created in Section 3.
ListOfDEResults <- list(DEResults_GSE113754, DEResults_GSE128770, DEResults_GSE132076, DEResults_GSE144957, DEResults_GSE144957, DEResults_GSE33491, DEResults_GSE6514, DEResults_GSE78215, DEResults_GSE93041)
AligningDatasets(ListOfDEResults)
# From there, specify some exclusion constraints for the meta-analysis function. NumberOfComparisons refers to how many comparisons are contained in the DEResults. CutOffForNAs refers to the maximum number of comparisons that can have "NA" data for a gene for that gene to be included in the analysis.
NumberOfComparisons=18
CutOffForNAs=5
# First, we run a basic meta analysis with no covariates and an FDR correction. This basic meta analysis function is our confirmatory analysis. This function takes a while to run. Set the wd to a specific file location.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results/BasicMetaAnalysis")
metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
FalseDiscoveryCorrection(metaOutput)
# Second, we can define covariates of interest. We began with duration of sleep deprivation and duration of recovery sleep, then tried duration of sleep deprivation with the existence recovery sleep as a categorical variable. The covariates are created by giving each comparison (in this case, all 18 comparisons) one value. This must be in the same order as the MetaAnalysis_FoldChanges object's columns, as so:
# colnames(MetaAnalysis_FoldChanges)
# [1] "x"                                  
# [2]"X5.h.SD.vs..GSE113754.ctrl"         
# [3]"X3.h.SD.vs..GSE128770.ctrl"        
# [4] "X6.h.SD.vs..GSE128770.ctrl"         
# [5]"X9.h.SD.vs..GSE128770.ctrl"         
# [6]"X12.h.SD.vs..GSE128770.ctrl"       
# [7] "X4.h.SD.vs..GSE132076.ctrl"         
# [8]"X6.h.SD.vs..GSE144957.ctrl"         
# [9]"X4.h.SD..2.h.RS.vs..GSE144957.ctrl"
# [10] "X6.h.SD.vs..GSE33491.ctrl"          
# [11]"X6.h.SD..18.h.RS.vs..GSE33491.ctrl" 
# [12]"X6.h.SD.vs..GSE6514.ctrl"          
# [13] "X9.h.SD.vs..GSE6514.ctrl"           
# [14]"X12.h.SD.vs..GSE6514.ctrl"          
# [15]"X6.h.SD..1.h.RS.vs..GSE78215.ctrl" 
# [16] "X6.h.SD..2.h.RS.vs..GSE78215.ctrl"  
# [17]"X5.h.SD..3.h.RS.vs..GSE78215.ctrl"  
# [18]"X5.h.SD..6.h.RS.vs..GSE78215.ctrl" 
# [19] "X12.h.SD.vs..GSE93041.ctrl" 
SD_uncentered<- c(5, 3, 6, 9, 12, 4, 6, 4, 6, 6, 6, 9, 12, 6, 6, 5, 5, 12)
DurationSD<- SD_uncentered - mean(SD_uncentered)
RS_uncentered<- c(0, 0, 0, 0, 0, 0, 0, 2, 0, 18, 0, 0, 0, 1, 2, 3, 6, 0)
DurationRS<- RS_uncentered - mean(RS_uncentered)
#ExistenceRS has 0 = no RS, 1 = yes RS
ExistenceRS<- as.factor(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0))
# Third, we run a meta analysis with covariates and its FDR correction as an exploratory analysis. This function also takes a while to run. This includes two meta-analyses run with covariates as specified above. Again, set the wd to a specific file location.
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results/MetaAnalysis_DurationSD vs. DurationRS")
Predictor1<-DurationSD
Predictor2<-DurationRS
metaOutputWTC<-RunMetaAnalysisWTC(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
FalseDiscoveryCorrectionWTC(metaOutputWTC)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results/MetaAnalysis_DurationSD vs. ExistenceRS")
Predictor1<-DurationSD
Predictor2<-ExistenceRS
metaOutputWTC<-RunMetaAnalysisWTC(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
FalseDiscoveryCorrectionWTC(metaOutputWTC)


### Section 5: Visualizing the Data and Results
cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs")
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-1]), use="pairwise.complete.obs"))
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results/ForestPlots")
MakeForestPlots("Slc12a6")
MakeForestPlots("Thg1l")
MakeForestPlots("B3gnt3")
MakeForestPlots("Vps26c")
MakeForestPlots("Fancf")
MakeForestPlots("Hspa12b")
MakeForestPlots("Nr3c1")
MakeForestPlots("Gmeb1")
MakeForestPlots("Ddx51")
MakeForestPlots("Nr2e1")
MakeForestPlots("Gm14285")
MakeForestPlots("Exosc4")
MakeForestPlots("Tmod3")
MakeForestPlots("Rimoc1")
MakeForestPlots("Cdc42ep3")
MakeForestPlots("4931429L15Rik")
MakeForestPlots("Cd7")
MakeForestPlots("Lrrc75a")
MakeForestPlots("Wdr61")

### Section 6: What Matters? Exploring the code to see how our analysis choices affect our results.

##finding the gene symbols that are repeated.
#GSE6514
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE6514/11418_GSE6514_diffExpAnalysis_92880")
GSE6514_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE6514/11418_GSE6514_diffExpAnalysis_92880/TempResultsJoined_NoNA.csv")
GSE6514_raw_sorted <- GSE6514_raw[order(GSE6514_raw$Gene_Symbol),]
write.csv(GSE6514_raw_sorted, "GSE6514_raw_sorted.csv")
GSE6514_unique<-unique(GSE6514_raw_sorted$Gene_Symbol[duplicated(GSE6514_raw_sorted$Gene_Symbol) | duplicated(GSE6514_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE33491
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE33491/13786_GSE33491_diffExpAnalysis_99952")
GSE33491_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE33491/13786_GSE33491_diffExpAnalysis_99952/TempResultsJoined_NoNA.csv")
GSE33491_raw_sorted <- GSE33491_raw[order(GSE33491_raw$Gene_Symbol),]
write.csv(GSE33491_raw_sorted, "GSE33491_raw_sorted.csv")
GSE33491_unique<-unique(GSE33491_raw_sorted$Gene_Symbol[duplicated(GSE33491_raw_sorted$Gene_Symbol) | duplicated(GSE33491_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE78215
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE78215/12969_GSE78215_diffExpAnalysis_99294")
GSE78215_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE78215/12969_GSE78215_diffExpAnalysis_99294/TempResultsJoined_NoNA.csv")
GSE78215_raw_sorted <- GSE78215_raw[order(GSE78215_raw$Gene_Symbol),]
write.csv(GSE78215_raw_sorted, "GSE78215_raw_sorted.csv")
GSE78215_unique<-unique(GSE78215_raw_sorted$Gene_Symbol[duplicated(GSE78215_raw_sorted$Gene_Symbol) | duplicated(GSE78215_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE93041
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE93041/13143_GSE93041_diffExpAnalysis_100030")
GSE93041_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE93041/13143_GSE93041_diffExpAnalysis_100030/TempResultsJoined_NoNA.csv")
GSE93041_raw_sorted <- GSE93041_raw[order(GSE93041_raw$Gene_Symbol),]
write.csv(GSE93041_raw_sorted, "GSE93041_raw_sorted.csv")
GSE93041_unique<-unique(GSE93041_raw_sorted$Gene_Symbol[duplicated(GSE93041_raw_sorted$Gene_Symbol) | duplicated(GSE93041_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE113754
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE113754/15159_GSE113754_diffExpAnalysis_164400")
GSE113754_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE113754/15159_GSE113754_diffExpAnalysis_164400/TempResultsJoined_NoNA.csv")
GSE113754_raw_sorted <- GSE113754_raw[order(GSE113754_raw$Gene_Symbol),]
write.csv(GSE113754_raw_sorted, "GSE113754_raw_sorted.csv")
GSE113754_unique<-unique(GSE113754_raw_sorted$Gene_Symbol[duplicated(GSE113754_raw_sorted$Gene_Symbol) | duplicated(GSE113754_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE128770
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE128770/15623_GSE128770_diffExpAnalysis_131464")
GSE128770_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE128770/15623_GSE128770_diffExpAnalysis_131464/TempResultsJoined_NoNA.csv")
GSE128770_raw_sorted <- GSE128770_raw[order(GSE128770_raw$Gene_Symbol),]
write.csv(GSE128770_raw_sorted, "GSE128770_raw_sorted.csv")
GSE128770_unique<-unique(GSE128770_raw_sorted$Gene_Symbol[duplicated(GSE128770_raw_sorted$Gene_Symbol) | duplicated(GSE128770_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE132076
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE132076/17237_GSE132076_diffExpAnalysis_176993")
GSE132076_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE132076/17237_GSE132076_diffExpAnalysis_176993/TempResultsJoined_NoNA.csv")
GSE132076_raw_sorted <- GSE132076_raw[order(GSE132076_raw$Gene_Symbol),]
write.csv(GSE132076_raw_sorted, "GSE132076_raw_sorted.csv")
GSE132076_unique<-unique(GSE132076_raw_sorted$Gene_Symbol[duplicated(GSE132076_raw_sorted$Gene_Symbol) | duplicated(GSE132076_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#GSE144957
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE144957/17613_GSE144957_diffExpAnalysis_168449")
GSE144957_raw<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/GSE144957/17613_GSE144957_diffExpAnalysis_168449/TempResultsJoined_NoNA.csv")
GSE144957_raw_sorted <- GSE144957_raw[order(GSE144957_raw$Gene_Symbol),]
write.csv(GSE144957_raw_sorted, "GSE144957_raw_sorted.csv")
GSE144957_unique<-unique(GSE144957_raw_sorted$Gene_Symbol[duplicated(GSE144957_raw_sorted$Gene_Symbol) | duplicated(GSE144957_raw_sorted$Gene_Symbol, fromLast=TRUE)])
#all_repeated_genesymbols is a list of all the gene symbols that are duplicated in any of the 7 datasets.
# note: it could be good to functionalize all of the code above. Ask Dr. H about that.
all_repeated_genesymbols<-unique(c(GSE132076_unique,GSE128770_unique,GSE113754_unique,GSE93041_unique,GSE78215_unique,GSE33491_unique,GSE6514_unique))
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results")
metaOutputFDR_dataright<-read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/Personal/Professional/UMich Summer Internship/R Coding for Meta Analysis/Results/BasicMetaAnalysis/metaOutputFDR.csv")
metaOutputFDR_sig<-filter(metaOutputFDR_dataright, metaOutputFDR_dataright$FDR < 0.1)
all_significant_genesymbols<-(metaOutputFDR_sig$X)
#print(all_repeated_genesymbols)
#print(all_significant_genesymbols)
y<-logical()
for (x in 1:length(all_repeated_genesymbols)) {
    y = c(y, all_repeated_genesymbols[[x]] %in% all_significant_genesymbols)
}
proportion_sig <- (sum(y, na.rm=TRUE))/(length(all_significant_genesymbols))
proportion_all <- (length(all_repeated_genesymbols))/(length(unique(metaOutputFDR_dataright$X)))
H <- c(proportion_sig,proportion_all)
M <- c("genes from significant results", "all genes")
barplot(H,ylab="Proportion of Genes that are Multimapped",names.arg=M, ylim= c(0,1))
print ("proportion of statistically significant genes that had more than one probe per gene symbol: ") 
print (proportion_sig)
print ("proportion of all genes that had more than one probe per gene symbol: ")
print (proportion_all)


