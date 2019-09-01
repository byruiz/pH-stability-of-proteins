# pH-stability-of-proteins
analysis of multiplexed yeast proteome dataset
---
title: "fresh analysis CPP_9Aug2019"
---

#Setup
```{r}

library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(GGally)
library(ggrepel)
library(viridis)
library(gplots)
library(RColorBrewer)
library(scales)
library(gifski)
library(png)
library(transformr)
library(plotly)
library(splines)
library(zoo)
```
#Read in data
```{r}
#set working directory
setwd("~/Documents/Villen Lab//CPP/pH selection 50C/")

#read the csv files into a dataframes
repA <- read.csv("hot_pH_repA.csv", header=TRUE)
repB <- read.csv("hot_pH_repB.csv", header=TRUE)

#this csv was downloaded from the Villen Lab team drive
GOtermsSlim <- read.csv("201710SGD_go_slim_mapping.csv", header=TRUE)
GOterms <- read.csv("Yeast_GO_terms_2017.csv", header=TRUE)

```
#Data wrangling (bothReps, tidy)
```{r}

#combine data from raw CSVs into one DF
bothReps <- rbind(repA, repB)

#tidy and filter
tidy <- bothReps %>%
  
  #tidy
  gather(
    key = TMT,
    value = Intensity, 
    -Sample.Name, -Retention.time, -Sequence, -Reference, -charge, -QScore) %>%
  
  #make new column specifying TMT labels
  mutate(
    #add column to specify reverse phase basic Fraction
    Fraction = str_sub(Sample.Name, start = 25, end = -11),
    #add columnn to specify biological Replicate A or B
    Replicate = str_sub(Sample.Name, start = 20, end = 20),
    #add column with cleaned up sequence, remove extra characters
    trimmedSeq = str_sub(Sequence, start = 3, end = -3),
    #add column specifying TMT channel
    TMT = str_replace_all(TMT, pattern = "Intensity.", replacement = ""))%>%
    
  #remove decoy peptides and 10th channel (only 9 channels in use in this experiment)
  filter(!grepl("DECOY", Reference),
         !grepl("131", TMT),
         !grepl("#", trimmedSeq),
         !grepl("[*]", trimmedSeq)) %>%
  
  #remove sample name column
  select(-Sample.Name)


###add pH to df

#make DF to left join, to assign pH condition to each TMT channel
TMT.v.pH <- data.frame("TMT" = c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C"), 
                       "pH" = c(5.616, 6.02, 6.462, 6.908, 7.28, 7.567, 7.845, 8.18, 8.635))

#replace TMT channel with respective pH condition
tidy <- left_join(tidy, TMT.v.pH, by = "TMT") %>%
  #just rearranging columns and removing extra column
  select(Replicate, Fraction, pH, Reference, trimmedSeq, Intensity, Retention.time, Sequence, charge, QScore, -TMT)

```
#Aggregating intensities
```{r}

#add measurements from each Fraction into weighted average for each peptide
tidyAggregated <- tidy %>%

  #combine measurements into weighted average per replicate (whole proteome)
  group_by(Replicate, pH) %>%
  mutate(proteomeSumInt = sum(Intensity)) %>%
  ungroup() %>%
  
  #combine measurements into weighted average per protein
  group_by(Replicate, pH, Reference) %>%
  mutate(protSumInt = sum(Intensity)) %>%
  ungroup() %>%
  
  #combine measurements into weighted average per peptide
  group_by(Replicate, pH, Reference, trimmedSeq) %>%
  mutate(pepSumInt = sum(Intensity)) %>%
  ungroup() %>%
  
  #remove peptides with duplicate measurements due to multiple charge states
  select(Replicate, Reference, trimmedSeq, pH, proteomeSumInt, protSumInt, pepSumInt) %>%
  distinct()

```
#Adding GO term annotation (tidyWithComps)
```{r}
#GOtermsSlim does not include reference numbers (referred to as Gene.secondaryIdentifier in GOterms), the next few lines add reference numbers from GOterms to GOtermsSlim
referenceIDs <- GOterms %>%
  select(Gene.primaryIdentifier, Gene.secondaryIdentifier) %>%
  #lots of duplicates in this set because proteins have multiple annotations
  #remove duplicates with distinct function
  distinct()

#leftjoin GOtermsSlim with referenceIDs to add Reference
goTermsComplete <- inner_join(GOtermsSlim, referenceIDs, by = "Gene.primaryIdentifier") %>%
  #rename columns
  rename(Reference = Gene.secondaryIdentifier.x,
         Symbol = Gene.symbol,
         goID = Gene.primaryIdentifier,
         term = Gene.ontologyAnnotations.ontologyTerm.name,
         termType = Gene.ontologyAnnotations.ontologyTerm.namespace,
         termID = Gene.ontologyAnnotations.ontologyTerm.identifier)%>%
  select(Reference, termType, term, Symbol) %>%
  mutate(termType = str_replace_all(termType, "C", "compartment"),
         termType = str_replace_all(termType, "F", "function"),
         termType = str_replace_all(termType, "P", "process"),
         term = as.character(term))
  
#these lines can be used to filter out proteins with more than one annotation
  #filter(termType == "C") %>%
  #group_by(Reference) %>%
  #mutate(NumTerms=n()) %>%
  #filter(NumTerms == 1) %>%
  #select(Reference, term)

goTermsComps <- goTermsComplete %>% filter(termType == "compartment")
tidyWithComps <- left_join(tidyAggregated, goTermsComps, by = "Reference")

goTermsFunctions <- goTermsComplete %>% filter(termType == "function")
tidyWithFunctions <- left_join(tidyAggregated, goTermsFunctions, by = "Reference")

goTermsProcess <- goTermsComplete %>% filter(termType == "process")
tidyWithProcess <- left_join(tidyAggregated, goTermsProcess, by = "Reference")
```
#spline
```{r}
#peptide level spline analysis

cs.fit <- function(ph, vals){
  spline_ys <- ns(ph, df = 3)  #curves
  model <- lm(vals~spline_ys) #weighted averages for each curve
  pred<- predict(model) #takes weighted average and uses them on the splines
  return(pred)
}

tidySpline <- tidyAggregated %>%
  
  #condense data down to one measurement per peptide per replicate
  distinct(Replicate, pH, Reference, trimmedSeq, pepSumInt) %>%
  
  #filter for peptides measured in both replicates
  group_by(Reference, trimmedSeq, pH)%>%
  mutate(numRepPeps = n())%>%
  filter(numRepPeps>1)%>%
  
  #change 0 values to lowest pepSumInt value after filtering
  mutate(pepSumInt = if_else(pepSumInt == 0, 	762.583, pepSumInt),
         
         #log pepSumInt values
         pepSumInt = log2(pepSumInt)) %>%
  
  #calculate spline predictions with cs.fit function (written above)
  ungroup()%>%
  group_by(Replicate, Reference, trimmedSeq) %>%
  mutate(splinePred = cs.fit(pH, pepSumInt)) %>%

  #this column is for ggplot grouping purposes
  ungroup() %>%
  mutate(trimSeq_Rep = paste(trimmedSeq, Replicate, sep="_"))


tidySplineError <- tidySpline %>%
  #normalizing to max spline
  group_by(Replicate, Reference, trimmedSeq) %>%
  mutate(normPredInt =splinePred/max(splinePred),
         normPepInt = pepSumInt/max(pepSumInt),
         #calculate error between spline predictions and actual values
         splineError = abs(normPredInt-normPepInt))

```
#plotting
```{r}
ggplot(data = tidySplineError %>% filter(Reference == "YML010W"), aes(x= pH, group = trimSeq_Rep, color = trimmedSeq))+
    #written like this because of the weird lineplot grouping situation i always run into...
  #peptide
  geom_point(aes(y = normPepInt), alpha = 0.7)+
  geom_line(aes(y = normPredInt))+
  theme_bw()+
  ylab("Normalized Intensity")

ggplot(data = tidySplineError %>% filter(Reference == "YLR109W"), aes(x= pH, group = trimSeq_Rep, color = trimmedSeq))+
    #written like this because of the weird lineplot grouping situation i always run into...
  #peptide
  geom_point(aes(y = normPepInt), alpha = 0.7)+
  geom_line(aes(y = normPredInt))+
  theme_bw()+
  ylab("Normalized Intensity")


```

#analysing AA composition
```{r}

tidyAAcomposition <- tidySplineError %>%
  mutate(numDorE = str_count(trimmedSeq, "D") + str_count(trimmedSeq, "E"),
         numRorKorH = str_count(trimmedSeq, "R") + str_count(trimmedSeq, "K")+ str_count(trimmedSeq, "H"))


ggplot(data = tidyAAcomposition%>% filter(pH == 5.616))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 6.020))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 6.462))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 6.908))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 7.280))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 7.567))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH))+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 7.845))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 8.180))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH), alpha = 0.7)+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")

ggplot(data = tidyAAcomposition%>% filter(pH == 8.635))+
  geom_density(aes(x = normPredInt, group = numRorKorH, color = numRorKorH))+
  scale_color_viridis()+
  theme_dark(base_size = 20)+
  ylab("Density")+
  xlab("Spline Normalized Intensity")


```
