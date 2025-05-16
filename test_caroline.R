###############################################################################/
#PARTIE DE CAROLINE 
###############################################################################/

#data herbicides on Myzus persicae clones
MyzHerbi<-read.table("data_Myzus_herbicides_20250228.txt",
                     header=TRUE,stringsAsFactors=TRUE,sep="\t")
#reorder factors
MyzHerbi$Dose<-factor(MyzHerbi$Dose,levels=c("NT","N/2","N","2N"))

#reorder the Active substances
MyzHerbi$Active_substance<-factor(MyzHerbi$Active_substance,
                                  levels=c("Clethodim","Cycloxydim",
                                           "Pinoxaden","Clodinafop",
                                           "Fluazifop","Quizalofop"))


temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]
temp<-drop.levels(temp,reorder=FALSE)

#temp$Dose <- trimws(toupper(data$Dose))

###############################################################################/
# CHI² TEST : Herbicide selection pressure on the susceptible M.persicae clone 
# - DOSE N 
###############################################################################/

#Filter data for clone '17-0153-0020' and doses 'NT' and 'N'
temp_filteredN<-temp%>%
  filter(Clone == "17-0153-0020", Dose %in% c("NT", "N"))

#Cleaning and standardisation of active substance names
temp$Active_substance <- trimws(toupper(temp$Active_substance))

#List of single active substances
substances <- unique(temp_filteredN$Active_substance)
print(substances)

#Loop to perform a Chi² test on each active substance
for (substance in substances) {
  
  #Filter data for the current substance
  temp_substances <- temp_filteredN %>%
    filter(Active_substance == substance) %>%
    group_by(Dose)%>%
    summarise(Total_death = sum(Total_death))
  
  #Check that there are two values (NT and N) for the test
  if (nrow(temp_substances) == 2) {
    
    #Building the contingency table
    cont_table <- matrix(temp_substances$Total_death, nrow = 2, byrow = TRUE)
    
    #Displaying results
    cat("\nSubstance :", substance, "\n")
    print(cont_table) 
    
    #Chi² test
    test_chi2 <- chisq.test(cont_table)
    print(test_chi2)
    
    
  } else {
    cat("\nSubstance :", substance, " - Données insuffisantes pour le test\n")
  }
}


###############################################################################/
# FISHER TEST : Herbicide selection pressure on the susceptible M.persicae clone 
# - DOSE N --> because the numbers are too small
###############################################################################/

#Filter data for clone '17-0153-0020' and doses 'NT' and 'N'
temp_filteredN1<-temp%>%
  filter(Clone == "17-0153-0020", Dose %in% c("NT", "N"))

#Cleaning and standardisation of active substance names
temp$Active_substance <- trimws(toupper(temp$Active_substance))

#List of single active substances
substances <- unique(temp_filteredN1$Active_substance)
print(substances)

#Loop to perform a Fisher test on each active substance
for (substance in substances) {
  
  #Filter data for the current substance
  temp_substancesN1 <- temp_filteredN1 %>%
    filter(Active_substance == substance) %>%
    group_by(Dose)%>%
    summarise(Total_death = sum(Total_death), Total=sum(Total))
  temp_substancesN1 <- temp_substancesN1 %>%
    mutate(Survivors = Total- Total_death)
  
  #Check that there are two values (NT and N) for the test
  if (nrow(temp_substancesN1) == 2) {
    
    #Building the contingency table
    cont_table <- matrix(c(temp_substancesN1$Total_death, temp_substancesN1$Survivors), nrow = 2, byrow = FALSE)
    
    #Displaying results
    cat("\nSubstance :", substance, "\n")
    print(cont_table) 
    
    #Fisher test
    test_fisher <- fisher.test(cont_table)
    print(test_fisher)
    
    
  } else {
    cat("\nSubstance :", substance, " - Données insuffisantes pour le test\n")
  }
}


###############################################################################/
# CHI² TEST : Herbicide selection pressure on the susceptible M.persicae clone 
#- DOSE 2N
###############################################################################/

#Filter data for clone '17-0153-0020' and doses 'NT' and '2N'
temp_filtered2N<-temp%>%
  filter(Clone == "17-0153-0020", Dose %in% c("NT", "2N"))

#Cleaning and standardisation of active substance names
temp$Active_substance <- trimws(toupper(temp$Active_substance)) 

#List of single active substances
substances2N <- unique(temp_filtered2N$Active_substance)
print(substances2N)


#Loop to perform a Chi² test on each active substance
for (substance in substances2N) {
  
  #Filter data for the current substance
  temp_substance2N <- temp_filtered2N %>%
    filter(Active_substance == substance) %>%
    group_by(Dose)%>%
    summarise(Total_death = sum(Total_death))
  
  #Check that there are two values (NT and 2N) for the test
  if (nrow(temp_substance2N) == 2) {
    
    #Building the contingency table
    cont_table <- matrix(temp_substance2N$Total_death, nrow = 2, byrow = TRUE)
    
    #Displaying results
    cat("\nSubstance :", substance, "\n")
    print(cont_table) 
    
    #Chi² test
    test_chi2 <- chisq.test(cont_table)
    print(test_chi2)
    
    
  } else {
    cat("\nSubstance :", substance, " - Données insuffisantes pour le test\n")
  }
}


###############################################################################/
# FISHER TEST : Herbicide selection pressure on the susceptible M.persicae clone 
# - DOSE 2N --> because the numbers are too small
###############################################################################/

#Filter data for clone '17-0153-0020' and doses 'NT' and '2N'
temp_filtered2N2<-temp%>%
  filter(Clone == "17-0153-0020", Dose %in% c("NT", "2N"))

#Cleaning and standardisation of active substance names
temp$Active_substance <- trimws(toupper(temp$Active_substance))

#List of single active substances
substances2N2 <- unique(temp_filtered2N2$Active_substance)
print(substances2N2)

#Loop to perform a Fisher test on each active substance
for (substance in substances2N2) {
  
  #Filter data for the current substance
  temp_substances2N2 <- temp_filtered2N2 %>%
    filter(Active_substance == substance) %>%
    group_by(Dose)%>%
    summarise(Total_death = sum(Total_death), Total=sum(Total))
  temp_substances2N2 <- temp_substances2N2 %>%
    mutate(Survivors = Total- Total_death)
  
  #Check that there are two values (NT and N) for the test
  if (nrow(temp_substances2N2) == 2) {
    
    #Building the contingency table
    cont_table <- matrix(c(temp_substances2N2$Total_death, temp_substances2N2$Survivors), nrow = 2, byrow = FALSE)
    
    #Displaying results
    cat("\nSubstance :", substance, "\n")
    print(cont_table) 
    
    #Fisher test
    test_fisher <- fisher.test(cont_table)
    print(test_fisher)
    
    
  } else {
    cat("\nSubstance :", substance, " - Données insuffisantes pour le test\n")
  }
}


###############################################################################/
# CHI² TEST : variability in the sensitivity of clones to AS
###############################################################################/

substances2 <- unique(temp$Active_substance)
doses2<- unique(temp$Dose)
clones2<- unique(temp$Clone)

if (length(clones2) !=2) {
  stop("STOP")
}

for (substance in substances2) {
  temp_substance2<-temp %>%
    filter(Active_substance == substance)
  
  for (dose in doses2) {
    temp_filtered2<-temp_substance2 %>%
      filter(Dose == dose) %>%
      group_by(Clone) %>%
      summarise(Total_death = sum(Total_death), .groups="drop") %>%
      arrange(Clone)
    
    cat("\nSubstance :", substance, "- Dose :", dose, "\n")
    print(temp_filtered2)
    
    if (nrow(temp_filtered2) == 2) {
      cont_table <- as.matrix(temp_filtered2$Total_death)
      cat("\nSubstance :", substance, "- Dose :", dose, "\n")
      print(cont_table)
      
      test_chi2<-chisq.test(cont_table)
      print(test_chi2)
      
    } else {
      cat("\nSubstance :", substance, "- Dose :", dose, "Données insuffisantes pour le test\n")
    }
  }
}

#########################################################################/
# FISHER TEST : variability in the sensitivity of clones to AS
#########################################################################/

# List of active substances, doses and clones
substances3 <- unique(temp$Active_substance)
doses3 <- unique(temp$Dose)
clones3 <- unique(temp$Clone)

if (length(clones3) != 2) {
  stop("STOP")
}

# Loop for each active substance
for (substance in substances3) {
  temp_substance3 <- temp %>%
    filter(Active_substance == substance)
  
  # Loop for each dose
  for (dose in doses3) {
    
    # Filter for current dose and group by clone
    temp_filtered3 <- temp_substance3 %>%
      filter(Dose == dose) %>%
      group_by(Clone) %>%
      summarise(Total_death = sum(Total_death), Live = sum(Live), .groups = "drop") %>%
      arrange(Clone)  # Ensuring a coherent order of clones
    
    # Check extracted values
    cat("\nSubstance :", substance, "- Dose :", dose, "\n")
    print(temp_filtered3)  
    
    # Check that you have 2 clones before testing
    if (nrow(temp_filtered3) == 2) {
      
      # Creating a well-ordered contingency table
      cont_table <- matrix(c(temp_filtered3$Total_death[1], temp_filtered3$Total_death[2], 
                             temp_filtered3$Live[1], temp_filtered3$Live[2]), 
                           nrow = 2, byrow = TRUE)
      
      
      # Display the table for checking
      cat("\nTable de contingence pour", substance, "- Dose :", dose, "\n")
      print(cont_table)
      
      
      #Fisher test
      test_fisher <- fisher.test(cont_table)
      print(test_fisher)
      
      
    } else {
      cat("\nSubstance :", substance, " - Données insuffisantes pour le test\n")
    }
  }
}

###############################################################################/

# Charger les bibliothèques nécessaires
library(dplyr)
library(ggplot2)
library(rstatix)
library(multcompView)

# Charger les données
data <- read.delim("data_Myzus_herbicides_20250228.txt", sep="\t")

# Filtrer pour le clone 17-0153-0020
clone_data <- data %>% filter(Clone == "17-0153-0020")

# Fonction pour analyser chaque substance active
test_toxicity <- function(substance_data) {
  cat("\nAnalyse pour la substance :", unique(substance_data$Active_substance), "\n")
  
  # Vérifier la normalité (Shapiro-Wilk) et l'homogénéité des variances (Levene)
  shapiro_result <- substance_data %>% group_by(Dose) %>% summarize(p_value = shapiro.test(Total_death)$p.value)
  print(shapiro_result)
  
  levene_result <- levene_test(substance_data, Total_death ~ Dose)
  print(levene_result)
  
  if (all(shapiro_result$p_value > 0.05) & levene_result$p > 0.05) {
    cat("ANOVA paramétrique\n")
    anova_result <- aov(Total_death ~ Dose, data = substance_data)
    tukey_result <- TukeyHSD(anova_result)
    print(tukey_result)
  } else {
    cat("Test de Kruskal-Wallis non paramétrique\n")
    kruskal_result <- kruskal_test(substance_data, Total_death ~ Dose)
    dunn_result <- dunn_test(substance_data, Total_death ~ Dose, p.adjust.method = "bonferroni")
    print(kruskal_result)
    print(dunn_result)
  }
  
  # Visualisation
  ggplot(substance_data, aes(x = Dose, y = Total_death)) +
    geom_boxplot(aes(fill = Dose)) +
    labs(title = paste("Effet des doses de", unique(substance_data$Active_substance)),
         x = "Dose", y = "Nombre total de morts") +
    theme_minimal() +
    theme(legend.position = "none")
}

# Appliquer l'analyse à chaque substance active
clone_data %>%
  group_by(Active_substance) %>%
  group_walk(~ test_toxicity(.x))




###########################################################################/
#BROUILLON
##########################################################################/
# TEST CAROLINE - Effect of repetition on Live and Total_death

library(lme4)

temp$Repetition<-factor(temp$Repetition,levels=c("1","2","3","4"))

MortRep<-glmer(cbind(Live,Total_death)~Clone + Active_substance + Dose + (1 | Clone/Dose) + (1 |Repetition),
               family =binomial, data=temp)
summary(MortRep)


# All factors
glmm_model1 <- glmer(cbind(Live, Total_death) ~ Clone * Active_substance * Dose * Repetition + 
                       (1 | Repetition) + (1 | Clone/Dose), 
                     family = binomial, data = temp)
summary(glmm_model1)

#Filter one clone, one AS, one dose

filtered_data<-temp %>%
  filter(Clone == "17-0153-0020", 
         Active_substance == "Clodinafop", 
         Dose == "N")
print(filtered_data)
glm_model<-glm(cbind(Live, Total_death)~ Repetition, 
               family = binomial, data = filtered_data)
summary(glm_model)

library(glmmTMB)

glmm_model <- glmmTMB(cbind(Live, Total_death) ~ Repetition + (1 | Repetition),
                      family = binomial, data = filtered_data)
summary(glmm_model)



#Firstly, do herbicides exert a selection pressure on Myzus persicae clone 17-0153-0020?


MyzHerbi$Active_substance<-factor(MyzHerbi$Active_substance,
                                  levels=c("Clodinafop","Cycloxydim",
                                           "Pinoxaden","Clethodim",
                                           "Fluazifop","Quizalofop"))
temp<-MyzHerbi[MyzHerbi$Dose!="N/2",]    

#Filter data for reference clone 17-0173-0020
clone_ref_temp <-temp%>% filter(Clone == "17-0153-0020")

#Calculate the mortality rate (proportion of deaths in relation to the total)
clone_ref_temp <- clone_ref_temp %>%
  mutate(mortality_rate = Total_death / Total)

#View the distribution of mortality rates according to active substances
ggplot(clone_ref_temp, aes(x = Active_substance, y = mortality_rate, fill = Active_substance)) +
  geom_boxplot() +
  labs(title = "Mortality by Active Substance (Clone 17-0173-0020)",
       x = "Active substance", y = "Mortality rate") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

model <- glm(cbind(Total_death, Total - Total_death) ~ Active_substance, 
             data = clone_ref_temp, 
             family = binomial(link = "logit"))
summary(model)
anova(model, test = "Chisq")
