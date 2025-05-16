#####code envoyé par Bertrand

library(gridExtra)
library(lme4)
library(DHARMa) # validité des modeles
library(car) #Anova
library(PROreg) # betaBinomial mixed model
library(multcomp)  # tuckey
library(effects) # plot model predicts
library(lsmeans)
library(glmmTMB)



#création variable y
Global_Campa_mod_all_c$y<-cbind(Global_Campa_mod_all_c$Nb_presence, 
                                Global_Campa_mod_all_c$Nb_absence)

table(Global_Campa_mod_all_c$f_Vegetation, useNA = "always")

glmerCampa_Sub_c<-glmmTMB(Nb_presence ~ # soit y (composee pres/abs), soit presence (0/1 sur rang ou parcelle) soit Nb_presence avec offset = log(Nb_inter)
                            Year + Season  + s_cumul_precipitation_15j + s_area +      # var globales
                            irrigation_simple + filet_type1 + enherbement_sur_rang + gestion_rang_buttage +  # var pratiques
                            pourcent_vergers_filet_monoparcelle_buffer_1000  + pourcent_tout_vergers_buffer_1000 + pourcent_tout_prairial_buffer_1000 + #s_pourcent_maraichage_buffer_1000 +
                            pourcent_vergers_filet_monoparcelle_buffer_1000:filet_type1 +
                            (1|Parcelle)+ (1|observateur_simplifie),           # effets aleatoires
                          family =nbinom1(),  # Binomiale  negbin1, negbin2, tweedie(), betabinomial(link = "logit")
                          offset = log(Nb_inter),
                          data = Global_Campa_mod_all_c,
                          na.action = na.omit)


AIC(glmerCampa_Sub_c)

overdisp_fun = function(model) {
  sum( residuals(model, type = "pearson")^2)/df.residual(model)
}

overdisp_fun(glmerCampa_Sub_c)


## petit script DHARma

rm(simulationOutput)
simulationOutput <- simulateResiduals(fittedModel =glmerCampa_Sub_c, n = 5000,use.u = T)
plot(simulationOutput)

#testZeroInflation(simulationOutput)
summary(glmerCampa_Sub_c)
Anova(glmerCampa_Sub_c)
vif(glmerCampa_Sub_c)


# plot des predit
p1<-predictorEffects(glmerCampa_Sub_c, ~ pourcent_vergers_filet_monoparcelle_buffer_1000)

plot(p1, axes=list(grid=TRUE,y=list(lab="Probabilite de présence")), layout= c(3,1))



# Tuckey tests LSMEANS = prise en compte de l'interaction

lsmeans(glmerCampa_Sub_c, pairwise ~ irrigation_simple)

lsmeans(glmerCampa_Sub_c, pairwise ~ gestion_rang_buttage)