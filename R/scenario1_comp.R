# Generation
# source("estimateurs.R")
# source("fonctions_simulations_competition.R")
# source("http://myweb.uiowa.edu/pbreheny/7210/f19/notes/fun.R")

#' Generer un echantillon selon un modele e risques competitifs
#'
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param p_cause2 probabilite pour la cause 1 e t_star
#' @param t_star fin de la fenêtre d'observation
#' @param nombre_obs taille d'echantillon
#' @param graine graine fixee pour la reproduction
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#'
#' @return un dataframe avec le statut, temps de toxicite et observation de toxicite des patients
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' p_cause2 <-1-p_cause1
#' t_star <- 6
#' nombre_obs <- 18
#' graine <- 133
#' type1 <- "constant"
#' type2 <- "constant"
#' generation_comp(p_cause1,p_cause2,t_star,nombre_obs,graine,type1,type2)
#'
generation_comp<-function(p_cause1,p_cause2,t_star,nombre_obs,graine,type1,type2){

  alpha1<-get_alpha(p_cause1,obswin=t_star,typ="weibull",typ_wb=type1)
  alpha2<-get_alpha(p_cause2,obswin=t_star,typ="weibull",typ_wb=type2)
  liste_dataset<-get_dataset0(n=nombre_obs,alpha1,alpha2,tstar=t_star,graine=graine,K=1,type="weibull")
  data<-liste_dataset$data_complete
  data<-as.data.frame(data)
  data$is_observed<-ifelse(data$status==0,0,1)
  data_estim<-data[,c("status","time","is_observed")]
  colnames(data_estim)<-c("status","tox_time","is_observed")
  return(data_estim)
}


# Estimateurs.

#' Generer un graphique de la proportion de censures dans les echantillons selon le modele
#' e risques competitifs
#'
#' @param N nombre d'echantillon
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param n taille d'echantillon
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param t_star fin de la fenêtre d'observation
#' @param graine graine fixee pour la reproduction
#'
#' @return Graphique de la proportion de censures
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' t_star <- 6
#' N <- 1000
#' n <- 18
#' type1 <- "constant"
#' type2 <- "constant"
#' prop_censure_alt(N,p_cause1,n,type1,type2,t_star,graine=133)
#'
prop_censure_alt <- function(N,p_cause1,n,type1,type2,t_star,graine=133){

  p <- 33
  p_cause2<-(1-p_cause1)

  res <- list()
  censures <- rep(NA,N)
  TDL <- rep(NA,N)
  gueris <- rep(NA,N)
  for(i in 1:N){

    df <- generation_comp(p_cause1=p_cause1, p_cause2=p_cause2,
                          t_star=t_star,nombre_obs=n,
                          type1=type1,type2=type2,graine=graine+i)


    res[[i]] <- df
    nb_status0 <- length(which(df$status==0))
    nb_status1 <- length(which(df$status==1))
    nb_status2 <- length(which(df$status==2))

    censures[i] <- (nb_status0/n)*100     # censures
    TDL[i] <- (nb_status1/n)*100    # non gueris (toxicit?)
    gueris[i] <- (nb_status2/n)*100    # gueris
  }
  # censures_mean <- mean(censures)
  # TDL_mean <- mean(TDL)
  # gueris_mean <- mean(gueris)

  dens_censure <- density(censures)
  dens_TDL <- density(TDL)

  par(mfrow = c(1,2))

  plot(x=dens_TDL$x, y=dens_TDL$y, main="TDL", type="l",xlab="TDL", ylab="densit?")
  abline(v=p, col="red")
  # Add legend
  legend("topright", # Position of the legend
         "cible TDL", # Labels
         col = "red",
         pch = c(19, NA), # Point markers (NA means no marker)
         lty = c(1, 1), # Line styles (1 means solid)
         lwd = c(2, 2),
         bty ="n",
         cex = 0.6) # Line widths
  plot(x=dens_censure$x, y=dens_censure$y, main="Censures", type="l", xlab="censure", ylab="densite")
  # M <- cbind(TDL, censures)
  # dnplot(M[,1])
  # dnplot(M[,1], pos=TRUE)
  # dnplot(M)
  # dnplot(M, labs=c('TDL', 'Censures'))
  # dnplot(M, labs=c('TDL', 'Censures'))

}

#' Calculer le biais des 3 estimateurs une fois selon le modele e risques competitifs
#'
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param n taille d'echantillon
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param t_star fin de la fenêtre d'observation
#' @param graine graine fixee pour la reproduction
#'
#' @return Biais des 3 estimateurs
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' t_star <- 6
#' n <- 18
#' type1 <- constant
#' type2 <- constant
#' fonction_estim_comp_once(N,p_cause1,n,type1,type2,t_star,graine=133)
#'
fonction_estim_comp_once<-function(p_cause1,n,type1,type2,t_star,graine=133){

  p_cause2<-(1-p_cause1)
  data<-generation_comp(p_cause1 = p_cause1,p_cause2=p_cause2,t_star=t_star,nombre_obs = n,type1=type1,type2=type2,graine = graine)
  data$tox_time<-ifelse(data$status==2,t_star+1,data$tox_time)
  data$is_observed<-ifelse(data$tox_time>t_star,0,1)
  indices_non_obs<-which(data$is_observed==0)
  if(length(indices_non_obs)==n){
    # tous pas observes donc censures
    estimateursurv<-0
    estimateurbern<-0
    estimateurcure<-fonction_cure(df=data,t_star)
    sous_liste<-list(estimateursurv, estimateurbern, estimateurcure)
    names(sous_liste)<-c("Survie","Bernoulli","Guerison")
    return(sous_liste)
  }
  estimateursurv<-fonction_KM(df=data,t_star)
  estimateurbern<-fonction_Bern(df=data)
  estimateurcure<-fonction_cure(df=data,t_star)
  sous_liste<-list(estimateursurv,estimateurbern,estimateurcure)
  names(sous_liste)<-c("Survie","Bernoulli","Guerison")
  return(sous_liste)
}

#' Calculer le biais des 3 estimateurs plusieurs fois selon le modele e risques competitifs
#'
#' @param K Nombre d'echantillon
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param n taille d'echantillon
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2  forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param t_star fin de la fenêtre d'observation
#' @param graine graine fixee pour la reproduction
#'
#' @return Biais des 3 estimateurs sur N echantillons
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' t_star <- 6
#' n <- 18
#' N <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' Simuler_estim_mult_times(N,p_cause1,n,type1,type2,t_star,graine=133)
#'
Simuler_estim_mult_times<-function(K,p_cause1,n,type1,type2,t_star,graine){

  graine_inf <- graine
  graine_sup <- graine + K-1
  ensemble_graine<-c(graine_inf:graine_sup)
  result<-cbind(sapply(ensemble_graine,fonction_estim_comp_once,p_cause1=p_cause1,type1=type1,type2=type2,t_star=t_star,n=n))
  result<-as.data.frame(t(result))
  colnames(result)<-c("Survie","Bernoulli","Guerison")
  return(colMeans(sapply(result,as.numeric)))
}

# Biais

#' Calcul le biais selon differentes tailles d'echantillon, selon le modele e risques comp.
#'
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param K Nombre d'echantillon
#' @param t_star fin de la fenêtre d'observation
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param graine graine fixee pour la reproduction
#'
#' @return Biais des 3 estimateurs
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' t_star <- 6
#' K <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' biais.selon.lambda_alt(p_cause1=p_cause1,K=K, type1=type1,type2=type2,t_star=t_star,graine=133)
#'
biais.selon.lambda_alt <-function(p_cause1,K,t_star,type1,type2,graine){

  results <- NULL
  n <- 20
  while(n<100){
    vec.biais <- Simuler_estim_mult_times(K=K,p_cause1=p_cause1,n=n,type1=type1,type2=type2,t_star=t_star,graine=graine)
    biais_surv<-vec.biais[[1]]-p_cause1
    biais.bern<-vec.biais[[2]]-p_cause1
    biais.cure<-vec.biais[[3]]-p_cause1
    results<-rbind(results,c(n,biais_surv,biais.cure,biais.bern))
    n <- n+5
  }
  return(results)
}


# EQM

#'  Calcul l'EQM selon differentes tailles d'echantillon, selon le modele e risques comp.
#'
#' @param p_cause1 probabilite pour la cause 1 e t_star
#' @param K Nombre d'echantillon
#' @param t_star fin de la fenêtre d'observation
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param graine graine fixee pour la reproduction
#'
#' @return EQM des 3 estimateurs
#' @export
#'
#' @examples
#' p_cause1 <- 0.5
#' t_star <- 6
#' K <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' eqm.selon.alpha(p_cause1=_cause1,K=K,type1=type1,type2=type2,t_star=t_star,graine=133)
#'
eqm.selon.alpha<-function(p_cause1,K,t_star,type1,type2,graine){

  results <- NULL
  n <- 20
  graine_inf <- graine
  graine_sup <- graine + K
  ensemble_graine<-c(graine_inf:graine_sup)
  while(n<200){
    liste_global<-as.data.frame(t(cbind.data.frame(sapply(ensemble_graine,fonction_estim_comp_once,
                                                          p_cause1=p_cause1,type1=type1,type2=type2,t_star=t_star,n=n))))
    liste_global$Guerison<-as.numeric(liste_global$Guerison)
    liste_global$Survie<-as.numeric(liste_global$Survie)
    liste_global$Bernoulli<-as.numeric(liste_global$Bernoulli)
    valeurs<-colMeans((liste_global-p_cause1)^2)
    results<-rbind(results,c(n,valeurs[1],valeurs[3],valeurs[2]))
    n <- n+10
  }
  return(results)
}


# Graphiques

#' Evolution du biais en fonction de la taille d'echantillon, selon le modele e risques comp.
#'
#' @param N Nombre d'echantillon
#' @param t_star fin de la fenêtre d'observation
#' @param p proportion de TDL
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param graine graine fixee pour la reproduction
#'
#' @return Graphique du biais en fonction de n
#' @export
#'
#' @examples
#' p <- 0.3
#' t_star <- 6
#' N <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' fonction_ggplot_evol_biais_alt(N=N,t_star=t_star, p=p,type1=type1,type2=type2,graine=133)
#'
fonction_ggplot_evol_biais_alt <- function(N,t_star, p,type1,type2,graine=133) {

  library(gridExtra)
  library(ggplot2)
  library(scales)
  # Generate the data
  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=p,K=N, t_star=t_star,type1,type2,graine=graine)
  RES0.3.3 <- data.frame(RES)
  colnames(RES0.3.3) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")
  borne_min <- min(RES0.3.3$mean.surv, RES0.3.3$mean.cure,RES0.3.3$mean.bernoulli)
  borne_max <- max(RES0.3.3$mean.surv, RES0.3.3$mean.cure,RES0.3.3$mean.bernoulli)

  gg1 <- ggplot(data =RES0.3.3, aes(x = n)) +
    geom_smooth(aes(y = mean.cure, col = "modele guerison"), size = 1, alpha = 0.5) +
    geom_smooth(aes(y = mean.bernoulli, col = "modele bernoulli"), size = 1, alpha = 0.5) +
    scale_color_manual(name = "Modeles", values=c("modele guerison"="red1","modele bernoulli"="blue1")) +
    ggtitle("Evolution du biais moyen en \n fonction de n") +
    xlab("Taille echantillon") + ylab("Biais moyen") +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24)
          , legend.text = element_text(family = "Helvetica", size = 20)
          ,text = element_text(size=rel(20))) +
    ylim(borne_min, borne_max)

  gg2 <- ggplot(data = RES0.3.3, aes(x = n)) +
    geom_smooth(aes(y = mean.cure, col = "modele guerison"), size = 1, alpha = 0.5) +
    geom_smooth(aes(y = mean.surv, col = "modele survie"), size = 1, alpha = 0.5) +
    scale_color_manual(name = "Modeles", values=c("modele guerison"="red1","modele survie"="darkgreen")) +
    ggtitle("Evolution du biais moyen en \n fonction de n") +
    xlab("Taille echantillon") + ylab("Biais moyen") +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24)
          , legend.text = element_text(family = "Helvetica", size = 20)
          ,text = element_text(size=rel(20))) +
    ylim(borne_min, borne_max)+
    labs(caption = sprintf("N = %s, n variant de %s a %s \n par pas de %s,type1=%s,type2=%s, p=%s" ,
                           as.character(N),
                           as.character(20),
                           as.character(100),
                           as.character(5),
                           as.character(type1),
                           as.character(type2),
                           as.character(p)))

  gg <- grid.arrange(gg1, gg2, ncol = 2, widths = c(8,8))
}


#' Boxplot des biais sur plusieurs echantillons selon le modele e risques comp.
#'
#' @param K Nombre d'echantillon
#' @param n taille d'echantillon
#' @param p proportion de TDL
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param t_star fin de la fenêtre d'observation
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param graine graine fixee pour la reproduction
#'
#' @return Boxplot des 3 estimateurs
#' @export
#'
#' @examples
#' p <- 0.3
#' t_star <- 6
#' K <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' n <- 18
#' plots_scenario_1_alt(K=K, n=n, t_star=t_star, p=p,type1=type1,type2=type2,graine=133)
#'
plots_scenario_1_alt <- function(K, n, p,type1,t_star,type2,graine=133){

  require(ggplot2)
  require(dplyr)
  require(tidyr)
  # df ? 3 colones (mod?le de gu?rison, mod?le de survie, mod?le de bernouilli)
  graine_liste<-graine+c(1:K)
  res <-as.data.frame(t(cbind.data.frame(sapply(graine_liste,fonction_estim_comp_once,n=n,p_cause1=p,type1=type1,type2=type2,t_star=t_star))))
  res$Survie<-as.numeric(res$Survie)
  res$Bernoulli<-as.numeric(res$Bernoulli)
  res$Guerison<-as.numeric(res$Guerison)
  res <- res - p
  # on renomme les colonnes

  # bornes
  borne_min <- min(res)
  borne_max <- max(res)
  # On tranforme les colonnes d?j? pr?sentes en une seule colonne (valeurs)
  # ensuite ajouter une nouvelle colonne modele qui servira a
  # distinguer les 2 mod?les
  df <- res %>% gather(key = "modele", value = "valeurs")

  # boxplot
  boxplot <- ggplot(df, aes(x = modele, y = valeurs, fill = modele)) +
    geom_violin(alpha = 0.8) +
    scale_fill_manual(values = c("#0072B2", "#E69F00","purple")) +
    # theme_classic()+
    ylim(borne_min, borne_max)

  # Add labels and title
  boxplot +
    labs(x = "Modeles", y = "Biais",
         title = "Comparaison du biais pour N simulations et n fixe",
         caption = sprintf("N = %s, p=%s,n=%s,type1=%s,type2=%s",as.character(K),as.character(p),as.character(n),type1,type2)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20)
          ,plot.subtitle = element_text(hjust = 0, size = 10)
          ,axis.text = element_text(size = 15)
          ,axis.title = element_text(size = 15)
          ,legend.text = element_text(size = 12)
          , legend.title = element_text(size = 15)
          , plot.caption = element_text(size = 12)
          # ,text = element_text(size=rel(8))
    )

}

# Utilisation des parametres


#' Evolution du biais pour differentes valeurs de p_cause1
#'
#' @param N Nombre d'echantillon
#' @param t_star fin de la fenêtre d'observation
#' @param vect_cause1 vecteur de differentes probabilite de cause 1
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param graine graine fixee pour la reproduction
#'
#' @return Graphique du biais selon les p_cause1
#' @export
#'
#' @examples
#' t_star <- 6
#' N <- 1000
#' type1 <- "constant"
#' type2 <- "constant"
#' fonction_compar_plotsn_lambda_alt_8p <- function(N,t_star,
#' vect_cause1=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8),type1,type2,graine=133)
#'
fonction_compar_plotsn_lambda_alt_8p <- function(N,
                                                 t_star,
                                                 vect_cause1=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8),
                                                 type1,type2,graine=133) {

  library(gridExtra)
  library(ggplot2)
  library(scales)
  # Generate the data
  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[1]],K=N, t_star=t_star,type1,type2,graine=graine)
  RES0.2 <- data.frame(RES)
  colnames(RES0.2) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")
  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[2]],K=N, t_star=t_star,type1,type2,graine=graine)
  RES0.3 <- data.frame(RES)
  colnames(RES0.3) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[3]],K=N,t_star=t_star,type1,type2,graine=graine)
  RES0.4 <- data.frame(RES)
  colnames(RES0.4) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[4]],K=N,t_star=t_star,type1,type2,graine=graine)
  RES0.5 <- data.frame(RES)
  colnames(RES0.5) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[5]],K=N,t_star=t_star,type1,type2,graine=graine)
  RES0.6 <- data.frame(RES)
  colnames(RES0.6) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[6]],K=N,t_star=t_star,type1,type2,graine=graine)
  RES0.7 <- data.frame(RES)
  colnames(RES0.7) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda_alt(p_cause1=vect_cause1[[7]],K=N,t_star=t_star,type1,type2,graine=graine)
  RES0.8 <- data.frame(RES)
  colnames(RES0.8) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")
  # Get the min and max bounds of each variable to be used in the plots
  # les bornes min et max du modele de survie
  borne_min <- min(
    min(RES0.2$mean.surv),
    min(RES0.4$mean.surv),
    min(RES0.5$mean.surv)
    ,min(RES0.6$mean.surv)
    ,min(RES0.7$mean.surv)
    ,min(RES0.8$mean.surv)
    ,min(RES0.7$mean.surv)
  )

  borne_max <- max(
    max(RES0.2$mean.surv),
    max(RES0.3$mean.surv),
    max(RES0.4$mean.surv)
    ,max(RES0.5$mean.surv)
    ,max(RES0.6$mean.surv)
    ,max(RES0.7$mean.surv)
    ,max(RES0.8$mean.surv)
  )
  # les bornes min et max du modele de guerison
  borne_min.c <- min(
    min(RES0.2$mean.cure),
    min(RES0.3$mean.cure),
    min(RES0.4$mean.cure)
    ,min(RES0.5$mean.cure)
    ,min(RES0.6$mean.cure)
    ,min(RES0.7$mean.cure)
    ,min(RES0.8$mean.cure)
  )

  borne_max.c <- max(
    max(RES0.2$mean.cure),
    max(RES0.3$mean.cure),
    max(RES0.4$mean.cure)
    ,max(RES0.5$mean.cure)
    ,max(RES0.6$mean.cure)
    ,max(RES0.7$mean.cure)
    ,max(RES0.8$mean.cure)
  )

  # les bornes min et max du modele de Bernoulli
  borne_min.b <- min(
    min(RES0.2$mean.bernoulli),
    min(RES0.3$mean.bernoulli),
    min(RES0.4$mean.bernoulli)
    ,min(RES0.5$mean.bernoulli)
    ,min(RES0.6$mean.bernoulli)
    ,min(RES0.7$mean.bernoulli)
    ,min(RES0.8$mean.bernoulli)
  )

  borne_max.b <- max(
    max(RES0.2$mean.bernoulli),
    max(RES0.3$mean.bernoulli),
    max(RES0.4$mean.bernoulli)
    ,max(RES0.5$mean.bernoulli)
    ,max(RES0.6$mean.bernoulli)
    ,max(RES0.7$mean.bernoulli)
    ,max(RES0.8$mean.bernoulli)
  )
  # Plot the data
  # le modele de survie
  gg1 <-  ggplot(RES0.2, aes(n, mean.surv)) +
    geom_line(aes(color = "0.2"), size = 0.6) +
    geom_line(data = RES0.3, aes(n, mean.surv, color = "0.3"), size = 0.6) +
    geom_line(data = RES0.4, aes(n, mean.surv, color = "0.4"), size = 0.6) +
    geom_line(data = RES0.5, aes(n, mean.surv, color = "0.5"), size = 0.6) +
    geom_line(data = RES0.6, aes(n, mean.surv, color = "0.6"), size = 0.6)+
    geom_line(data = RES0.7, aes(n, mean.surv, color = "0.7"), size = 0.6) +
    geom_line(data = RES0.8, aes(n, mean.surv, color = "0.8"), size = 0.6) +
    scale_color_manual(name = "p1", values = c("#0072B2", "red", "#009E73", "#F0E442",
                                               "purple", "#D55E00","blue")) +
    # scale_colour_colorblind() +
    ylim(borne_min -0.04, borne_max+0.04)+
    labs(
      title = "Modele de survie",
      x = "n",
      y = "biais moyen",
      color = "n")+
    theme_bw()
  # le modele de guerison
  gg2 <-  ggplot(RES0.2, aes(n, mean.cure)) +
    geom_line(aes(color = "0.2"), size = 0.6) +
    geom_line(data = RES0.3, aes(n, mean.cure, color = "0.3"), size = 0.6)+
    geom_line(data = RES0.4, aes(n, mean.cure, color = "0.4"), size = 0.6) +
    geom_line(data = RES0.5, aes(n, mean.cure, color = "0.5"), size = 0.6) +
    geom_line(data = RES0.6, aes(n, mean.cure, color = "0.6"), size = 0.6) +
    geom_line(data = RES0.7, aes(n, mean.cure, color = "0.7"), size = 0.6) +
    geom_line(data = RES0.8, aes(n, mean.cure, color = "0.8"), size = 0.6) +
    scale_color_manual(name = "p1", values = c("#0072B2", "red", "#009E73", "#F0E442",
                                               "purple", "#D55E00","blue")) +
    # scale_colour_colorblind() +
    ylim(borne_min.c -0.04, borne_max.c+0.04)+
    labs(
      title = "Modele de guerison",
      x = "n",
      y = "biais moyen",
      color = "n")+
    theme_bw()
  #le modele de Bernoulli
  gg3 <-  ggplot(RES0.2, aes(n, mean.bernoulli)) +
    geom_line(aes(color = "0.2"), size = 0.6) +
    geom_line(data = RES0.3, aes(n, mean.bernoulli, color = "0.3"), size = 0.6) +
    geom_line(data = RES0.4, aes(n, mean.bernoulli, color = "0.4"), size = 0.6) +
    geom_line(data = RES0.5, aes(n, mean.bernoulli, color = "0.5"), size = 0.6) +
    geom_line(data = RES0.6, aes(n, mean.bernoulli, color = "0.6"), size = 0.6) +
    geom_line(data = RES0.7, aes(n, mean.bernoulli, color = "0.7"), size = 0.6) +
    geom_line(data = RES0.8, aes(n, mean.bernoulli, color = "0.8"), size = 0.6) +
    scale_color_manual(name = "p1", values = c("#0072B2", "red", "#009E73", "#F0E442",
                                               "purple", "#D55E00","blue")) +
    # scale_colour_colorblind() +
    ylim(borne_min.b -0.04, borne_max.b+0.04)+
    labs(
      title = "Modele de Bernoulli",
      x = "n",
      y = "biais moyen",
      color = "n")+
    theme_bw()
  # on remet tout dans un seul graphique
  g <- grid.arrange(gg1, gg2, gg3, top = sprintf("Influence de n et de p1 pour un alpha de type %s", type1)
                    ,bottom = sprintf("genere avec N = %s pour chaque taille n", N),nrow=2)
  return(g)

}




### Evolution ####

######### EQM ##########
#' Evolution de l'EQM pour plusieurs echantillons selon le modele e risques competitifs
#'
#' @param K Nombre d'echantillon
#' @param type1 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantane (constant, increasing or decreasing)
#' @param p proportion de TDL
#' @param graine graine fixee pour la reproduction
#' @param t_star fin de la fenêtre d'observation
#'
#' @return Graphiques de l'EQM
#' @export
#'
#' @examples
#' K <- 1000
#' type1 <- "constant"
#' p <- 0.3
#' graine <- 133
#' t_star <- 6
#' eqm.selon.taille_echantillon_alt(K=K, type1=type1, p=p,graine=graine,t_star=t_star,type2="constant")
#'
eqm.selon.taille_echantillon_alt<-function(K, type1, p,graine,t_star,type2){

  require(ggplot2)
  require(gridExtra)
  # On fixe un n de d?part ? 10 individus et on incr?ment par 5 jusqu'a 100
  # On calcule le biais pour K simulations et n-?chantillons
  result_final<-as.data.frame(eqm.selon.alpha(K=K, type1=type1, p_cause1=p,graine=graine,type2=type2,t_star=t_star))
  colnames(result_final) <- c("taille_echantillon","modele_survie","modele_guerison", "modele_bernoulli")
  # plot
  borne_min <- min(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  borne_max <- max(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)

  gg1 <- ggplot(data = result_final, aes(x = taille_echantillon)) +
    geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
    geom_smooth(aes(y = modele_bernoulli, col = "bernoulli"), size = 1.2, alpha = 0.5) +
    scale_color_manual(name = "Mod?les", values = c("guerison" = "red1", "bernoulli" = "blue1")) +
    ggtitle("Evolution de l'eqm en \nfonction de n") +
    xlab("Taille echantillon") + ylab("EQM") +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=10),
          axis.title=element_text(family = "Helvetica", size=12),
          plot.title = element_text(family = "Helvetica", size = 10)) +
    # ylim(borne_min, borne_max)
    ylim(0.01, 0.03)

  gg2 <- ggplot(data = result_final, aes(x = taille_echantillon)) +
    geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
    geom_smooth(aes(y = modele_survie, col = "survie"), size = 1.2, alpha = 0.5) +
    scale_color_manual(name = "Mod?les", values = c("guerison" = "red1", "survie" = "darkgreen")) +
    ggtitle("Evolution de l'eqm en \nfonction de n") +
    xlab("Taille echantillon") + ylab("EQM") +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=10),
          axis.title=element_text(family = "Helvetica", size=12),
          plot.title = element_text(family = "Helvetica", size = 10)) +
    # ylim(borne_min, borne_max)+
    ylim(0.01, 0.03)+
    labs(caption = sprintf("N = %s, p=%s,type=%s,type2=%s" ,
                           as.character(K),
                           as.character(p),
                           as.character(type1),
                           as.character(type2)))

  gg <- grid.arrange(gg1, gg2, ncol = 2, widths = c(7,7))
}
