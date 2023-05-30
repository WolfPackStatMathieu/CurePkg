function_estim_doses_comp<-function(n,probabilite_a_priori,t_star,type1,type2,graine=133){
  nb_doses<-length(probabilite_a_priori)
  require(dfcrm)
  require(flexsurvcure)
  require(cuRe)
  require(npcure)
  df<-matrix(NA,n,4)
  df<-as.data.frame(df)
  data_returns<-as.data.frame(matrix(NA,nb_doses,4))
  colnames(data_returns)<-c("estimateur_bernoulli","estimateur_survie","estimateur_guerison","p")
  colnames(df)<-c("dose","statut","tox_time","is_observed")
  df$dose<-sample(c(1:nb_doses),n,replace=TRUE)
  for (k in c(1:nb_doses)){
    index_dosek<-which(df$dose==k)
    # print("passé par là")
    p<-probabilite_a_priori[k]
    n_k<-length(index_dosek)
    df[index_dosek,]<-cbind(rep(k,n_k),generation_comp(p_cause1=p,p_cause2=1-p,t_star,nombre_obs=n_k,type1=type1,type2=type2,graine=graine))
    df$tox_time<-ifelse(df$statut==2,t_star+1,df$tox_time)
    df$is_observed<-ifelse(df$tox_time>t_star,0,1)
    data_returns[k,"estimateur_bernoulli"]<-fonction_Bern(df[index_dosek,])
    data_returns[k,"p"]<-p
  }

#' Calcul le biais pour chaque dose pour taille d'échantillon fixée
#'
#' @param K Nombre d'échantillon
#' @param n taille d'échantillon
#' @param probabilite_a_priori probabilité à priori
#' @param t_star fin de la fenêtre d'observation
#' @param type1 forme de la fonction de risque instantané (constant, increasing or decreasing)
#' @param type2 forme de la fonction de risque instantané (constant, increasing or decreasing)
#' @param graine_depart graine fixée pour la reproduction
#'
#' @return Biais pour chaque dose
#' @export
#'
#' @examples
#' K <- 1000
#' t_star <- 6
#' n <- 18
#' graine_depart <- 133
#' type1 <- constant
#' type2 <- constant
#' probabilite_a_priori <- c(0,5,0.75,0.2, 0.3)
#' generation_comp_mean(K,n,probabilite_a_priori,t_star,type1,type2,graine_depart)
#'
generation_comp_mean<-function(K,n,probabilite_a_priori,t_star,type1,type2,graine_depart){
    require(ggplot2)
    require(gridExtra)
    graine_debut<-graine_depart+1
    graine_fin<-graine_depart+K
    ensemble_graine<-c(graine_depart:graine_fin)
    result<-lapply(ensemble_graine,function_estim_doses_comp,n=n,probabilite_a_priori=probabilite_a_priori,t_star=t_star,type1=type1,type2=type2)
    nb_doses<-length(probabilite_a_priori)
    matrice<-as.data.frame(matrix(NA,nb_doses,5))
    colnames(matrice)<-c("numero_dose","modele_bernoulli","modele_survie","modele_guerison","p")
    matrice$numero_dose<-c(1:nb_doses)
    for(j in c(1:nb_doses)){
      ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
      ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
      ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
      ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
      ensemble_obs_dosek$estimateur_survie<-as.numeric(unlist(ensemble_obs_dosek$estimateur_survie))
      ensemble_obs_dosek$p<-as.numeric(ensemble_obs_dosek$p)
      matrice[j,c("modele_bernoulli","modele_survie","modele_guerison","p")]<-colMeans(ensemble_obs_dosek)
    }
    return(matrice)
  }
  fonction_surv<-Surv(as.numeric(df$tox_time),event=df$is_observed)
  indice_cens<-which(df$is_observed==0)
  df$factdose<-as.factor(df$dose)
  if(length(indice_cens)==0){
    df2<-df[,c("dose","factdose","is_observed","tox_time")]
    estimateur_surv<-rep(1,nb_doses)
    #Prob_whole_cure<-probcure(x=factdose,t=tox_time,dataset = df,d=is_observed,x0=dose_recalibree,h=c(1,1.5,2),local=FALSE)
    #estimateur_cure<-1-Prob_whole_cure[,2]
    #Prob_whole_curefx<-result<-flexsurvcure(Surv(tox_time,is_observed)~factdose+0, data=df, dist="weibullPH",
    #  anc=list(scale=~factdose+0))
    Prob_whole_cure<-fit.cure.model(Surv(tox_time,is_observed) ~ factdose,formula.surv=list(~1,~factdose),
                                    data =df2,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2,nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimateur_surv,estimateur_cure)
  }
  else{
    # print("passé par là 2")
    fit_surv <- survfit(fonction_surv ~factdose, data = df)
    # print("passé par là 3")
    df2<-df[,c("dose","factdose","is_observed","tox_time")]
    Prob_whole_cure<-fit.cure.model(Surv(tox_time,is_observed) ~ factdose, formula.surv=list(~1,~factdose),
                                    data =df2,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    #print(estimation_cure)
    #print("----------")
    estimation_surv<-rep(NA,nb_doses)
    sommaire<-summary(fit_surv,t_star,extend=TRUE)$surv
    for (j in c(1:nb_doses)){
      indice<-which(((df2$dose==j) & (df2$is_observed==1)))
      somme<-length(indice)
      if(somme>0){
        estimation_surv[j]<-1-sommaire[j]}
      else{estimation_surv[j]<-0}}
  }
  data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimation_surv,estimation_cure)
  return(data_returns)
}
########## Biais #####

evol_n_par_dose<-function(results,n,i,K=K,type1,type2){
  longueur_resultats<-c(1:length(n))
  function_intermed<-function(x,results,i){
    return(unlist(results[[x]][i,]))
  }
  result_final<-as.data.frame(t(cbind(sapply(longueur_resultats,function_intermed,results=results,i=i))))
  result_final$taille_echantillon<-n
  result_final$modele_guerison<-result_final$modele_guerison-result_final$p
  result_final$modele_bernoulli<-result_final$modele_bernoulli-result_final$p
  result_final$modele_survie<-result_final$modele_survie-result_final$p
  borne_min <- min(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  borne_max <- max(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  gg1 <- {ggplot(data = result_final, aes(x = taille_echantillon)) +
      geom_smooth(aes(y = modele_guerison, col = "modele guerison"), size = 1, alpha = 0.5) +
      geom_smooth(aes(y = modele_bernoulli, col = "modele bernoulli"), size = 1, alpha = 0.5) +
      scale_color_manual(name = "Modeles", values = c("modele guerison"="red","modele bernoulli"="blue"))+
      ggtitle("Evolution du biais  \n  en fonction de la taille d'echantillon") +
      xlab("Taille echantillon") + ylab("Biais moyen") +
      theme_classic() +
      theme(legend.title=element_blank(),
            axis.text=element_text(family = "Helvetica", size=10),
            axis.title=element_text(family = "Helvetica", size=12),
            plot.title = element_text(family = "Helvetica", size = 10)) +
      ylim(borne_min, borne_max) }
  gg2 <- {ggplot(data = result_final, aes(x = taille_echantillon)) +
      geom_smooth(aes(y = modele_guerison, col = "modele guerison"), size = 1, alpha = 0.5) +
      geom_smooth(aes(y = modele_survie, col = "modele survie"), size = 1, alpha = 0.5) +
      scale_color_manual(name = "Modeles", values = c("modele guerison"="red","modele survie"="darkgreen")) +
      ggtitle("Evolution du biais \n en fonction de la taille") +

      xlab("Taille echantillon") + ylab("Biais moyen") +
      theme_classic() +
      theme(legend.title=element_blank(),
            axis.text=element_text(family = "Helvetica", size=10),
            axis.title=element_text(family = "Helvetica", size=12),
            plot.title = element_text(family = "Helvetica", size = 10)) +
      ylim(borne_min, borne_max)+
      labs(caption = sprintf("p=%s,N=%s,type1=%s,type2=%s",
                             as.character(result_final$p),
                             as.character(K),as.character(type1),as.character(type2)))}

  gg <- {grid.arrange(gg1, gg2, ncol = 2, widths = c(8,8))}
  return(gg)
}
evol_biais_comp<-function(K,probabilite_a_priori,t_star,type1,type2,graine_depart){
  debut <- 20
  fin <- 100
  pas <- 5
  n <- seq(debut,fin , pas)
  results<-lapply(n,generation_comp_mean,K=K,probabilite_a_priori=probabilite_a_priori,t_star=t_star,type1=type1,type2=type2,graine_depart=graine_depart)
  ensemble_ggplots_par_dose<-lapply(c(1:length(probabilite_a_priori)),evol_n_par_dose,results=results,n=n,K=K,type1,type2)
  return(ensemble_ggplots_par_dose)
}


#### EQM ####
evol_eqm_comp<-function(K,probabilite_a_priori,t_star,type1,graine_depart,type2){
  debut <- 20
  fin <- 100
  pas <- 5
  n <- seq(debut,fin , pas)
  results<-lapply(n,generation_comp_eqm,K=K,probabilite_a_priori=probabilite_a_priori,t_star=t_star,type1=type1,graine_depart=graine_depart,type2=type2)
  ensemble_ggplots_par_dose<-lapply(c(1:length(probabilite_a_priori)),evol_n_par_dose_eqm,results=results,n=n,K=K,type1=type1,type2=type2)
  return(ensemble_ggplots_par_dose)
}
evol_n_par_dose_eqm<-function(results,n,i,K=K,type1,type2){
  require(ggplot2)
  require(gridExtra)
  longueur_resultats<-c(1:length(n))
  function_intermed<-function(x,results,i){
    return(unlist(results[[x]][i,]))
  }
  result_final<-as.data.frame(t(cbind(sapply(longueur_resultats,function_intermed,results=results,i=i))))
  result_final$taille_echantillon<-n
  result_final$modele_guerison<-result_final$modele_guerison
  result_final$modele_bernoulli<-result_final$modele_bernoulli
  result_final$modele_survie<-result_final$modele_survie
  borne_min <- min(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  borne_max <- max(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)

  gg1 <- {ggplot(data = result_final, aes(x = taille_echantillon)) +
      geom_smooth(aes(y = modele_guerison, col = "modele guerison"), size = 1, alpha = 0.5) +
      geom_smooth(aes(y = modele_survie, col = "modele survie"), size = 1, alpha = 0.5) +
      scale_color_manual(name = "Modeles", values =  c("modele guerison"="red", "modele survie"="darkgreen" )) +
      # ggtitle("Evolution de l'EQM en fonction de la \ntaille d'echantillon") +
      xlab("Taille echantillon") + ylab("EQM") +
      theme_classic() +

      ylim(borne_min, borne_max) +
      labs(caption = sprintf("p=%s, N=%s, type1=%s, type2=%s",
                             as.character(result_final$p),
                             as.character(K),as.character(type1),as.character(type2))
      )
  }+
    theme(legend.title=element_blank(),
          axis.text = element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24)
          , legend.text = element_text(family = "Helvetica", size = 20)
          ,text = element_text(size=rel(20))
    )

  gg2 <- {ggplot(data = result_final, aes(x = taille_echantillon)) +
      geom_smooth(aes(y = modele_guerison, col = "modele guerison"), size = 1, alpha = 0.5) +
      geom_smooth(aes(y = modele_bernoulli, col = "modele bernoulli"), size = 1, alpha = 0.5) +
      scale_color_manual(name = "Modeles", values = c("modele guerison"="red", "modele bernoulli"="blue" )) +
      # ggtitle("Evolution de l'EQM en fonction de la \ntaille d'echantillon") +
      xlab("Taille echantillon") + ylab("EQM") +
      theme_classic() +
      ylim(borne_min, borne_max)+
      labs(caption = "")}+
    theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24)
          , legend.text = element_text(family = "Helvetica", size = 20)
          ,text = element_text(size=rel(20))
    )

  gg <- {grid.arrange(gg2, gg1, ncol = 2, widths = c(8,8)
                      ,top =textGrob(paste("Evolution de l'EQM en fonction de la taille d'echantillon, dose ",i),gp=gpar(fontsize=24,font=3))
  )}
  return(gg)
}
generation_comp_eqm<-function(K,n,probabilite_a_priori,t_star,type1,graine_depart,type2){
  require(ggplot2)
  require(gridExtra)
  require(grid)
  graine_ensemble<-graine_depart+c(1:K)
  result<-lapply(graine_ensemble,function_estim_doses_comp,n=n,probabilite_a_priori=probabilite_a_priori,t_star=t_star,type1=type1,type2=type2)
  nb_doses<-length(probabilite_a_priori)
  matrice<-as.data.frame(matrix(NA,nb_doses,5))
  # print("chien")
  colnames(matrice)<-c("numero_dose","modele_bernoulli","modele_survie","modele_guerison","p")
  # print("chat")
  matrice$numero_dose<-c(1:nb_doses)
  for(j in c(1:nb_doses)){
    ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
    ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
    ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
    ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
    ensemble_obs_dosek$estimateur_survie<-as.numeric(unlist(ensemble_obs_dosek$estimateur_survie))
    ensemble_obs_dosek$p<-probabilite_a_priori[j]
    matrice[j,c("modele_bernoulli","modele_survie","modele_guerison")]<-colMeans((ensemble_obs_dosek-ensemble_obs_dosek$p)^2)[1:3]
    matrice[j,"p"]<-probabilite_a_priori[j]
  }
  return(matrice)
}
