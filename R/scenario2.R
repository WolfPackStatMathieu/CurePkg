######## Premiere methode de generation ########
#' Calcul du biais pour chaque dose pour chaque estimateur avec un seul échantillon.
#'
#' @param n nombre d'individus considérés.
#' @param nb_doses Nombre de doses. Ce nombre correspond à la longueur de liste_params.
#' @param liste_params liste contenant plusieurs sous-listes. Chaque sous-liste contient le lambda et le k de chaque dose.
#' @param t_star fin de la fenetre d'observation
#' @return
#' @export
#'
#' @examples
#' ######Test ######
#' n<-100
#' k1<-1
#' lambda1<-3
#' sous_liste1<-list(k1,lambda1,0.33)
#' names(sous_liste1)<-c("k","lambda","p")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2,0.5)
#' names(sous_liste2)<-c("k","lambda","p")
#' ls<-list(sous_liste1,sous_liste2)
#' result<-function_estim_doses(n=100,liste_params=ls,nb_doses=length(ls),t_star=6)
function_estim_doses<-function(n,liste_params,nb_doses,t_star){
  require(dfcrm)
  df<-matrix(NA,n,4)
  df<-as.data.frame(df)
  data_returns<-as.data.frame(matrix(NA,nb_doses,4))
  colnames(data_returns)<-c("estimateur_bernoulli","estimateur_survie","estimateur_guerison","p")
  colnames(df)<-c("dose","sensible","tox_time","is_observed")
  df$dose<-sample(c(1:nb_doses),n,replace=TRUE)
  for (k in c(1:nb_doses)){
    index_dosek<-which(df$dose==k)
    sous_liste<-liste_params[[k]]
    n_k<-length(index_dosek)
    df[index_dosek,]<-cbind(rep(k,n_k),Generation_un_ech(n=n_k,lambda=sous_liste[["lambda"]],t_star=t_star,p=sous_liste[["p"]],k=sous_liste[["k"]]))
    data_returns[k,"estimateur_bernoulli"]<-fonction_Bern(df[index_dosek,])
    data_returns[k,"p"]<-sous_liste[["p"]]
  }
  fonction_surv<-survival::Surv(as.numeric(df$tox_time),event=df$is_observed)
  indice_cens<-which(df$is_observed==0)
  df$factdose<-as.factor(df$dose)
  if(length(indice_cens)==0){

    estimateur_surv<-rep(1,nb_doses)
    Prob_whole_cure<-cuRe::fit.cure.model(survival::Surv(tox_time,is_observed) ~ factdose, formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimateur_surv,estimateur_cure)
  }
  if(length(indice_cens)==nrow(df)){
    estimateur_surv<-rep(0,nb_doses)
    Prob_whole_cure<-cuRe::fit.cure.model(survival::Surv(tox_time,is_observed) ~ factdose,  formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    estimation_surv<-rep(0,nb_doses)
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimateur_surv,estimateur_cure)
  }
  else{
    fit_surv <- survival::survfit(fonction_surv ~factdose, data = df)
    Prob_whole_cure<-cuRe::fit.cure.model(survival::Surv(tox_time,is_observed) ~ factdose,  formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    estimation_surv<-rep(NA,nb_doses)
    sommaire<-summary(fit_surv,t_star,extend=TRUE)$surv
    for (j in c(1:nb_doses)){
      indice<-which(((df$dose==j) & (df$is_observed==1)))
      somme<-length(indice)
      if(somme>0){
        estimation_surv[j]<-1-sommaire[j]}
      else{estimation_surv[j]<-0}
    }
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimation_surv,estimation_cure)
  }

  return(data_returns)
}
############# BIAIS####
#' Calcul du biais pour chaque dose pour chaque estimateur pour K échantillons.
#'
#' @param n nombre d'individus considérés.
#' @param K nombre d'échantillons générés.
#' @param nb_doses Nombre de doses. Ce nombre correspond à la longueur de liste_params.
#' @param liste_params liste contenant plusieurs sous-listes. Chaque sous-liste contient le lambda et le k de chaque dose.
#' @param t_star fin de la fenetre d'observation
#' @return Liste. Valeur du biais pour chaque échantillon pour chaque estimateur.
#' @export
#'
#' @examples
#' ######Test ######
#' K<-10
#' n<-100
#' k1<-1
#' lambda1<-3
#' sous_liste1<-list(k1,lambda1,0.33)
#' names(sous_liste1)<-c("k","lambda","p")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2,0.5)
#' names(sous_liste2)<-c("k","lambda","p")
#' ls<-list(sous_liste1,sous_liste2)
#' biais_K<-Realisations_estim_cas_mult(K=K,n=n,liste_params=ls,nb_doses=2,t_star=6)
Realisations_estim_cas_mult<-function(K,n,liste_params,nb_doses,t_star){
  result<-lapply(rep(n,K),function_estim_doses,liste_params=liste_params,nb_doses=nb_doses,t_star=t_star)
  matrice<-list(rep(NA,nb_doses))
  for(j in c(1:nb_doses)){
    ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
    ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
    ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
    ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
    ensemble_obs_dosek$estimateur_survie<-as.numeric(ensemble_obs_dosek$estimateur_survie)
    ensemble_obs_dosek$p<-as.numeric(ensemble_obs_dosek$p)
    matrice[[j]]<-ensemble_obs_dosek}
  return(matrice)
}
#' Calculer la densité du biais dans le cas à plusieurs doses.
#'
#' @param n nombre d'individus considérés.
#' @param nb_doses Nombre de doses. Ce nombre correspond à la longueur de liste_params.
#' @param liste_params liste contenant plusieurs sous-listes. Chaque sous-liste contient le lambda et le k de chaque dose.
#' @param t_star fin de la fenetre d'observation
#' @param K nombre de répétitions de l'expérience.
#' @return
#' @export
#'
#' @examples
#' ######Test ######
#' p<-0.33
#' lambda_test<-0.33
#' t_star<-6
#' k1<-1
#' liste_parameter<-list(lambda_test,t_star,p,k1)
#' names(liste_parameter)<-c("lambda","t_star","p","k")
#' lb_test2<-0.2
#' t_star2<-7
#' p2<-0.5
#' k2<-1
#' liste_2<-list(lb_test2,t_star2,p2,k2)
#' names(liste_2)<-c("lambda","t_star","p","k")
#'vecteur_param<-list(liste_parameter,liste_2)
#'test<-plots_scenario_mult(K=10,n=100,liste_params = vecteur_param,t_star=6,nb_doses=2)
plots_scenario_mult <- function(K, n,liste_params, t_star,nb_doses){
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  res <- Realisations_estim_cas_mult(K,n,liste_params,nb_doses,t_star)
  result<-list(rep(NA,nb_doses))
  # on renomme les colonnes
  for (j in c(1:nb_doses)){
    p<-liste_params[[j]][["p"]]
    lambda<-liste_params[[j]][["lambda"]]
    k<-liste_params[[j]][["k"]]
    res_j<-as.data.frame(res[[j]][,c("estimateur_bernoulli","estimateur_guerison","estimateur_survie")])-res[[j]]$p
    colnames(res_j) <- c("Bernoulli", "Guerison", "Survie")

    # bornes
    borne_min <- min(res_j)
    borne_max <- max(res_j)


    # On tranforme les colonnes d?j? pr?sentes en une seule colonne (valeurs)
    # ensuite ajouter une nouvelle colonne modele qui servira a
    # distinguer les 2 mod?les
    df <- res_j %>% gather(key = "modele", value = "valeurs")

    # boxplot
    boxplot <- ggplot2::ggplot(df, ggplot2::aes(x = modele, y = valeurs, fill = modele)) +
      ggplot2::geom_violin(alpha = 0.8) +
      ggplot2::scale_fill_manual(values = c("#0072B2", "#E69F00","Sky blue")) +
      ggplot2::theme_classic()+
      ggplot2::ylim(borne_min, borne_max)

    # Add labels and title
    result[[j]]<-{boxplot +
        ggplot2::labs(x = "Mod?les", y = "Biais moyen",
             title = "Comparaison du biais moyen pour K n-echantillons",
             caption = sprintf("K = %s, lambda = %s, k = %s, n = %s,p=%s" ,
                               as.character(K),
                               as.character(lambda),
                               as.character(k),
                               as.character(n),
                               as.character(p))) +
        ggplot2::theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 12, face = "bold"))}
  }
  return(result)
}
#calcul du biais moyen pour k echantillons, size fixee.
calcul_mean_size_Ktimes<-function(size,vecteur_param,nb_doses,K,t_star){

  matrice_mean_doses<-as.data.frame(matrix(NA,nb_doses,5))
  colnames(matrice_mean_doses)<-c("mean_Bernoulli","mean_guerison","mean_survie","p","n")
  result<-lapply(rep(size,K),function_estim_doses,liste_params=vecteur_param,nb_doses=nb_doses,t_star=t_star)
  for(j in c(1:nb_doses)){
    ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
    ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
    p<-vecteur_param[[j]][["p"]]
    ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
    ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
    ensemble_obs_dosek$estimateur_modele_survie<-as.numeric(ensemble_obs_dosek$estimateur_survie)
    ensemble_obs_dosek$p<-as.numeric(ensemble_obs_dosek$p)
    mean_bern<-mean(ensemble_obs_dosek$estimateur_bernoulli)-p
    mean_surv<-mean(ensemble_obs_dosek$estimateur_modele_survie)-p
    mean_cure<-mean(ensemble_obs_dosek$estimateur_guerison)-p
    vecteur<-c(mean_bern,mean_cure,mean_surv,mean(ensemble_obs_dosek$p),size)
    matrice_mean_doses[j,]<-vecteur}
  return(matrice_mean_doses)
}

#calcul du biais pour plusieurs tailles d'echantillon.
fonction_simul_doses_mean<-function(vector_size,vecteur_parametres,K){
  vector_size<-vector_size[order(vector_size)]
  nb_doses<-length(vecteur_parametres)
  liste_gg<-list(rep(NA,nb_doses))
  t_star<-vecteur_parametres[[1]][["t_star"]]
  resultat_all_sizes<-lapply(vector_size,calcul_mean_size_Ktimes,nb_doses=nb_doses,K=K,vecteur_param=vecteur_parametres,t_star=t_star)
  result_by_dose<-list(c(1:nb_doses))
  for (j in c(1:nb_doses)){
    result_by_dose[[j]]<-t(cbind(sapply(resultat_all_sizes,function(x,indice){return(x[indice,])},indice=j)))
    result_by_dose[[j]]<-as.data.frame(result_by_dose[[j]])
  }
  return(result_by_dose)
}
#plot des résultats.
print_mean_mult_doses<-function(N,liste_parameter,limit_inf,limit_sup)
{
  require(gridExtra)
  require(ggplot2)
  vector_size<-seq.int(limit_inf,limit_sup,10)
  vector_size<-vector_size[order(vector_size)]
  nombre_doses<- length(liste_parameter)
  MEAN<-fonction_simul_doses_mean(vector_size=vector_size,
                                  vecteur_parametres=liste_parameter,K=N)
  vecteur_gg<-rep(NA,nombre_doses)
  result<-list(rep(NA,nombre_doses))
  for (j in c(1:nombre_doses)){
    data<-MEAN[[j]]
    data <- as.data.frame(lapply(data, unlist))
    minimum<-min(min(data$mean_guerison),min(data$mean_Bernoulli),min(data$mean_survie))
    maximum<-max(max(data$mean_guerison),max(data$mean_Bernoulli),max(data$mean_survie))
    k<-liste_parameter[[j]][["k"]]
    lambda<-liste_parameter[[j]][["lambda"]]
    p<-liste_parameter[[j]][["p"]]
    gg1<-ggplot(data=data,aes(x=n,y=mean_guerison,col="guerison"))+
      geom_line()+
      ggtitle(paste("Evolution du biais moyen en \n fonction de n, dose",as.character(j)))+
      geom_line(data=data,aes(x=n,y=mean_Bernoulli,col="Bernoulli"))+
      scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "Bernoulli" = "blue")) +
      ylim(minimum,maximum)+
      xlab("Taille echantillon") + ylab("Moyenne du biais")

    gg2<-ggplot(data=data,aes(x=n,y=mean_guerison,col="guerison"))+
      geom_line()+
      ggtitle(paste("Evolution du biais moyen en \n fonction de n, dose",as.character(j)))+
      geom_line(data=data,aes(x=n,y=mean_survie,col="survie"))+
      scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "survie" = "darkgreen")) +
      ylim(minimum,maximum)+
      xlab("Taille echantillon") + ylab("Moyenne du biais")+
      labs(caption = sprintf("lambda = %s, alpha= %s, p=%s,N=%s" ,
                             as.character(round(lambda,2)),
                             as.character(k),
                             as.character(p),
                             as.character(N)))
    gg <- grid.arrange(gg1, gg2, ncol = 2, widths = c(7,7))
    result[[j]]<-gg
  }
  return(result)
}
########## EQM ###########
#calcul pour K echantillons de même taille.
#' Calcul de l'eqm pour chaque dose pour chaque estimateur pour K échantillons.
#'
#' @param size nombre d'individus considérés.
#' @param K nombre d'échantillons générés.
#' @param nb_doses Nombre de doses. Ce nombre correspond à la longueur de liste_params.
#' @param vecteur_param liste contenant plusieurs sous-listes. Chaque sous-liste contient le lambda et le k de chaque dose.
#' @param t_star fin de la fenetre d'observation
#' @return Liste. Valeur du biais pour chaque échantillon pour chaque estimateur.
#' @export
#'
#' @examples
#' ######Test ######
#' K<-10
#' n<-100
#' k1<-1
#' lambda1<-3
#' sous_liste1<-list(k1,lambda1,0.33)
#' names(sous_liste1)<-c("k","lambda","p")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2,0.5)
#' names(sous_liste2)<-c("k","lambda","p")
#' ls<-list(sous_liste1,sous_liste2)
#' eqm<-calcul_eqm_size_Ktimes(K=K,size=n,vecteur_param=ls,nb_doses=2,t_star=6)
calcul_eqm_size_Ktimes<-function(size,vecteur_param,nb_doses,K,t_star){

  matrice_eqm_doses<-as.data.frame(matrix(NA,nb_doses,5))
  colnames(matrice_eqm_doses)<-c("eqm_Bernoulli","eqm_guerison","eqm_survie","p","n")
  result<-lapply(rep(size,K),function_estim_doses,liste_params=vecteur_param,nb_doses=nb_doses,t_star=t_star)
  for(j in c(1:nb_doses)){
    ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
    ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
    p<-vecteur_param[[j]][["p"]]
    ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
    ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
    ensemble_obs_dosek$estimateur_modele_survie<-as.numeric(ensemble_obs_dosek$estimateur_survie)
    ensemble_obs_dosek$p<-as.numeric(ensemble_obs_dosek$p)
    eqm_bern<-mean((ensemble_obs_dosek$estimateur_bernoulli-p)^(2))
    eqm_surv<-mean((ensemble_obs_dosek$estimateur_modele_survie-p)^(2))
    eqm_cure<-mean((ensemble_obs_dosek$estimateur_guerison-p)^(2))
    matrice_eqm_doses[j,c("eqm_Bernoulli","eqm_guerison","eqm_survie","p","n")]<-c(eqm_bern,
                                                                                   eqm_cure,
                                                                                   eqm_surv,
                                                                                   mean(ensemble_obs_dosek$p),
                                                                                   size)
  }
  return(matrice_eqm_doses)
}
#### calcul pour plusieurs tailles.######
fonction_simul_doses_eqm<-function(vector_size,vecteur_parametres,K){
  vector_size<-vector_size[order(vector_size)]
  nb_doses<-length(vecteur_parametres)
  liste_gg<-list(rep(NA,nb_doses))
  t_star<-vecteur_parametres[[1]][["t_star"]]
  result_by_dose<-list(c(1:nb_doses))
  resultat_all_sizes<-lapply(vector_size,calcul_eqm_size_Ktimes,nb_doses=nb_doses,K=K,vecteur_param=vecteur_parametres,t_star=t_star)
  for (j in c(1:nb_doses)){
    result_by_dose[[j]]<-t(cbind(sapply(resultat_all_sizes,function(x,indice){return(x[indice,])},indice=j)))
    result_by_dose[[j]]<-as.data.frame(result_by_dose[[j]])
  }
  return(result_by_dose)
}
### plot du résultat.#####
print_eqm_mult_doses<-function(N,liste_parameter,limit_inf,limit_sup,nombre_doses)
{
  require(gridExtra)
  vector_size<-seq.int(limit_inf,limit_sup,10)
  vector_size<-vector_size[order(vector_size)]
  EQM<-fonction_simul_doses_eqm(vector_size=vector_size,
                                vecteur_parametres=liste_parameter,K=N)
  vecteur_gg<-rep(NA,nombre_doses)
  result<-list(rep(NA,nombre_doses))
  for (j in c(1:nombre_doses)){
    data<-as.data.frame(EQM[[j]])
    data <- as.data.frame(lapply(data, unlist))
    minimum<-min(min(data$eqm_guerison),min(data$eqm_Bernoulli),min(data$eqm_survie))
    maximum<-max(max(data$eqm_guerison),max(data$eqm_Bernoulli),max(data$eqm_survie))
    k<-liste_parameter[[j]][["k"]]
    lambda<-liste_parameter[[j]][["lambda"]]
    p<-liste_parameter[[j]][["p"]]
    gg1<-ggplot(data=data,aes(x=n,y=eqm_guerison,col="guerison"))+
      geom_line()+
      ggtitle(paste("Evolution de l'EQM, dose",as.character(j)))+
      geom_line(data=data,aes(x=n,y=eqm_Bernoulli,col="Bernoulli"))+
      scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "Bernoulli" = "blue")) +
      ylim(minimum,maximum)+
      xlab("Taille echantillon") + ylab("EQM")

    gg2<-ggplot(data=data,aes(x=n,y=eqm_guerison,col="guerison"))+
      geom_line()+
      geom_line(data=data,aes(x=n,y=eqm_survie,col="survie"))+
      ggtitle(paste("Evolution de l'EQM, dose",as.character(j)))+
      scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "survie" = "darkgreen")) +
      ylim(minimum,maximum)+
      xlab("Taille echantillon") + ylab("EQM")+
      labs(caption = sprintf("lambda = %s, alpha = %s, p=%s,N=%s" ,
                             as.character(round(lambda,2)),
                             as.character(k),
                             as.character(p),
                             as.character(N)))
    gg <- grid.arrange(gg1, gg2, ncol = 2, widths = c(7,7))
    result[[j]]<-gg
  }
  return(result)
}



