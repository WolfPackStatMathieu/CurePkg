######## premiere methode de generation ########
#' Calcul du biais pour chaque dose pour chaque estimateur avec un seul échantillon.
#'
#' @param n nombre d'individus considérés.
#' @param nb_doses Nombre de doses. Ce nombre correspond à la longueur de liste_params.
#' @param liste_params liste contenant plusieurs sous-listes. Chaque sous-liste contient le lambda et le k de chaque dose.
#' @param t_star fin de la fenetre d'observation
#' @param k paramètre de la loi Weibull (shape).
#' @return Ggplot. Valeur du biais moyen pour les trois estimateurs pour n allant de 20 à 100.
#' @export
#'
#' @examples
#' ######Test ######
#' K<-10
#' n<-100
#' k1<-1
#' lambda1<-3
#' sous_liste1<-list(k1,lambda1)
#' names(sous_liste1)<-c("k","lambda")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2)
#' names(sous_liste2)<-c("k","lambda")
#' ls<-list(sous_liste1,sous_liste2)
#' result<-function_estim_doses(n=100,liste_params=ls,nb_doses=len(ls),t_star=6)
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
  fonction_surv<-Surv(as.numeric(df$tox_time),event=df$is_observed)
  indice_cens<-which(df$is_observed==0)
  df$factdose<-as.factor(df$dose)
  if(length(indice_cens)==0){

    estimateur_surv<-rep(1,nb_doses)
    Prob_whole_cure<-fit.cure.model(Surv(tox_time,is_observed) ~ factdose, formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimateur_surv,estimateur_cure)
  }
  if(length(indice_cens)==nrow(df)){
    estimateur_surv<-rep(0,nb_doses)
    Prob_whole_cure<-fit.cure.model(Surv(tox_time,is_observed) ~ factdose,  formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
    beta0<-as.numeric(Prob_whole_cure$coefs[1]$'1')[1]
    reste_beta<-as.numeric(Prob_whole_cure$coefs[1]$'1')[c(2:nb_doses)]
    coeffs<-beta0+c(0,reste_beta)
    estimation_cure<-1-plogis(coeffs)
    estimation_surv<-rep(0,nb_doses)
    data_returns[,c("estimateur_survie","estimateur_guerison")]<-c(estimateur_surv,estimateur_cure)
  }
  else{
    fit_surv <- survfit(fonction_surv ~factdose, data = df)
    Prob_whole_cure<-fit.cure.model(Surv(tox_time,is_observed) ~ factdose,  formula.surv=list(~1,~factdose),data =df,dist="weibull",link="logit")
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
#' sous_liste1<-list(k1,lambda1)
#' names(sous_liste1)<-c("k","lambda")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2)
#' names(sous_liste2)<-c("k","lambda")
#' ls<-list(sous_liste1,sous_liste2)
#' biais_K<-Realisations_estim_cas_mult(K=K,n=n,liste_params=ls,nb_doses=2,t_star=6)
Realisations_estim_cas_mult<-function(K,n,liste_params,nb_doses,t_star){
  ### Genere la moyenne des estimateurs pour la taille n
  result<-lapply(rep(n,K),function_estim_doses,liste_params=liste_params,nb_doses=nb_doses,t_star=t_star)
  matrice<-list(rep(NA,nb_doses))
  for(j in c(1:nb_doses)){
    ensemble_obs_dosek<-t(cbind.data.frame(sapply(result,function(x,indice){return(x[indice,])},indice=j)))
    ensemble_obs_dosek<-as.data.frame(ensemble_obs_dosek)
    ensemble_obs_dosek$estimateur_bernoulli<-as.numeric(ensemble_obs_dosek$estimateur_bernoulli)
    ensemble_obs_dosek$estimateur_guerison<-as.numeric(ensemble_obs_dosek$estimateur_guerison)
    ensemble_obs_dosek$estimateur_modele_survie<-as.numeric(ensemble_obs_dosek$estimateur_survie)
    ensemble_obs_dosek$p<-as.numeric(ensemble_obs_dosek$p)
    matrice[[j]]<-ensemble_obs_dosek}
  return(matrice)
}
#' Calcul de l'eqm pour chaque dose pour chaque estimateur pour K échantillons.
#'
#' @param size nombre d'individus considérés.
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
#' sous_liste1<-list(k1,lambda1)
#' names(sous_liste1)<-c("k","lambda")
#' k2<-2
#' lambda2<-0.5
#' sous_liste2<-list(k2,lambda2)
#' names(sous_liste2)<-c("k","lambda")
#' ls<-list(sous_liste1,sous_liste2)
#' eqm<-Realisations_estim_cas_mult(K=K,size=n,vecteur_param=ls,nb_doses=2,t_star=6)
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
