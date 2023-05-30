
#' Simuler un échantillon de taille n selon la loi exponentiel de paramètre lambda.
#'
#' @param n : taille de l'échantillon
#' @param lambda : paramètre de l'exponentiel.
#'
#' @return Un échantillon de taille n.
#' @export
#'
#' @examples
#' simul_exp(10,0.3)
simul_exp<-function(n,lambda){
  ###generer un echantillon de taille n suivant une loi exponentielle de param?tre lambda.
  return(rexp(n,lambda))
}
#' Retourner la probabilité de manifester la toxicité de 0 à t_star.
#'
#' @param n : taille de l'échantillon
#' @param lambda : paramètre de la loi exponentielle.
#' @param t_star : int, fin de la fenêtre d'observation.
#'
#' @return le quantile de la fonction de survie en t_star.
#' @export
#'
#' @examples
#'set.seed(133)
#'n<-100
#'test_surv<-simul_survie(n,0.5,6)
simul_survie<-function(n,lambda,t_star){
  #### Calculer la probabilite que la toxicite apparaisse de 0 ? t_star.
  ### Se base sur la simulation de temps via la fonction simul_exp.
  ### KAPLAN-MEIER.
  donnees<-simul_exp(n,lambda)
  donnees_censure_tstar<-ifelse(donnees<t_star,donnees,t_star)
  donnees_indicatrice_observee<-ifelse(donnees<t_star,1,0)
  donnees_ensemble<-cbind.data.frame(donnees_censure_tstar,donnees_indicatrice_observee)
  colnames(donnees_ensemble)<-c("tox_time","isobserved")
  surv_object<-survival::Surv(donnees_ensemble$tox_time,event=donnees_ensemble$isobserved)
  fit <- survival::survfit(surv_object ~1, data = donnees_ensemble)
  # on cherche a recuperer les donnees au temps T=6
  #afin de pouvoir tracer la droite Toxicite =f(dose)
  quantile <-quantile(fit)
  quantile$quantile
  centiles <- quantile(fit, 1:100/100)
  cent <-centiles$quantile
  m<-1
  individu <- cent[m]
  # on touche la proportion de tstar au premier NA.
  while (is.na(individu)==FALSE) {
    individu <- cent[m]
    m<- m+1
  }
  transformation <- m-1
  transformation <- transformation / 100
  return(transformation)
}
#' Generer N fois un echantillon de taille n
#'
#' @param N nombre de fois generation echantillon
#' @param n taille echantillonn
#' @param lambda parametre loi exponentielle.
#' @param t_star fin de la fenetre.
#'
#' @return Generer N fois un echantillon de taille n
#' @export
#'
#' @examples
#' lambda_test<-1/3
#' n<-100
#' t_star<-6
#' N<-100
#' test_simul_total<-Simuler_Nfois_n_echantillons(N,n,lambda_test,t_star)

Simuler_Nfois_n_echantillons<-function(N,n,lambda,t_star){
  vecteur_biais<-rep(NA,N)
  vecteur_taille<-rep(n,N)
  vecteur_biais<-sapply(vecteur_taille,fonction_biais_survie,lambda=lambda,t_star=t_star)
  return(vecteur_biais)
}
#' Biais.
#'
#' @param n taille d'échantillon
#' @param lambda paramètre d'échelle de la loi de weibull
#' @param t_star fin de la fenêtre d'observation
#'
#' @return Calcul du biais de la probabilite de toxicite estimee par la fonction simul_survie.
#' @export
#'
#' @examples
fonction_biais_survie<-function(n,lambda,t_star){
  #### Calcul du biais de la probabilite de toxicite estimee par la fonction simul_survie.
  ### Comparaison avec la fonction de repartition d'une exp(lambda) en t_star.
  estimateur<-simul_survie(n,lambda,t_star)
  valeur_theorique<-pexp(t_star,rate=lambda)
  return(estimateur-valeur_theorique)
}
#' Generer les estimateurs pour un echantillon de taille n de la quantite p.
#'
#' @param n : taille de l'echantillon
#' @param k : parametre de la loi weibull.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : fin de la fenetre d'observation
#' @param p : la proportion de non gueris.
#' @return un vecteur des estimateurs.
#' @export
#'
#' @examples
#' test_retour<-Simuler_biais_un_n_ech(n=10,lambda=0.5,t_star=6,0.33,2)
Simuler_biais_un_n_ech<-function(n,lambda,t_star,p,k){
  vecteur_censure<-simul_bernoulli(n,p)
  vecteur_temp<-rep(NA,n)
  estimateur_cure<-mean(vecteur_censure)
  df<-cbind.data.frame(vecteur_censure,vecteur_temp)
  colnames(df)<-c("censure","temps")
  id_censures<-which(df$censure==1)
  id_obs<-which(df$censure==0)
  df[id_censures,2]<-t_star
  df[id_obs,2]<-simul_weibull(length(id_obs),lambda,k)
  df$temps<-ifelse(df$temps<t_star,df$temps,t_star)
  colnames(df)<-c("isobserved","tox_time")
  df$isobserved<-ifelse(df$tox_time<t_star,1,0)
  surv_object<-survival::Surv(df$tox_time,event=df$isobserved)
  fit <- survival::survfit(surv_object ~1, data = df)
  # on cherche a recuperer les donnees au temps T=6
  #afin de pouvoir tracer la droite Toxicite =f(dose)
  quantile <-quantile(fit)
  quantile$quantile
  centiles <- quantile(fit, 1:100/100)
  cent <-centiles$quantile
  m<-1
  individu <- cent[m]
  # on touche la proportion de tstar au premier NA
  while (is.na(individu)==FALSE) {
    individu <- cent[m]
    m<- m+1
  }
  transformation <- m-1
  estimateur_survie <- transformation / 100
  liste_biais<-list(estimateur_cure,estimateur_survie)
  names(liste_biais)<-c("Modele_guerison","Modele_survie")
  return(liste_biais)
}
#' Generer les estimateurs pour plusieurs echantillons de taille n de la quantite p.
#'
#' @param n : taille de l'echantillon
#' @param k : parametre de la loi weibull.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : fin de la fenetre d'observation
#' @param p : la proportion de non gueris.
#' @param K : nombre de simulations des echantillons.
#' @return un vecteur des estimateurs.
#' @export
#'
#' @examples
#' test_several_times<-Simuler_biais_taillen(n=10,lambda=0.5,t_star=6,p=0.33,k=2,K=10)
Simuler_biais_taillen<-function(K,n,lambda,t_star,p,k){
  df_biases<-as.data.frame(t(cbind.data.frame(sapply(rep(n,K),Simuler_biais_un_n_ech,lambda=lambda,t_star=t_star,p=p,k=k))))
  df_biases$Modele_guerison<-as.numeric(df_biases$Modele_guerison)
  df_biases$Modele_survie<-as.numeric(df_biases$Modele_survie)
  return(df_biases)
}
#' Calculer le biais pour les deux modeles par rapport a p.
#'
#' @param n : taille de l'echantillon
#' @param k : parametre de la loi weibull.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : fin de la fenetre d'observation
#' @param p : la proportion de non gueris.
#' @param K : nombre de simulations des echantillons.
#' @return un vecteur des biais.
#' @export
#'
#' @examples
#' test_biais_moy<-Calcul_biais_moyen_taillen(n=10,lambda=0.5,t_star=6,p=0.33,k=2,K=10)
Calcul_biais_moyen_taillen<-function(K,n,lambda,t_star,p,k){
  data<-Simuler_biais_taillen(K,n,lambda,t_star,p,k)
  result<-rep(NA,2)
  result<-colMeans(data)
  result[1]<-result[1]-p
  result[2]<-result[2]-p
  return(result)
}


