#' Renvoie une simulation d'une Weibull de parametres k et lambda.
#'
#' @param n Nombre d'individus a generer.
#' @param lambda parametre correspondant a l'echelle (scale) de la loi Weibull generee.
#' @param k parametre correspondant a la forme (shape) de la loi Weibull generee.
#'
#' @return Renvoie une simulation d'une Weibull de parametres k et lambda.
#' @export
#'
#' @examples
#' test_simul<-simul_weibull(100,0.1,1)
simul_weibull<-function(n,lambda,k){
 return(stats::rweibull(n,shape=k,scale=1/lambda))
}
#' Trouver le quantile a t_star de la fonction de survie (Kaplan-Meier).
#'
#' @param n Taille de l'echantillon.
#' @param lambda parametre de la fonction Weibull (scale).
#' @param k parametre de la fonction Weibull (shape).
#' @param t_star fin de la fenetre d'observation.
#'
#' @return Trouver le quantile a t_star de la fonction de survie (Kaplan-Meier).
#' @export
#'
#' @examples
#' lambda_test<-1/3
#' k<-3
#' n<-100
#' t_star<-6
simul_survie_weibull<-function(n,lambda,k,t_star){
  donnees<-simul_weibull(n,lambda,k)
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
  # on touche la proportion de tstar au premier NA
  while (is.na(individu)==FALSE) {
    individu <- cent[m]
    m<- m+1
  }
  transformation <- m-1
  transformation <- transformation / 100
  return(transformation)
}
#' Calculer le biais d'une Weibull k,lambda.
#'
#' @param n taille d'échantillon
#' @param lambda paramètre d'échelle
#' @param k paramètre de forme
#' @param t_star fin de la fenêtre d'observation
#'
#' @return
#' @export
#'
#' @examples
#' N<-20
#' lambda_test<-1/3
#' k<-3
#' n<-100
#' t_star<-6
#' test_biais_weibull<-fonction_biais_survie_weibull(n,lambda=lambda_test,k,t_star)
fonction_biais_survie_weibull<-function(n,lambda,k,t_star){
  estimateur<-simul_survie_weibull(n,lambda,k,t_star)
  valeur_theorique<-pweibull(t_star,scale=lambda,shape=k)
  return(estimateur-valeur_theorique)
}
#' Generer N fois un echantillon de taille n et trouver le biais.
#'
#' @param N nombre de simulation
#' @param n taille de l'echantillon
#' @param lambda parametre de la Weibull, scale.
#' @param k parametre de la weibull, shape.
#' @param t_star fin de la fenetre d'observation.
#'
#' @return valeur du biais pour les N simulations.
#' @export
#'
#' @examples
#' lambda_test<-1/3
#' N<-20
#' k<-3
#' n<-100
#' t_star<-6
#' vecteur<-Simuler_Nfois_n_weibull(N,n,lambda=lambda_test,k,t_star)
Simuler_Nfois_n_weibull<-function(N,n,lambda,k,t_star){
  vecteur_biais<-rep(NA,N)
  vecteur_taille<-rep(n,N)
  vecteur_biais<-sapply(vecteur_taille,fonction_biais_survie_weibull,lambda=lambda,k=k,t_star=t_star)
  return(vecteur_biais)
}
#' Generer le biais en fonction d'un nombre de k (shape) differents.
#'
#' @param n taille de l'echantillon.
#' @param lim_moins valeur minimale de la fenetre des k.
#' @param lim_plus valeur maximale de la fenetre des k.
#' @param lambda valeur de la loi Weibull (scale).
#' @param t_star fin de la fenetre d'observation.
#' @param number_k Nombre de k differents etudies.
#'
#' @return Valeur du biais pour tous les k differents etudies.
#' @export
#'
#' @examples
#' number_trials<-10
#' n <-100
#' l_plus<-5
#' l_moins<-0.1
#' t_star <- 6
#' lambda_test <- 1
#' vecteur_k_bias<-function_influence_rate(n,lim_moin=l_moins,lim_plus=l_plus,
#' lambda=lambda_test,number_k=number_trials,t_star=t_star)
function_influence_rate<-function(n,lim_moins,lim_plus,lambda,t_star,number_k){

  vector_k<-as.vector(seq.int(lim_moins,lim_plus,length.out = number_k))
  curve_k<-sapply(vector_k,fonction_biais_survie_weibull,lambda=lambda,t_star=t_star,n=n)
  return(curve_k)
}
