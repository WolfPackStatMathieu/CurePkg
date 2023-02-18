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
  return(rweibull(n,shape=k,scale=lambda))
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
#' test_biais_weibull<-fonction_biais_survie_weibull(n,lambda=lambda_test,k,t_star)
simul_survie_weibull<-function(n,lambda,k,t_star){

  donnees<-simul_weibull(n,lambda,k)
  donnees_censure_tstar<-ifelse(donnees<t_star,donnees,t_star)
  donnees_indicatrice_observee<-ifelse(donnees<t_star,1,0)
  donnees_ensemble<-cbind.data.frame(donnees_censure_tstar,donnees_indicatrice_observee)
  colnames(donnees_ensemble)<-c("tox_time","isobserved")
  surv_object<-Surv(donnees_ensemble$tox_time,event=donnees_ensemble$isobserved)
  fit <- survfit(surv_object ~1, data = donnees_ensemble)
  summary(fit)
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
