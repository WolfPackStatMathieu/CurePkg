
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
  surv_object<-Surv(donnees_ensemble$tox_time,event=donnees_ensemble$isobserved)
  fit <- survfit(surv_object ~1, data = donnees_ensemble)
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

