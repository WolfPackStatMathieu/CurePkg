
# source("surv.R")
# source("bernoulli.R")
# library("ggplot2")

#' Calcule le biais à partir d'un vecteur de lambdas
#'
#' @param lambda_vec : vecteur de lambdas
#' @param list_params : liste de paramètres nécessaires (ici n : taille de l'échantillon, t_star : fenêtre d'obs)
#'
#' @return retourne un vecteur de biais
#' @export
#'
#' @examples
#' n<-100
#' t_star<-6
#' lambda_vec <- c(0.8,0.41, 0.28, 0.25,0.42, 0.04)
#' liste_parameter<-list(n,t_star)
#' names(liste_parameter)<-c("n","t_star")
#' test <-lambda_func(lambda_vec=lambda_vec,list_params =  liste_parameter)
#'
#' #plot du résultat
#' #----------------
#' ggplot2::ggplot(data.frame(x = lambda_vec, y = test), ggplot2::aes(x = x, y = y)) + ggplot2::geom_line(color = "blue", size = 1.5) + ggplot2::ggtitle("Bias with different lambdas") + ggplot2::xlab("lambda") + ggplot2::ylab("Biais") + ggplot2::theme_minimal()
lambda_func <- function(lambda_vec, list_params){
  vec_biais_surv <- sapply(lambda_vec,fonction_biais_survie,n=list_params[["n"]],t_star=list_params[["t_star"]])
  return(vec_biais_surv)
}


#' Génère le biais à partir d'un vecteur de probabilité p
#'
#' @param p_vec : vecteur de probabilités
#' @param n : taille de l'échantillon
#'
#' @return retourne un vecteur de biais
#' @export
#'
#' @examples
#' vec_p <- c(0.2,0.33,0.24,0.12,0.18)
#' n <- 100
#' test1 <- pbinom_func(vec_p, n)
#'
#' #plot du résultat
#' #----------------
#' ggplot2::ggplot(data.frame(x = vec_p, y = test1), ggplot2::aes(x = x, y = y)) + ggplot2::geom_line(color = "blue", size = 1.5) + ggplot2::ggtitle("Bias with different probabilities") + ggplot2::xlab("p") + ggplot2::ylab("Bias") + ggplot2::theme_minimal()
pbinom_func <- function(p_vec, n){
  vec_biais_binom <- sapply(p_vec,biais_pi,n=n)
  return(vec_biais_binom)
}




