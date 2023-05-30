#source("surv.R")
#source("bernoulli.R")

#' Permet de renvoyer le biais pour des tailles differentes.
#'
#' @param vector_size vecteur des différentes tailles.
#' @param model renverra a surv ou a bern. Correspondra aux simulations differentes.
#' @param liste_parameter Correspondra a p(bernoulli) ou a lambda et t_star pour exp.
#'
#' @return Vecteur des biais pour chaque taille.
#' @export
#'
#' @examples
#' #################TEST EXP #####################
#' vecteur_size<-sample(c(1:100),10)
#' lambda_test<-0.33
#' t_star<-6
#' liste_parameter<-list(lambda_test,t_star)
#' names(liste_parameter)<-c("lambda","t_star")
#' modele<-"surv"
#' test_generation_taillediff_exp<-fonction_generation_taille_differente(vector_size=vecteur_size,model=modele,liste_parameter = liste_parameter)
#'
#' #################TEST BERNOULLI ###################
#' prop<-0.33
#' list_param<-list(prop)
#' names(list_param)<-c("p")
#' modele2<-"bernoulli"
#' test_generation_taillediff_bern<-fonction_generation_taille_differente(vector_size = vecteur_size,model=modele2,liste_parameter = list_param)

fonction_generation_taille_differente<-function(vector_size,model,liste_parameter){

  if (model=="bernoulli"){
    vecteur_realisation<-sapply(vector_size[order(vector_size)],biais_pi,liste_parameter[["p"]])
    return(vecteur_realisation)
  }
  else{if(model=="surv"){
    vecteur_realisation<-sapply(vector_size[order(vector_size)],fonction_biais_survie,lambda=liste_parameter[["lambda"]],t_star=liste_parameter[["t_star"]])
    return(vecteur_realisation)
  }}
}

#' Renvoie la valeur du biais et la taille des echantillons dans la meme table.
#'
#' @param vector_size vecteur correspondant aux differentes tailles d'echantillons.
#' @param model chaine de caracteres.
#' @param liste_parameter : liste de parametres. Vaudra p (pour Bernoulli) ou lambda et t_star pour survie
#'
#' @return Renvoie la valeur du biais et la taille des echantillons dans la meme table.
#' @export
#'
#' @examples
#' ######1) exp#####
#' N<-50
#' vecteur_size<-sample(c(10:1000),N)
#' lambda_test<-0.33
#' t_star<-6
#' liste_parameter<-list(lambda_test,t_star)
#' names(liste_parameter)<-c("lambda","t_star")
#' modele<-"surv"
#' test_graph_exp<-fonction_graph_fonc_size(vector_size = vecteur_size,model=modele,liste_parameter=liste_parameter)
fonction_graph_fonc_size<-function(vector_size,model,liste_parameter){

  vector_size<-vector_size[order(vector_size)]
  results_vector<-fonction_generation_taille_differente(vector_size = vector_size,model=model,liste_parameter = liste_parameter)
  bias_size_data<-cbind.data.frame(vector_size,results_vector)
  colnames(bias_size_data)<-c("Size","Bias")
  return(bias_size_data)
}

