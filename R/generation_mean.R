# source("surv.R")
# source("bernoulli.R")
#' Génére n fois le modèle
#'
#' @param vector_size : un vecteur de taille d'échantillon
#' @param liste_parameter : liste des paramètres
#' @param K : nombre de fois que l'on génère un échantillon de taille n
#'
#' @return : Génére n fois le modèle (bernouilli ou survie)
#' @export
#'
#' @examples
#' #test avec le modèle "survie"
#' #---------------------------
#' N<-10
#' vecteur_size<-sample(c(1:100),N)
#' lambda_test<-3
#' t_star<-6
#' p<-0.33
#' k<-1
#' liste_parameter<-list(lambda_test,t_star,p=p,k)
#' names(liste_parameter)<-c("lambda","t_star","p","k")
#' K2<-20
#' test_exp_taillemoy<-fonction_generation_taille_mean(vector_size=vecteur_size,K=K2,liste_parameter = liste_parameter)
#' #################### Plot des r?sultats en fonction de la taille.########
#' donnees_taille_biaismoyen<-cbind.data.frame(vecteur_size[order(vecteur_size)],test_exp_taillemoy)
#' colnames(donnees_taille_biaismoyen)<-c("Size","Mean_Bias_Cure","Mean_Bias_Surv")
#' plot(donnees_taille_biaismoyen$Size,donnees_taille_biaismoyen$Mean_Bias_Cure,main="The mean bias according to the size with Survival function")
#' points(x=donnees_taille_biaismoyen$Size,y=donnees_taille_biaismoyen$Mean_Bias_Surv,col="red")
fonction_generation_taille_mean<-function(vector_size,liste_parameter,K){
  vector_size<-vector_size[order(vector_size)]
  Value_bias<-lapply(vector_size,Simuler_biais_taillen,K=K,lambda=liste_parameter[['lambda']],t_star=liste_parameter[["t_star"]],
                     p=liste_parameter[["p"]],k=liste_parameter[["k"]])
  value_means<-as.data.frame(t(sapply(Value_bias,colMeans)))
  value_means$Modele_survie<-value_means$Modele_survie-liste_parameter[["p"]]
  value_means$Modele_guerison<-value_means$Modele_guerison-liste_parameter[["p"]]
  return(value_means)
}

fonction_sapply<-function(x){
  return(sapply(x,var))
}
#' Genere l'ecart_quadratique_moyen pour un ensemble de taille.
#'
#' @param vector_size : un vecteur de taille d'échantillon
#' @param liste_parameter : liste des paramètres
#' @param K : nombre de fois que l'on genere un échantillon de taille n
#'
#' @return : Génére n fois le modèle (bernouilli ou survie)
#' @export
#'
#' @examples
#' N<-10
#' vecteur_size<-sample(c(1:100),N)
#' lambda_test<-3
#' t_star<-6
#' p<-0.33
#' k<-1
#' liste_parameter<-list(lambda_test,t_star,p=p,k)
#' names(liste_parameter)<-c("lambda","t_star","p","k")
#' K2<-20
#' test_exp_eqm<-fonction_generation_eqm(vector_size=vecteur_size,K=K2,liste_parameter = liste_parameter)

fonction_generation_eqm<-function(vector_size,liste_parameter,K){
  ### renvoie la generation avec des tailles differentes avec un lambda,k,t_star,p.
  vector_size<-vector_size[order(vector_size)]
  Value_bias<-lapply(vector_size,Simuler_biais_taillen,K=K,lambda=liste_parameter[['lambda']],t_star=liste_parameter[["t_star"]],
                     p=liste_parameter[["p"]],k=liste_parameter[["k"]])
  value_means<-as.data.frame(t(sapply(Value_bias,colMeans)))
  value_variance<-as.data.frame(t(sapply(Value_bias,fonction_sapply)))
  value_eqm<-(value_means)^(2)+value_variance
  return(value_eqm)
}
