# source("surv.R")
# source("bernoulli.R")
#' Génére n fois le modèle
#'
#' @param vector_size : un vecteur de taille d'échantillon
#' @param model : modèle bernouilli ou surv
#' @param liste_parameter : liste des paramètres
#' @param K : nombre de fois que l'on génère un échantillon de taille n
#'
#' @return : Génére n fois le modèle (bernouilli ou survie)
#' @export
#'
#' @examples
#' test avec le modèle "survie"
#' ---------------------------
#' N<-100
#' vecteur_size<-sample(c(1:1000),N)
#' lamdba_test<-0.33
#' t_star<-6
#' liste_parameter<-list(lambda_test,t_star)
#' names(liste_parameter)<-c("lambda","t_star")
#' modele<-"surv"
#' k<-20
#' res <- fonction_generation_taille_mean(vector_size=vecteur_size,model=modele,K=k,liste_parameter = liste_parameter)
#'
#' plot des résultats
#' ------------------
#' donnees_taille_biaismoyen<-cbind.data.frame(vecteur_size[order(vecteur_size)],res)
#' colnames(donnees_taille_biaismoyen)<-c("Size","Mean_Bias")
#' plot(donnees_taille_biaismoyen,main="The mean bias according to the size with Survival function")
#' lines(x=donnees_taille_biaismoyen$Size,y=estimation_ymoy,col="red")
#'
#' test avec le modèle bernouilli
#' ------------------------------
#' prop<-0.33
#' liste_param<-list(prop)
#' names(liste_param)<-c("p")
#' modele1<-"bernoulli"
#' test_bern_taillemoy<-fonction_generation_taille_mean(vector_size=vecteur_size,model=modele1,K=k,liste_parameter = liste_param)
#'
#' plot des resultats
#' ------------------
#' donnees_bern_biaismoyen<-cbind.data.frame(vecteur_size[order(vecteur_size)],test_bern_taillemoy)
#' plot(donnees_bern_biaismoyen,main="The mean bias according to the size with Bernoulli")
#' lines(x=donnees_bern_biaismoyen$Size,y=estimation_ymoy,col="blue")
#'
fonction_generation_taille_mean<-function(vector_size,model,liste_parameter,K){
  vector_size<-vector_size[order(vector_size)]
  if (model=="bernoulli"){
    vecteur_realisation<-sapply(vector_size,Simuler_Nfois_n_echantillons_bern,N=K,p=liste_parameter[["p"]])
    return(colMeans(vecteur_realisation))
  }
  else{if(model=="surv"){
    vecteur_realisation<-sapply(vector_size,Simuler_Nfois_n_echantillons,N=K,lambda=liste_parameter[["lambda"]],t_star=liste_parameter[["t_star"]])
    return(colMeans(vecteur_realisation))
  }
    if(model=="weibull"){
      vecteur_realisation<-sapply(vector_size,Simuler_Nfois_n_weibull,N=K,lambda=liste_parameter[["lambda"]],t_star=liste_parameter[["t_star"]],
                                  k=liste_parameter[["k"]])
      return(colMeans(vecteur_realisation))
    }}
}



