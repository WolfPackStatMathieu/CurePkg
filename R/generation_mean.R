# source("surv.R")
# source("bernoulli.R")
#' G√©n√©re n fois le mod√®le
#'
#' @param vector_size : un vecteur de taille d'√©chantillon
#' @param model : mod√®le bernouilli ou surv
#' @param liste_parameter : liste des param√®tres
#' @param K : nombre de fois que l'on g√©n√®re un √©chantillon de taille n
#'
#' @return : G√©n√©re n fois le mod√®le (bernouilli ou survie)
#' @export
#'
#' @examples
#' #test avec le mod√®le "survie"
#' #---------------------------
#' N<-10
#' vecteur_size<-sample(c(1:100),N)
#' lamdba_test<-3
#' t_star<-6
#' p<-0.33
#' k<-1
#' liste_parameter<-list(lambda_test,t_star,p=p,k)
#' names(liste_parameter)<-c("lambda","t_star","p","k")
#' K2<-20
#' test_exp_taillemoy<-fonction_generation_taille_mean(vector_size=vecteur_size,K=K2,liste_parameter = liste_parameter)
#' #################### Plot des rÈsultats en fonction de la taille.########
#' donnees_taille_biaismoyen<-cbind.data.frame(vecteur_size[order(vecteur_size)],test_exp_taillemoy)
#' colnames(donnees_taille_biaismoyen)<-c("Size","Mean_Bias_Cure","Mean_Bias_Surv")
#' plot(donnees_taille_biaismoyen$Size,donnees_taille_biaismoyen$Mean_Bias_Cure,main="The mean bias according to the size with Survival function")
#' points(x=donnees_taille_biaismoyen$Size,y=donnees_taille_biaismoyen$Mean_Bias_Surv,col="red")
fonction_generation_taille_mean<-function(vector_size,liste_parameter,K){
  ### renvoie la gÈnÈration avec des tailles diffÈrentes du modËle model (string) ayant comme paramËtre la liste_parameter, 
  ### liste de paramËtres avec le modËles. 1 seul comme bernoulli et 2 pour exp() (lambda et t_star).
  vector_size<-vector_size[order(vector_size)]
  ##### idÈe. 
  Value_bias<-lapply(vector_size,Simuler_biais_taillen,K=K,lambda=liste_parameter[['lambda']],t_star=liste_parameter[["t_star"]],
                     p=liste_parameter[["p"]],k=liste_parameter[["k"]])
  value_means<-as.data.frame(t(sapply(Value_bias,colMeans)))
  value_means$Modele_survie<-value_means$Modele_survie-p
  value_means$Modele_guerison<-value_means$Modele_guerison-p
  return(value_means)
}


