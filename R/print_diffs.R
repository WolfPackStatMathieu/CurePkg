#' Plot des valeurs des biais moyens selon la taille des echantillons.
#'
#' @param limit_inf taille minimale des echantillons
#' @param limit_sup taille maximale des echantillons
#' @param N nombre de tailles d'echantillon differents.
#' @param p probabilite selon le modele de guerison de toxicite.
#' @param lambda parametre de la fonction de survie.
#' @param t_star fin de la fenetre d'observation
#' @param K nombre de simulations de la meme taille d'echantillon.
#'
#' @return Plot des valeurs des biais moyens selon la taille des echantillons.
#' @export
#'
#' @examples
#' ######Test ######
#' p2<-0.3
#' k<-50
#' lambda7<-0.33
#' t_star<-6
#' lmoins<-1
#' l_plus<-1000
#' N<-50
#' test_plot<-fonction_compar_plots(limit_sup = l_plus,limit_inf = lmoins,N=N,p=p2,lambda=lambda7,t_star=t_star,K=k)

fonction_compar_plots<-function(limit_inf,limit_sup,N,p,lambda,t_star,K){
  vector_size<-sample(c(limit_inf:limit_sup),N)
  vector_size<-vector_size[order(vector_size)]
  liste_param1<-list(p)
  names(liste_param1)<-c("p")
  modele_bern<-"bernoulli"
  result1<-fonction_generation_taille_mean(vector_size,modele_bern,liste_param1,K)
  liste_param2<-list(lambda,t_star)
  names(liste_param2)<-c("lambda","t_star")
  modele_exp<-"surv"
  result2<-fonction_generation_taille_mean(vector_size,modele_exp,liste_param2,K)
  whole_data_expbern<-cbind.data.frame(vector_size,result1,result2)
  colnames(whole_data_expbern)<-c("Size","Mean_Bias_Bern","Mean_Bias_Surv")
  ####plot
  gg1<-ggplot2::ggplot(data=whole_data_expbern,ggplot2::aes(x=Size,y=Mean_Bias_Bern))+
    ggplot2::geom_smooth(colour="red")+
    ggplot2::labs(y="Mean Bias with Bern model")

  gg2<-ggplot2::ggplot(data=whole_data_expbern,ggplot2::aes(x=Size,y=Mean_Bias_Surv))+
    ggplot2::geom_smooth(colour="blue")+
    ggplot2::labs(y="Mean Bias with Surv model")

  whole_g<-gridExtra::grid.arrange(gg1,gg2,ncol=2,top="Comparison of the two methods")
  return(whole_g)
}
