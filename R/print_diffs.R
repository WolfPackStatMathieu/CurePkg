#' Plot des valeurs des biais moyens selon la taille des echantillons.
#'
#' @param limit_inf taille minimale des echantillons
#' @param limit_sup taille maximale des echantillons
#' @param N nombre de tailles d'echantillon differents.
#' @param p probabilite selon le modele de guerison de toxicite.
#' @param lambda parametre d'échelle (scale) de la loi weibull.
#' @param sh paramètre de la loi weibull.
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
#' shape<-2
#' test_plot<-fonction_compar_plots(limit_sup = l_plus,limit_inf = lmoins,
#' N=N,p=p2,lambda=lambda7,t_star=t_star,K=k,sh=shape)
fonction_compar_plots<-function(limit_inf,limit_sup,N,p,lambda,t_star,K,sh){
  #### N corresponds to the number of sizes. K correspond to the number of samples for each size.
  require(gridExtra)
  require(ggplot2)
  vector_size<-sample(c(limit_inf:limit_sup),N)
  vector_size<-vector_size[order(vector_size)]
  liste_param<-list(lambda,t_star,sh,p)
  names(liste_param)<-c("lambda","t_star","k","p")
  result<-fonction_generation_taille_mean(vector_size,liste_param,K)
  result$taille<-vector_size
  colnames(result)<-c("Mean_Bias_Cure","Mean_Bias_Surv","Size")
  ####plot
  borne_min<-min(min(result$Mean_Bias_Cure),min(result$Mean_Bias_Surv))
  borne_max<-max(max(result$Mean_Bias_Cure),max(result$Mean_Bias_Surv))
  gg1<-ggplot2::ggplot(data=result,ggplot2::aes(x=Size,y=Mean_Bias_Cure))+
    ggplot2::geom_point(colour="red")+
    ggplot2::labs(y="Mean Bias with cure model")+ggplot2::ylim(borne_min,borne_max)

  gg2<-ggplot2::ggplot(data=result,ggplot2::aes(x=Size,y=Mean_Bias_Surv))+
    ggplot2::geom_point(colour="blue")+
    ggplot2::labs(y="Mean Bias with Surv model")+ggplot2::ylim(borne_min,borne_max)

  whole_g<-gridExtra::grid.arrange(gg1,gg2,ncol=2,top="Comparison of the two methods")
  return(whole_g)
}
NSimulations.selon.n<-function(N,lambda,t_star){
  #' Matrice composee des biais moyens associes a la taille de l'echantillon de n=20 a n=200 par saut de 20.
  #'
  #' @param N nombre de tailles d'echantillon differents.
  #' @param lambda parametre de la loi exponentielle.
  #' @param t_star fin de la fenetre d'observation
  #'
  #' @return Valeur du biais moyen selon n dans l'intervalle (20,200).
  #' @export
  #'
  #' @examples
  #' ######Test ######
  #' t_star<-3
  #' N<-10
  #' lambda<-c(0.2)
  #' result<-NSimulations.selon.n(N,lambda,t_star)

  results<- NULL
  n<- 20
  while (n<200)
  {
    vecteur_biais<-rep(NA,N)
    biais<-  Simuler_Nfois_n_echantillons(N,n,lambda,t_star)
    results<-rbind(results,c(n,mean(biais)  ))
    n<- n+20
  }
  return(results)
}

fonction_compar_plotsn_lambda<-function(N,window_lambda,t_star){
  #' Plot des valeurs des biais moyens selon la taille des echantillons et du lambda.
  #'
  #' @param N nombre de tailles d'echantillon differents.
  #' @param window_lambda vecteur de valeurs pour lambda
  #' @param t_star fin de la fenetre d'observation
  #'
  #' @return Plot des valeurs des biais moyens en fonction du lambda et de la taille des echantillons.
  #' @export
  #'
  #' @examples
  #' ######Test ######
  #' t_star<-3
  #' N<-10
  #' window_lambda<-c(0.2,0.5,0.1)
  #' result<-fonction_compar_plotsn_lambda(N,window_lambda,t_star)

  set.seed(12345)
  RES<- NULL
  RES<- NSimulations.selon.n(N,window_lambda[1],t_star)
  RES0.2.3<-data.frame(RES)
  colnames( RES0.2.3)<- c("n","mean.bias")
  set.seed(12345)
  RES<- NULL
  RES<- NSimulations.selon.n(N,window_lambda[2],t_star)
  RES0.5.3<-data.frame(RES)
  colnames( RES0.5.3)<- c("n","mean.bias")

  set.seed(12345)
  RES<- NULL
  RES<- NSimulations.selon.n(N,window_lambda[3],t_star)
  RES0.1.3<-data.frame(RES)
  colnames( RES0.1.3)<- c("n","mean.bias")

  plot(RES0.2.3$n,RES0.2.3$mean.bias,title=paste("Influence of n"),
       ylim=c(0,0.1),type='b',bty="n",xlab="nbre sujets",ylab="biais moyen")
  title("Influence de n et lambda")
  lines(RES0.5.3$n,RES0.5.3$mean.bias,type="b",col="blue")
  lines(RES0.1.3$n,RES0.1.3$mean.bias,type="b",col="red")
  legend("topright",c("0.1","0.2","0.5"),col=c("red","black","blue"),lty=1,bty="n")
}


