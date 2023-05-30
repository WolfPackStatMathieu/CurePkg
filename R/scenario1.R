#source("estimateurs.R")
#' Generer un echantillon de taille n selon le modèle de guérison à mélange.
#'
#' @param n : taille de l'echantillon permettant d'obtenir un estimateur.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : Fin de la fenetre d'observation.
#' @param k : paramètre de la loi weibull.
#'
#' @return : Un n-echantillon.
#' @export
#'
#' @examples
#' ####test#####
#' k<-2
#' lambda<-0.2
#' ech<-Generation_un_ech(n=n,lambda=lambda,t_star=6,p=0.33,k=2)
Generation_un_ech<-function(n,lambda,t_star,p,k){
  vecteur_censure<-rbinom(n,1,p)
  vecteur_temp<-rep(NA,n) # cree un vecteur des temps associ?s
  # l'estimateur du modele de Bernoulli est la moyenne des 1 du vecteur_censure
  # on cree un dataframe qui acolle la DLT au vecteur_temps
  df<-cbind.data.frame(vecteur_censure,vecteur_temp)
  #on les renomme
  colnames(df)<-c("sensible","temps")
  #recuperation des numeros de lignes des individus censures
  id_non_sensibles<-which(df$sensible==0)
  #recuperation des numeros de lignes des individus a risque de DLT
  id_sensibles<-which(df$sensible==1)
  #tous les individu censure se voient attribues comme temps la limite de la
  #fenetre d observation
  df[id_non_sensibles,2]<-t_star+1
  # les autres individus se voient attribuer un temps simule a partir d une
  # loi de Weibull (qui peut etre une loi exponentielle si k=1)
  df[id_sensibles,2]<-simul_weibull(length(id_sensibles),lambda,k)
  # bien sur, si le temps observe est superieur a la fenetre d observation, alors
  # on le remplace par la fin de fenetre d observation
  # on renomme les colonnes pour une meilleure interpretation
  # on remplit la colonne isobserved avec des 1 si on observe une toxicite avant
  #la fin de la fenetre d observation, sinon on met des 0
  colnames(df)<-c("sensible","tox_time")
  df$is_observed<-ifelse(df$tox_time<t_star,1,0)
  return(df)
}
#' Calculer la proportion de censure.
#'
#' @param n : taille de l'echantillon permettant d'obtenir un estimateur.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : Fin de la fenetre d'observation.
#' @param k : paramètre de la loi weibull.
#' @param p : proportion de non-guéris.
#' @param N : nombre de répétitions de l'expérience.
#' @return : Un violin_plot de la censure due au temps ou au fait d'être guéri.
#' @export
#'
#' @examples
#' ####test#####
#' graph<-calcule_prop_censure(N=100,n=100,lambda=0.2,t_star=6,p=0.33,k=1)
calcule_prop_censure<-function(N, n, lambda, t_star, p, k){
  #initialisation du dataframe
  result_censures<- as.data.frame(matrix(NA, nrow = length(N), 3))
  colnames(result_censures)<-c("prop_censure_totale", "prop_censure_gueris", "prop_censure_additionnelle")
  j<-1 #compteur de ligne
  for (i in 1:N){
    ech <-Generation_un_ech(n,lambda,t_star,p,k) #genere un echantillon
    #pourcentage de censure totale dans l echantillon
    result_censures[j,"prop_censure_totale"] <-as.numeric(table(ech$is_observed == 0)["TRUE"])/n
    # censure due au pourcentage de gueris, au debut de la generation de l echantillon
    result_censures[j,"prop_censure_gueris"] <- as.numeric(table(ech$sensible == 0)["TRUE"])/n
    # calcul de la censure additionnelle= censure_totale - censure_gueris
    result_censures[j,"prop_censure_additionnelle"]<-result_censures[j,"prop_censure_totale"] - result_censures[j,"prop_censure_gueris"]

    j<- j+1
  }

  # boxplot_plot_censures<-ggplot(result_censures) +
  violin_plot <- result_censures %>%
    gather(key="Type_de_censure", value="Val") %>%
    ggplot2::ggplot( aes(x=Type_de_censure, y=Val, fill=Type_de_censure) ) +
    ggplot2::geom_violin() +
    ggplot2::ggtitle("Distribution des types de censure, modele de generation 1")+
    ggplot2::ylab("Pourcentage de censure") + xlab("Type de censure")+
    ggplot2::theme(axis.text=element_text(family = "Helvetica", size=18),
          axis.title=element_text(family = "Helvetica", size=18),
          plot.title = element_text(family = "Helvetica", size = 25)) +
    ggplot2::labs(caption = sprintf("N = %s, n = %s, lambda= %s,t_star= %s, p= %s, alpha= %s", as.character(N),as.character(n), as.character(lambda), as.character(t_star),as.character(p), as.character(k)))
  print(violin_plot)
  return(result_censures)
}

#################Simulation en plusieurs fois.##########
#'Calculer les trois estimateurs pour un n-échantillon.
#'
#' @param n : taille de l'echantillon permettant d'obtenir un estimateur.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : Fin de la fenetre d'observation.
#' @param k : paramètre de la loi weibull.
#' @param p : proportion réelle de non-guéris.
#' @return : Vecteur,valeur des trois estimateurs pour un n-échantillon.
#' @export
#'
#' @examples
#' ####test#####
#' vecteur_estimateurs<-Simuler_biais_un_n_ech(n=100,lambda=2,t_tsar=6,p=0.33,k=1)
Simuler_biais_un_n_ech<-function(n,lambda,t_star,p,k){
  database<-Generation_un_ech(n=n,lambda=lambda,t_star=t_star,p=p,k=k)
  estimateur_bern<-fonction_Bern(df=database)
  estimateur_surv<-fonction_KM(df=database)
  estimateur_cure<-fonction_cure(df=database,t_star=t_star)
  # on prepare une liste avec les deux estimateurs calcules
  liste_biais<-list(estimateur_bern,estimateur_surv,estimateur_cure)
  names(liste_biais)<-c("Modele_bernoulli","Modele_survie","Modele_guerison")
  return(liste_biais)
}
#'Calculer la moyenne des trois estimateurs pour un K n-échantillon.
#'
#' @param n : taille de l'echantillon permettant d'obtenir un estimateur.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : Fin de la fenetre d'observation.
#' @param k : paramètre de la loi weibull.
#' @param p : proportion réelle de non-guéris.
#' @param K : nombre de répétitions de l'expérience.
#' @return : Vecteur. Moyenne des trois estimateurs.
#' @export
#'
#' @examples
#' ####test#####
#' vecteur_estims<-Simuler_biais_taillen(n=100,lambda=2,t_tsar=6,p=0.33,k=1,K=100)
Simuler_biais_taillen<-function(K,n,lambda,t_star,p,k){
  # Simuler_biais_un_n_ech retourne le biais du modele de guerison
  # et le biais du modele de survie
  # on cr?e un dataframe de K lignes de ces deux biais, pour des ?chantillons de
  #taille n.
  df_biases<-as.data.frame(t(cbind.data.frame(sapply(rep(n,K),Simuler_biais_un_n_ech,lambda=lambda,t_star=t_star,p=p,k=k))))
  # on donne un nom aux deux colonnes
  df_biases$Modele_bernoulli<-as.numeric(df_biases$Modele_bernoulli)
  df_biases$Modele_survie<-as.numeric(df_biases$Modele_survie)
  df_biases$Modele_guerison<-as.numeric(df_biases$Modele_guerison)
  return(df_biases)
}

#############Biais##############
########influence des params.#####
#'Calculer le biais moyen des trois estimateurs pour un K n-échantillon pour plusieurs valeurs de k.
#'
#' @param n : taille de l'echantillon permettant d'obtenir un estimateur.
#' @param lambda : parametre de la loi weibull.
#' @param t_star : Fin de la fenetre d'observation.
#' @param p : proportion réelle de non-guéris.
#' @param K : nombre de répétitions de l'expérience.
#' @return : Dataframe. Valeur du biais moyen pour chaque valeur du paramètre k.
#' @export
#'
#' @examples
#' ####test#####
#' liste_biais<-biais.selon.k(n=100,lambda=2,t_tsar=6,p=0.33,K=100)
biais.selon.k <-function(K, n, lambda, t_star,p){
  k <- seq(0.8, 5, by = 0.1)
  results <- NULL

  for(i in c(1:length(k))){
    vec.biais <- Simuler_biais_taillen(K, n, lambda, t_star, p, k[i])
    biais.surv <- vec.biais$Modele_survie
    biais.cure <- vec.biais$Modele_guerison
    biais.bernoulli <- vec.biais$Modele_bernoulli
    results <- rbind(results, c(k[i], mean(biais.surv) - p, mean(biais.cure) -p, mean(biais.bernoulli) -p ))
  }
  return(results)
}

#' Ggplot des valeurs des biais moyens selon la taille des echantillons et du lambda.
#'
#' @param N  nombre de tailles d'echantillon differents.
#' @param window_lambda  ensemble des valeurs de lambda considérées.
#' @param t_star fin de la fenetre d'observation
#' @param n nombre d'individus considérés.
#' @param t_tsar fin de la fenêtre d'observation.
#'
#' @return Plot des valeurs des biais moyens en fonction du lambda et de la taille des echantillons.
#' @export
#'
#' @examples
#' ######Test ######
#' plot<-fnct_compar_plt_biais.selon.k1(N=100,n=100,window_lambda=c(1,2,4),t_tsar=6,p=0.33)
fnct_compar_plt_biais.selon.k1 <- function(N, n, window_lambda, t_star, p) {
  library(gridExtra)
  library(ggplot2)
  # Generate the data
  set.seed(12345)
  RES <- biais.selon.k(N, n, window_lambda[1], t_star, p = p)
  RES0.1.3 <- data.frame(RES)
  colnames(RES0.1.3) <- c("k", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.k(N, n, window_lambda[2], t_star, p = p)
  RES0.2.3 <- data.frame(RES)
  colnames(RES0.2.3) <- c("k", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.k(N, n, window_lambda[3], t_star, p = p)
  RES0.5.3 <- data.frame(RES)
  colnames(RES0.5.3) <- c("k", "mean.surv", "mean.cure", "mean.bernoulli")

  # Get the min and max bounds of each variable to be used in the plots
  borne_min <- min(
    min(RES0.5.3$mean.surv),
    min(RES0.1.3$mean.surv),
    min(RES0.2.3$mean.surv)
  )

  borne_max <- max(
    max(RES0.5.3$mean.surv),
    max(RES0.1.3$mean.surv),
    max(RES0.2.3$mean.surv)
  )

  borne_min.c <- min(
    min(RES0.5.3$mean.cure),
    min(RES0.1.3$mean.cure),
    min(RES0.2.3$mean.cure)
  )

  borne_max.c <- max(
    max(RES0.5.3$mean.cure),
    max(RES0.1.3$mean.cure),
    max(RES0.2.3$mean.cure)
  )

  borne_min.b <- min(
    min(RES0.5.3$mean.bernoulli),
    min(RES0.1.3$mean.bernoulli),
    min(RES0.2.3$mean.bernoulli)
  )

  borne_max.b <- max(
    max(RES0.5.3$mean.bernoulli),
    max(RES0.1.3$mean.bernoulli),
    max(RES0.2.3$mean.bernoulli)
  )

  # Plot the data

  gg1 <-  ggplot2::ggplot(RES0.2.3, aes(k, mean.bernoulli)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(k, mean.bernoulli, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(k, mean.bernoulli, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min.b -0.1, borne_max.b+0.1)+
    ggplot2::labs(
      title = "Mod?le de Bernoulli",
      x = expression(alpha),
      y = "biais moyen",
      color = expression(alpha))+
    theme_bw()


  gg2 <-  ggplot2::ggplot(RES0.2.3, aes(k, mean.surv)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(k, mean.surv, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(k, mean.surv, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min -0.1, borne_max+0.1)+
    ggplot2::labs(
      title = "Mod?le de Survie",
      x = expression(alpha),
      y = "biais moyen",
      color = expression(alpha))+
    ggplot2::theme_bw()

  gg3 <-  ggplot2::ggplot(RES0.2.3, aes(k, mean.cure)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(k, mean.cure, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(k, mean.cure, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min.c -0.1, borne_max.c+0.1)+
    ggplot2::labs(
      title = "Mod?le de Gu?rison",
      x = expression(alpha),
      y = "biais moyen",
      color = expression(alpha))+
    ggplot2::theme_bw()
  g <- gridExtra::grid.arrange(gg1, gg2, gg3, top = paste("influence de", expression(alpha)))

  return(g)

}
#' Evolution du biais en fonction de la taille de l'échantillon.
#'
#' @param K  nombre d'échantillons créées par taille d'échantillon.
#' @param lambda  paramètre de la loi weibull.
#' @param t_star fin de la fenetre d'observation
#' @param p proportion de non-guéris.
#' @return Dataframe. Pour chaque estimateur, valeur du biais moyen selon la taille de l'échantillon.
#' @export
#'
#' @examples
#' ######Test ######
#' df<-biais.selon.lambda(K=100,n=100,k=1,lambda=2,t_tsar=6,p=0.33)
biais.selon.lambda <-function(K, lambda, t_star,p, k){
  results <- NULL

  n <- 20
  while(n<200){
    vec.biais <- Simuler_biais_taillen(K, n, lambda, t_star, p, k)
    biais.surv <- vec.biais$Modele_survie
    biais.cure <- vec.biais$Modele_guerison
    biais.bernoulli <- vec.biais$Modele_bernoulli
    results <- rbind(results, c(n, mean(biais.surv) - p, mean(biais.cure) -p, mean(biais.bernoulli) -p ))
    n <- n+20
  }
  return(results)
}



#' Ggplot des valeurs des biais moyens selon la taille des echantillons et du lambda.
#'
#' @param N nombre de tailles d'echantillon differents.
#' @param window_lambda un vecteur de valeurs pour lambda
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
#' p <- 0.3
#' k <- 1
#' result<-fonction_compar_plotsn_lambda1(N=N,window_lambda=window_lambda,t_star=t_star, p=p, k=k)
fonction_compar_plotsn_lambda1 <- function(N, window_lambda, t_star, p, k) {
  library(gridExtra)
  library(ggplot2)
  # Generate the data
  set.seed(12345)
  RES <- biais.selon.lambda(K=N, lambda=window_lambda[1], t_star=t_star, p = p, k)
  RES0.1.3 <- data.frame(RES)
  colnames(RES0.1.3) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda(K=N, lambda=window_lambda[2], t_star=t_star, p = p, k)
  RES0.2.3 <- data.frame(RES)
  colnames(RES0.2.3) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  set.seed(12345)
  RES <- biais.selon.lambda(K=N, lambda=window_lambda[3], t_star, p = p, k)
  RES0.5.3 <- data.frame(RES)
  colnames(RES0.5.3) <- c("n", "mean.surv", "mean.cure", "mean.bernoulli")

  # Get the min and max bounds of each variable to be used in the plots
  borne_min <- min(
    min(RES0.5.3$mean.surv),
    min(RES0.1.3$mean.surv),
    min(RES0.2.3$mean.surv)
  )

  borne_max <- max(
    max(RES0.5.3$mean.surv),
    max(RES0.1.3$mean.surv),
    max(RES0.2.3$mean.surv)
  )

  borne_min.c <- min(
    min(RES0.5.3$mean.cure),
    min(RES0.1.3$mean.cure),
    min(RES0.2.3$mean.cure)
  )

  borne_max.c <- max(
    max(RES0.5.3$mean.cure),
    max(RES0.1.3$mean.cure),
    max(RES0.2.3$mean.cure)
  )

  borne_min.b <- min(
    min(RES0.5.3$mean.bernoulli),
    min(RES0.1.3$mean.bernoulli),
    min(RES0.2.3$mean.bernoulli)
  )

  borne_max.b <- max(
    max(RES0.5.3$mean.bernoulli),
    max(RES0.1.3$mean.bernoulli),
    max(RES0.2.3$mean.bernoulli)
  )

  # Plot the data

  gg1 <-  ggplot2::ggplot(RES0.2.3, aes(n, mean.bernoulli)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(n, mean.bernoulli, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(n, mean.bernoulli, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min.b -0.1, borne_max.b+0.1)+
    ggplot2::labs(
      title = "Mod?le de Bernoulli",
      x = "n",
      y = "biais moyen",
      color = "n" )+
    ggplot2::theme_bw()
  gg2 <-  ggplot2::ggplot(RES0.2.3, aes(n, mean.surv)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(n, mean.surv, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(n, mean.surv, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min -0.1, borne_max+0.1)+
    ggplot2::labs(
      title = "Mod?le de Survie",
      x = "n",
      y = "biais moyen",
      color = "n")+
    ggplot2::theme_bw()

  gg3 <-  ggplot2::ggplot(RES0.2.3, aes(n, mean.cure)) +
    ggplot2::geom_line(aes(color = "0.2"), size = 0.6) +
    ggplot2::geom_line(data = RES0.5.3, aes(n, mean.cure, color = "0.5"), size = 0.6) +
    ggplot2::geom_line(data = RES0.1.3, aes(n, mean.cure, color = "0.1"), size = 0.6) +
    ggplot2::scale_color_manual(name = expression(lambda), values = c("red", "black", "blue")) +
    ggplot2::ylim(borne_min.c -0.1, borne_max.c+0.1)+
    ggplot2::labs(
      title = "Mod?le de Gu?rison",
      x = "n",
      y = "biais moyen",
      color = "n")+
    ggplot2::theme_bw()


  g <- gridExtra::grid.arrange(gg1, gg2, gg3, top = "influence de n et lambda" )
  return(g)

}

########## en fixant les parametres.
#' Calcul du biais moyen des trois estimateurs pour K n-échantillons.
#'
#' @param K nombre d'échantillons.
#' @param lambda paramètre de la loi Weibull (scale).
#' @param t_star fin de la fenetre d'observation
#' @param k paramètre de la loi Weibull (shape).
#' @return Liste. Valeur du biais moyen pour les trois estimateurs.
#' @export
#'
#' @examples
#' ######Test ######
#' t_star1<-6
#' K<-10
#' n<-100
#' result<-Calcul_biais_moyen_taillen(K=K,n=n,t_star=t_star1,lambda=0.5,p=0.33,k=2)
Calcul_biais_moyen_taillen<-function(K,n,lambda,t_star,p,k){
  # on effectue la simulation des biais pour K ?chantillons de taille n selon
  # les deux mod?les (de gu?rison, de survie)
  data<-Simuler_biais_taillen(K,n,lambda,t_star,p,k)
  # on va calculer le biais moyen. On pr?pare donc un vecteur pour stocker les
  # deux biais moyens
  result<-rep(NA,3)
  #on calcule les biais moyens
  result<-colMeans(data)
  # Rappel : biais = estimateur - valeur th?orique
  result[1]<-result[1]-p
  result[2]<-result[2]-p
  result[3]<-result[3]-p
  return(result)
}
###########Evolution du biais#############
#' Evolution du biais moyen des trois estimateurs pour K n-échantillons.
#'
#' @param K nombre d'échantillons.
#' @param lambda paramètre de la loi Weibull (scale).
#' @param t_star fin de la fenetre d'observation
#' @param k paramètre de la loi Weibull (shape).
#' @return Ggplot. Valeur du biais moyen pour les trois estimateurs pour n allant de 20 à 100.
#' @export
#'
#' @examples
#' ######Test ######
#' t_star1<-6
#' K<-10
#' n<-100
#' result<-biais.selon.taille_echantillon(K=K,t_star=t_star1,lambda=0.5,p=0.33,k=2)
biais.selon.taille_echantillon <- function(K, lambda, t_star, p, k){
  require(ggplot2)
  require(gridExtra)
  # On fixe un n de d?part ? 10 individus et on incr?ment par 5 jusqu'a 100
  debut <- 20
  fin <- 100
  pas <- 5
  n <- seq(debut,fin , pas)

  # On calcule le biais pour K simulations et n-?chantillons
  liste_parameter <- list(lambda, t_star, p, k)
  names(liste_parameter)<-c("lambda","t_star","p","k")
  result_final <- fonction_generation_taille_mean(vector_size = n, liste_parameter = liste_parameter, K=K)
  result_final$n <- n

  colnames(result_final) <- c("modele_bernoulli","modele_survie", "modele_guerison", "taille_echantillon")
  # plot
  borne_min <- min(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  borne_max <- max(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)

  borne_min <- min(result_final$modele_bernoulli, result_final$modele_guerison, result_final$modele_survie)
  borne_max <- max(result_final$modele_bernoulli, result_final$modele_guerison, result_final$modele_survie)

  gg1 <- ggplot2::ggplot(data = result_final, aes(x = taille_echantillon)) +
    ggplot2::geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
    ggplot2::geom_smooth(aes(y = modele_bernoulli, col = "bernoulli"), size = 1.2, alpha = 0.5) +
    ggplot2::scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "bernoulli" = "blue")) +
    ggplot2::ggtitle("Evolution du biais en \nfonction de n") +
    ggplot2::xlab("Taille echantillon") + ggplot2::ylab("Biais") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24)
          ,legend.text = element_text(size = 20)
          # , legend.title = element_text(size = 22)
          , plot.caption = element_text(size = 20)) +
    ggplot2::ylim(borne_min, borne_max)

  gg2 <- ggplot2::ggplot(data = result_final, aes(x = taille_echantillon)) +
    ggplot2::geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
    ggplot2::geom_smooth(aes(y = modele_survie, col = "survie"), size = 1.2, alpha = 0.5) +
    ggplot2::scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "survie" = "darkgreen")) +
    ggplot2::ggtitle("Evolution du biais moyen en \n fonction de n") +
    ggplot2::xlab("Taille echantillon") + ylab("Biais") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title=element_blank(),
          axis.text=element_text(family = "Helvetica", size=20),
          axis.title=element_text(family = "Helvetica", size=20),
          plot.title = element_text(family = "Helvetica", size = 24),
          legend.text = element_text(size = 20)
          # , legend.title = element_text(size = 22)
          , plot.caption = element_text(size = 20)) +
    ggplot2::ylim(borne_min, borne_max)+
    ggplot2::labs(caption = sprintf("N = %s, p=%s,lambda=%s,alpha=%s" ,
                           as.character(K),
                           as.character(p),
                           as.character(round(lambda,2)),
                           as.character(alpha)))

  gg <- gridExtra::grid.arrange(gg1, gg2, ncol = 2, widths = c(7,7))

}

######### EQM ################
#' Evolution de l'EQM moyen des trois estimateurs pour K n-échantillons.
#'
#' @param K nombre d'échantillons.
#' @param lambda paramètre de la loi Weibull (scale).
#' @param t_star fin de la fenetre d'observation
#' @param k paramètre de la loi Weibull (shape).
#' @return Ggplot. Valeur de l'eqm moyen pour les trois estimateurs pour n allant de 20 à 100.
#' @export
#'
#' @examples
#' ######Test ######
#' t_star1<-6
#' K<-10
#' n<-100
#' result<-eqm.selon.taille_echantillon(K=K,t_star=t_star1,lambda=0.5,p=0.33,k=2)
eqm.selon.taille_echantillon <- function(K, lambda, t_star, p, k){
  require(ggplot2)
  require(gridExtra)
  # On fixe un n de d?part ? 10 individus et on incr?ment par 5 jusqu'a 100
  debut <- 20
  fin <- 100
  pas <- 5
  n <- seq(debut,fin , pas)

  # On calcule le biais pour K simulations et n-?chantillons
  liste_parameter <- list(lambda, t_star, p, k)
  names(liste_parameter)<-c("lambda","t_star","p","k")
  result_final <- fonction_generation_taille_eqm(vector_size = n, liste_parameter = liste_parameter, K=K)
  result_final$n <- n
  colnames(result_final) <- c("modele_bernoulli","modele_survie", "modele_guerison", "taille_echantillon")
  # plot
  borne_min <- min(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)
  borne_max <- max(result_final$modele_guerison, result_final$modele_survie,result_final$modele_bernoulli)

  # define color palette
  palette <- c("#0072B2", "#D55E00", "#E69F00")

  gg1 <- {ggplot2::ggplot(data = result_final, aes(x = taille_echantillon)) +
      ggplot2::geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
      ggplot2::geom_smooth(aes(y = modele_bernoulli, col = "bernoulli"), size = 1.2, alpha = 0.5) +
      ggplot2::scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "bernoulli" = "blue1")) +
      ggplot2::xlab("Taille echantillon") + ggplot2::ylab("EQM") +
      #theme_classic() +
      ggplot2::theme(legend.title=element_blank(),
            axis.text=element_text(family = "Helvetica", size=20),
            axis.title=element_text(family = "Helvetica", size=20),
            plot.title = element_text(family = "Helvetica", size = 24)
            , legend.text = element_text(family = "Helvetica", size = 20)
            ,text = element_text(size=rel(20))) +
      ggplot2::ylim(borne_min, borne_max)}

  gg2 <- {ggplot2::ggplot(data = result_final, aes(x = taille_echantillon)) +
      ggplot2::geom_smooth(aes(y = modele_guerison, col = "guerison"), size = 1.2, alpha = 0.5) +
      ggplot2::geom_smooth(aes(y = modele_survie, col = "survie"), size = 1.2, alpha = 0.5) +
      ggplot2::scale_color_manual(name = "Modeles", values = c("guerison" = "red1", "survie" = "darkgreen")) +
      ggplot2::xlab("Taille echantillon") + ggplot2::ylab("EQM") +
      # theme_classic() +
      ggplot2::theme(legend.title=element_blank(),
            axis.text=element_text(family = "Helvetica", size=20),
            axis.title=element_text(family = "Helvetica", size=20),
            plot.title = element_text(family = "Helvetica", size = 24)
            , legend.text = element_text(family = "Helvetica", size = 20)
            ,text = element_text(size=rel(20))) +
      ggplot2::ylim(borne_min, borne_max)+
      ggplot2::labs(caption = sprintf("N = %s, p=%s,lambda=%s,alpha=%s" ,
                             as.character(K),
                             as.character(p),
                             as.character(round(lambda,2)),
                             as.character(alpha)))}

  gg <- {gridExtra::grid.arrange(gg1, gg2, ncol = 2, widths = c(7,7)
                      ,top =textGrob("Evolution de l'EQM en fonction de la taille d'echantillon n",gp=gpar(fontsize=24,font=3)))}

}
