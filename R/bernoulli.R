#install.packages("Rlab")                         # Install Rlab package
# library("Rlab")
#'  simule un n_echantillon d'une Bernoulli(p).
#'
#' @param n taille de l'echantillon
#' @param p parametre de la loi de Bernoulli
#'
#' @return un vecteur de taille n representant un n-echantillon de Bernoulli
#' @export
#'
#' @examples
#' n<-10
#' test<-simul_bernoulli(100,0.1)
#' test
simul_bernoulli<-function(n,p){
  ### simuler un n_echantillon d'une Bernoulli(p).
  return(rbern(n,p))
}


#' calcule le biais de la probabilite estimee de toxicite selon une Bernoulli(p)
#' Compare avec la valeur theorique.
#'
#' @param n taille de l'echantillon
#' @param p parametre de la loi de Bernoulli utilisee
#'
#' @return valeur du biais retourné
#' @export
#'
#' @examples
#' test_biais<-biais_pi(18,0.33)
biais_pi<-function(n,p){
  ### calculer le biais de la probabilite estimee de toxicite selon une Bernoulli(p).
  ### Comparaison avec la valeur theorique.
  simulation<-simul_bernoulli(n,p)
  biais<-abs(mean(simulation)-p)
  return(biais)
}


#' retourne un vecteur de taille N constitue des biais issu d'une simulation de
#' N echantillons de taille n pour une loi de Bernoulli
#'
#' @param N le nombre d'echantillons generes
#' @param n la taille de chacun des N echantillons generes
#' @param p parametre de la loi de Bernoulli
#'
#' @return un vecteur de taille N
#' @export
#'
#' @examples
#' n<-10
#' N<-100
#' p<-0.33
#' test_simul_bern_total<-Simuler_Nfois_n_echantillons_bern(N,n,p)
Simuler_Nfois_n_echantillons_bern<-function(N,n,p){
# boxplot(test_simul_bern_total,main="Distribution du biais pour le modele de guerison",col="red")
  #### Simuler
  vecteur_biais<-rep(NA,N)
  vecteur_taille<-rep(n,N)
  ##On utilise le p donn? en entr?e.
  vecteur_biais<-sapply(vecteur_taille,biais_pi,p=p)
  return(vecteur_biais)
}

# boxplot(test_simul_bern_total,main="Distribution du biais pour le mod?le de guerison",col="red")
