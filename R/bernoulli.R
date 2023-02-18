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
# n<-10
# test<-simul_bernoulli(100,0.1)

#' calcule le biais de la probabilite estimee de toxicite selon une Bernoulli(p)
#' Compare avec la valeur theorique.
#'
#' @param n taille de l'echantillon
#' @param p parametre de la loi de Bernoulli utilisee
#'
#' @return valeur du biais
#' @export
#'
#' @examples
biais_pi<-function(n,p){
  ### calculer le biais de la probabilite estimee de toxicite selon une Bernoulli(p).
  ### Comparaison avec la valeur theorique.
  simulation<-simul_bernoulli(n,p)
  biais<-abs(mean(simulation)-p)
  return(biais)
}
test_biais<-biais_pi(18,0.33)

Simuler_Nfois_n_echantillons_bern<-function(N,n,p){
  #### Simuler
  vecteur_biais<-rep(NA,N)
  vecteur_taille<-rep(n,N)
  ##On utilise le p donn? en entr?e.
  vecteur_biais<-sapply(vecteur_taille,biais_pi,p=p)
  return(vecteur_biais)
}
N<-100
p<-0.33
test_simul_bern_total<-Simuler_Nfois_n_echantillons_bern(N,n,p)
boxplot(test_simul_bern_total,main="Distribution du biais pour le mod?le de guerison",col="red")
