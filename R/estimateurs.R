fonction_KM<-function(df,t_star){
  df<-df[,c("tox_time","is_observed")]
  indice_cens<-which(df$is_observed==0)
  if(length(indice_cens)==0){return(1)}
  surv_object<-survival::Surv(as.numeric(df$tox_time),event=df$is_observed)
  fit <- survival::survfit(surv_object ~1, data = df)
  estimateur_survie<-1-tp.surv(fit,6)[3][[1]]
  return(estimateur_survie)
}
fonction_Bern<-function(df){
  return(mean(df$is_observed))
}
fonction_cure<-function(df,t_star){
  # retourne la probabilite de ne pas avoir fait de DLT a T_star
  indice_observed<-which(df$is_observed==1)
  indice_censored<-which(df$is_observed==0)
  df$covar<-rep(1,nrow(df))
  df2<-df[,c("covar","is_observed","tox_time")]
  require(cuRe)
  fit.wei <- cuRe::fit.cure.model(survival::Surv(tox_time,is_observed) ~ 1, formula.surv=list(~1),data =df2,
                            dist="weibull",link="logit")
  prob_cure<-plogis(as.numeric(fit.wei$coefs[1]))
  estimateur_tox<-1-prob_cure
  return(estimateur_tox)
}
