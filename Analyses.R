

library(mgcv)

library(ggplot2)


###########################################################################


dat<-read.csv("Post_disturbance_recovery_dataset.csv",stringsAsFactors=F,header=T,
              sep=";",dec=",")
names(dat)

#> names(dat)
#[1] "site"                          "plot"                          "Pre.disturbance.management"    "Post.disturbance.management"  
#[5] "Years.since.disturbance"       "Disturbance.agent"             "Elevation"                     "Aridity.Change"               
#[9] "HLI"                           "Late.successional.WDensity"    "Early.successional.WDensity"   "Pioneer.WDensity"             
#[13] "Total.WDensity"                "Pioneer.RawDensity"            "Early.successional.RawDensity" "Late.successional.RawDensity" 
#[17] "Total.RawDensity"              "Pre_Density"                   "Pre_Sp1"                       "Pre_Sp2"                      
#[21] "WDensity_Sp1"                  "WDensity_Sp2"


###########################################################################

# plots and summaries

table(dat$Pre.disturbance.management,dat$Post.disturbance.management)
table(dat$Disturbance.agent,paste(dat$Pre.disturbance.management,dat$Post.disturbance.management,sep="-"))
pairs(dat[,c("Years.since.disturbance","Elevation",
             "HLI","Aridity.Change")],pch=".")
pairs(dat[,c("Late.successional.WDensity","Early.successional.WDensity",
             "Pioneer.WDensity","Total.WDensity")],pch=".")
pairs(dat[,c("Pioneer.RawDensity","Early.successional.RawDensity",
             "Late.successional.RawDensity","Total.RawDensity")],pch=".")
plot(dat$Late.successional.RawDensity,dat$Late.successional.WDensity)
plot(dat$Early.successional.RawDensity,dat$Early.successional.WDensity)
plot(dat$Pioneer.RawDensity,dat$Pioneer.WDensity)
plot(dat$Total.RawDensity,dat$Total.WDensity)
#
pairs(log(dat[,c("Years.since.disturbance","Elevation",
             "HLI","Aridity.Change")]),pch=".")
pairs(log(dat[,c("Late.successional.WDensity","Early.successional.WDensity",
             "Pioneer.WDensity","Total.WDensity")]),pch=".")
pairs(log(dat[,c("Pioneer.RawDensity","Early.successional.RawDensity",
             "Late.successional.RawDensity","Total.RawDensity")]),pch=".")

plot(log(dat$Late.successional.RawDensity),log(dat$Late.successional.WDensity))
plot(log(dat$Early.successional.RawDensity),log(dat$Early.successional.WDensity))
plot(log(dat$Pioneer.RawDensity),log(dat$Pioneer.WDensity))
plot(log(dat$Total.RawDensity),log(dat$Total.WDensity))
hist(dat$Total.RawDensity)
hist(dat$Total.WDensity)
hist(log(dat$Total.RawDensity))
hist(log(dat$Total.WDensity))



##################################################################################

# analyses (only wind-disturbed sites)


smaz<-dat[dat$Disturbance.agent=="Wind",]
smaz$f.site<-factor(smaz$site)
smaz$f.pre.dist<-factor(smaz$Pre.disturbance.management,levels=c("N","L","Y"))
smaz$f.post.dist<-factor(smaz$Post.disturbance.management,levels=c("None","Logging only","Intensive"))


### Total weighted density

v<-gam(log(Total.WDensity+1)~s(f.site,bs="re")+
            s(Elevation)+
            f.pre.dist+f.post.dist+
            Years.since.disturbance+
            Aridity.Change+
            HLI,
          data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                  Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                  HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.total.reg.density.csv",row.names=F)


### Pioneer Weighted Density

v<-gam(log(Pioneer.WDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.pioneer.reg.density.csv",row.names=F)


### Early-successional Weighted Density

v<-gam(log(Early.successional.WDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.early.reg.density.csv",row.names=F)

### Late-successional Weighted Density

v<-gam(log(Late.successional.WDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.late.reg.density.csv",row.names=F)


### Total Raw Density


v<-gam(log(Total.RawDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.total.raw.density.csv",row.names=F)


### Pioneer Raw Density


v<-gam(log(Pioneer.RawDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.pioneer.raw.density.csv",row.names=F)


### Early-Successional Raw Density


v<-gam(log(Early.successional.RawDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.early.raw.density.csv",row.names=F)



### Late-Successional Raw Density


v<-gam(log(Late.successional.RawDensity+1)~s(f.site,bs="re")+
         s(Elevation)+
         f.pre.dist+f.post.dist+
         Years.since.disturbance+
         Aridity.Change+
         HLI,
       data=smaz)
anova(v)
summary(v)
AIC(v)
plot(v,scale=0)
plot(v,scale=0,select=2,xlab="Elevation",ylab="effect",lwd=2)
gam.check(v)
concurvity(v)
#
x<-seq(from=min(smaz$Elevation,na.rm=T),to=max(smaz$Elevation,na.rm=T),by=5)
n<-length(x)
pom<-data.frame(f.site=rep(smaz$f.site[1],n),
                Elevation=x,
                f.pre.dist=rep(smaz$f.pre.dist[1],n),
                f.post.dist=rep(smaz$f.post.dist[1],n),
                Years.since.disturbance=rep(mean(smaz$Years.since.disturbance,na.rm=T),n),
                Aridity.Change=rep(mean(smaz$Aridity.Change,na.rm=T),n),
                HLI=rep(mean(smaz$HLI,na.rm=T),n))
a<-predict(v,newdata=pom,se.fit=T,type="terms")
vysl<-data.frame(Elevation=x,effect.estimate=a$fit[,"s(Elevation)"])
vysl$L<-vysl$effect.estimate-1.96*a$se.fit[,"s(Elevation)"]
vysl$U<-vysl$effect.estimate+1.96*a$se.fit[,"s(Elevation)"]
plot(vysl$Elevation,vysl$effect.estimate,type="l",ylim=range(c(vysl$L,vysl$U)))
lines(vysl$Elevation,vysl$L,type="l",lty=2)
lines(vysl$Elevation,vysl$U,type="l",lty=2)
write.csv(vysl,"log.late.raw.density.csv",row.names=F)

