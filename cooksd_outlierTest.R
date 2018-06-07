library(ggplot2)
library(ggrepel)
library(Cairo)
library(car)
my_data <- read.csv("ns.chr1_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr1_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot1<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 1')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12, face = "bold"), axis.title.y=element_text(size=12, face = "bold"))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr1.pdf", plot1, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr2_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr2_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot2<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 2')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12, face = "bold"), axis.title.y=element_text(size=12, face = "bold"))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr2.pdf", plot2, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr3_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr3_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot3<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 3')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr3.pdf", plot3, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr4_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr4_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot4<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 4')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr4.pdf", plot4, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr5_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr5_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot5<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 5')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr5.pdf", plot5, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))



my_data <- read.csv("ns.chr6_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr6_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot6<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 6')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr6.pdf", plot6, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))



my_data <- read.csv("ns.chr7_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr7_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot7<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 7')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr7.pdf", plot7, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr8_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr8_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot8<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 8')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr8.pdf", plot8, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))

my_data <- read.csv("ns.chr9_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr9_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot9<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 9')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr9.pdf", plot9, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr10_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr10_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot10<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 10')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr10.pdf", plot10, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr11_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr11_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot11<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 11')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr11.pdf", plot11, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr12_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr12_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot12<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 12')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr12.pdf", plot12, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr13_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr13_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot13<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 13')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr13.pdf", plot13, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))

my_data <- read.csv("ns.chr14_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr14_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot14<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 14')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr14.pdf", plot14, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr15_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr15_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot15<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 15')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr15.pdf", plot15, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))



my_data <- read.csv("ns.chr16_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr16_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot16<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 16')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr16.pdf", plot16, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr17_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr17_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot17<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 17')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr17.pdf", plot17, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr18_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr18_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot18<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 18')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr18.pdf", plot18, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))


my_data <- read.csv("ns.chr19_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr19_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot19<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 19')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr19.pdf", plot19, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))



my_data <- read.csv("ns.chr20_50k.csv", header=TRUE)
mod <- lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data)
cooksd <- cooks.distance(mod)
labels <- cooksd>20*mean(cooksd, na.rm=T)
T<-which(labels==TRUE)
my_data[T,]
write.csv(my_data[T,], "ns_chr20_50k_cooksdOut_20.csv")
ggda<-data.frame(x=seq(length(cooksd)),cooksd,labels)
p<-ggplot(data=ggda)
plot20<-p+geom_point(aes(x=x,y=cooksd,col=labels))+
labs(x="Window's index", y="Cook's distance", title='Chromosome 20')+
geom_hline(yintercept=20*mean(cooksd, na.rm=T), col="blue")+
geom_text_repel(aes(x=x,y=cooksd),label=ifelse(cooksd>20*mean(cooksd, na.rm=T),names(cooksd),""),size=3.5)+
theme(legend.position=c(.85,.92), legend.key=element_blank(), legend.background = element_blank(), legend.title=element_blank())+
theme(plot.title = element_text(size=16,colour = "blue",face = "bold"),axis.title.x =element_text(size=12), axis.title.y=element_text(size=12))+
theme(legend.text = element_text(size = 12, hjust = 3, vjust = 3, face = 'bold'))+
scale_colour_discrete(labels = c('Normal window','Outlier window'))
ggsave("chr20.pdf", plot20, width = 15, height = 8) 
outlierTest(lm(my_data$CE1+my_data$CE2+my_data$CR1+my_data$SM1+my_data$SM2+my_data$TM1+my_data$XH1+my_data$PM1~my_data$end, data=my_data))









