library(ggplot2)
library(gridExtra)

#.txt files contain effect size infor for features with ES greater than |1|
#If more than 10 are significant, only the most 10 divergent features were shown (5 top 5 bottom)

tax<-read.table("aldex/plot_tax.txt", sep="\t", quote="", check.names=F, header=T, comment.char="")
tax<-data.frame(tax)
enz<-read.table("aldex/plot_enz.txt", sep="\t", quote="", check.names=F, header=T, comment.char="")
enz<-data.frame(enz)
paths<-read.table("aldex/plot_paths.txt", sep="\t", quote="", check.names=F, header=T, comment.char="")
paths<-data.frame(paths)

t1<-tax[grep("1",tax$Donor), ]
t2<-tax[grep("2",tax$Donor), ]
e1<-enz[grep("1",enz$Donor), ]
e2<-enz[grep("2",enz$Donor), ]
p1<-paths[grep("1",paths$Donor), ]
p2<-paths[grep("2",paths$Donor), ]

# lock in factor level order
t1$Feature <- factor(t1$Feature, levels = t1$Feature)

t1<-ggplot(data=t1, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

# lock in factor level order
t2$Feature <- factor(t2$Feature, levels = t2$Feature)

t2<-ggplot(data=t2, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

# lock in factor level order
# so that ggplot doesn't order all the features alphabetically
e1$Feature <- factor(e1$Feature, levels = e1$Feature)

e1<-ggplot(data=e1, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

# lock in factor level order
e2$Feature <- factor(e2$Feature, levels = e2$Feature)

e2<-ggplot(data=e2, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

# lock in factor level order
p1$Feature <- factor(p1$Feature, levels = p1$Feature)

p1<-ggplot(data=p1, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

# lock in factor level order
p2$Feature <- factor(p2$Feature, levels = p2$Feature)

p2<-ggplot(data=p2, aes(x=Feature, y=EffectSize))+
geom_bar(stat="identity")+
coord_flip()+
theme_minimal()

#from gridExtra package, arrange the plots into a panel
grid.arrange(t1, e1, p1, t2, e2, p2, nrow=2)
