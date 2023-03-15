library(ggplot2)
# Note we convert the cyl variable to a factor here in order to fill by cylinder
ggplot(mpg) + 
  geom_bar(aes(x = class, fill = factor(cyl)), position = position_dodge(preserve = 'single'))


download.file("http://r-bio.github.io/data/surveys_complete.csv",
              "data/surveys_complete.csv")

surveys_complete <- read.csv(file = "surveys_complete.csv")

ggplot(subset(surveys_complete, species_id %in% c("DO", "DM", "DS") & sex %in% c("F", "M")),
       aes(x = sex, y = weight,  colour = interaction(sex, species_id))) + facet_wrap( ~ species_id) +
  geom_point(alpha = 0.3, position = "jitter") +
  geom_boxplot(alpha = 0, colour = "black")

###### Bulk Soil ######
z = read.csv("7-levels_samples_stat_2_for_plotting.csv", header = T)
pdf("Bulksoil_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("CO-CL-BO", "CO-CL-YO","CO-CY-BU","CO-CY-YO",
                              "CO-SC-HE", "CO-SC-SH", "CO-SL-AN", 
                              "CO-SL-BE", "CO-SL-SH")),
       aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Sugarbeet ######
pdf("Sugarbeet_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("SU-CL-BO",
                              "SU-CL-YO",
                              "SU-CY-BU",
                              "SU-CY-YO",
                              "SU-SC-HE",
                              "SU-SC-SH",
                              "SU-SL-AN",
                              "SU-SL-BE",
                              "SU-SL-SH"
                              )),
       aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()


###### Barley ######
pdf("Barley_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("BA-CL-BO",
                              "BA-CL-YO",
                              "BA-CY-BU",
                              "BA-CY-YO",
                              "BA-SC-HE",
                              "BA-SC-SH",
                              "BA-SL-AN",
                              "BA-SL-BE",
                              "BA-SL-SH"
                              
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Wheat ######
pdf("Wheat_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("WH-CL-BO",
                              "WH-CL-YO",
                              "WH-CY-BU",
                              "WH-CY-YO",
                              "WH-SC-HE",
                              "WH-SC-SH",
                              "WH-SL-AN",
                              "WH-SL-BE",
                              "WH-SL-SH"
                              
                              
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()


###### Oilseed Rape ######
pdf("OSR_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("OR-CL-BO",
                              "OR-CL-YO",
                              "OR-CY-BU",
                              "OR-CY-YO",
                              "OR-SC-HE",
                              "OR-SC-SH",
                              "OR-SL-AN",
                              "OR-SL-BE",
                              "OR-SL-SH"
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Beans ######
pdf("Beans_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("BE-CL-BO",
                              "BE-CL-YO",
                              "BE-CY-BU",
                              "BE-CY-YO",
                              "BE-SC-HE",
                              "BE-SC-SH",
                              "BE-SL-AN",
                              "BE-SL-BE",
                              "BE-SL-SH"
                              
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Beans ######
pdf("Beans_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("BE-CL-BO",
                              "BE-CL-YO",
                              "BE-CY-BU",
                              "BE-CY-YO",
                              "BE-SC-HE",
                              "BE-SC-SH",
                              "BE-SL-AN",
                              "BE-SL-BE",
                              "BE-SL-SH"
                              
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Oats ######
pdf("Oats_ASV.pdf",6.5 ,4)
ggplot(subset(z, Group %in% c("OA-CL-BO",
                              "OA-CL-YO",
                              "OA-CY-BU",
                              "OA-CY-YO",
                              "OA-SC-HE",
                              "OA-SC-SH",
                              "OA-SL-AN",
                              "OA-SL-BE",
                              "OA-SL-SH"
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()

###### Controls ######
pdf("controls_ASV.pdf",3,4)
ggplot(subset(z, Group %in% c("positive", "negative"
)),
aes(x = Group, y = ASV,  colour = interaction(Group, Type))) +
  geom_point(alpha = 1, position = "jitter", size = 4) + 
  geom_boxplot(alpha = 0, colour = "black", size = 0.8)+ 
  theme_classic() +   theme(text = element_text(size=15, colour = "black"), 
                            axis.ticks = element_line(colour = "black", size = 1.25),
                            axis.line = element_line(colour = 'black', size = 1.25),
                            axis.text.x = element_text(angle=45, hjust=1, colour = "black", size = 13),
                            axis.text.y = element_text(angle=0, hjust=0.5, colour = "black",size = 13),
                            axis.title.y = element_text(color="black", size=15,face="bold"), legend.position = "none")
dev.off()
