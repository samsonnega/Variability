CV <- read.csv("cv trait means.csv")

CVplot <- ggplot(CV, aes(x=trait, y=cv)) + 
  geom_point(col="burlywood4", size=5) +   # Draw points
  geom_errorbar(width=.1, aes(ymin=cv-ci, ymax=cv+ci)) +
  labs(y = "Mean (squared) standardized variation", x = "") +  
  coord_flip() +
  theme_classic(base_size = 20)
  
  

ggsave(CVplot, 
       file = "standard var.jpg", 
       width = 12, height = 6, unit = "in", dpi = 300)
