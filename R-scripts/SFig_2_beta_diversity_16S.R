
## bARLEY VS bULKSOIL
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Barley", "Bulksoil"))
my_palette <- c("darkgoldenrod", "limegreen")
NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)


plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))


## Beans vs Bulksoil
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Beans", "Bulksoil"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))

## Oatsvs Bulk Soil
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Oats", "Bulksoil"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))


## OilSeedRapre vs Bulksoil
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("OilseedRape", "Bulksoil"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))


### Wheat vs Bulksoil
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Wheat", "Bulksoil"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))


physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Sugarbeet", "Bulksoil"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))

### Beans vs Wheat vs Barley vs Oats
physeq.SU <- physeq.norm %>% subset_samples(Type %in% c("Beans", "Wheat","Barley" ,"Oats"))

NMDS <- ordinate(physeq = physeq.SU, 
                 method = "NMDS", 
                 distance = "bray"
)

my_palette <- c("darkgoldenrod", "limegreen", "brown", "grey")

plot_ordination(physeq = physeq.SU,
                ordination = NMDS,
                color = "Type",
                shape = "Type"
) +
  theme_classic() + 
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  geom_point(aes(color = Type), alpha = 1, size = 3) +
  theme(
    text = element_text(size = 18, colour = "black"), 
    axis.ticks = element_line(colour = "black", size = 1.1),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, 
                               size = 13, face = "bold"),
    axis.title.y = element_text(color = "black", size = 20, face = "bold"), 
    axis.title.x = element_text(color = "black", size = 20, face = "bold")
  ) + stat_ellipse(geom="polygon", type="norm", alpha=0.10, aes(fill=Type)) +
  scale_shape_manual(values=c(15,17,3,4,16,18,21,22,23))
