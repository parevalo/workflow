library(ggplot2)
library(dplyr)
library(tidyr)
library(extrafont)

#font_import() # Only needs to be run once
loadfonts() #(device="postscript")

me_biannual = read.csv("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_sampling/results/post_katelyn/tables/margin_error_original_lut_buffered_3B.csv")
me_single = read.csv("/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/original_sampling/results/original_buffer3B/margin_error_step2_original_lut.csv")

# df's created manually bc the reshape2 package is currently unavailable for me
time = seq(2002, 2014, 2)
label = rep("Biannual", 7)
me = me_biannual$Forest.to.pasture
biannual_df = cbind(time, me, label)

label = rep("Single", 7)
me = me_single$Forest.to.pasture
single_df = cbind(time, me, label)

df = as_tibble(rbind(biannual_df, single_df))
df$me = as.numeric(df$me) * 100

df_plot = ggplot(data=df, aes(x=time, y=me, fill=label, order=label)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Time") + ylab("Margin of error [%]") +
  labs(fill='Sampling type') + scale_fill_brewer(palette="Set1") +
  theme(text=element_text(size=8, family="Times New Roman"))

outfile = "/media/paulo/785044BD504483BA/OneDrive/Lab/area_calculation/biannual_sampling/results/post_katelyn/figures/me_compare.pdf" 
ggsave(outfile, device="pdf", width = 90, height = 70, units='mm')
embed_fonts(outfile)
