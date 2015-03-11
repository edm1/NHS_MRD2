# Make a bubble plot from a .tab file
# 
# It can be ran from the command line or interactively. To run from cmd line:
#    $ Rscript make-bubbleplot.R <tab file>

library("ggplot2")
library("plyr")
library("tools")

# Options if being ran interactively
in.tab = "/Users/ed/Work/jack-bubble-plot-150305/MA_CSF_top1000.tab"

# Get command line args
isinteractive = interactive()
if (!isinteractive) {
  args = commandArgs(trailingOnly = TRUE)
  in.tab = args[1]
  # Set plotting device to pdf
  out.prefix = file_path_sans_ext(in.tab)
  out.name = paste0(out.prefix, "_bubbleplot.pdf")
  pdf(file=out.name)
}

# Load unprocessed and get columns of interest
df = read.table(in.tab, sep="\t", header=F, stringsAsFactors=F)
data = df[, c("V5", "V7", "V3")]
colnames(data) = c("v_full", "j_full", "reads")
head(data)

# Split v by "-" and take first part
get.v.type = function(x) { unlist(strsplit(x, "*", fixed=T))[1] }
data["v_type"] = sapply(data$v_full, get.v.type)
head(data) # Show table
unique(data$v_type) # Show unique V parts

# Add up read totals for V J pairs
data.totals = ddply(data, .(v_type, j_full), summarize, total=sum(reads))
head(data.totals)

# Get the V family
get.v.fam = function(x) { unlist(strsplit(x, "-", fixed=T))[1] }
data.totals["v_fam"] = factor(sapply(data.totals$v_type, get.v.fam))
head(data.totals)

# Make plot title
title.label = basename(file_path_sans_ext(in.tab))
title.label = paste0("IgH Clones - ", title.label)

# Make bubble plot
ggplot(data.totals, aes(y=v_type, x=j_full, size=total, colour=v_fam)) + 
  geom_point() +
  scale_size(range = c(2, 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(colour="V family", x="J type", y="V type", size="Total reads") +
  ggtitle(title.label)


# ~~~~~~~~~~~~~~~~~~~old coce ~~~~~~~~~~~~~~~~~~~~~~~~#

#   ggtitle('IgH Clones REH')+
#   labs(x="V gene", y="J gene")+
#   theme(plot.title = element_text(size=20, face="bold", vjust=2))+
#   theme(axis.title.x = element_text(color="forestgreen", vjust=-0.35))+
#   theme(axis.title.y = element_text(color="cadetblue" , vjust=0.35)) +
#   theme(axis.text.x=element_text(angle=90, size=9, vjust=0.5))+
#   theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))+
#   scale_colour_tableau()+
#   scale_color_discrete(name="Vh gene")