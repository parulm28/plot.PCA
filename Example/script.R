library(ggplot2)

ifrm <- function(obj, env = globalenv()) 
{
    obj <- deparse(substitute(obj))
    if(exists(obj, envir = env)) 
    {
        rm(list = obj, envir = env)
    }
}

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) 
{
  stop("2 input files required: \n1. Gene intensity file\n2. Metadata", call.=FALSE)
}

Colors <- c("dodgerblue2","#E31A1C","green4","#6A3D9A",
            "#FF7F00", "yellow","skyblue2","#FB9A99","palegreen2","red4")

# Reading gene expression data
df = read.csv(args[1], header = TRUE, row.names = 1)
if (colnames(df)[1] != 'symbol')
{
  stop ("The second column of gene intensity file should be gene symbols with a header symbol.\n", call.=FALSE)
}

ifrm(a)

# identify rows containing junk data
for (i in 2:ncol(df))
{
  if (exists("a") == F)
  {
    a = which(!grepl('^[0-9]',df[,i]))
  }
  else
  {
    a = c(a, which(!grepl('^[0-9]',df[,i])))
  }
}
a = sort.int(unique(a))
if (length(a) > 0.5*(nrow(df)))
{
	cat ("WARNING: More than 50% genes contain junk data.", "\n", "These genes will be removed from the analysis.", "\n")
}

if (length(a) > 0.8*(nrow(df)))
{
	stop ("WARNING: More than 80% genes contain junk data.\nThese genes will be removed from the analysis.", call.=FALSE)
}

#removing gene expression data for which values are junk

df2 = df[-c(a), ]
df3 = df2
df3[2:ncol(df2)] = as.data.frame(sapply(df2[2: ncol(df2)], as.numeric))
df3 = df3[!duplicated(df3[2:ncol(df3)]), ]
rm(df, df2, i, a)

#transpose gene intensities for pca

data = t(df3[2:ncol(df3)])
colnames(data) = df3$symbol

# removing sample S2 and S13
row.names.remove <- c("S2", "S13")
data = data[!(row.names(data) %in% row.names.remove), ]

#performing pca using prcomp.
data.pca = prcomp(log(data), center = T, scale = T)

df_out <- as.data.frame(data.pca$x)
df_out$sIdx <- sapply( strsplit(as.character(row.names(data)), "_"), "[[", 1 )

#Reading the second input file
data1 = read.csv(args[2], header=T)

if (is.element('FALSE', unique(row.names(data)) %in% unique(data1$sIdx) )) 
{
  stop("Metadata does not match the Gene intensity file. Make sure that samples name are similar.\n The sample names should have a header sIdx in the Metadata file", call.=FALSE)
}

if (!is.element('sIdx', colnames(data1)))
{
  stop ("The second column of Metadata file should be sample names with a header sIdx.\n", call.=FALSE)
}

#adding metadata to principal components
d2 = merge(x = df_out, y = data1, by='sIdx')
rm(data1, df_out)
  
data.pca_r <- as.data.frame(data.pca$rotation)
data.pca_r = data.pca_r[order(-abs(data.pca_r$PC1)),]
data.pca_r$feature <- row.names(data.pca_r)
PCAloadings = head(data.pca_r[1:2], n = 3)
PCAloadings$symbol = row.names(PCAloadings)



xmax = max(d2$PC1)
ymax = max(d2$PC2)
x_r = (xmax-0.1*xmax)/PCAloadings$PC1[1]
y_r = (ymax-0.1*ymax)/PCAloadings$PC2[1]
x_r1 = (xmax-0.2*xmax)/PCAloadings$PC1[1]
y_r1 = (ymax-0.2*ymax)/PCAloadings$PC2[1]
x_r2 = (xmax-0.03*xmax)/PCAloadings$PC1[1]
y_r2 = (ymax-0.03*ymax)/PCAloadings$PC2[1]

tiff(filename = "PCA.tiff", width=1400, height=1400, compression="lzw", res = 300)

ggplot(d2,aes(x=PC1,y=PC2))  + geom_segment(data = PCAloadings, aes(x = PC1*x_r1, y = PC2*y_r1, xend = 
  PC1*x_r, yend = PC2*y_r), arrow = arrow(length = unit(0.5, "cm"), type = "closed"), color = "black") +
  geom_point(aes(color=factor(Time)), size = 4.2, pch = 1) +
  geom_point(size = 4, aes(color=factor(Time), alpha=0.4)) + 
  scale_color_manual(values= Colors) + guides(alpha=FALSE) +
  theme(legend.position = c(0.1, 0.67), panel.background = element_rect(fill = "white", 
  colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.ticks.length = unit(.25, "cm"), axis.text = element_text(colour = "black"), 
  legend.key = element_rect(colour = NA, fill = NA), legend.title=element_blank()) +
  annotate("text", x = PCAloadings$PC1*x_r2, y = PCAloadings$PC2*y_r2,
  label = PCAloadings$symbol)

dev.off()

