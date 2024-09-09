# Vietnam
Here we present the NGS metabarcoding characterization of prokaryotic extremely halophilic communities, thriving in the salt crystallizer ponds of the Hon Khoi solar saltwork fields (HKsf), South Vietnam.

# Enigma-lake

## DADA2- https://benjjneb.github.io/dada2/tutorial.html

	library(dada2)
	library(microbiome) # data analysis and visualisation
	library(phyloseq) # also the basis of data object. Data analysis and visualisation
	library(RColorBrewer) # nice color options
	library(dplyr) # data handling
	library(network) # networks
	library(intergraph) # networks
	library(ggnet)  # network plotting with ggplot
	library(igraph) # networks
	library(phyloseq) # ASV ecological analysis package
	library(ggplot2) # plotting library
	library(gridExtra) # gridding plots
	library(ape) # importing and handling phylogenetic trees
	library(ggthemes) # additional themes fro ggplot2
	library(magrittr)
	library(rioja) # plotting poackages for tabular bubbleplots
	library(ggpubr)
	library(ggtern)
	library(plyr)
	library(coda.base)
	library(vegan)
	library(propr)
	library(msa)
	library(phangorn)

	path <- "//your path/" # directory containing the fastq files after unzipping.
	list.files(path)

### Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
	fnFs <- sort(list.files(path, pattern="_R1_trimmed.fastq", full.names = TRUE))
	fnRs <- sort(list.files(path, pattern="_R2_trimmed.fastq", full.names = TRUE))

### Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
	sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) 
	ls()

### Inspect read quality profiles
	plotQualityProfile(fnFs[1:2])
	plotQualityProfile(fnRs[1:2])

### Place filtered files in filtered/ subdirectory
	filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
	filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
	names(filtFs) <- sample.names
	names(filtRs) <- sample.names

	out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, , #change truncLen value depending on qprofile c(155,145)
                     maxN=0, truncQ=0, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
	head(out)
	
    				reads.in	reads.out
	# H1_R1_trimmed.fastq		  31218		31202
	# H10_R1_trimmed.fastq		  55941		55917
	# H11_R1_trimmed.fastq		  54938		54907
	# H14_R1_trimmed.fastq		  44480		44447
	# H15_R1_trimmed.fastq		  50093		50066
	# H16_R1_trimmed.fastq		  147518	146693

### Learn the Error Rates

	errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
	errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)

	plotErrors(errF, nominalQ=TRUE)

	derepFs <- derepFastq(filtFs, verbose=TRUE)
	derepRs <- derepFastq(filtRs, verbose=TRUE)

	names(derepFs) <- sample.names
	names(derepRs) <- sample.names

### Sample Inference
	dadaFs <- dada(derepFs, err=errF, pool="pseudo", multithread=TRUE)
	dadaRs <- dada(derepRs, err=errR, pool="pseudo", multithread=TRUE)

### Merge paired reads
	mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

### Inspect the merger data.frame from the first sample
	head(mergers[[1]])

### Construct sequence table
	seqtab <- makeSequenceTable(mergers)
	dim(seqtab)
### Inspect distribution of sequence lengths
	table(nchar(getSequences(seqtab)))

### Remove chimeras
	seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
	dim(seqtab.nochim)

	sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline
	getN <- function(x) sum(getUniques(x))
	track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
	colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
	rownames(track) <- sample.names
	head(track)
	write.csv(track, "summary_chim_clean.csv")
	head(track)

	#	input		filtered	denoisedF	denoisedR	merged		nonchim
	#H1	31218		31202		30327		30476		20968		14139
	#H10	55941		55917		54348		54858		46377		34852

### Assign taxonomy
	taxa <- assignTaxonomy(seqtab.nochim, "//home/ your path/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
	taxa <- addSpecies(taxa, "//your path/silva_species_assignment_v138.1.fa.gz")

	taxa.print <- taxa # Removing sequence rownames for display only
	rownames(taxa.print) <- NULL
	head(taxa.print)

	save.image()

### Write table in .csv format
	write.table(t(seqtab.nochim), "seqtab-nochim2.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
	uniquesToFasta(seqtab.nochim, fout='rep-seqs2.fna', ids=colnames(seqtab.nochim))
	write.csv(taxa.print, "taxa_print2.csv")

### Make the phylogenetic tree
	seqs <- getSequences(seqtab.nochim)
	names(seqs) <- seqs # This propagates to the tip labels of the tree
	mult <- msa(seqs, method="ClustalW", type="dna", order="input")
	save.image()


	library(phangorn)
	phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
	dm <- dist.ml(phang.align)
	treeNJ <- NJ(dm) # Note, tip order != sequence order
	fit = pml(treeNJ, data=phang.align)
	fitGTR <- update(fit, k=4, inv=0.2)
	fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
	fit
	save.image()

### Construct a phyloseq object
	prok_sample <- read.csv("dataset_env_2_mod.csv", header=T, sep=",", row.names=1)
	summary(prok_sample)
	head(prok_sample) 
	prok_data <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), 
                      phy_tree(treeNJ), 
                      tax_table(taxa), 
                      sample_data(prok_sample))
	prok_data
	readcount(prok_data)
	write_phyloseq(prok_data, type = "all")
	save.image()


### Eliminate contaminates
	Viet_no_cont <- subset_taxa(prok_data,  (Genus != "Acidovorax") | is.na(Genus))
 	Viet_no_cont <- subset_taxa(Viet_no_cont,  (Genus != "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium") | is.na(Genus))  


### Remove sample from phyloseq object
	vietnam2024 = subset_samples(Viet_no_cont, ID != "38") #change name of phyloseq object

# Alpha diversity estimates

### Write alpha diversity table
	test<-estimate_richness(vietnam2024, split = TRUE, measures = NULL)
	write.csv(test, "diversity table.csv")#

### Alpha diversity estimates
	plot_richness(vietnam2024, measures=c("Observed", "Chao1", "Shannon", "Simpson", "Fisher"), x="ID", color="s_type")
	
# Beta diversity
### transform OTU table in presence/absence table

	library(vegan)
	library(metagMisc)
	library(microViz)
	presenceAbsencevietnam2024 <- vietnam2024  %>% tax_transform("binary")

### Transform  abundance to relative abundances for plotting
	bac_ra = transform_sample_counts(vietnam2024, function(x){(x / sum(x))*100})
	
 or
 
	bac_ra = transform_sample_counts(vietnam2024, function(x){x / sum(x)})

### PCoA unifrac weighted with relative counts
	prok_data_w <- ordinate(bac_ra, method = "PCoA", distance = "unifrac", weighted=T)
	evals_w <- prok_data_w $values$Eigenvalues
	plot_ordination(bac_ra, prok_data_w, type = "sample", color = "s_type", label="ID", title="PCoA weighted Unifrac") +
  	labs(col = "Type of sample") +
  	coord_fixed(sqrt(evals_w[2] / evals_w[1]))

### PCoA unifrac unweighted with relative counts
	prok_data_un <- ordinate(bac_ra, method = "PCoA", distance = "unifrac", weighted=F)
	evals <- prok_data_un$values$Eigenvalues
	plot_ordination(bac_ra, prok_data_un, type = "sample", color = "s_type", label="ID", title="PCoA unweighted Unifrac") +
	  labs(col = "Type of sample") +
	  coord_fixed(sqrt(evals[2] / evals[1]))


### nMDS with Jaccard and Bray-Curtis distance
	bac_nmds_j <- ordinate(bac_ra, method = "NMDS", distance = "jaccard", weighted=T, trymax=100)
	stressplot(bac_nmds_j)
	bac_nmds_bc <- ordinate(bac_ra, method = "NMDS", distance = "bray", weighted=T, trymax=100)
	stressplot(bac_nmds_bc)
	plot_ordination(bac_ra, bac_nmds_j, color="Samples", label="Samples", title="nMDS Jaccard diversity colored by Area") +
 	 theme_bw()
	plot_ordination(bac_ra, bac_nmds_bc, color="Samples", label="Samples", title="nMDS Bray-Curtis diversity colored by Area") +
 	 theme_bw()

### Plot PCoA and nMDS together

	grid.arrange(nrow = 2, ncol=2,
             	plot_ordination(prok_data2, bac_pcoa_w, type = "Samples", color = "area", label="Samples", title="PCoA weighted Unifrac") +
              	 theme(legend.position = "none") +
              	 coord_fixed(sqrt(evals_w[2] / evals_w[1])),
             	plot_ordination(prok_data2, bac_pcoa_un, type = "Samples", color = "area", label="Samples",title="PCoA unweighted Unifrac") +
              	 theme(legend.position = "none")+
               	coord_fixed(sqrt(evals[2] / evals[1])),
            	 plot_ordination(bac_ra, bac_nmds_j, color="area", label="Samples", title="nMDS Jaccard distances") +  theme(legend.position = "none"),
             	plot_ordination(bac_ra, bac_nmds_bc, color="area", label="Samples", title="nMDS Bray-Curtis distances") +   theme(legend.position = "none")
	)

## Mantel's test
### https://github.com/Hy4m/linkET/blob/master/README.Rmd
### https://cloud.tencent.com/developer/article/2136590



	library(dplyr)
	library(linkET)
	library(readxl) 

	Phy<- read_excel("~/...your path.../Phylum.xlsx")#file with your phylum list
	environmental<- read_excel("...your path.../environmental.xlsx")#file with enviromental parameters

	#The first column becomes the name of the rows
		Phy <- Phy %>%
  	tibble::column_to_rownames("sample") 

	environmental <- environmental %>%
  	tibble::column_to_rownames("samples_id")  


### Data processing

## matrix_data
	matrix_data(list(environmental=environmental))

## md_tbl
	matrix_data(list(environmental = environmental)) %>% 
	as_md_tbl()

## as method
	as_matrix_data(environmental)
	as_md_tbl(environmental)

## special function for correlation matrix
	correlate(environmental) %>% 
		as_matrix_data()

	correlate(environmental) %>% 
  	as_md_tbl()

### Heatmap
	library(ggplot2)
	matrix_data(list(environmental = environmental)) %>% 
  	hyplot(aes(fill = environmental)) +
  	geom_tile()

	as_md_tbl(environmental) %>% 
 	 hyplot(aes(size = environmental)) +
  	geom_point(shape = 21, fill = NA)

	correlate(environmental) %>% 
	as_md_tbl() %>% 
	qcorrplot() +
 	geom_square()

	library(vegan)

	correlate(Phy[1:37], Phy) %>% 
  	qcorrplot() +
  	geom_square() +
  	scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

	qcorrplot(Phy[1:37], type = "lower") +
 		geom_square() +
		scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

## you can set your style
	set_corrplot_style()
	qcorrplot(environmental) + geom_square()

## reset to default style
	set_default_style()

## mantel test
	mantel <- mantel_test(environmental, Phy,
                      	spec_select = list(pH=1,
                                        	T=2,
                                        	S=3,
                                         	Redox=4)) %>% 
  	mutate(rd = cut(r, breaks = c(-Inf, -0.3, 0.3, Inf),
                  	labels = c("<= -0.3", ">= -0.3 - <= 0.3", ">= 0.3")),
         	pd = cut(p, breaks = c(-Inf, 0.011, 0.05, Inf),
                  	labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))



	qcorrplot(correlate(Phy), type = "lower", diag = FALSE) +
  	geom_square() +
  	geom_couple(aes(colour = pd, size = rd), 
              	data = mantel, 
              	curvature = nice_curvature()) +
  	scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  	scale_size_manual(values = c(0.5, 1, 2)) +
  	scale_colour_manual(values = color_pal(3)) +
  	guides(size = guide_legend(title = "Mantel's r",
                             	override.aes = list(colour = "grey35"), 
                             	order = 2),
         	colour = guide_legend(title = "Mantel's p", 
					override.aes = list(size = 3), 
                              order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))

 ## False Discovery Rate (FDR)

	a <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)#nomi campioni
	p_values <- 		c(0.001,0.002,0.006,0.009,0.011,0.011,0.011,0.024,0.036,0.038,0.046,0.048,0.040)
	df <- data.frame(a, p_values)
	df$FDR <- p.adjust(df$p_values, method = "BH")


  # Bar plot dendrogram

	library(phyloseq)
	library(ggtree)
	library(ggplot2)
	library(reshape)
	library(readxl)
	library(tidyverse)
	library(phangorn)
	library(vegan)
	library(ggdendro)
	library(dendextend)
	library(ggsci)
	library(cowplot)


### Import data from excel 

	Kingdom<- read_excel("~/your path")
	PhylumB<- read_excel("~/your path")
	PhylumA<- read_excel("~/your path")

### First column becomes row names
	Kingdom <- Kingdom %>%
  	tibble::column_to_rownames("Kingdom") 

	PhylumB <- PhylumB %>%
  	tibble::column_to_rownames("Phylum") 

	PhylumA <- PhylumA %>%
  	tibble::column_to_rownames("Phylum") 

### Generate example data
	set.seed(500)
	combined_matrix <- data.frame(PhylumA)
	row.names(combined_matrix) <- paste0("s", seq(1,12)) #replaces the names with S1, S2, ...etc)



### Dendrogram
	UniFrac(your phyloseq file, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
	(y <- UniFrac(your phyloseq file, TRUE))

### Create a distance object
	unifrac_dist <- as.matrix(y)

### Cluster the distance matrix
	hc <- hclust(as.dist(unifrac_dist))
	dendro <- as.dendrogram(hc)

### Plot the dendrogram
	plot(dendro)
	plot(dendro, horiz = T)#plot horizontal

### calculate UPGMA tree with phangorn::upgma() and convert to dendrogram
	dendUPGMA <- dendro # if you want to associate the dendrogram of the DADA2 data

	dendUPGMA <- as.dendrogram(upgma(dm))#if you want to associate other data

	plot_dendro_bars_v <- function(df, dend, taxonomy) {
 	 #convert dendrogram to segment data
 	 dend_data <- dendro_data(dend, type="rectangle")
 	 segment_data <- dend_data[["segments"]]
 	 #sample positions df
 	 sample_pos_table <- with(dend_data$labels, 
                          	 data.frame(x_center = x, sample = as.character(label), width = 0.9))
	  #prepare input data
 	 ptdf <- rownames_to_column(df, var = "sample") %>%
   	 pivot_longer(-sample, names_to = taxonomy, values_to = "Frequency") %>%
   	 group_by(sample) %>%
   	 mutate(Frequency = Frequency/100,
          	 ymax = cumsum(Frequency/sum(Frequency)),
          	 ymin = ymax - Frequency/sum(Frequency),
          	 y_center = ymax-(Frequency/2)) %>%
   	 left_join(sample_pos_table) %>%
  	  mutate(xmin = x_center-width/2,
          	 xmax = x_center+width/2)
	  #plot stacked bars
 	 axis_limits <- with(sample_pos_table, 
                      	c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))) + 
   	 0.1 * c(-1, 1) # extra spacing: 0.1
 	 plt_hbars <- ggplot(ptdf, 
                      	aes_string(x = "x_center", y = "y_center", fill = taxonomy, xmin = "xmin", xmax = "xmax",
                                	 height = "Frequency", width = "width")) + 
    	geom_tile() +
   	 geom_rect(ymin = 0, ymax = 1, color = "black", fill = "transparent") +
   	 scale_fill_rickandmorty() +
  	  scale_fill_manual(values=c("darkviolet", "deepskyblue4", "gold", "orange", "springgreen4", "hotpink1", "darkturquoise", "green", "blue4", "gold3", "darkolivegreen1", "pink", "cornflowerblue", "beige", "chartreuse3", "grey", "magenta", "red4", "red", "black", "yellow3", "darkorange2", "plum4", "lightgoldenrod1", "brown2","cornsilk2", "coral4", "chocolate2","chartreuse4", "cadetblue", "burlywood", "bisque4", "azure3", "aquamarine4", "dodgerblue", "gray45", "lemonchiffon2", "mediumorchid", "khaki3", "palegreen2", "seashell3", "violetred3", "indianred1", "lightslategray", "orange2", "thistle", "olivedrab", "moccasin", "mistyrose2", "grey70", "firebrick", "cyan3", "gold3", "lightpink2", "forestgreen", "ghostwhite", "slateblue1", "peru", "gold4", "orange1", "gray63", "mediumspringgreen", "goldenrod2", "lightcyan1", "purple4", "grey92", "firebrick", "plum2", "snow2", "steelblue", "darkorchid", "blue1", "darkolivegreen", "yellow4", "tan2", "rosybrown", "gray77", "darkorange4", "lightpink", "steelblue4", "darkkhaki", "mediumturquoise", "peachpuff3", "gray46", "antiquewhite", "cyan4", "gold3", "grey80", "lightblue", "mediumorchid1", "springgreen2", "lightgoldenrod", "midnightblue", "grey67", "purple4", "goldenrod", "azure2", "darkolivegreen3", "lightpink3", "navajowhite3", "grey50", "gray2", "darkorchid1", "darkorange4", "oldlace", "thistle", "tomato3", "lightblue4", "lightyellow2", "seagreen", "violetred2", "chocolate1", "gray32", "limegreen", "linen", "sienna1", "yellowgreen", "grey63", "palevioletred", "lawngreen", "gray39", "deepskyblue", "yellow3")) +
   	 scale_y_continuous(expand = c(0, 0)) + 
  	  # For the y axis, alternatively set the labels as: gene_position_table$gene
   	 scale_x_continuous(breaks = sample_pos_table[, "x_center"], 
                     	  labels = sample_pos_table$sample,
                      	 limits = axis_limits, 
                      	 expand = c(0, 0)) + 
   	 labs(x = "", y = "Frequency") +
   	 theme_bw() +
   	 theme(# margin: top, right, bottom, and left
    	  plot.margin = unit(c(-0.9, 0.2, 1, 0.2), "cm"), 
     	 panel.grid.minor = element_blank())
	  #plot dendrogram
	  plt_dendr <- ggplot(segment_data) + 
  	  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  	  scale_y_continuous(expand = c(0, 0.05)) + 
  	  scale_x_continuous(breaks = sample_pos_table$x_center, 
                   	    labels = rep("", nrow(sample_pos_table)), 
                   	    limits = axis_limits, 
                    	   expand = c(0, 0)) + 
  	  labs(x = "", y = "Distance", colour = "", size = "") +
 	   theme_bw() + 
  	  theme(panel.grid.minor = element_blank(),
        	  panel.grid.major = element_blank())
	  #combine plots
	  comb <- plot_grid(plt_dendr, plt_hbars, align = 'v', ncol = 1, axis = "lr", rel_heights = c(0.3, 1))
	  comb
	}

	plot_dendro_bars_v(df = combined_matrix, dend = dendUPGMA, taxonomy = "Phylum")

