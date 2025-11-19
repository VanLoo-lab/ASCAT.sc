sc_plot_heatmap <- function(res,
                            outdir,
                            projectname)
{
    Nchr = length(res$chr)
    ## ##########################################################
    require(dplyr)
    require(tidyverse)
    require(ComplexHeatmap)
    require(magrittr)
    require(fastcluster)
    ## ##########################################################
    ## ##########################################################
    ## plotting functions total CN by Yixiao Cheng
    ## inspired by Haixi Yan's allele-specific CN functions below
    ## ##########################################################
    sc_totCN_heatmap <- function(CN_mtx,
                                 km = NULL,
                                 row_split = NULL,
                                 probes,
                                 row_ann = NULL, max_CN = 10,
                                 column_title = "Single Cell Copy Number Heatmap",
                                 title = "CN",
                                 colour_scheme = c("#5CA4A9", "#9BC1BC", "#FFFFFF",
                                                   "#EDEECE", "#F4F1BB", "#F1AE8B",
                                                   "#EF8C73", "#D78C96", "#BC69A0",
                                                   "#B24C9E", "#A72F9B"))
    {
        show_raw_dend = T
        if (!is.null(km))
        {
            myhclust = T
            show_raw_dend = F
        } #Cluster each kmeans cluster but don't show dendrogram
        CN_mtx[CN_mtx<0] <- 0 #set min CN value
        CN_mtx[CN_mtx>max_CN] <- max_CN #set max CN value
        names(colour_scheme) <- 0:max_CN #Colours
                                        #Annotation
                                        # Create a vector for each bin: a for even chr; b for odd chr
        pattern_vector <- character(ncol(CN_mtx))  # Create an empty character vector
        repeat_pattern <- c("a", "b")
        current_index <- 1
        for (i in seq_len(nrow(tenX_chr_probes)-1)) {
            repeat_length <- tenX_chr_probes$total.probes[-1][i]
            repeat_value <- rep(repeat_pattern[i %% 2 + 1], repeat_length)
            pattern_vector[current_index:(current_index + repeat_length - 1)] <- repeat_value
            current_index <- current_index + repeat_length
        }
        colnames(CN_mtx) <- pattern_vector
        ha_column <- HeatmapAnnotation(
            chr = colnames(CN_mtx),
            col = list(chr = c("a" = "#5F618E", "b" = "#E5D4ED")),
            chr1 = anno_empty(border = F),
            show_legend = F)
                                        #Chromosome annotation at bottom
        myhclust = function(x) fastcluster::hclust(dist(x,method = "manhattan"))
        print(ComplexHeatmap::Heatmap(matrix = CN_mtx, cluster_rows = myhclust,
                                      row_km = km, row_split = row_split,
                                      cluster_row_slices = FALSE, show_row_dend = show_raw_dend,
                                      row_title_rot = 0, cluster_columns = F,
                                      show_row_names = F, show_column_names = F,
                                      bottom_annotation = ha_column, left_annotation = row_ann,
                                      col = colour_scheme, column_title = column_title, name = title))
                                        #Add chr lines
        probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
        decorate_heatmap_body(heatmap = title, {
            for (k in 1:(length(probes))) {
                grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="grey", lty = 1, lwd = 1))
            }
        }) #Chr lines
        if (!is.null(row_split))
        {
            for (i in 2:length(unique(row_split)))
            {
                decorate_heatmap_body(heatmap = title, row_slice = i,
                {
                    for (k in 1:(length(probes)))
                    {
                        grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="grey", lty = 1, lwd = 1.5))
                    }
                }) #Chr lines
            }
        }
        if (!is.null(km))
        {
            for (i in 2:K)
            {
                decorate_heatmap_body(heatmap = title,
                                      row_slice = i,
                                      {
                                          for (k in 1:(length(probes)))
                                          {
                                              grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
                                          }
                                      }) #Chr lines
            }
        }
        decorate_annotation(annotation = "chr1",
        {
            for (k in 2:(length(probes)))
            {
                {
                    grid.text(names(probes)[k],
                              x=(probes[k-1]/ncol(CN_mtx)) + 0.001, y=0.5,
                              just = "left",
                              gp=gpar(frontsize = 2))
                }
            }
        }) #Chr numbers
############################################################    ​
    }
    ## ##########################################################
    tenX_probes <- (do.call(rbind,
                            lapply(res$allTracks.processed[[1]][["lSegs"]], function(s) return(s[["output"]]))) %>%
                    select(chrom, num.mark, loc.start, loc.end) %>%
                    dplyr::rename(n.probes = num.mark, startpos = loc.start, endpos = loc.end) %>%
                    filter(n.probes>1)) #Remove segments of just 1/2 probes(these are NAs in ascat.sc outputs)
############################################################​
    tenX_chr_probes <- (rbind(data.frame(chrom = 0, total.probes = 0),
                              tenX_probes[,1:2] %>%
                              group_by(chrom) %>%
                              summarise(total.probes = sum(n.probes))) %>%
                        mutate(cum.probes = cumsum(total.probes)))
############################################################​
    scDNA_CN_mtx <- do.call(rbind, lapply(1:length(res$allProfiles), function(b)
    {
        rep(na.omit(as.numeric(res$allProfiles[[b]][,"total_copy_number"])), tenX_probes$n.probes)
    }))
    ## ##########################################################​
    ha_row <- NULL
    ## ##########################################################​
    pdf(paste0(outdir,"/",projectname,"_heatmap_tot.pdf"))
    (sc_totCN_heatmap(CN_mtx = scDNA_CN_mtx,
                     probes = tenX_chr_probes$cum.probes[-1] %>%
                         set_names(nm = c(1:22, "X","Y")[1:Nchr]),
                     title = "Total CN"))
    dev.off()
############################################################
    ## ##########################################################
    if(any(grepl("allProfiles_AS",names(res))))
    {
        ## ##########################################################
        ## plotting functions allele-specific CN by Haixi Yan
        ## ##########################################################
        scDNA_asCN_mtx <- do.call(rbind, lapply(which(lengths(res[["allProfiles_AS"]])!=0), function(b)
        {
            prof <- res$allProfiles_AS[[b]]$nprof.fixed[!is.na(res$allProfiles_AS[[b]]$nprof.fixed$total_copy_number),] #Using fixed
            prof$nA[is.na(prof$nA)] <- 0
            prof$nB[is.na(prof$nB)] <- prof$total_copy_number[which(is.na(prof$nB))]
            asCN_df <- data.frame(allele1 = prof$nA,
                                  allele2 = prof$nB) %>% #flip alleles so A is larger
                mutate(alleleA = ifelse(allele1>=allele2, allele1, allele2),
                       alleleB = ifelse(allele1>=allele2, allele2, allele1)) %>%
                mutate(alleleT = alleleA + alleleB) %>% #get total CN
                mutate(plotT = alleleT*10+alleleB) %>% select(-c(1:2))#Convert to matrix format
            rep(na.omit(as.numeric(asCN_df$plotT)), tenX_probes$n.probes)#asCN_probes$n.probes) #NAs are omitted...
        })) %>% set_rownames(names(res[["allProfiles_AS"]]))
        ## ##########################################################
        sc_asCN_heatmap <- function(CN_mtx,
                                    hclust = T,
                                    km = NULL,
                                    row_split = NULL,
                                    probes,
                                    row_ann = NULL,
                                    column_title = "Single Cell Copy Number Heatmap", title = "CN",
                                    colour_scheme = c("royalblue3", "skyblue2", "grey80", "white",
                                                      "gold1", "khaki1", "darkorange3", "darkorange1",
                                                      "orange", "red4", "red", "orangered2", "purple4"))
        {
            show_raw_dend = T
            if (!is.null(km))
            {
                hclust = T
                show_raw_dend = F
            } #Cluster each kmeans cluster but don't show dendrogram
            CN_mtx[CN_mtx<0] <- 0 #set min CN value
            CN_mtx[CN_mtx>60] <- 60 #set max CN value
            names(colour_scheme) <- c(0, 10, 20, 21, 30, 31, 40, 41, 42, 50, 51, 52, 60)#Colours
            CN_labels <- c("0+0", "1+0", "2+0",
                           "1+1", "3+0", "2+1",
                           "4+0", "3+1", "2+2",
                           "5+0", "4+1", "3+2", ">6")
            names(CN_labels)=names(colour_scheme)
            CN_labels=CN_labels[which(names(colour_scheme) %in% sort(unique(as.numeric(CN_mtx))))] #Get labels in legend
            ha_column = HeatmapAnnotation(odd = anno_empty(border = F),
                                          even = anno_empty(border = F)) #Chromosome annotation at bottom
            print(ComplexHeatmap::Heatmap(matrix = CN_mtx,
                                          cluster_rows = hclust,
                                          row_km = km,
                                          row_split = row_split,
                                          cluster_row_slices = FALSE,
                                          show_row_dend = show_raw_dend,
                                          row_title_rot = 0,
                                          cluster_columns = F,
                                          show_row_names = F,
                                          show_column_names = F,
                                          heatmap_legend_param = list(at=as.numeric(names(CN_labels)), labels = CN_labels),
                                          bottom_annotation = ha_column,
                                          left_annotation = row_ann,
                                          col = colour_scheme,
                                          column_title = column_title, name = title))
            probes <- c(0, probes) #Adds posiiton zero for start of first chromosome
            decorate_heatmap_body(heatmap = title, {
                for (k in 1:(length(probes))) {
                    grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1), gp=gpar(col="black", lty = 1, lwd = 1.5))
                }
            }) #Chr lines
            if (!is.null(row_split)) {
                for (i in 2:length(unique(row_split))) {
                    decorate_heatmap_body(heatmap = title, row_slice = i, {
                        for (k in 1:(length(probes))) {
                            grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
                        }
                    }) #Chr lines
                }
            }
            if (!is.null(km)) {
                for (i in 2:K) {
                    decorate_heatmap_body(heatmap = title, row_slice = i, {
                        for (k in 1:(length(probes))) {
                            grid.lines(x=probes[k]/ncol(CN_mtx), y=c(0,1.1), gp=gpar(col="black", lty = 1, lwd = 1.5))
                        }
                    }) #Chr lines
                }
            }
            decorate_annotation(annotation = "odd", {
                for (k in 2:(length(probes))) {
                    if (k%%2 == 0) {grid.text(names(probes)[k], x=(probes[k-1]/ncol(CN_mtx)) + 0.001,
                                              y=0.5, just = "left", gp=gpar(cex = 2))}
                }
            }) #Chr numbers
            decorate_annotation(annotation = "even", {
                for (k in 2:(length(probes))) {
                    if (k%%2 == 1) {grid.text(names(probes)[k], x=(probes[k-1]/ncol(CN_mtx)) + 0.001,
                                              y=0.5, just = "left", gp=gpar(cex = 2))}
                }
            }) #Chr numbers
        }
        ##Function to calculate allele specific distance
        ascn_dist <- function(cna1, cna2)
        {
            ## ind1 <- is.na(cna1[,"nA"]) | is.na(cna2[,"nA"])
            ## ind2 <- !ind1
            ind1 <- (is.na(cna1[,"nA"]) | is.na(cna2[,"nB"])) #Edited as one segment has NA
            ind2 <- !(is.na(cna1[,"nA"]) | is.na(cna2[,"nB"]))
            sizes <- cna1[,"endpos"]-cna1[,"startpos"]
            sizes <- sizes/500000/sum(sizes/500000)
            dists <- sum(pmin(2,abs(cna1[ind2,"nA"]-cna2[ind2,"nA"])+abs(cna1[ind2,"nB"]-cna2[ind2,"nB"]))*sizes[ind2])
            ## dists <- dists+sum(pmin(2,abs(cna1[ind1,"fitted"]-cna2[ind1,"fitted"])*sizes[ind1]))
            dists
        }
        print("Caculating all pairwise sample distances in AS space for heatmap")
        alldists <- matrix(NA,length(res[["allProfiles_AS"]]),length(res[["allProfiles_AS"]]))
        for(i in 1:(length(res[["allProfiles_AS"]])-1))
        {
            cat(".")
            for(j in (i+1):length(res[["allProfiles_AS"]]))
                alldists[i,j] <- ascn_dist(res[["allProfiles_AS"]][[i]][["nprof.fixed"]],
                                           res[["allProfiles_AS"]][[j]][["nprof.fixed"]])
        }
        rownames(alldists) <- names(res[["allProfiles_AS"]])
        colnames(alldists) <- names(res[["allProfiles_AS"]])
        alldists.dist <- as.dist(t(alldists))
        asCN_hclust_ward <- hclust(alldists.dist, method = "ward.D2")
        ##Plot CN profiles
        ha_row <- NULL
        png(filename = paste0(outdir,"/",
                              projectname,"_all_asfixed_asCN_heatmap.png"),
            width = 4000,
            height = 4000,
            res = 200)
        print(sc_asCN_heatmap(CN_mtx = scDNA_asCN_mtx,
                              hclust = asCN_hclust_ward,
                              row_ann = NULL, #ha_row,
                              probes = tenX_chr_probes$cum.probes[-1] %>% set_names(nm = c(1:22, "X","Y")[1:Nchr])))
        dev.off()
    }
}
