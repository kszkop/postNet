######
goDotplot <- function(ptn,
                      #goIn,
                      category,
                      pool = TRUE,
                      termSel = NULL,
                      #colours = NULL,
                      nCategories=10,
                      size = 'Count',
                      pdfName=NULL){
  #
  check_ptn(ptn)
  if(!check_logical(pool)){
    stop("'pool' must be logical: TRUE of FALSE ")
  }
  if(!check_number(nCategories)) {
    stop("Please provide a postive integer for 'nCategories'")
  }
  check_size(size)
  #
  check_category(category)
  #
  for(sel in category){
    if(is.null(slot(ptn@analysis@GO, sel))){
      stop("No results found for the selected category. Please run the goAnalysis function 
           first for all of the selected categories.")
    } else {
      goIn <- slot(ptn@analysis@GO, sel)
    }
    #
    if(isTRUE(pool)){
      nameOut<- ifelse(is.null(pdfName), paste0("pooled_GOdotplot", "_", sel,".pdf"), 
                       paste(pdfName, paste0("pooled_GOdotplot", "_", sel,".pdf"), sep = "_"))
      #
      goDf <- data.table::rbindlist(lapply(goIn, function(x) x@result),use.names=TRUE, idcol=TRUE)
      colnames(goDf)[1] <- 'regulation'
      #
      if(!is.null(termSel)){
        goDf <- goDf[goDf$ID %in% termSel,]
      }
      #
      idx <- order(goDf$p.adjust, decreasing = FALSE)
      #
      goDf <- goDf[idx,]
      if(nCategories<nrow(goDf)){
        goDf <- goDf[1:nCategories,]
      }
      
      #to plot
      goDf$log10fdr <- -log10(goDf$p.adjust)
      
      #
      if(size=='geneRatio'){
        goDf$scale <- goDf$Count/goDf$Size
      } else {
        goDf$scale <- goDf$Count
      }
      
      #
      colOut <- colPlot(ptn)[-1]
      names(colOut) <- names(ptn@dataIn@geneList)
      #
      if (all(goDf$scale == floor(goDf$scale))) {
        # Integer scale values
        size_breaks <- pretty(range(goDf$scale, na.rm = TRUE), n = 4)
        size_breaks <- size_breaks[size_breaks == floor(size_breaks)]
        size_limits <- range(size_breaks, na.rm = TRUE)
      } else {
        # Decimal scale values between 0 and 1
        rng <- range(goDf$scale, na.rm = TRUE)
        step <- 0.05  # finer step for better resolution in small ranges
        
        min_val <- floor(rng[1] / step) * step
        max_val <- ceiling(rng[2] / step) * step
        
        size_breaks <- seq(min_val, max_val, by = step)
        size_limits <- c(min_val, max_val)
      }
      #
      pdf(nameOut, width = 8, height = 8, useDingbats = F)
      par(mar = c(5, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
      pOut <- ggplot2::ggplot(goDf, ggplot2::aes(x=log10fdr, y=reorder(Description,log10fdr), size=scale,colour = regulation)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = colOut) +  
        ggplot2::scale_size_continuous(name = size,
                                       breaks = size_breaks,
                                       limits = size_limits,
                                       range = c(2, 8)) + 
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 'dashed', linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.5, 'cm')) +   
        ggplot2::xlab('-log10 FDR') +
        ggplot2::scale_x_continuous(limits = c(0, NA)) +
        ggplot2::ylab(" ") 
      #+
      #ggplot2::ggtitle(paste(names(goIn)[i],sel, sep='_'))
      plot(pOut)
      dev.off()
    } else{
      for(i in 1:length(goIn)){
        nameOut <- ifelse(is.null(pdfName), paste(names(goIn)[i], paste0("GOdotplot", "_", sel,".pdf"),sep='_'),
                          paste(pdfName, names(goIn)[i], paste0("GOdotplot", "_", sel,".pdf"), sep = "_"))
        
        goDf <- goIn[[i]]@result
        if(nrow(goDf) == 0){
          message(paste('For the geneListL: ', names(goIn)[i], ' there are no categories to plot.', sep=''))
        } else {
          if(!is.null(termSel)){
            goDf <- goDf[goDf$ID %in% termSel,]
          }
          #
          idx <- order(goDf$p.adjust, decreasing = FALSE)
          #
          goDf <- goDf[idx,]
          if(nCategories<nrow(goDf)){
            goDf <- goDf[1:nCategories,]
          }
          #to plot
          goDf$log10fdr <- -log10(goDf$p.adjust)
        
          # 
          if(size=='geneRatio'){
            goDf$scale <- goDf$Count/goDf$Size
          } else {
            goDf$scale <- goDf$Count
          }
        
          #
          colOut <- colPlot(ptn)[-1]
          names(colOut) <- names(ptn_geneList(ptn))
          #
          # Determine breaks and limits for the size scale
          if (all(goDf$scale == floor(goDf$scale))) {
            # Integer scale values
            size_breaks <- pretty(range(goDf$scale, na.rm = TRUE), n = 4)
            size_breaks <- size_breaks[size_breaks == floor(size_breaks)]
            size_limits <- range(size_breaks, na.rm = TRUE)
          } else {
            # Decimal scale values between 0 and 1
            rng <- range(goDf$scale, na.rm = TRUE)
            step <- 0.05  # finer step for better resolution in small ranges
            
            min_val <- floor(rng[1] / step) * step
            max_val <- ceiling(rng[2] / step) * step
            
            size_breaks <- seq(min_val, max_val, by = step)
            size_limits <- c(min_val, max_val)
          }
          #
          pdf(nameOut, width = 8, height = 8, useDingbats = F)
          par(mar = c(5, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
          pOut <- ggplot2::ggplot(goDf, ggplot2::aes(x=log10fdr, y=reorder(Description,log10fdr), size=scale)) +
            ggplot2::geom_point(color= colOut[names(goIn)[i]]) +
          #  ggplot2::scale_color_manual(values = colOut[names(goIn)[i]]) +   
            ggplot2::theme_bw() +
            ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 'dashed', linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.5, 'cm')) +   
            ggplot2::xlab('-log10 FDR') +
            ggplot2::scale_size_continuous(name = size,
                                           breaks = size_breaks,
                                           limits = size_limits,
                                           range = c(2, 8)) + 
            ggplot2::ylab(" ") +
            ggplot2::ggtitle(paste(names(goIn)[i],sel, sep='_')) 
          print(pOut)
          dev.off()
        }
      }
    }
  }
}
