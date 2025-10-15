plotSignatures <- function(ptn,
                           signatureList,
                           signature_colours=NULL,
                           dataName, 
                           generalName,
                           xlim = c(-2,2), 
                           tableCex = 0.7,
                           pdfName = NULL){
  #
  check_ptn(ptn)
  check_geneList(signatureList)

  signNames <- names(signatureList)
  if(!is.null(signature_colours)){
    if (!is.character(signature_colours) || !length(signature_colours)== length(signatureList)) {
      stop("'signature_colours' should be a character vector of the same length as number of signatures in signatureList. These colours will be used for plotting.")
    }
  } else {
    signature_colours <- paste0("#", sprintf("%06X", sample(0:16777215, length(signatureList))))
  }
  if(is.null(dataName) | is.null(generalName)){
    stop('Please provide name describing signatures (generalName) and descibing data used (dataName)')
  }

  if(!check_number(tableCex)){
    stop("please provide number for tableCex to scale size of table with statistics")
  }
  
  if(!is_numeric_vector(xlim) | length(xlim) != 2){
    stop("please provide numeric vector for xlim. Exactly too numbers")
  }
  #
  effIn <- ptn_effect(ptn)
  regData <- data.frame(geneSymb = names(effIn))
  regData$effIn <- as.numeric(effIn)
  
  if(any(duplicated(unlist(signatureList)))){
    cat("There are some genes that overlap between signatures. Separate background for each will be used")
    
    ##collect signatures
    for(i in 1:length(signatureList)){
      regData[,2+i] <- 'bkg'
      regData[,2+i][regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
    }
  } else {
    regData$signature <- 'bkg'
    for(i in 1:length(signatureList)){
      regData$signature[regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
    }
  }
  
  ##
  tableOut <- matrix(NA, nrow= length(signNames), ncol= 5)
  colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
  tableOut[,1] <- signNames
  
  if(any(duplicated(unlist(signatureList)))){
    tmpBgOut <- list()
    for(i in 1:length(signatureList)){
      #
      tmpBg <- sort(as.numeric(regData[regData[,(2+i)]=='bkg','effIn']))
      ecdfBg <- 1:length(tmpBg)/length(tmpBg)
      bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
      bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
      bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
      #
      tmpBgOut[[i]] <- tmpBg
      #
      tableOut[i,2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData[,(2+i)]==signNames[i],'effIn']),as.numeric(regData[regData[,(2+i)]=='bkg','effIn']),alternative='two.sided')[3]),scientific = TRUE,digits=2)
      
      #C
      tmpSign <- sort(as.numeric(regData[regData[,(2+i)]==signNames[i],'effIn']))
      ecdfSign <- 1:length(tmpSign)/length(tmpSign)
      
      tableOut[i,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
      tableOut[i,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
      tableOut[i,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
      
    }   #
  } else {
    #
    tmpBg <- sort(as.numeric(regData[regData$signature=='bkg','effIn']))
    ecdfBg <- 1:length(tmpBg)/length(tmpBg)
    bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
    bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
    bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
    
    #
    for(i in 1:length(signatureList)){
      tableOut[i,2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData$signature==signNames[i],'effIn']),as.numeric(regData[regData$signature=='bkg','effIn']),alternative='two.sided')[3]),scientific = TRUE,digits=2)
      
      #
      tmpSign <- sort(as.numeric(regData[regData$signature==signNames[i],'effIn']))
      ecdfSign <- 1:length(tmpSign)/length(tmpSign)
      
      tableOut[i,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
      tableOut[i,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
      tableOut[i,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
    }
  }
  #
  pdf(ifelse(is.null(pdfName),
             paste0('data_', dataName, '_signature_', generalName, '.pdf'),
             paste0(pdfName, '_data_', dataName, '_signature_', generalName, '.pdf')),
      width = 8, height = 8, useDingbats = FALSE)
  
  layout(matrix(c(1,2), ncol=1), heights=c(25,100))
  
  # 
  if (!is.null(xlim)) {
    xmin <- xlim[1]
    xmax <- xlim[2]
  } else {
    IQRv <- stats::quantile(regData$effIn, 0.75) - stats::quantile(regData$effIn, 0.25)
    xmin <- as.numeric(postNet::roundNice(stats::quantile(regData$effIn,0.25) - 1.5*IQRv, direction='up'))
    xmax <- as.numeric(postNet::roundNice(stats::quantile(regData$effIn,0.75) + 1.5*IQRv, direction='up'))
  }
  
  xmid <- (xmin + xmax)/2

  n_rows <- nrow(tableOut)
  n_cols <- ncol(tableOut)
  
  #
  par(mar=c(0,4,2,2))
  plot.new()
  plot.window(xlim=c(xmin,xmax), ylim=c(0,1))
  
  #
  tableCex <- max(0.8, min(1.5, 1.5 - n_rows/25))
  xpad <- 0.3
  ypad <- 0.6
  
  #
  if(n_cols > 10) {
    tableCex <- tableCex * 10/n_cols
    xpad <- xpad * 10/n_cols
  }
  
  # 
  if(n_rows > 25) {
    ypad <- ypad * 25/n_rows
    tableCex <- tableCex * 25/n_rows
  }
  
  #
  approx_char_width <- max(nchar(unlist(tableOut)), na.rm=TRUE)*0.05
  table_width <- approx_char_width * n_cols
  scale_factor <- (xmax - xmin)/table_width * 1.2
  tableCex <- tableCex * scale_factor
  xpad <- xpad * scale_factor
  
  #
  plotrix::addtable2plot(
    x = xmid, y = 0.1,
    table = tableOut,
    display.rownames = FALSE,
    hlines = TRUE,
    vlines = TRUE,
    bg = signature_colours,
    xpad = xpad,
    ypad = ypad,
    cex = tableCex,
    bty = "o",
    xjust = 0.5
  )
  
  #
  par(mar=c(5,5,0,5))
  plot(stats::ecdf(as.numeric(regData[,2])),
                 col='white', main='', xlab='effect',
                 verticals=TRUE, do.p=FALSE, lwd=3,
                 xlim=c(xmin,xmax))
  legend(xmin,0.95,fill='grey55',border='grey55','Background',bty='n',cex=1.3)
  
  # 
  if(any(duplicated(unlist(signatureList)))) {
    for(i in 1:length(signatureList)) {
      regData[, 2+i] <- 'bkg'
      regData[, 2+i][regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
      
      lines(stats::ecdf(as.numeric(regData[regData[,2+i]=='bkg','effIn'])),
                      col='grey55', verticals=TRUE, do.p=FALSE, lwd=3)
      lines(stats::ecdf(as.numeric(regData[regData[,2+i]==names(signatureList)[i],'effIn'])),
                      col=signature_colours[i], verticals=TRUE, do.p=FALSE, lwd=3)
    }
  } else {
    regData$signature <- 'bkg'
    for(i in 1:length(signatureList)) {
      regData$signature[regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
    }
    lines(stats::ecdf(as.numeric(regData[regData$signature=='bkg','effIn'])),
                    col='grey55', verticals=TRUE, do.p=FALSE, lwd=3)
    for(i in 1:length(signatureList)) {
      lines(stats::ecdf(as.numeric(regData[regData$signature==names(signatureList)[i],'effIn'])),
                      col=signature_colours[i], verticals=TRUE, do.p=FALSE, lwd=3)
    }
  }
  
  dev.off()
  
  
  
}
    