#' Kullback-Leibler for zero-mean Gaussian data
#' 
#' Compute the Kullback-Leibler for zero-mean Gaussian distributed data. Can be used to measure the distance between the true underlying covariance (or precision) matrix and the estimated one.
#' @param Ktr the true matrix. Can be both covariance or precision
#' @param K the estimated matrix. Can be both covariance or precision
#' @return scalar, the distance between N(0, Ktr) and N(0, K)
#'
#' @export
KL_dist = function(Ktr,K){
  p = dim(K)[1]
  inv_Ktr = solve(Ktr)
  A = sum( diag(inv_Ktr%*%K) ) 
  B = log( det(K)/det(Ktr) )
  return (0.5*(A - p - B))
} 

#' Upper to complete matrix
#' 
#' Maps upper triangular matrix into symmetric matrix. Keeps the same diagonal.
#' @param A an upper triangular matrix. Lower part is not used and overwritten
#' @return symmetric matrix
#' 
#' @export
upper2complete = function(A){
  diag_A = diag(A)
  A = A + t(A)
  diag(A) = diag_A
  return(A)
}

#' Compute number of errors
#' 
#' Takes two matrices containing 0s and 1s and compute the number of elements that are different. The diagonal may or may not be checked.
#' @param G1 binary matrix, only the upper part is used.
#' @param G2 binary matrix, only the upper part is used.
#' @param diag boolean, if the diagonal has to be checked or not
#' @return the number of elements that are different
#' 
#' @export
Graph_distance = function(G1,G2,diag = F){
  
  G1_vett = G1[upper.tri(G1,diag = diag)]
  G2_vett = G2[upper.tri(G1,diag = diag)]
  Table   = table(G1_vett, G2_vett)
  
  return ( Table[1,2] + Table[2,1] )
}


#' Makes colors trasparent
#'
#' Get a color and returns its transparent version. New color can be named.
#' @param color the color name
#' @param percent percentage of transparecy. Low values stands for more transparecy.
#' @param name an optional name for the new color 
#' @return the new color, which is also saved.
#'
#' @export
t_col <- function(color, percent = 30, name = NULL) {
  #    color = color name
  #    percent = % transparency
  #    name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
} #funzione per rendere i colori trasparenti




#' Heatmap of a matrix. 
#'
#' It plots the heatmap of a given matrix. It is usefull because it allows the user to use custom, non-symmetric palettes, centered around a desired value.
#' @param Mat the matrix to be plotted.
#' @param col.upper the color for the highest value in the matrix.
#' @param col.center the color for the centered value (if desired) or for the middle value in the matrix (oterhwise).
#' @param col.lower the color for the lowest value in the matrix.
#' @param center_value the value around which the palette has to be centered. Set NULL for symmetric palette.
#' @param col.n_breaks has to be odd. The refinement of the palette.
#' @param use_x11_device boolean, if x11() device has to be activated or not.
#' @param title the title of the plot.
#' @param x_label the name of the x-axis in the plot.
#' @param y_label the name of the y-axis in the plot.
#' @param remove_diag boolean, if the diagonal has to be removed or not.
#' @param horizontal boolean, if the heatmap bar has to be horizontal or not.
#' @return this function does not return anything
#'
#' @export
ACheatmap = function(Mat, center_value = 0, col.upper = "darkred", col.center = "grey95", col.lower = "#3B9AB2", 
                     col.n_breaks = 59, use_x11_device = TRUE, remove_diag = FALSE, title = "Title", x_label = "x_axis", y_label = "y_axis", horizontal=TRUE ) {

  library(fields)
  #Check for NA
  if(any(is.na(Mat))){
    cat('\n NA values have been removed from Matrix  \n')
  }

  #Create color palette
  if(col.n_breaks %% 2 == 0){
    warning( 'col.n_breaks is even but it has to be odd. Adding 1' )
    col.n_breaks = col.n_breaks + 1
  }
  colorTable = designer.colors(col.n_breaks, c( col.lower, col.center, col.upper) )
  col_length = (col.n_breaks + 1) / 2

  #Check Matrix
  if(remove_diag){
    diag(Mat) = NA
  }
  min_val = min(Mat,na.rm=T)
  max_val = max(Mat,na.rm=T)
  p_row = dim(Mat)[1]
  p_col = dim(Mat)[2]


  #Plot
  if(!is.null(center_value)){

    if(!(min_val < center_value & max_val > center_value)){
      stop('\n The lowest value has to be smaller than center_value and the highest value has to be larger. \n')
    }
    brks = c(seq( min_val, center_value-0.0001,l=col_length), seq( center_value+0.0001, max_val, l=col_length))

  }else{
    brks = seq( min_val, max_val, l=2*col_length)
  }
    
  colnames(Mat) = 1:p_col
  rownames(Mat) = 1:p_row

  if(use_x11_device){
    x11()
  }

  par(mar=c(5.1, 4.1, 4.1, 2.1),mgp=c(3,1,0))
  if(horizontal){
    
    fields::image.plot(Mat, axes=F, horizontal=T, main = title,
                       col=colorTable,breaks=brks,xlab=x_label,ylab=y_label)
  }else{

    fields::image.plot(Mat, axes=F, horizontal=FALSE, main = title,
                     col=colorTable,breaks=brks,xlab=x_label,ylab=y_label)
  }
  box()

}



#' Heatmap of a matrix with hcl palettes. 
#'
#' It plots the heatmap of a given matrix using one of the hcl palettes. It is a quicker version of the more complete function \code{\link{ACheatmap}}
#' @param Mat the matrix to be plotted.
#' @param palette.name the name of the hcl palette to be used. See here for possible names https://developer.r-project.org/Blog/public/2019/04/01/hcl-based-color-palettes-in-grdevices/  
#' @param col.n_breaks the refinement of the palette.
#' @inheritParams ACheatmap
#' @return this function does not return anything
#'
#' @export
ACheatmap_hcl = function(Mat, palette.name = "RdGy", col.n_breaks = 64, use_x11_device = TRUE, remove_diag = FALSE, title = 'title', x_label = "x_axis", y_label = "y_axis", horizontal = TRUE ) {

  library(fields)

  #Check for NA
  if(any(is.na(Mat))){
    cat('\n NA values have been removed from Matrix  \n')
  }

  #Check Matrix
  if(remove_diag){
    diag(Mat) = NA
  }

  #Plot 
  if(use_x11_device){
    x11()
  }

  par(mar=c(5.1, 4.1, 4.1, 2.1),mgp=c(3,1,0))
  if(horizontal){
    fields::image.plot(Mat, col= hcl.colors(col.n_breaks,palette.name), main = title, xlab=x_label,ylab=y_label, horizontal = TRUE)
  }else{
    fields::image.plot(Mat, col= hcl.colors(col.n_breaks,palette.name), main = title, xlab=x_label,ylab=y_label, horizontal = FALSE)
  }

  

}


#' Plot boxplot from matrix
#'
#' Given a (nxp)-matrix, it plots the boxplots of the columns.
#' @param data an (nxp)-matrix.
#' @param use_ggplot if boxplots should be done using ggplot or not.
#' @inheritParams ACheatmap
#' @return this function does not return anything
#'
#' @export
boxplot_from_matrix = function(data, use_x11_device = TRUE, use_ggplot = TRUE, title = "Title", x_label = "x_axis", y_label = "y_axis"){

  library(reshape)

  # Turn into data.frame
  df = suppressWarnings ( reshape::melt(data, varnames = c("row.index", "variable")) )#creates a datafram of p*n observations and 3 columns.

  #Plot 
  if(use_x11_device){
    x11()
  }
  if(use_ggplot){
    library(ggplot2)
    title_theme = labs(title=title, x=x_label, y=x_label) #FA SCHIFO
    ggplot(data = df, aes(x = factor(variable), y = value, fill = factor(variable)) ) + geom_boxplot() + title_theme + theme(panel.background = element_blank())
  }else{
    boxplot(value~variable, data =  df, main = title, xlab = x_label, ylab = y_label ,col = 'red')
  }

}


#' Plot boxplot from vectors
#'
#' Takes an arbitrary number of vectors and plot the boxplots. When calling this function, any arguments after ... must be fully named. 
#' @param ... an arbitrary number of vectors. The function does not checks if they are really vectors, so be careful.
#' @param names the names of the groups of the different boxplots. Not defaulted for safety, can be set to 1:length(...).
#' @inheritParams ACheatmap
#' @return this function does not return anything
#' @examples 
#' boxplot_from_vectors(v1,v2,v3, names = 1:3)
#' boxplot_from_vectors(v1,v2,v3, names = c("one","two", "three"), title = "Nice Boxplot")
#'
#' @export
boxplot_from_vectors = function(..., names, use_x11_device = TRUE, use_ggplot = TRUE, title = "Title", x_label = "x_axis", y_label = "y_axis"){

  #read ...
  l = list(...)
  L = length(l)
  
  #check length and consistency with names
  if(L < 1)
    stop("Number of vectors passed in ... has to be positive")
  
  if(L != length(names)){
    stop('Number of vector in ... is not consistent with the length of names')
  }
  
  #create df for the plot
  df = data.frame(  value=l[[1]], variable = rep(names[1], length(l[[1]]))  )
  if(L>1){
    for(i in 2:L){
      df = rbind(df, data.frame(  value=l[[i]], variable = rep(names[i], length(l[[i]]))  ) )
    }
  }
  
  #Plot 
  if(use_x11_device){
    x11()
  }
  if(use_ggplot){
    library(ggplot2)
    title_theme = labs(title=title, x=x_label, y=x_label) #FA SCHIFO
    ggplot(data = df, aes(x = factor(variable), y = value, fill = factor(variable)) ) + geom_boxplot() + title_theme + theme(panel.background = element_blank())
  }else{
    boxplot(value~variable, data =  df, main = title, xlab = x_label, ylab = y_label ,col = 'red')
  }


}