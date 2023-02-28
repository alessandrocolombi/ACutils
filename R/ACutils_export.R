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
#' @param main the title of the plot.
#' @param x_label the name of the x-axis in the plot.
#' @param y_label the name of the y-axis in the plot.
#' @param remove_diag boolean, if the diagonal has to be removed or not.
#' @param horizontal boolean, if the heatmap bar has to be horizontal or not.
#' @return this function does not return anything
#'
#' @export
ACheatmap = function(Mat, center_value = 0, col.upper = "darkred", col.center = "grey95", col.lower = "#3B9AB2",
                     col.n_breaks = 59, use_x11_device = TRUE, remove_diag = FALSE, main = " ", x_label = " ",
                     y_label = " ", horizontal=TRUE )
{

  #library(fields)
  #Check for NA
  if(any(is.na(Mat))){
    cat('\n NA values have been removed from Matrix  \n')
  }

  #Create color palette
  if(col.n_breaks %% 2 == 0){
    warning( 'col.n_breaks is even but it has to be odd. Adding 1' )
    col.n_breaks = col.n_breaks + 1
  }
  colorTable = fields::designer.colors(col.n_breaks, c( col.lower, col.center, col.upper) )
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

    fields::image.plot(Mat, axes=F, horizontal=T, main = main,
                       col=colorTable,breaks=brks,xlab=x_label,ylab=y_label)
  }else{

    fields::image.plot(Mat, axes=F, horizontal=FALSE, main = main,
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
ACheatmap_hcl = function(Mat, palette.name = "RdGy", col.n_breaks = 64, use_x11_device = TRUE, remove_diag = FALSE,
                         main = ' ', x_label = " ", y_label = " ", horizontal = TRUE )
{

  #library(fields)

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
boxplot_from_matrix = function(data, use_x11_device = TRUE, use_ggplot = TRUE, main = " ", x_label = " ",
                               y_label = " ")
{

  #library(reshape)

  # Turn into data.frame
  df = suppressWarnings ( reshape::melt(data, varnames = c("row.index", "variable")) )#creates a datafram of p*n observations and 3 columns.

  #Plot
  if(use_x11_device){
    x11()
  }
  if(use_ggplot){
    library(ggplot2)
    title_theme = labs(title=main, x=x_label, y=x_label) #FA SCHIFO
    ggplot(data = df, aes(x = factor(variable), y = value, fill = factor(variable)) ) + geom_boxplot() + title_theme + theme(panel.background = element_blank())
  }else{
    boxplot(value~variable, data =  df, main = main, xlab = x_label, ylab = y_label ,col = 'red')
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
boxplot_from_vectors = function(..., names, use_x11_device = TRUE, use_ggplot = TRUE, main = " ",
                                x_label = " ", y_label = " ")
{

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
    title_theme = labs(title=main, x=x_label, y=x_label) #FA SCHIFO
    ggplot(data = df, aes(x = factor(variable), y = value, fill = factor(variable)) ) + geom_boxplot() + title_theme + theme(panel.background = element_blank())
  }else{
    boxplot(value~variable, data =  df, main = main, xlab = x_label, ylab = y_label ,col = 'red')
  }
}


#' Plot curves
#'
#' \loadmathjax This functions gets one or two dataset representig functional data and plot them. It does not smooth the curves, indeed it requires as input the data, not
#' the regression coefficients. Use \code{\link{smooth_curves}} function for that.
#' @param data1 matrix of dimension \mjseqn{n\_curves \times n\_grid\_points} representing the first functional dataset to be plotted.
#' @param data2 matrix of dimension \mjseqn{n\_curves \times n\_grid\_points} representing the second dataset to be plotted, if needed.
#' @param range the range where the curves has to be plotted. Not needed if \code{n_plot} is 0.
#' @param n_plot the number of curves to be plotted. Set 0 for no plot.
#' @param grid_points vector of size \mjseqn{n\_grid\_points} with the points where the splines are evaluated. If defaulted they are uniformly generated. Not needed if \code{n_plot} is 0.
#' @param internal_knots vector with the internal knots used to construct the splines. Default is null.
#' If provided, \code{n_basis} are displayed in the plot. The \code{k}-th interval represents the segment where the \code{k}-th spline dominates the others.
#' @param highlight_band1 a vector that states if a particular band of the plot has to be highlighted. It has to be a vector within the range \mjseqn{\[1,n\_basis\]}.
#' @param highlight_band2 a vector that states if a particular band of the plot has to be highlighted. It has to be a vector within the range \mjseqn{\[1,n\_basis\]}.
#' @param main the title of the plot.
#' @param xtitle the title of the x-axis.
#' @param ytitle the title of the x-axis.
#' @param legend_name1 the name for \code{data1} to be printed in the legend.
#' @param legend_name2 the name for \code{data2} to be printed in the legend. Used only is two datasets are actually plotted.
#'
#' @return No values are returned.
#' @export
ACplot_curves = function( data1, data2 = NULL, range, n_plot = 1, grid_points = NULL,
                          internal_knots = NULL, highlight_band1 = NULL, highlight_band2 = NULL,
                          main = "Curves", xtitle = " ", ytitle = " ", legend_name1 = "data1", legend_name2 = "data2")
{
  #Plot
  if(n_plot > 0) #Plot n_plot curves
  {
    #Check dimensions
    if(!(length(range)==2 && range[1] < range[2]))
      stop("Invalid range, it has to be a vector of length 2 containing first the lower bound of the interval and then the upper bound.")
    #Computes grid_points
    if(!is.null(grid_points)){
      if(length(grid_points) != r)
        stop("The number of points provided in grid_points is not equal to the size of BaseMat.")
      X = grid_points;
    }else{
      X = seq(range[1], range[2], length.out = r)
    }
    #Classical plot, does not depend on the size of the graph
    if(is.null(internal_knots)){
      if(n_plot > 1){
        x11(height=4)
        matplot( x = X, t(data1[1:n_plot,]), type = 'l', lty = 1,
                col = c('darkolivegreen','darkgreen','darkolivegreen4','forestgreen'),
                lwd = 3, ylim = c(min(data1[1:n_plot,]),max(data1[1:n_plot,])), axes = T,
                main = title_plot, xlab = xtitle,
                ylab = ytitle)
        if(!is.null(data2)){
            matplot( x = X, t(data2[1:n_plot,]), type = 'l', lty = 1,
            col = c('steelblue','steelblue1','skyblue3', 'lightsteelblue1'), lwd = 3, add = T)
            legend("topright", legend=c(legend_name1, legend_name2), col = c('darkolivegreen','steelblue'), lty = c(1,1), lwd = 3)
        }
      }
      else if(n_plot == 1){
        x11(height=4)
        matplot( x = X, (data1[1:n_plot,]), type = 'l', lty = 1,
                col = c('darkolivegreen','darkgreen','darkolivegreen4','forestgreen'),
                lwd = 3, ylim = c(min(data1[1:n_plot,]),max(data1[1:n_plot,])), axes = T,
                main = title_plot, xlab = xtitle,
                ylab = ytitle)
       if(!is.null(data2)){
           matplot( x = X, (data2[1:n_plot,]), type = 'l', lty = 1,
           col = c('steelblue','steelblue1','skyblue3', 'lightsteelblue1'), lwd = 3, add = T)
           legend("topright", legend=c(legend_name1, legend_name2), col = c('darkolivegreen','steelblue'), lty = c(1,1), lwd = 3)
       }
      }
    }
    else{ #Plot with bands representing the domanin of the spline
        knots <- c(range[1],
                 range[1] + (internal_knots[1]-range[1])/2,
                 internal_knots,
                 range[2] - (range[2]-internal_knots[length(internal_knots)])/2,
                 range[2] )
        knots_name = round(knots, digits = 2)
        names <- rep("", length(knots_name))
        for (i in 1:length(knots_name)) {
          names[i] <- paste0(knots_name[i])
        }
        names_y = round(seq(min(data1[1:n_plot,]), max(data1[1:n_plot,]), length.out = 10), digits = 2)
        if(n_plot > 1){
          x11(height=4)
            matplot(x = X, t(data1[1:n_plot,]), type = 'l', lty = 1,
                    col = c('darkolivegreen','darkgreen','darkolivegreen4','forestgreen'),
                    lwd = 3, ylim = c(min(data1[1:n_plot,]),max(data1[1:n_plot,])), axes = F,
                    main = title_plot, xlab = xtitle,
                    ylab = ytitle)
            if(!is.null(data2)){
                matplot( x = X, t(data2[1:n_plot,]), type = 'l', lty = 1,
                col = c('steelblue','steelblue1','skyblue3', 'lightsteelblue1'), lwd = 3, add = T)
                legend("topright", legend=c(legend_name1, legend_name2), col = c('darkolivegreen','steelblue'), lty = c(1,1), lwd = 3)
            }
            abline(v = knots, lty = 2, col = 'black')
                if(!is.null(highlight_band1)){
                  for(i in c(highlight_band1, highlight_band1[length(highlight_band1)]+1) )
                    abline(v = knots[i], lty = 2, col = 'red')
                }
                if(!is.null(highlight_band2)){
                  for(i in c(highlight_band2[1]-1,highlight_band2) )
                    abline(v = knots[i], lty = 2, col = 'red')
                }
            mtext(text = names, side=1, line=0.3, at = knots , las=2, cex=0.7)
            mtext(text = names_y, side=2, line=0.3, at=names_y, las=1, cex=0.9)
        }
        else if(n_plot == 1){
            x11(height=4)
              matplot(x = X, (data1[1:n_plot,]), type = 'l', lty = 1,
                      col = c('darkolivegreen','darkgreen','darkolivegreen4','forestgreen'),
                      lwd = 3, ylim = c(min(data1[1:n_plot,]),max(data1[1:n_plot,])), axes = F,
                      main = title_plot, xlab = xtitle,
                      ylab = ytitle)
              if(!is.null(data2)){
                  matplot( x = X, (data2[1:n_plot,]), type = 'l', lty = 1,
                  col = c('steelblue','steelblue1','skyblue3', 'lightsteelblue1'), lwd = 3, add = T)
                  legend("topright", legend=c(legend_name1, legend_name2), col = c('darkolivegreen','steelblue'), lty = c(1,1), lwd = 3)
              }
              abline(v = knots, lty = 2, col = 'black')
                  if(!is.null(highlight_band1)){
                    for(i in c(highlight_band1, highlight_band1[length(highlight_band1)]+1) )
                      abline(v = knots[i], lty = 2, col = 'red')
                  }
                  if(!is.null(highlight_band2)){
                    for(i in c(highlight_band2[1]-1,highlight_band2) )
                      abline(v = knots[i], lty = 2, col = 'red')
                  }
              mtext(text = names, side=1, line=0.3, at = knots , las=2, cex=0.7)
              mtext(text = names_y, side=2, line=0.3, at=names_y, las=1, cex=0.9)
        }
    }
  }
  else
    stop("The number of curves to be plotted has to be positive.")
}


#' Plot undirected graph from matrix or vector
#'
#' \loadmathjax Given a (pxp)-matrix or a n_elem-dimensional vector, it plots the corresponding undirected graph.
#' @param G it may be (pxp)-matrix, where p is the number or nodes or a vector containing the elements of the adiacency
#' matrix. It is a vector, it may contain or not the diagonal elements. Its length should be equal to
#' \code{choose(p,2)} if diagonal elements are not included or \code{choose(p,2) + p}
#' if they are not included. Values should be 1 or 0.
#' @param col.links the color of the links.
#' @param col.zeros the color of the absent links, the zeros in the adiacency matrix.
#' @param graph.size the number of nodes. Used only is \code{G} is a vector.
#' @param byrow boolean, used only is \code{G} is a vector. Set TRUE if the elements of G are taken by row from its
#' adiacency matrix. Set FALSE if they are taken by column.
#' @inheritParams ACheatmap
#' @return this function does not return anything
#'
#' @export
ACplot_graph = function(G, col.links = "#3B9AB2", col.zeros = "grey95", graph.size = NULL, byrow = TRUE,
                        main =' ', x_label = " ", y_label = " ", use_x11_device = TRUE)
{
  if(!is.matrix(G) & !is.vector(G))
    stop('Graph G is not recognized as vector or matrix')
  if(is.matrix(G)){
    if(dim(G)[1]!=dim(G)[2])
      stop('G is recognized as a matrix but it is not squared')
    p = dim(G)[1]
  }
  if(is.vector(G)){
    if(is.null(graph.size))
      stop('G is recognized as a vector but parameter graph.size is null. Please insert the number of nodes')

    p = graph.size
    n_elem = choose(p,2)
    n_elem_diag = n_elem + p
    if(length(G)==n_elem){
      diag = FALSE
    }else if( length(G) == n_elem_diag){
      diag = TRUE
    }else{
      stop("G is recognized as a vector but its length is not compatible with graph.size. length(G) should be choose
        (graph.size,2) if the diagonal is included or choose(graph.size,2) - graph.size if the diagonal is not included
        ")
    }

    Gplot = matrix(0,p,p)
    if(byrow == TRUE){
      Gplot[lower.tri(Gplot, diag = diag)] = G
    }else{
      Gplot[upper.tri(Gplot, diag = diag)] = G
    }
    G = Gplot
  }

  if(!isSymmetric(G)){
    G = G + t(G)
  }
  diag(G) = rep(1,p)

  ACheatmap(Mat = G, center_value = 0.5, col.lower = col.zeros, col.upper = col.links, col.n_breaks = 63,
            use_x11_device = use_x11_device, main = main, x_label = x_label, y_label = y_label, horizontal=FALSE )
}


#' Plot undirected network from matrix or vector
#'
#' @export
ACplot_network = function(G, labels.name = NULL, col.nodes = NULL, col.links = NULL, col.labels = "white",
                          node.size = 12 ,label.size = 3.5, use_x11_device = T)
{
  p = dim(G)[1]
  if(is.null(labels.name))
    labels.name = 1:p

  if(is.null(col.nodes))
    col.nodes = rgb(32,102,164,max=255)
  if(is.null(col.links))
    col.links = rgb(212,167,48,max=255)

  net = network(G, directed = FALSE)
  if(use_x11_device){
    x11()
  }
  GGally::ggnet2(net, label=TRUE, size = node.size, label.size = label.size, color=col.nodes, label.color = 'white',
         edge.color=col.links)
}


#' Non Central Student t - Density
#'
#' \loadmathjax This function evaluates the density of a Non Central Student t in a given point. Such a distribution is defined as follow:
#' \mjsdeqn{X~\sim~nct(n_{0},\mu_{0},\gamma_{0})} if
#' \mjtdeqn{$$\begin{eqnarray*}f(X|n_{0},\mu_{0},\gamma_{0}) = \frac{\Gamma(\frac{n_{0}+1}{2})}{\Gamma(\frac{n_{0}}{2})}\frac{1}{\sqrt{\pi\gamma_{0}n_{0}}}\left(1+\frac{1}{n_{0}}(\frac{x-\mu_{0}}{\gamma_{0}})^{2}\right)^{-\frac{n_{0}+1}{2}} {eqnarray*}$$}{\begin{eqnarray*}f(X|n_{0},\mu_{0},\gamma_{0}) = \frac{\Gamma(\frac{n_{0}+1}{2})}{\Gamma(\frac{n_{0}}{2})}\frac{1}{\sqrt{\pi\gamma_{0}n_{0}}}\left(1+\frac{1}{n_{0}}(\frac{x-\mu_{0}}{\gamma_{0}})^{2}\right)^{-\frac{n_{0}+1}{2}}  \end{eqnarray*}}{\begin{eqnarray*} f(X|n_{0},\mu_{0},\gamma_{0}) = \frac{\Gamma(\frac{n_{0}+1}{2})}{\Gamma(\frac{n_{0}}{2})}\frac{1}{\sqrt{\pi\gamma_{0}n_{0}}}\left(1+\frac{1}{n_{0}}(\frac{x-\mu_{0}}{\gamma_{0}})^{2}\right)^{-\frac{n_{0}+1}{2}}  \end{eqnarray*}}
#' where \mjseqn{n_{0}} are the degree of freedom, \mjseqn{\mu_{0}} is the location parameter and \mjseqn{\gamma_{0}} is the scale parameter. The usual t density can be recovered by placing \mjseqn{\mu_{0}=0} and \mjseqn{\gamma_{0} = 1}.
#' @param x the point where to evaluate the density.
#' @param n0 the degree of fredoom.
#' @param mu0 the location parameter.
#' @param gamma0 the scale parameter.
#' @return scalar representing the evaluation of the density at x
#' @export
dnct = function( x, n0, mu0, gamma0 )
{
  if(gamma0 <= 0)
    stop("The scale parameter has to be strictly positive.")
  return( 1/sqrt(gamma0) * dt(x = (x-mu0)/gamma0, df = n0 ) )
}

#' Non Central Student t - Random Number Generator
#'
#' This function generates a random number distributed as Non Central Student t. See \code{\link{dnct}} for the definition of such a distribution.
#' @param n the number of points to be generated.
#' @inheritParams ACheatmap
#' @return scalar generated from Non Central Student t.
#' @export
rnct = function( n = 1, n0, mu0, gamma0 )
{
  if(gamma0 <= 0)
    stop("The scale parameter has to be strictly positive.")
  return(  gamma0*rt(n = n, df = n0) + mu0  )
}


#' Plot the desity function
#'
#' \loadmathjax This function plots the density function of some know distributions. Current version, does not allow to overlap different plots
#' if use_ggplot is active. Though, it is possible when using standard plots. In such a case, the legend can be added a posteriori.
#' @param ... the list of parameters needed by the distribution specified in dist.
#' @param dist string with the name of the density to be plotted. Possibilities are,
#' \code{"norm"} (requires mean and sd), \code{"gamma"} (requires shape and rate), \code{"inv-gamma"} (requires shape and scale, see below), \code{"chi-sq"} (requires dof), \code{"exp"} (requires rate),
#' \code{"t"} (requires dop), \code{"Fisher"} (requires dof1 and dof2), \code{"nct"} (requires degree of freedom, location and scale. See \code{\link{dnct}} for the definition of such a distribution.).
#' \code{"laplace"} (requires location and scale).
#' The \code{"inv-gamma"} distribution is defined as follow:
#' \mjsdeqn{X~\sim~inv-gamma(shape = \alpha, scale = \beta)} if
#' \mjtdeqn{$$\begin{eqnarray*}f(X|\alpha,\beta) = \frac{\beta^{\alpha}}{\Gamma(\alpha)}\left(\frac{1}{x}\right)^{\alpha+1}e^{- \frac{\beta}{x}}{eqnarray*}$$}{\begin{eqnarray*}f(X|\alpha,\beta) = \frac{\beta^{\alpha}}{\Gamma(\alpha)}\left(\frac{1}{x}\right)^{\alpha+1}e^{- \frac{\beta}{x}}\end{eqnarray*}}{\begin{eqnarray*} f(X|\alpha,\beta) = \frac{\beta^{\alpha}}{\Gamma(\alpha)}\left(\frac{1}{x}\right)^{\alpha+1}e^{- \frac{\beta}{x}}\end{eqnarray*}}
#' @param xlim set the range of x-axis. Not used if add is TRUE.
#' @param ylim set the range of y-axis. Not used if add is TRUE.
#' @param col set color of density line.
#' @param add boolean, if current plot has to be added to previous plot or not.
#' @inheritParams ACheatmap
#' @return this function does not return anything
#' @export
plot_density = function( ..., dist, main = " ", x_label = " ", y_label = " ", xlim = NULL, ylim = NULL,
                         col = "gray48", use_x11_device = TRUE, use_ggplot = FALSE, add = FALSE )
{
  #check density name
  if(!(dist == "norm" || dist == "gamma" || dist == "inv-gamma" || dist == "chi-sq" || dist == "exp" || dist == "t" || dist == "Fisher" || dist == "nct" || dist == "laplace"))
    stop("Unknown distribution requested, the only one availables are norm, gamma, chi-sq, exp, t, Fisher, nct, laplace")

  #read ...
  l = list(...)
  L = length(l)

  npoints = 100000

  if(dist == "norm"){
    if(L < 2){
      stop("norm density requires 2 parameters, mean and sd.")
    }else if(L > 2){
      warning("norm density requires 2 parameters, mean and sd but more has been provided. Use only the first two.")
    }else{
      y = rnorm(n = npoints, mean = l[[1]], sd = l[[2]])
    }

  }else if(dist == "gamma"){
    if(L < 2){
      stop("gamma density requires 2 parameters, shape and rate.")
    }else if(L > 2){
      warning("gamma density requires 2 parameters, shape and rate but more has been provided. Use only the first two.")
    }else{
      y = rgamma(n = npoints, shape = l[[1]], rate = l[[2]])
    }

  }else if(dist == "inv-gamma"){
    if(L < 2){
      stop("gamma density requires 2 parameters, shape and scale")
    }else if(L > 2){
      warning("gamma density requires 2 parameters, shape and scale but more has been provided. Use only the first two.")
    }else{
      y = 1 / ( rgamma(n = npoints, shape = l[[1]], rate = l[[2]]) )
    }

  }else if(dist == "chi-sq"){
    if(L < 1){
      stop("gamma density requires 1 parameter, the degree of freedom.")
    }else if(L > 1){
      warning("gamma density requires 1 parameter, the degree of freedom but more has been provided. Use only the first one.")
    }else{
      y = rchisq(n = npoints, df = l[[1]])
    }

  }else if(dist == "exp"){
    if(L < 1){
      stop("exp density requires 1 parameter, rate.")
    }else if(L > 1){
      warning("gamma density requires 1 parameter, rate but more has been provided. Use only the first one.")
    }else{
      y = rexp(n = npoints, rate = l[[1]])
    }

  }else if(dist == "t"){
    if(L < 1){
      stop("student t density requires 1 parameter, the degree of freedom.")
    }else if(L > 1){
      warning("gamma density requires 1 parameter, the degree of freedom but more has been provided. Use only the first one.")
    }else{
      y = rt(n = npoints, df = l[[1]])
    }

  }else if(dist == "Fisher"){
    if(L < 2){
      stop("Fisher density requires 2 parameters, df1 and df2.")
    }else if(L > 2){
      warning("Fisher density requires 2 parameters, df1 and df2 but more has been provided. Use only the first two.")
    }else{
      y = rf(n = npoints, df1 = l[[1]], df2 = l[[2]])
    }

  }else if(dist == "nct"){
    if(L < 3){
      stop("nct density requires 3 parameters, degree of freedom, location and scale.")
    }else if(L > 3){
      warning("nct density requires 2 parameters, degree of freedom, location and scale but more has been provided. Use only the first three.")
    }else{
      y = ACutils:::rnct(n = npoints, n0 = l[[1]], mu0 = l[[2]], gamma0 = l[[3]])
    }

  }else if(dist == "laplace"){
    if(L < 2){
      stop("laplace density requires 2 parameters, location and scale.")
    }else if(L > 2){
      warning("laplace density requires 2 parameters, location and scale but more has been provided. Use only the first two.")
    }else{
      y = ExtDist::rLaplace(n = npoints, mu = l[[1]], b = l[[2]])
    }

  }else{
    stop("Unknown distribution requested, the only one availables are norm, gamma, chi-sq, exp, t, Fisher, nct")
  }

  #Density estimation
  density = density(y)

  #check and set xlim and ylim
  if(!is.null(xlim)){
    if(length(xlim) != 2){
      warning("xlim has not length equal to 2")
      xlim = NULL
    }
    if(xlim[2]<xlim[1]){
      warning("xlim[2] can not be smaller than xlim[1]")
      xlim = NULL
    }
  }
  if(is.null(xlim)){
    xlim = c(0,0)
    xlim[1] = min(density$x)
    xlim[2] = max(density$x)
    range = xlim[2] - xlim[1]
    xlim[1] = xlim[1] - range/10
    xlim[2] = xlim[2] + range/10
  }

  if(!is.null(ylim)){
    if(length(ylim) != 2){
      warning("ylim has not length equal to 2")
      ylim = NULL
    }
    if(ylim[2]<ylim[1]){
      warning("ylim[2] can not be smaller than ylim[1]")
      ylim = NULL
    }
  }
  if(is.null(ylim)){
    ylim = c(0,0)
    ylim[1] = min(density$y)
    ylim[2] = max(density$y)
    range = ylim[2] - ylim[1]
    ylim[1] = ylim[1] - range/10
    ylim[2] = ylim[2] + range/10
  }


  #Plot
  if(use_x11_device){
    x11()
  }
  if(use_ggplot){
    #library(ggplot2)
    df = data.frame(  value=y  )
    title_theme = ggplot2::labs(title=main, x=x_label, y=y_label) #FA SCHIFO
    bg = ggplot2::theme_bw() #ggplot2::theme(panel.background = ggplot2::element_blank())
    ggplot_obj = ggplot2::ggplot(data = df, ggplot2::aes(x = value), y = value ) + ggplot2::geom_density( color = col ) + title_theme + bg

    # non so come aggiungere i plot in questo caso!
  }else{
    if(!add){ #new plot
      plot(x = density$x, y = density$y, main = main, xlab = x_label, ylab = y_label, col = col, type = 'l', lwd = 2, xlim = xlim, ylim = ylim)
    }else{
      points(x = density$x, y = density$y, main = main, xlab = x_label, ylab = y_label, col = col, type = 'l', lwd = 2)
    }

  }
}

#' sum_check
#' check if elements in a vector sum to an expected value. Elements of the vector
#'             are supposed to be computed as expected_sum*p. If expectations are not respected
#'             return a corrected version of vec, according to p
#'
#' @param vec vector to integer value to be checked
#' @param expected_sum expected value for sum(vec)
#' @param p theoretical weights applied to divide expected_sum in vec's elements
#' @return corrected version of vec (if needed)
#' @import dplyr
#' @export
sum_check <- function(vec, expected_sum, p){
  if( !dplyr::near( sum(p), 1) )
    stop("weights must sum to 1!")

  round_diff = vec - expected_sum*p # large if value in vec is OVER-estimate

  if(sum(vec) > expected_sum ){
    i_max = which.max(round_diff)
    vec[i_max] = vec[i_max] - 1
    vec = sum_check(vec, expected_sum, p)
  }
  if(sum(vec) < expected_sum){
    i_min = which.min(round_diff)
    vec[i_min] = vec[i_min] + 1
    vec = sum_check(vec, expected_sum, p)
  }

  return(vec)
}

#' dmix
#'
#' density of a mixture of gaussian distributions
#' @param x point or vector of values where I want to evaluate the density
#' @param w vector of weights of the
#' @param mu_vec vector of mean for the gaussian distribution
#' @param sigma_vec vector of standard deviation for the gaussian distribution
#' @return named matrix with (Inf, estimate, Sup)
#' @import dplyr
#' @export
dmix <- function(x, w_j, mu_vec, sigma_vec){
  if( !dplyr::near(sum(w_j), 1) )
    stop("weigths don't sum to one")

  if(length(w_j) != length(mu_vec) || length(w_j) != length(mu_vec) )
    stop("length of w_j, mu_vec and sigma_vec differs")

  K = length(w_j)

  n = length(x)
  mix_density = numeric(n)

  for(k in 1:K){
    mix_density = mix_density + w_j[k]*dnorm(x, mu_vec[k], sigma_vec[k])
  }

  return(mix_density)
}


#' rmix
#'
#' function to extract a random sample of dimension n from a Gaussian mixture
#' defined by the 3 vectors p, mu, sigma (that must have same length!)
#'
#' @param n dimension of the random sample that has to be extracted
#' @param p vector of weights for the components of the mixture
#' @param mu vector of means for the components
#' @param sigma vector of standard deviations for the components
#' @return random sample derived my the specified mixture
#' @import dplyr
#' @export
rmix <- function(n, p, mu, sigma){
  if( !dplyr::near(sum(p),1) )
    stop("weights for the components must sum to 1")

  if( length(p) != length(mu) || length(p) != length(sigma) )
    stop("number of weights for the components must agree with the dimensions of
         mu's vector and sigma's vector AND *mu* and *sigma* must have
         length")
  # rnormm cannot take 0 weighted components ==> eliminate them
  ind_0 = which( p == 0 )

  if(length(ind_0) == 0){
    sample_p = p
    sample_mu = mu
    sample_sigma = sigma
  }else{
    sample_p = p[-ind_0]
    sample_mu = mu[-ind_0]
    sample_sigma = sigma[-ind_0]
  }

  n_p = sapply(sample_p, function(x){round(x*n)})
  n_p = sum_check(n_p, n, sample_p)
  n_p = as.vector(n_p, mode = "integer")

  mix_sample = numeric(n)
  i = 1
  M = length(sample_p)

  for(m in 1:M){
    ind_m = seq(i, i+n_p[m]-1)
    mix_sample[ind_m] = rnorm(n_p[m], sample_mu[m], sample_sigma[m])
    i = i + n_p[m]
  }

  # shuffle sampled values
  mix_sample = mix_sample[sample(n)]

  return( mix_sample )
}



#' myboxplot2
#'
#' Useful function to produce boxplot within for loops.
#' @param mydata [tibble] the dataset
#' @param myexposure [string] the name of the numerical variable whose boxplot must the produced.
#' @param myoutcome [string] the name of the categorical variable that defines the different levels. One boxplot for each level.
#'
#' @export
#'
#' @examples myboxplot2(mydata = data, myexposure =  "treatment", myoutcome =  "age")
myboxplot2 <- function(mydata, myexposure = NULL, myoutcome){
  if(!is.null(myexposure)){
    bp <- ggplot(mydata, aes_(as.name(myexposure), as.name(myoutcome), fill = as.name(myexposure) )) +
      geom_boxplot() +
      labs(title= myoutcome ) + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))
    print(bp)
  }else{
    bp <- ggplot(mydata, aes_(y=as.name(myoutcome) )) +
      geom_boxplot() +
      labs(title= myoutcome ) + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))
    print(bp)
  }
}

#' mybarplot2
#'
#' Useful function to produce barplot within for loops.
#' @param mydata [tibble] the dataset
#' @param myoutcome [string] the name of the categorical variable whose counts must the shown in a barplot.
#' @param myexposure [string] the name of the categorical variable that defines the different levels.
#' If it is not NULL, different bar levels are produced.
#'
#' @return One bar for each level of level of the outcome variable. If there is exposure, an additional bar for each
#' level of the exposure is produced.
#' @export
#' @examples mybarplot2(mydata = data2, myexposure =  "treatment", myoutcome =  "race")
mybarplot2 <- function(mydata, myexposure = NULL, myoutcome){
  if(is.null(myexposure)){
    bp = mydata %>% ggplot( aes_(x=as.name(myoutcome) ) ) +
      geom_bar(stat="count", position=position_dodge()) +
      labs(title= myoutcome ) + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))
    print(bp)
  }
  else{
    bp = mydata %>% ggplot( aes_(x=as.name(myoutcome), fill=as.name(myexposure) ) ) +
      geom_bar(stat="count", position=position_dodge()) +
      labs(title= myoutcome ) + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))
    print(bp)
  }

}

