.okTOM = function (tom) {
  if (!is.matrix(tom) | nrow(tom)<2) {
    stop("tom must be a matrix with more than one row")
  }
  if (!isSymmetric(tom[1:2,1:2])) {
    stop("tom must be symmetric")
  }
  if (any(tom<0) | any(tom>1) | any(is.infinite(tom))) {
    stop('tom must contain values between 0 and 1')
  }
  if (any(is.na(tom))) {
    stop('tom cannot have missing values')
  }
  return(TRUE)
}



.getNNeighbors = 
function (gene, tom, n) {
  if (.okTOM(tom)) {
    return(sort(tom[rownames(tom)==gene, ],
                decreasing=TRUE,
                na.last=TRUE)
           [2:(n+1)])
  }
}


.grey2red = 
function (n, base, gamma) {
  red = seq(from=base^gamma, to=1, length.out = n)^(1/gamma)
  green = blue = seq(from = base^gamma, to=0, length.out = n)^(1/gamma);
  col = rgb(red, green, blue, maxColorValue = 1); 
}

.numbers2colors = function (x, signed = NULL, centered = signed, lim = NULL, commonLim = FALSE, 
          colors = if (signed) blueWhiteRed(100) else blueWhiteRed(100)[51:100], 
          naColor = "grey") 
{
  x = as.matrix(x)
  if (!is.numeric(x)) 
    stop("'x' must be numeric. For a factor, please use as.numeric(x) in the call.")
  if (is.null(signed)) {
    if (any(x < 0, na.rm = TRUE) & any(x > 0, na.rm = TRUE)) {
      signed = TRUE
    }
    else signed = FALSE
  }
  if (is.null(centered)) 
    centered = signed
  if (is.null(lim)) {
    if (signed & centered) {
      max = apply(abs(x), 2, max, na.rm = TRUE)
      lim = as.matrix(cbind(-max, max))
    }
    else {
      lim = as.matrix(cbind(apply(x, 2, min, na.rm = TRUE), 
                            apply(x, 2, max, na.rm = TRUE)))
    }
    if (commonLim) 
      lim = c(min(lim[, 1], na.rm = TRUE), max(lim[, 2], 
                                               na.rm = TRUE))
  }
  if (is.null(dim(lim))) {
    if (length(lim) != 2) 
      stop("'lim' must be a vector of length 2 or a matrix with 2 columns.")
    if (!is.numeric(lim)) 
      stop("'lim' must be numeric")
    if (sum(is.finite(lim)) != 2) 
      stop("'lim' must be finite.")
    lim = t(as.matrix(lim))
  }
  else {
    if (ncol(x) != nrow(lim)) 
      stop("Incompatible numbers of columns in 'x' and rows in 'lim'.")
    if (!is.numeric(lim)) 
      stop("'lim' must be numeric")
    if (sum(is.finite(lim)) != length(lim)) 
      stop("'lim' must be finite.")
  }
  xMin = matrix(lim[, 1], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  xMax = matrix(lim[, 2], nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
  if (sum(xMin == xMax) > 0) 
    warning("(some columns in) 'x' are constant. Their color will be the color of NA.")
  xx = x
  xx[is.na(xx)] = ((xMin + xMax)[is.na(xx)])/2
  if (sum(x < xMin, na.rm = TRUE) > 0) {
    warning("Some values of 'x' are below given minimum and will be truncated to the minimum.")
    x[xx < xMin] = xMin[xx < xMin]
  }
  if (sum(x > xMax, na.rm = TRUE) > 0) {
    warning("Some values of 'x' are above given maximum and will be truncated to the maximum.")
    x[xx > xMax] = xMax[xx > xMax]
  }
  mmEq = xMin == xMax
  nColors = length(colors)
  xCol = array(naColor, dim = dim(x))
  xInd = (x - xMin)/(xMax - xMin)
  xInd[xInd == 1] = 1 - 0.5/nColors
  xCol[!mmEq] = colors[as.integer(xInd[!mmEq] * nColors) + 
                         1]
  xCol[is.na(xCol)] = naColor
  xCol
}
.checkAdjMat = 
function (adjMat, min = 0, max = 1){
  dim = dim(adjMat)
  if (is.null(dim) || length(dim) != 2) 
    stop("adjacency is not two-dimensional")
  if (!is.numeric(adjMat)) 
    stop("adjacency is not numeric")
  if (dim[1] != dim[2]) 
    stop("adjacency is not square")
  if (max(abs(adjMat - t(adjMat)), na.rm = TRUE) > 1e-12) 
    stop("adjacency is not symmetric")
  if (min(adjMat, na.rm = TRUE) < min || max(adjMat, na.rm = TRUE) > 
      max) 
    stop("some entries are not between", min, "and", max)
}


.circlePlot = 
function (adjacency, labels, order,
          startNewPlot=TRUE, plotBox=c(-1,1,-1,1),
          center=c(0,0), radii=c(0.8, 0.8), startAngle=0,
          variable.cex.labels=TRUE,
          min.cex.labels=1, max.cex.labels=1.5,
          variable.cex.points=TRUE,
          min.cex.points=1, max.cex.points=3,
          variable.line.width=TRUE,
          min.line.width=1, max.line.width=5,
          lineColors=.grey2red(50,0.6,1),
          pch=21, pointBg="black",
          labelColors="black", pointColors="black",
          xMargin=1-radii[1], yMargin=1-radii[2],
          xLabelOffset=0.01, yLabelOffset=0.01,
          variableLabelAngle=TRUE,
          ...) {
  
  if (startNewPlot) {
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);
  }
  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  .checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);
  
  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = apply(abs(adjx), 2, sum)-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);
  
  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);
  
  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];
  
  actualCexPts = rep(0, n);
  for (node in 1:n) {
    cex = min.cex.points;
    if (variable.cex.points)
      cex = min.cex.points + (max.cex.points - min.cex.points) * 
        (oConn[node] - minConn)/(maxConn - minConn)
    actualCexPts[node] = cex
  }
  
  diag(oAdj) = 0;
  maxA = max(abs(oAdj));
  if (sum(oAdj < 0) > 0) {
    adjCol = .numbers2colors(oAdj, signed = TRUE, lim = c(-maxA, maxA), colors = lineColors); #colors = lineColors added by ATH
  } else {
    adjCol = .numbers2colors(oAdj, signed = FALSE, lim = c(0, maxA), colors = lineColors); #colors = lineColors added by ATH
  }
  
  ltA = oAdj;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;
  
  adjOrder = order(c(abs(ltA)))
  rows = row(oAdj)[adjOrder];
  cols = col(oAdj)[adjOrder];
  
  nLines = n*(n-1)/2;
  for (line in 1:nLines) {
    n1 = rows[line];
    n2 = cols[line];
    a = oAdj[n1, n2];
    normA = abs(a)/maxA;
    
    w = min.line.width;
    if (variable.line.width) {
      w = min.line.width + (max.line.width - min.line.width) * normA;
    }
    #pRadius1 = par("cxy") * actualCexPts[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * actualCexPts[n2]/35;
    lineLen = sqrt( (x[n1] - x[n2])^2 + (y[n1] - y[n2])^2);
    x1 = x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen
    
    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }
  
  for (node in 1:n) {
    points(x[node], y[node], pch = oPch[node], cex = actualCexPts[node], bg = oPBg[node], col = oPColors[node]);
  }
  for (node in 1:n) {
    cex = min.cex.labels;
    if (variable.cex.labels)
      cex = min.cex.labels + (max.cex.labels - min.cex.labels) *
        (oConn[node] - minConn)/(maxConn - minConn)
    textWidth = strwidth(oLabs[node], cex = cex);
    textHeight = strheight(oLabs[node], cex = cex);
    if (variableLabelAngle) {
      ang = angles[node]/pi * 180;
      if (ang < 180) {
        dir = 1;
      } else {
        dir = -1;
        ang = ang - 180;
      }
      ang = (90 - ang)/2
      xDir = 1;
      yDir = 1;
      cosAng = cos(ang/180*pi);
      sinAng = sin(ang/180*pi);
    } else {
      ang = 0;
      xDir = x[node];
      yDir = y[node];
      cosAng = 1;
      sinAng = 1;
      dir = 1;
    }
    angRad = ang/180*pi;
    pRadius = par("cxy") * actualCexPts[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(sum(c(cosAng^2, sinAng^2) * pRadius^2));
    rotMat = matrix( c(cosAng, sinAng, -sinAng, cosAng), 2, 2);
    labelShift = rotMat %*% as.matrix(c(textWidth, textHeight));
    text(x[node] + dir * xDir * (labelShift[1]/2 + cosAng * effPointRadius + xLabelOffset[node]), 
         y[node] + dir * yDir * (labelShift[2]/2 + sinAng * effPointRadius + yLabelOffset[node]), 
         labels = oLabs[node], adj = c(0.5, 0.5), 
         cex = cex, col = oLColors[node], srt = ang, xpd = TRUE);
  }
  
}