# ---------------------------------------------------------------------------------------------------------------
# --
# --Depends on: .computekME
#	              .getModuleGenes
#	              .getModGenesRankedBykME
# --
# --
.getNetworkBasics = function (net, data, suffix='.1') {
  assign(paste('dendro', suffix, sep=''), net$dendrograms, envir=.GlobalEnv)
  assign(paste('blockGenes', suffix, sep=''), net$blockGenes, envir=.GlobalEnv)
  assign(paste('colors', suffix, sep=''), net$colors, envir=.GlobalEnv)
  assign(paste('MEs', suffix, sep=''), net$MEs, envir=.GlobalEnv)
  assign(paste('kME', suffix, sep=''), .computekME(data, net$MEs)$all, envir=.GlobalEnv)
  assign(paste('modGenes', suffix, sep=''), .getModuleGenes(data, net$colors), envir=.GlobalEnv)
  assign(paste('modkMEs', suffix, sep=''), 
         .getModGenesRankedBykME(names(table(net$colors)),
                                 net$colors,
                                 .computekME(data, net$MEs)$all), 
         envir=.GlobalEnv)
}

# ---------------------------------------------------------------------------------------------------------------
# --Compute eigengene based connectivity (kME) for all genes/transcripts in all modules
# --Depends on:	WGCNA::orderMEs()
#				        WGCNA::corPvalueStudent()
# --Arguments
#     DATA: data frame of numeric expression values, rows=samples, cols=genes/transcripts
#     MEs: data frame of module eigengenes, as computed by WGCNA::blockwiseModules() or WGCNA::moduleEigengenes()
#     combine:
# --Value
#     list containing: $k = data frame of kME for each network member (rows) in each module (columns)
#                      $pval = data frame of pvalues corresponding to values in $k
.computekME = function(DATA, MEs, modules=F, colors=NULL) {
  MEs = orderMEs(MEs)
  kME = as.data.frame(cor(DATA, MEs, use='p'))
  names(kME) = paste('k', names(kME), sep='')
  kMEpval = as.data.frame(corPvalueStudent(as.matrix(kME), nrow(DATA)))
  names(kMEpval) = paste('p.', names(kMEpval), sep='')
  
  out = as.data.frame(matrix(nrow=nrow(kME), ncol=(2*ncol(kME))))
  kcols = seq(1, ncol(out), 2) 
  pcols = kcols+1
  out[, kcols] = kME
  names(out)[kcols] = names(kME)
  out[, pcols] = kMEpval
  names(out)[pcols] = names(kMEpval)
  rownames(out) = rownames(kME)
  
  modList = NULL
  if (modules) {
    modList = list()
    modNames = names(table(colors))
    for (mod in 1:length(modNames)) {
      modList[[mod]] = out[colors==modNames[mod], ]
      names(modList)[mod] = modNames[mod]
    }
  }
  
  return(list(all=out, mods=modList))
}

# ---------------------------------------------------------------------------------------------------------------
# --Get names of genes in each module
# --Depends on: NA
# --Arguments
#     DATA: data frame of numeric expression values, rows=samples, cols=genes/transcripts
#     colors: character vector of module colors, length(colors) must equal ncol(DATA)		
# --Value
#     modGenes:	list of character vectors containing module gene ids	
.getModuleGenes = function(DATA, colors) {
  modNames = names(table(colors))
  modGenes = list()
  for (mod in modNames) {
    modGenes[[mod]] = names(DATA)[colors==mod]
  }
  return(modGenes)
}

# ---------------------------------------------------------------------------------------------------------------
# --
# --Depends on: .getModulekME
# --
# --
.getModGenesRankedBykME = function (module_names, colors, kME) {
  outList = list()
  for (m in 1:length(module_names)) {
    outList[[m]] = .getModulekME(module_names[m],colors,kME)
  }
  names(outList) = module_names
  return(outList)
}

# ---------------------------------------------------------------------------------------------------------------
# --
# --Depends on: NA
# --
# --
.getModulekME = function(module, colors, kMEtable, kMEpvals=NULL, orderByRank=T) {
  modGenes = colors==module
  if (is.null(kMEpvals)) {
    #modCols = grep(module, names(kMEtable))#######
    modCols = which(gsub('ME|kME|p.kME', '', names(kMEtable)) == module)
    out = kMEtable[modGenes, modCols]
  } else {
    modCol = match(module, gsub('kME', '', names(kMEtable)))
    out = cbind(kMEtable[modGenes, modCol], kMEpvals[modGenes, modCol])
    dimnames(out) = list(rownames(kMEtable)[modGenes], 
                         c(names(kMEtable)[modCol], names(kMEpvals)[modCol]))
  }
  if (orderByRank) {
    out = out[order(out[, 1], decreasing=T), ]
  }
  return(out)
}