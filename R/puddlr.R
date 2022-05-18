# puddlr.R
#     Copyright (C) 2022  Rahul Dhodapkar
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

################################################################################
##
################################################################################

#' Create a 'puddlr' object.
#'
#' @param response vector of observed values for the response variable.
#'                 the response variable type must be compatible with the
#'                 `glm` family expected.  For example, for binomial (logistic)
#'                 regression, response should be a factor, while for
#'                 gaussian (ordinary least squares), response should be a
#'                 numeric vector.
#' @param predictors matrix of predictors, formatted as observations (rows)
#'                   by variables (columns)
#'
#' @return a Puddlr object
#'
#' @rdname CreatePuddlrObject
#' @export CreatePuddlrObject
#'
CreatePuddlrObject <- function(response, predictors) {
    # ***TODO*** add guards to do input validation
    return (
        structure(
            .Data = list(
                raw.response=response,
                raw.predictors=predictors
            ), class='puddlr'
        )
    )
}

#' Center and normalize predictors for linear dimensionality reduction.
#' Also prunes all observations where NA values are identified for either
#' the response or any of the predictors.
#'
#' @param puddlr A puddlr object to normalize in-place
#' @param scale whether to scale predictor matrix
#' @param center whether to center predictor matrix
#' @param max.na.ratio exclude predictors where the percent of NA observations
#'                     is \code{> max.na.ratio}. Default = \code{0.5}.
#'
#' @return A normalized puddlr object
#'
#' @rdname NormalizePredictors
#' @export NormalizePredictors
#'
NormalizePredictors <- function(puddlr,
                                scale=TRUE,
                                center=TRUE,
                                max.na.ratio=0.5) {

    na.columns <- ((
            colSums(is.na(puddlr$raw.predictors)) / nrow(puddlr$raw.predictors)
        ) >= max.na.ratio)

    puddlr$predictors <- puddlr$raw.predictors[,!na.columns]

    observations.with.na.values <- ( 
        (is.na(puddlr$raw.response))
        | (rowSums(is.na(puddlr$predictors)) > 0)
    )

    puddlr$response <- puddlr$raw.response[!observations.with.na.values]
    puddlr$predictors <- puddlr$predictors[!observations.with.na.values,]

    # scale and center predictors
    puddlr$predictors <- scale(puddlr$predictors, scale=scale, center=center)

    # remove predictor columns where there is 0 variance (all values equal)
    nan.columns <- colSums(is.nan(puddlr$predictors)) > 0
    puddlr$predictors <- puddlr$predictors[,!nan.columns]

    return(puddlr)
}

#' Calculates dispersion of predictors, using supplied measures.  Defaults to
#' using median absolute deviation (MAD).
#'
#' @param puddlr A puddlr object to calculate dispersion of predictors
#' @param dispersion.measure whether to scale predictor matrix. Default
#'               value = 'MAD' (Median Absolute Deviation).
#'
#' @return A puddlr object with dispersion calculated and stored in the
#'         'predictor.meta' attribute.
#'
#' @importFrom stats mad
#'
#' @rdname CalculatePredictorDispersion
#' @export CalculatePredictorDispersion
#'
CalculatePredictorDispersion <- function(puddlr,
                                         dispersion.measure='MAD') {

    if (!"predictors" %in% attributes(puddlr)$names) {
        stop('ERROR: cannot run "CalculatePredictorDispersion" without first running "NormalizePredictors".')
    }

    if (dispersion.measure == 'MAD') {
        dispersion.estimates <- apply(X=puddlr$predictors, MARGIN=2, FUN=mad, na.rm=TRUE)
    } else {
        stop(paste0('ERROR: invalid "dispersion.measure [', dispersion.measure, '] provided.'))
    }

    if (!"predictor.meta" %in% attributes(puddlr)$names) {
        puddlr$predictor.meta <- data.frame(
            predictor.name = colnames(puddlr$predictors)
        )
    }
    puddlr$predictor.meta[,dispersion.measure] <- dispersion.estimates

    return(puddlr)
}

#' Masks normalized predictors from downstream analysis, either by a character
#' vector of predictor names, or a boolean vector where TRUE indicates that a
#' predictor should be masked.  Masking predictors will remove all
#' dimensionality reductions from a puddlr object, and these will need to be
#' re-computed if required.
#'
#' @param puddlr A puddlr object to calculate dispersion of predictors
#' @param predictors.to.keep either a character vevtor or boolean vector containing
#'                           which predictors to keep. Indexed from the
#'                           "predictors" attribute
#'
#' @return A puddlr object with predictors masked. The "raw.predictors"
#'         attribute will remain unchanged.
#'
#' @rdname MaskPredictors
#' @export MaskPredictors
#'
MaskPredictors <- function(puddlr, predictors.to.keep) {

    if (!"predictors" %in% attributes(puddlr)$names) {
        stop('ERROR: cannot run "CalculatePredictorDispersion" without first running "NormalizePredictors".')
    }

    # remove predictor-dependent attributes
    attr(puddlr, 'predictor.meta') <- NULL
    attr(puddlr, 'reductions') <- NULL

    # perform predictor selection
    puddlr$predictors <- puddlr$predictors[,predictors.to.keep]

    return(puddlr)
}

#' Runs PCA for linear dimensionality reduction on the predictors of a 
#' normalized puddlr object
#'
#' @param puddlr    A puddlr object
#' @param center    whether to zero-center
#' @param scale     whether to center data to unit variance
#'
#' @return puddlr object with pca reduction
#'
#' @importFrom stats prcomp
#'
#' @rdname RunPCA
#' @export RunPCA
#'
RunPCA <- function(puddlr, center=TRUE, scale=TRUE) {
    pca <- prcomp(x=puddlr$predictors, center = center, scale = scale)

    if (is.null(puddlr$reductions)) {
        puddlr$reductions <- list()
    }
    puddlr$reductions$pca <- list()

    puddlr$reductions$pca$embedding <- pca$x
    puddlr$reductions$pca$rotation <- pca$rotation

    puddlr$reductions$pca$misc <- list()
    puddlr$reductions$pca$misc$sdev <- pca$sdev

    return(puddlr)
}

#' Runs ICA for linear dimensionality reduction on the predictors of a 
#' normalized puddlr object
#'
#' @param puddlr    A puddlr object
#' @param nc        number of components to extract. Default = 20
#'
#' @return puddlr object with pca reduction
#'
#' @importFrom ica icafast
#'
#' @rdname RunICA
#' @export RunICA
#'
RunICA <- function(puddlr, nc=20) {
    ica <- icafast(X=puddlr$predictors, nc=nc)

    t(ica$Q) %*% ica$R

    if (is.null(puddlr$reductions)) {
        puddlr$reductions <- list()
    }
    puddlr$reductions$ica <- list()

    puddlr$reductions$ica$embedding <- ica$S
    puddlr$reductions$ica$rotation <- t(ica$Q) %*% ica$R

    return(puddlr)
}

#' Runs generalized linear model on a puddlr object.  Intended to be run after
#' normalization and removal of incomplete observations.
#'
#' @param puddlr    A puddlr object
#' @param formula   an object of class "formula" (or one that can be coerced to
#'                  that class): a symbolic description of the model to be 
#'                  fitted. The details of model specification can be found 
#'                  under `glm`.  The name of the response variable in the
#'                  formula will be used to name the response variable.
#' @param family    a description of the error distribution and link function
#'                  to be used in the model.  Passed to the argument of the
#'                  same name under `glm`.
#' @param reduction a string specifying the linear dimensionality reduction to
#'                  use. Valid options are \code{'pca'}.  Default = \code{'pca'}
#' @param n.components an integer specifying the total number of reduction
#'                  components to use in the GLM.
#' @param adj.rsq   boolean flag specifying whether to adjust the pseudo-R^2
#'                  goodness-of-fit calculated value.
#'
#' @return puddlr object with glm
#'
#' @importFrom DescTools PseudoR2
#' @importFrom stats as.formula glm pnorm vcov
#'
#' @rdname RunGLM
#' @export RunGLM
#'
RunGLM <- function(puddlr,
                   formula,
                   family,
                   reduction,
                   n.components,
                   adj.rsq=FALSE) {

    response.var.name <- all.vars(formula)[1]

    data.df <- data.frame(puddlr$reductions[[reduction]]$embedding[,1:n.components])
    data.df[,response.var.name] <- puddlr$response

    puddlr$model <- list()

    puddlr$model$obj <- glm(
      formula = as.formula(paste0(response.var.name, '~.')),
      family = family,
      data = data.df
    )

    puddlr$model$rsq <- PseudoR2(puddlr$model$obj, which='McFadden')

    # project back to original feature space
    sum.feature.estimate <- (
        rowSums(t(
            t(puddlr$reductions[[reduction]]$rotation[,1:n.components])
            * puddlr$model$obj$coefficients[(1:n.components)+1]
        ))
    )

    sum.feature.var <- c()
    for (i in 1:ncol(puddlr$predictors)) {
        vconv.unscaled <- vcov(puddlr$model$obj)[
                                1:n.components+1,
                                1:n.components+1]
        v <- puddlr$reductions[[reduction]]$rotation[i,1:n.components]
        var.i <- sum(t(vconv.unscaled * v) * v)
        sum.feature.var <- c(sum.feature.var, var.i)
    }
    names(sum.feature.var) <- colnames(puddlr$predictors)

    stats.df <- data.frame(
        VariableName = names(sum.feature.var),
        GLMEstimate = sum.feature.estimate,
        StdErr = sqrt( sum.feature.var )
    )
    stats.df$ZScore <- stats.df$GLMEstimate / stats.df$StdErr
    stats.df$PValue <- 2*pnorm(-abs(stats.df$ZScore))

    puddlr$model$df <- stats.df

    return(puddlr)
}

#' Iteratively runs GLM fitting on a set number of components to identify the
#' appropriate number of dimensions to restrict the model to.  Performs k-fold
#' cross validation to assess for overfitting with averaged root mean squared
#' error, and pseudo-R squared (McFadden's) to assess for goodness of fit.
#' Data for the search are saved in the puddlr object for later plotting and
#' visual inspection
#'
#' @param puddlr    A puddlr object
#' @param n.to.scan.vec Vector of possible values for n.components to scan for
#'                  GLM fitting.
#' @param k.cross   number of train-test splits. (corresponds to k-fold cross
#'                  validation).  Default = 7.
#' @param rand.seed random seed for cross validation split, default=42.
#' @param formula   an object of class "formula" (or one that can be coerced to
#'                  that class): a symbolic description of the model to be 
#'                  fitted. The details of model specification can be found 
#'                  under `glm`.  The name of the response variable in the
#'                  formula will be used to name the response variable.
#' @param family    a description of the error distribution and link function
#'                  to be used in the model.  Passed to the argument of the
#'                  same name under `glm`.
#' @param reduction a string specifying the linear dimensionality reduction to
#'                  use. Valid options are \code{'pca'}.  Default = \code{'pca'}
#' @param adj.rsq   boolean flag specifying whether to adjust the pseudo-R^2
#'                  goodness-of-fit calculated value.
#'
#' @return puddlr object with scanned n components by rsq plot.
#'
#' @importFrom caret trainControl
#' @importFrom stats predict
#'
#' @rdname ScanComponentsSubset
#' @export ScanComponentsSubset
#'
ScanComponentsSubset <- function(puddlr,
                                 n.to.scan.vec,
                                 k.cross = 7,
                                 rand.seed = 42,
                                 formula,
                                 family,
                                 reduction,
                                 adj.rsq) {
    set.seed(rand.seed)

    avg.train.rsq <- c()
    avg.rmse <- c()
    for(n.components in n.to.scan.vec) {
        print(paste0("Testing with [", 
                         as.character(n.components) ,"] components"))
        #Create equally size folds
        folds <- cut(sample(1:length(puddlr$response)),breaks=k.cross,labels=FALSE)
        #Perform k-fold cross validation
        rmse.values <- c()
        rsq.values <- c()
        for(k.i in 1:k.cross){
            #Segement data by fold
            test.ixs <- which(folds==k.i,arr.ind=TRUE)
            train.ixs <- which(folds!=k.i,arr.ind=TRUE)

            # train data
            train.puddlr <- puddlr
            train.puddlr$reductions[[reduction]]$rotation <- (
                train.puddlr$reductions[[reduction]]$rotation)
            train.puddlr$reductions[[reduction]]$embedding <- (
                train.puddlr$reductions[[reduction]]$embedding[train.ixs,])

            train.puddlr$predictors <- (train.puddlr$predictors[train.ixs,])
            train.puddlr$response <- (train.puddlr$response[train.ixs])

            # test data
            test.puddlr <- puddlr
            test.puddlr$reductions[[reduction]]$rotation <- (
                test.puddlr$reductions[[reduction]]$rotation)
            test.puddlr$reductions[[reduction]]$embedding <- (
                test.puddlr$reductions[[reduction]]$embedding[test.ixs,])

            test.puddlr$predictors <- (test.puddlr$predictors[test.ixs,])
            test.puddlr$response <- (test.puddlr$response[test.ixs])

            train.puddlr <- RunGLM(train.puddlr, 
                 formula=y ~ .,
                 family=family,
                 reduction=reduction,
                 n.components=n.components,
                 adj.rsq=FALSE)
            y.pred <- predict(train.puddlr$model$obj, 
                data.frame(test.puddlr$reductions[[reduction]]$embedding))
            rmse.values <- c(rmse.values, 
                sqrt(mean((test.puddlr$response - y.pred)^2)))
            rsq.values <- c(rsq.values, train.puddlr$model$rsq)
        }
        avg.rmse <- c(avg.rmse, mean(rmse.values))
        avg.train.rsq <- c(avg.train.rsq, mean(rsq.values))
    }

    if (is.null(puddlr$optimization)) {
        puddlr$optimization <- list()
    }

    local.optim.df <- data.frame(
        n.components = n.to.scan.vec,
        family = family$family,
        reduction = reduction,
        k.folds = k.cross,
        rmse = avg.rmse,
        avg.train.rsq = avg.train.rsq
    )

    if (is.null(puddlr$optimization$n.components.scan.df)) {
        puddlr$optimization$n.components.scan.df <- local.optim.df
    } else {
        for (i in 1:nrow(local.optim.df)) {
            matching.ix <- which(
                (puddlr$optimization$n.components.scan.df$n.components
                    == local.optim.df$n.components[[i]])
                & (puddlr$optimization$n.components.scan.df$family
                    == local.optim.df$family[[i]])
                & (puddlr$optimization$n.components.scan.df$reduction
                    == local.optim.df$reduction[[i]])
                & (puddlr$optimization$n.components.scan.df$k.folds
                    == local.optim.df$k.folds[[i]])
            )
            if (length(matching.ix) == 0) {
                puddlr$optimization$n.components.scan.df <- rbind(
                    puddlr$optimization$n.components.scan.df,
                    local.optim.df[i,]
                )
            } else if (length(matching.ix) == 1) {
                puddlr$optimization$n.components.scan.df[matching.ix,] <- (
                    local.optim.df[i,])
            } else {
                warning('WARN:ScanComponentsSubset:: encountered illegal df contents.')
            }
        }
    }

    return(puddlr)
}

