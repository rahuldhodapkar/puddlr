# puddlr.R
#     Copyright (C) 2021  Rahul Dhodapkar
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


#' Create Puddlr object
#'
#' @param response vector of observed values for the response variable
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
#' @importFrom stats na.omit
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
        ) > max.na.ratio)

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

################################################################################
## Methods
################################################################################

#' Pretty print details of puddlr object
#'
#' @param puddlr A puddlr object to normalize in-place
#'
#' @return NULL
#'
#' @rdname print.puddlr
#' @export print.puddlr
#'
print.puddlr <- function(puddlr) {
    print("A Puddlr Object")
}

