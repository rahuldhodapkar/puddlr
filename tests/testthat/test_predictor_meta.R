# test_glm.R
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


context("Test Predictor Meta-information Extraction")
library(puddlr)
library(stringr)


test_that("Test Median Absolute Dispersion (MAD) Calculation", {
    response.test <- as.factor(c(1,0,1,1,0,0))
    predictors.test <- matrix(c( 0, 0, 0, 0, 0, 1,
                                 1, 3, 2, 1, 3, 2,
                                 3, 3, 1, 0, 2, 1,
                                 2, 1, 2, 1, 1, 2,
                                 1, 6, 5, 1, 3, 5,
                                 40, 20, 30, 40, 30, 60),
                              nrow=6, ncol=6)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- CalculatePredictorDispersion(puddlr.obj, dispersion.measure='MAD')

    expect_equal(puddlr.obj$predictor.meta$MAD[1], 0)
})


test_that("Test Predictor Masking by Boolean Vector", {
    response.test <- as.factor(c(1,0,1,1,0,0))
    predictors.test <- matrix(c( 0, 0, 0, 0, 0, 1,
                                 1, 3, 2, 1, 3, 2,
                                 3, 3, 1, 0, 2, 1,
                                 2, 1, 2, 1, 1, 2,
                                 1, 6, 5, 1, 3, 5,
                                 40, 20, 30, 40, 30, 60),
                              nrow=6, ncol=6)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- MaskPredictors(puddlr.obj, colnames(puddlr.obj$predictors) %in% c('c1', 'c2'))

    expect_equal(colnames(puddlr.obj$predictors), c('c1', 'c2'))
})


test_that("Test Predictor Masking by Character Vector", {
    response.test <- as.factor(c(1,0,1,1,0,0))
    predictors.test <- matrix(c( 0, 0, 0, 0, 0, 1,
                                 1, 3, 2, 1, 3, 2,
                                 3, 3, 1, 0, 2, 1,
                                 2, 1, 2, 1, 1, 2,
                                 1, 6, 5, 1, 3, 5,
                                 40, 20, 30, 40, 30, 60),
                              nrow=6, ncol=6)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- MaskPredictors(puddlr.obj, c('c1', 'c4'))

    expect_equal(colnames(puddlr.obj$predictors), c('c1', 'c4'))
})


