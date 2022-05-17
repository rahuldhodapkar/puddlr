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


context("Test GLM")
library(puddlr)
library(stringr)

test_that("Test PCA", {
    response.test <- as.factor(c(1,0,1,1))
    predictors.test <- matrix(c( 1, 2, 3, 4, 
                                 4, 3, 2, 1, 
                                 0, 3, 1, 0, 
                                 0.5, 1, 2, 1,
                                 5, 6, 5, 1,
                                 10, 20, 30, 40),
                              nrow=4, ncol=6)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- RunPCA(puddlr.obj)
    puddlr.obj <- RunGLM(puddlr.obj, 
                         formula=response ~ .,
                         family=binomial(link='logit'),
                         reduction='pca',
                         n.components=3)

    expect_equal(class(puddlr.obj), 'puddlr')
})

test_that("Test 'gaussian' GLM", {
    response.test <- c(0.3, 0.1, 2.2, 4.3)
    predictors.test <- matrix(c( 1, 2, 3, 4, 
                                 4, 3, 2, 1, 
                                 0, 3, 1, 0, 
                                 0.5, 1, 2, 1,
                                 5, 6, 5, 1,
                                 10, 20, 30, 40),
                              nrow=4, ncol=6)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- RunPCA(puddlr.obj)
    puddlr.obj <- RunGLM(puddlr.obj, 
                         formula=response ~ .,
                         family=gaussian(link='identity'),
                         reduction='pca',
                         n.components=3)

    expect_equal(class(puddlr.obj), 'puddlr')
})

test_that("Test 'ScanComponentsSubset' with 'binomial' GLM", {
    set.seed(42)

    response.test <- sample(c(0,1), size=100, replace=TRUE)
    predictors.test <- matrix(runif(100 * 1000),
                              nrow=100, ncol=1000)

    rownames(predictors.test) <- paste0('r', 1:nrow(predictors.test))
    colnames(predictors.test) <- paste0('c', 1:ncol(predictors.test))

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj)
    puddlr.obj <- RunPCA(puddlr.obj)
    puddlr.obj <- RunGLM(puddlr.obj, 
                         formula=y ~ .,
                         family=binomial(link='logit'),
                         reduction='pca',
                         n.components=3,
                         adj.rsq=FALSE)
    puddlr.obj <- ScanComponentsSubset(puddlr.obj,
                         n.to.scan.vec=3:5,
                         k.cross=7,
                         formula=y ~ .,
                         family=binomial(link='logit'),
                         reduction='pca',
                         adj.rsq=FALSE)
    expect_equal(nrow(puddlr.obj$optimization$n.components.scan.df), 3)

    puddlr.obj <- ScanComponentsSubset(puddlr.obj,
                         n.to.scan.vec=3:7,
                         k.cross=7,
                         formula=y ~ .,
                         family=binomial(link='logit'),
                         reduction='pca',
                         adj.rsq=FALSE)
    expect_equal(nrow(puddlr.obj$optimization$n.components.scan.df), 5)

    expect_equal(class(puddlr.obj), 'puddlr')
})


