# test_preprocessing.R
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


context("Test data transformations")
library(puddlr)

test_that("Test 'puddlr' constructor", {
    response.test <- c(1,0,1,1)
    predictors.test <- matrix(c( 1, 2, 3, 4, 
                                 4, 3, 2, 1, 
                                 0, 0, 0, 0, 
                                 NA, NA, 2, 1,
                                 5, 5, 5, 5,
                                 10, 20, 30, 40),
                              nrow=4, ncol=6)
    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    expect_equal(class(puddlr.obj), 'puddlr')
})

test_that("Test 'NormalizePredictors'", {
    response.test <- c(1,0,1,1)
    predictors.test <- matrix(c( 1, 2, 3, 4, 
                                 4, 3, 2, 1, 
                                 0, 0, 0, 0, 
                                 NA, NA, 2, 1,
                                 5, 5, 5, 5,
                                 10, 20, 30, 40),
                              nrow=4, ncol=6)
    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj,
                                scale=TRUE,
                                center=TRUE,
                                max.na.ratio=0.4)

    expect_equal(ncol(puddlr.obj$predictors), 3)
    expect_equal(puddlr.obj$predictors[,1], puddlr.obj$predictors[,3])

    puddlr.obj <- CreatePuddlrObject(
        response = response.test,
        predictors = predictors.test
    )
    puddlr.obj <- NormalizePredictors(puddlr.obj,
                                scale=TRUE,
                                center=TRUE,
                                max.na.ratio=0.9)
    expect_equal(ncol(puddlr.obj$predictors), 4)
})
