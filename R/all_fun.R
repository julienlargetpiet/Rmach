#' @import stringr
#' @title Rmach

#' poly_model
#'
#' @description Take a datasets of x and y values and a function tha could fit all the data with the missing coefficients, and returns a list containing the coefficients that fit the best the data for a given function, as a vector for the first index, and  at the second index, the actual sum of difference between each data point and the function at the same x values.
#'
#' @param inpt_datf is the input data as a dataframe, first column is the x values and the second is the y values
#' @param degree is how many coefficients will be involved (each coefficient multiplies either an x to the power of something, an exponential of something or a base something logarithm for a something value)
#' @param twk_val is the value used for finding the best coefficients, it is directly linked to the accuracy of the coefficients, see the description for more information. Defaults to (max(yval) - min(yval)) / n
#' @param sensi_val is the value from which two variations of a coefficient brings a so small accuracy contribution that the algorythm does not continue to find better coefficients. For example, if i set sensi_val = 0.001, so if coefficients alpha1 and beta1 brings a total difference between the function and the actual data of 10.8073 and then the algorythm find alpha2 and beta1 that brings a total difference equal to 10.8066, so the algorythm will stop running. But the coefficients returned will still be the best, that is alpha2 and beta1
#' @param coeff_v is a vector containing the original coefficients for the function, so the closest those are from the best one, the fastest the algorythm will compute the best coefficients. The first value of coeff is always the constant.
#' @param mth_symb is a vector containing the elemnts linked to the coefficients from the second element. It can be x, e (exp(x)) or log-X (log(x)-base), and their reverse like 1/x. If the numerator varies the element should be entered like tis list/x, list/e or list/log-base. See numrtr_v for the values related to list  
#' @param powers is a vector containing the exponent, or related value to mth_symb. powers can be a vector if those values are constants or it could be a list of vectors the length of observed individuals, if those values varies like in the examples. Notthat if you use variables in powers (list), each values of a vector from this list has to be at the exact same x coordinates of each observed individuals in the input dataframe. Ex: datf <- data.frame("x"=c(4, 4, 3, 2, 1, 1), "y"=c(1:6)), so vector(s) from powers that contain varying value must be of length 4. Also, the values are not ascendly sorted, don't worry values are ascendly sorted under the hood, so fill your powers vectors in the intuitive ascendly way 
#' @param numrtr_v is a vector containing the values for the numerator related to mth_symb if on element is like this: list/x or list/e
#' @examples
#'
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -5), powers=c(1), mth_symb=c("x"),
#' 
#'                  numrtr_v=NA))
#' 
#' [[1]]
#' [1] 33.234375 -4.265625
#' 
#' [[2]]
#' [1] 74.78275
#' 
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -3), powers=c(1), mth_symb=c("x"),
#' 
#'                  numrtr_v=NA))
#' 
#' [[1]]
#' [1] 31.765625 -3.734375
#' 
#' [[2]]
#' [1] 80.36228
#' 
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -3), powers=list(c(1:length(mtcars$wt))), mth_symb=c("x"),
#' 
#'                  numrtr_v=NA))
#' 
#' [[1]]
#' [1] 32.5 -3.0
#' 
#' [[2]]
#' [1] 1.067436e+24
#' 
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -3), powers=list(c(1:length(mtcars$wt))), mth_symb=c("list/x"),
#' 
#'                  numrtr_v=list(c(length(mtcars$wt):1))))
#' [[1]]
#' [1] 19.28125 -0.06250
#' 
#' [[2]]
#' [1] 35839.44
#'
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -3), powers=list(c(1:length(mtcars$wt))), mth_symb=c("1/x"),
#' 
#'                  numrtr_v=NA))
#' 
#' [[1]]
#' [1] 27.359375 -8.140625
#' 
#' [[2]]
#' [1] 160.2263
#' 
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=1, coeff_v=c(32.5), powers=NA, mth_symb=NA,
#' 
#'                  numrtr_v=NA))
#'
#' [[1]]
#' [1] 19.28125
#' 
#' [[2]]
#' [1] 148.7625
#'
#' print(poly_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), degree=2, coeff_v=c(32.5, -3, 2), powers=list(c(1:length(mtcars$wt)), 2), mth_symb=c("1/x", "x"),
#'  
#'                   numrtr_v=NA))
#' 
#' [[1]]
#' [1]  0.921875 -5.203125  2.000000
#' 
#' [[2]]
#' [1] 455.6017
#'
#' @export

poly_model <- function(inpt_datf, degree, twk_val=NA, sensi_val=twk_val, coeff_v=NA, 

                        powers=NA, mth_symb=c("x"), numrtr_v=NA){ 

   calc_fun <- function(x, coeff_v, powers){

           calc_fun2 <- function(inpt_val, mth_symb, x, numrtr_idx){

              if (!(str_detect(string=mth_symb, pattern="/"))){

                      return(c(inpt_val, numrtr_idx))

              }else{

                      numrtr <- unlist(strsplit(x=mth_symb, split="/"))[1]

                      if (numrtr == "list"){

                          numrtr_idx = numrtr_idx + 1

                          return(unlist(numrtr_v[numrtr_idx])[match(x=x, table=inpt_datf[, 1])]/inpt_val)

                      }else{

                          return(c(1/inpt_val, numrtr_idx))

                      }

              }

           }

           calc_fun = coeff_v[1]

           numrtr_idx = 0

           if (length(coeff_v) > 1){

                   if (typeof(powers) != "list"){

                      for ( i in 2:length(coeff_v)){

                        if (str_detect(string=mth_symb[(i-1)], pattern="x")){

                           cur_calc <- calc_fun2(inpt_val=powers[(i-1)], mth_symb=mth_symb[(i-1)], 
                           x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * x ** cur_calc[1]

                           numrtr_idx = cur_calc[2]

                       }else if (str_detect(string=mth_symb[(i-1)], pattern="e")){

                           cur_calc <- calc_fun2(inpt_val=exp(powers[(i-1)]), mth_symb=mth_symb[(i-1)],
                           x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                           numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="log")){

                           cur_calc <- calc_fun2(inpt_val=log(x=powers[(i-1)], base=unlist(strsplit(mth_symb[(i-1)], split="-"))[2]), mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                           numrtr_idx = cur_calc[2]

                        }

                      }

                   }else{

                      for ( i in 2:length(coeff_v)){

                        cur_powers <- unlist(powers[(i-1)]) 

                        if (length(cur_powers) > 1){

                                cur_powers_idx <- match(x=x, table=inpt_datf[, 1])

                        }else{

                                cur_powers_idx <- 1

                        }

                        if (str_detect(string=mth_symb[(i-1)], pattern="x")){

                            cur_calc <- calc_fun2(inpt_val=cur_powers[cur_powers_idx], mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[i] * x ** cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="e")){

                            cur_calc <- calc_fun2(inpt_val=exp(cur_powers[cur_powers_idx]), mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="log")){

                            cur_calc <- calc_fun2(inpt_val=log(x=cur_powers[cur_powers_idx], base=unlist(strsplit(x=mth_symb[(i-1)], split="-"))[2]), mth_symb=mth_symb[(i-1), x=x], numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[I] * cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }

                      }

                   }

           }

        return(calc_fun)

   }

   if (all(is.na(powers))){

           powers <- seq(from=1, to=degree, by=1)

   }

   if (all(is.na(coeff_v))){

           coeff_v <- c()

           for (i in 1:(degree+1)){

                   coeff_v <- c(coeff_v, 1)

           }

   }

   if (is.na(twk_val)){

        twk_val <- (max(inpt_datf[, 2]) - min(inpt_datf[, 2])) / nrow(inpt_datf)

   }

   inpt_datf[, c(1, 2)] <- inpt_datf[match(x=sort(x=inpt_datf[, 1], decreasing=FALSE), table=inpt_datf[, 1]), c(1, 2)]

   lst_diff  = 0

   for (pt in 1:nrow(inpt_datf)){

        lst_diff = lst_diff + abs((calc_fun(x=inpt_datf[, 1][pt], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

   }

   cur_diff = lst_diff - sensi_val

   no_stop = TRUE

   while (no_stop & (lst_diff - cur_diff) >= sensi_val){

           snapshot_coeff_v <- coeff_v

           lst_diff = cur_diff

           cur_diff = 0

           coeff_v[1] = coeff_v[1] + twk_val

           for (pt in 1:nrow(inpt_datf)){

                cur_diff = cur_diff + abs((calc_fun(x=inpt_datf[pt, 1], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

           }

           cur_diff2 = 0

           coeff_v[1] = coeff_v[1] - 2 * twk_val

           for (pt in 1:nrow(inpt_datf)){

                cur_diff2 = cur_diff2 + abs((calc_fun(x=inpt_datf[pt, 1], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

           }

           if (cur_diff < cur_diff2){

                coeff_v[1] = coeff_v[1] + 2 * twk_val

           }else if (cur_diff2 < cur_diff){

                   cur_diff <- cur_diff2

           }

           if (degree > 1){

                   for (I in 2:degree){

                        cur_diff2b = 0

                        coeff_v[I] = coeff_v[I] + twk_val

                        for (pt in 1:nrow(inpt_datf)){

                                cur_diff2b = cur_diff2b + abs((calc_fun(x=inpt_datf[pt, 1], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

                        }

                        cur_diff2 = 0

                        coeff_v[I] = coeff_v[I] - 2 * twk_val

                        for (pt in 1:nrow(inpt_datf)){

                                cur_diff2 = cur_diff2 + abs((calc_fun(x=inpt_datf[pt, 1], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

                        }

                        if (cur_diff2b < cur_diff2){

                            coeff_v[I] = coeff_v[I] + 2 * twk_val

                            cur_diff2 = cur_diff2b

                        } 

                        if (cur_diff2 < cur_diff){

                                cur_diff <- cur_diff2

                        }

                   }

           }

       if (cur_diff > lst_diff){

               no_stop = FALSE

               coeff_v <- snapshot_coeff_v

       }

   }

   return(list(coeff_v, lst_diff))

}

#' best_model
#' 
#' Returns the best input models. The coefficient of the best model can be found with the poly_model function  
#'
#' @param inpt_datf is the input dataframe, first column for the x values and second column for the y values 
#' @param Degree is a vector containing all the degrees. Each degree represents how many coefficients the model has.
#' @param Coeff_v is a list containing the vector containing the coefficients for each model. The first value of each coefficient vector is always the constant, so it is not linked to any math symbol 
#' @param Powers is a list containing all the values associated with the math symbols of mth_symb list for each model. Because you can have multiple models in the function, so Powers is separated with the "-" separator between the different powers values for each model like in the examples
#' @param Mth_symb is a list containing the vector of the different math symbols linked to the coefficients from the second value
#' @param Numrtr_v is a list containing the different numerator values for each math symbol for each model, see examples
#' @examples
#'
#' print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", 32.5, -3, "-", 32.5, -5, "-"), Powers=c("-", 1, "-", 1, "-"), Mth_symb=c("-", "x", "-", "x", "-"), Numrtr_v=NA))
#'
#' [1] 2
#'
#' print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3), "-"), Powers=c("-", c(1), "-", c(1), "-"), Mth_symb=c("-", c("x"), "-", c("1/x"), "-"), Numrtr_v=NA))
#'
#' [1] 1
#'
#' print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3), "-"), Powers=c("-", c(1), "-", list(c(1:length(mtcars$wt))), "-"), Mth_symb=c("-", c("x"), "-", c("1/x"), "-"), Numrtr_v=NA))
#'
#' [1] 1
#' 
#' print(best_model(inpt_datf=data.frame(mtcars$wt, mtcars$mpg), Degree=c(2, 2), Coeff_v=c("-", c(32.5, -5), "-", c(32.5, -3, 2), "-"), Powers=c("-", c(1), "-", list(c(1:length(mtcars$wt)), 2), "-"), Mth_symb=c("-", c("x"), "-", c("list/x", "x"), "-"), Numrtr_v=c("-", list(c(1:length(mtcars$wt))), "-")))
#'
#' #' [1] 1
#'
#' @export

best_model <- function(inpt_datf, Degree, Coeff_v=NA, 

                        Powers=NA, Mth_symb, Numrtr_v=NA){ 

   calc_fun <- function(x, coeff_v, powers){

           calc_fun2 <- function(inpt_val, mth_symb, x, numrtr_idx){

              if (!(str_detect(string=mth_symb, pattern="/"))){

                      return(c(inpt_val, numrtr_idx))

              }else{

                      numrtr <- unlist(strsplit(x=mth_symb, split="/"))[1]

                      if (numrtr == "list"){

                          numrtr_idx = numrtr_idx + 1

                          return(unlist(numrtr_v[numrtr_idx])[match(x=x, table=inpt_datf[, 1])]/inpt_val)

                      }else{

                          return(c(1/inpt_val, numrtr_idx))

                      }

              }

           }

           calc_fun = coeff_v[1]

           numrtr_idx = 0

           if (length(coeff_v) > 1){

                   if (typeof(powers) != "list"){

                      for ( i in 2:length(coeff_v)){

                        if (str_detect(string=mth_symb[(i-1)], pattern="x")){

                           cur_calc <- calc_fun2(inpt_val=powers[(i-1)], mth_symb=mth_symb[(i-1)], 
                           x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * x ** cur_calc[1]

                           numrtr_idx = cur_calc[2]

                       }else if (str_detect(string=mth_symb[(i-1)], pattern="e")){

                           cur_calc <- calc_fun2(inpt_val=exp(powers[(i-1)]), mth_symb=mth_symb[(i-1)],
                           x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                           numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="log")){

                           cur_calc <- calc_fun2(inpt_val=log(x=powers[(i-1)], base=unlist(strsplit(mth_symb[(i-1)], split="-"))[2]), mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                           calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                           numrtr_idx = cur_calc[2]

                        }

                      }

                   }else{

                      for ( i in 2:length(coeff_v)){

                        cur_powers <- unlist(powers[(i-1)]) 

                        if (length(cur_powers) > 1){

                                cur_powers_idx <- match(x=x, table=inpt_datf[, 1])

                        }else{

                                cur_powers_idx <- 1

                        }

                        if (str_detect(string=mth_symb[(i-1)], pattern="x")){

                            cur_calc <- calc_fun2(inpt_val=cur_powers[cur_powers_idx], mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[i] * x ** cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="e")){

                            cur_calc <- calc_fun2(inpt_val=exp(cur_powers[cur_powers_idx]), mth_symb=mth_symb[(i-1)], x=x, numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[i] * cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }else if (str_detect(string=mth_symb[(i-1)], pattern="log")){

                            cur_calc <- calc_fun2(inpt_val=log(x=cur_powers[cur_powers_idx], base=unlist(strsplit(x=mth_symb[(i-1)], split="-"))[2]), mth_symb=mth_symb[(i-1), x=x], numrtr_idx=numrtr_idx)

                            calc_fun = calc_fun + coeff_v[I] * cur_calc[1]

                            numrtr_idx = cur_calc[2]

                        }

                      }

                   }

           }

        return(calc_fun)

   }

   cnt = match(table=Coeff_v, x="-") + 1

   Coeff_v[(cnt-1)] <- "?"

   coeff_v <- c()

   while (Coeff_v[cnt] != "-"){

           coeff_v <- append(x=coeff_v, values=as.numeric(unlist(Coeff_v[cnt])))

           cnt = cnt + 1

   }

   cnt = match(table=Powers, x="-") + 1

   Powers[(cnt-1)] <- "?"

   pre_mtc <- match(table=Powers, x="-")

   if ((pre_mtc - cnt) > 1){

           powers <- list()

           while (Powers[cnt] != "-"){

                   powers <- append(x=powers, values=Powers[cnt]) 

                   cnt = cnt + 1

           }

   }else{

           powers <- c()

           while (Powers[cnt] != "-"){

                   powers <- append(x=powers, values=as.numeric(unlist(Powers[cnt]))) 

                   cnt = cnt + 1

           }

   }

   cnt = match(table=Mth_symb, x="-") + 1

   Mth_symb[(cnt-1)] <- "?"

   mth_symb <- c()

   while (Mth_symb[cnt] != "-"){

           mth_symb <- append(x=mth_symb, values=Mth_symb[cnt])

           cnt = cnt + 1

   }

   if (!(all(is.na(Numrtr_v))) & length(grep(x=mth_symb, pattern="list")) > 0){

           cnt = match(table=Numrtr_v, x="-") + 1

           Numrtr_v[(cnt-1)] <- "?"

           numrtr_v <- list()

           while (Numrtr_v[cnt] != "-"){

                   numrtr_v <- append(x=numrtr_v, values=Numrtr_v[cnt]) 

                   cnt = cnt + 1

           }

   }

   degree <- Degree[1]

   if (length(powers) == 0){

           powers <- seq(from=1, to=degree, by=1)

   }

   if (length(coeff_v) == 0){

           coeff_v <- c()

           for (i in 1:(degree+1)){

                   coeff_v <- c(coeff_v, 1)

           }

   }

   inpt_datf[, c(1, 2)] <- inpt_datf[match(x=sort(x=inpt_datf[, 1], decreasing=FALSE), table=inpt_datf[, 1]), c(1, 2)]

   pre_lst_diff  = 0

   for (pt in 1:nrow(inpt_datf)){

        pre_lst_diff = pre_lst_diff + abs((calc_fun(x=inpt_datf[, 1][pt], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

   }

   coeff_rtn <- 1

   if (length(Degree) > 1){

           for (I in 2:length(Degree)){

                   cnt = match(table=Coeff_v, x="-") + 1

                   Coeff_v[(cnt-1)] <- "?"

                   coeff_v <- c()

                   while (Coeff_v[cnt] != "-"){

                           coeff_v <- append(x=coeff_v, values=as.numeric(unlist(Coeff_v[cnt]))) 

                           cnt = cnt + 1

                   }

                   cnt = match(table=Powers, x="-") + 1

                   Powers[(cnt-1)] <- "?"

                   pre_mtc <- match(table=Powers, x="-")
                   
                   if ((pre_mtc - cnt) > 1){

                           powers <- list()

                           while (Powers[cnt] != "-"){

                                   powers <- append(x=powers, values=Powers[cnt]) 

                                   cnt = cnt + 1

                           }

                   }else{

                           powers <- c()

                           while (Powers[cnt] != "-"){

                                   powers <- append(x=powers, values=as.numeric(unlist(Powers[cnt]))) 

                                   cnt = cnt + 1

                           }

                   }

                   cnt = match(table=Mth_symb, x="-") + 1

                   Mth_symb[(cnt-1)] <- "?"

                   mth_symb <- c()

                   while (Mth_symb[cnt] != "-"){

                           mth_symb <- append(x=mth_symb, values=Mth_symb[cnt]) 

                           cnt = cnt + 1

                   }

                   if (!(all(is.na(Numrtr_v))) & length(grep(x=mth_symb, pattern="list")) > 0){

                           cnt = match(table=Numrtr_v, x="-")[1] + 1

                           Numrtr_v[(cnt-1)] <- "?"

                           numrtr_v <- list()

                           while (Numrtr_v[cnt] != "-"){

                                   numrtr_v <- append(x=numrtr_v, values=Numrtr_v[cnt]) 

                                   cnt = cnt + 1

                           }

                   }

                   degree <- Degree[I]

                   if (length(powers) == 0){

                           powers <- seq(from=1, to=degree, by=1)

                   }

                   if (length(coeff_v) == 0){

                           coeff_v <- c()

                           for (i in 1:(degree+1)){

                                   coeff_v <- c(coeff_v, 1)

                           }

                   }

                   lst_diff  = 0

                   for (pt in 1:nrow(inpt_datf)){

                        lst_diff = lst_diff + abs((calc_fun(x=inpt_datf[, 1][pt], coeff_v=coeff_v, powers=powers)) - inpt_datf[pt, 2])

                   }

                   if (lst_diff < pre_lst_diff){

                       coeff_rtn <- I 

                   }

                   pre_lst_diff <- lst_diff

           }

   }

   return(coeff_rtn)

}

#' calcall
#'
#' Takes a formula as a character as an input and makes the calculation. Accepts also variables, in this case the part of the formula that contains the variable wont be calculated, but the others part will be as usual.
#' @param inpt is the input formula as a character
#' @examples
#'
#' print(calcall(inpt="ze+(yu*((fgf)-(-12+8-Y+4-T+4+97+a)+tt))"))
#' 
#' [1] "ze+(yu*(fgf-(-4-Y+4-T+101+a)+tt))"
#' 
#' print(calcall(inpt="ze+(yu*((fgf)-(-12+8-7+3-67+4+97+1)+tt))"))
#'
#' [1] "ze+(yu*(fgf-27+tt))"
#' 
#' print(calcall(inpt="ze+(yu*((fgf)+(12*3/2+4)+tt))"))
#'
#' [1] "ze+(yu*(fgf+22+tt))"
#' 
#' print(calcall(inpt="1+3*2+(-2/-3*-3*((fgf)-(--12-6)+2))+5-3*5"))
#'
#' [1] "7+(-2*(fgf-4))+20"
#' 
#' print(calcall(inpt="1+3*2+(-2/-3*-3*((fgf)-(--12-6)+2))+(-log_e_1_e_2+t+2^3)+m-log_e_1_e_2+2^3-m-6*2+(e_ii-2-6+log_im_4-67)+-6+2+(y-5+7)"))
#'
#' [1] "7+(-2*(fgf-4))+(-2+t+8)+m+6-m-12+(e_ii-8+log_im_4-67)-4+(y+2)"
#'
#' print(calcall("(6+4*-(4-5))+3/3"))
#'
#' [1] "11" 
#' 
#' print(calcall(inpt="1+3*2+(-2/-3*-3*((fgf)-(--12-6)+2))+(-log_e_1_e_2+t+2^3)+m-log_e_1_e_2+2^3-m-6*2+-6+2"))
#'
#' [1] "7+(-2*(fgf-4))+(-2+t+8)+m+6-m-16"
#' 
#' print(calcall(inpt="(log_5_Z-2-6+5)+-6+2"))
#'
#' [1] "(log_5_Z-3)-4"
#'
#' print(calcall(inpt="m--2+-5"))
#'
#' [1] "m-3"
#'
#' print(calcall(inpt="(-2-6)+-6+2"))
#'
#' [1] "-12" 
#' 
#' print(calcall(inpt="m-6"))
#'
#' [1] "m-6"
#'
#' print(calcall(inpt="--6"))
#'
#' [1] "6"
#' 
#' @export

calcall <- function(inpt){
        
        can_be_num <- function(x){

                regex_spe_detect <- function(inpt){

                        fillr <- function(inpt_v, ptrn_fill="\\.\\.\\.\\d"){
                          
                          ptrn <- grep(ptrn_fill, inpt_v)

                          while (length(ptrn) > 0){
                           
                            ptrn <- grep(ptrn_fill, inpt_v)

                            idx <- ptrn[1] 
                            
                            untl <- as.numeric(c(unlist(strsplit(inpt_v[idx], split="\\.")))[4]) - 1
                           
                            pre_val <- inpt_v[(idx - 1)]

                            inpt_v[idx] <- pre_val

                            if (untl > 0){
                            
                              for (i in 1:untl){
                                
                                inpt_v <- append(inpt_v, pre_val, idx)
                                
                              }
                              
                            }

                          ptrn <- grep(ptrn_fill, inpt_v)
                            
                          }
                          
                          return(inpt_v)
                          
                        }

                   inpt <- unlist(strsplit(x=inpt, split=""))

                   may_be_v <- c("[", "]", "{", "}", "-", "_", ".", "(", ")", "/", "%", "*", "^", "?", "$")

                   pre_idx <- unique(match(x=inpt, table=may_be_v))

                   pre_idx <- pre_idx[!(is.na(pre_idx))]

                   for (el in may_be_v[pre_idx]){

                           for (i in grep(pattern=paste("\\", el, sep=""), x=inpt)){

                                   inpt <- append(x=inpt, values="\\", after=(i-1))

                           }

                   }

                
                   return(paste(inpt, collapse=""))

            }
    
            if (typeof(x) == "double"){

                    return(TRUE)

            }else{

                vec_bool <- c()

                v_ref <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "0", ".")    

                v_wrk <- unlist(str_split(x, ""))

                alrd <- TRUE

                for (i in 1:length(v_wrk)){ 

                        if (v_wrk[i] == "." & alrd){ 

                                vec_bool <- append(vec_bool, 1) 

                                alrd <- FALSE

                        }else{

                                vec_bool <- append(vec_bool, sum(grepl(pattern=regex_spe_detect(v_wrk[i]), x=v_ref))) 

                        }

                }

                if (sum(vec_bool) == length(vec_bool)){

                        return(TRUE)

                }else{

                        return(FALSE)

                }

            }

        }

        lst <- unlist(strsplit(x=inpt, split=""))

        lst_par <- c()

        lst_par_calc <- c()

        lst_pos <- c()

        paires = 1

        pre_paires = 1

        pre_paires2 = 1

        if ((length(grep(x=lst, pattern="\\(")) * 2) > 0){

                for (i in 1:(length(grep(x=lst, pattern="\\(")) * 2)){ 

                        lst_par <- c(lst_par, 0)

                        lst_par_calc <- c(lst_par_calc, 0)

                        lst_pos <- c(lst_pos, 0)


                }

        }

        vec_ret <- c()

        par_ = 1

        lvl_par = 0

        for (el in 1:length(lst)){

           if (lst[el] == "("){

                   if (!(is.null(vec_ret))){

                           lst_par_calc[pre_paires2:pre_paires][-vec_ret] <- lst_par_calc[pre_paires2:pre_paires][-vec_ret] + 1

                   }else{

                           lst_par_calc[pre_paires2:pre_paires] <- lst_par_calc[pre_paires2:pre_paires] + 1

                   }

                   pre_paires = pre_paires + 1

                   pre_cls <- TRUE

                   lst_pos[par_] <- el

                   par_ = par_ + 1

                   lvl_par = lvl_par + 1

           }

           if (lst[el] == ")"){

                   lvl_par = lvl_par - 1

                   if (!(is.null(vec_ret))){

                        lst_par_calc[c(pre_paires2:pre_paires)][-vec_ret] <- lst_par_calc[pre_paires2:pre_paires][-vec_ret] - 1

                        pre_val <- lst_par_calc[pre_paires2:pre_paires][vec_ret]

                        lst_par_calc[pre_paires2:pre_paires][vec_ret] <- (-2)
                   
                   }else{

                        lst_par_calc[c(pre_paires2:pre_paires)] <- lst_par_calc[pre_paires2:pre_paires] - 1

                   }

                   if (!(is.null(vec_ret))){ 

                           pre_mtch <- match(x=c(0, -1), table=lst_par_calc[pre_paires2:pre_paires][-vec_ret])

                           lst_par_calc[pre_paires2:pre_paires][vec_ret] <- pre_val 

                   }else{

                           pre_mtch <- match(x=c(0, -1), table=lst_par_calc[pre_paires2:pre_paires])

                   }

                   cnt_par = 1

                   cnt2 = 0

                   if (!(is.null(vec_ret))){

                           vec_ret <- sort(vec_ret)

                           if (pre_mtch[1] >= min(vec_ret)){

                                cnt2 = 2

                                while (pre_mtch[1] > cnt_par & cnt2 <= length(vec_ret)){

                                        if ((vec_ret[cnt2] - vec_ret[(cnt2 - 1)]) > 1){

                                                cnt_par = cnt_par + (vec_ret[cnt2] - vec_ret[(cnt2 - 1)]) - 1

                                        }

                                        cnt2 = cnt2 + 1

                                }

                                if (pre_mtch[1] > cnt_par){

                                        cnt_par = length(vec_ret) / 2 + 1

                                }

                                cnt2 = cnt2 - 1

                           }

                   }

                   lst_par[pre_mtch[1] + (pre_paires2 - 1) + ifelse(cnt2 %% 2 == 0, cnt2, (cnt2 - 1))] <- paires 

                   lst_par[pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret)] <- paires 

                   if ((pre_mtch[1] + (pre_paires2 - 1)) == 1){

                        pre_paires2 = pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret) + 1

                        vec_ret <- c()

                        cnt_par = 0

                   } else if (lst_par_calc[(pre_mtch[1] + (pre_paires2 - 1) - 1)] == -1 & ifelse(is.null(vec_ret), TRUE, 
                                is.na(match(x=-1, table=lst_par_calc[pre_paires2:pre_paires][-vec_ret])))){

                        pre_paires2 = pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret) + 1

                        vec_ret <- c()

                        cnt_par = 0

                   } else{

                        vec_ret <- c(vec_ret, (pre_mtch[1]) + ifelse(cnt2 %% 2 == 0, cnt2, (cnt2 - 1)), 
                                     (pre_mtch[2] + length(vec_ret)))

                   }

                   paires = paires + 1

                   pre_paires = pre_paires + 1

                   pre_cls <- FALSE

                   lst_pos[par_] <- el

                   par_ = par_ + 1

           }

        }

        el = 2

        vec_btwn <- c()

        while (length(lst_pos) > 0 & el <= length(lst_pos)){

                if (lst_par[el] == lst_par[(el - 1)]){

                        if (lst_pos[(el-1)] > 1){

                                pre_lst <- lst[1:(lst_pos[(el - 1)] - 1)]

                        }else{

                                pre_lst <- c()

                        }

                        cur_lst <- lst[lst_pos[(el - 1)]:lst_pos[el]]

                        if (lst_pos[el] < length(lst)){

                                post_lst <- lst[(lst_pos[el]+1):length(lst)] 

                        }else{

                                post_lst <- c()

                        }

                        if (cur_lst[2] == "-" & cur_lst[3] == "-"){

                                cur_lst <- cur_lst[-c(2, 3)]

                        }else if (cur_lst[2] == "+" & cur_lst[3] == "-"){

                                cur_lst <- cur_lst[-2]

                        }

                        i2 = 3

                        while (i2 < length(cur_lst)){

                                if (cur_lst[i2] == "-" & cur_lst[(i2 - 1)] == "-"){
                                        
                                        cur_lst <- cur_lst[-i2]

                                        cur_lst[(i2 - 1)] <- "+"

                                }

                                if (cur_lst[i2] == "-" & cur_lst[(i2 - 1)] == "+"){
                                        
                                        cur_lst <- cur_lst[-(i2 - 1)]

                                }

                                i2 = i2 + 1

                        }

                        pre_lngthf <- length(lst)

                        pre_mtch1 <- match(x="o", table=cur_lst)

                        pre_mtch2 <- match(x="e", table=cur_lst)

                        pre_mtch3 <- match(x="!", table=cur_lst)

                        pre_mtch4 <- match(x="^", table=cur_lst)

                        var_idx1 <- c()

                        var_idx2 <- c()

                        var_idx3 <- c()

                        var_idx4 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2)) | !(is.na(pre_mtch3)) | !(is.na(pre_mtch4))){

                                cur_pre_lngth <- length(cur_lst)

                                if (!(is.na(pre_mtch1)) & pre_mtch1 > max(c(ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        base_ <- c()

                                        cnt = 0

                                        while (can_be_num(cur_lst[pre_mtch1 + 3 + cnt])){

                                                base_ <- c(base_, cur_lst[pre_mtch1 + 3 + cnt])

                                                cnt = cnt + 1

                                        }

                                        cur_nb <- c()

                                        cnt = 4 + length(base_)

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch1 + cnt])

                                                cnt = cnt + 1

                                        }

                                        if (all(!(is.null(cur_nb)), !(is.null(base_)))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch1-2)], 
                                                             abs(log(x=as.numeric(paste(cur_nb, collapse="")), 
                                                base=as.numeric(paste(base_, collapse="")))), 
                                                cur_lst[(pre_mtch1 + cnt):length(cur_lst)])

                                                if (log(x=as.numeric(paste(cur_nb, collapse="")), 
                                                base=as.numeric(paste(base_, collapse=""))) < 0){

                                                        if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                cur_lst[pre_mtch1+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                cur_lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }

                                if (!(is.na(pre_mtch2)) & pre_mtch2 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        cur_nb <- c()

                                        cnt = 2

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt + 1

                                        }

                                        if (!(is.null(cur_nb))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2-1)], 
                                                             exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                             cur_lst[(pre_mtch2+cnt):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx2 <- c(var_idx2, pre_mtch2)

                                        }

                                }

                                if (!(is.na(pre_mtch3)) & pre_mtch3 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        cur_nb <- c()

                                        cnt = -1

                                        while (can_be_num(cur_lst[pre_mtch3 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch3 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(cur_nb))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch3+cnt)], 
                                                factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                cur_lst[(pre_mtch3 + 1):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx3 <- c(var_idx3, pre_mtch3)

                                        }

                                }

                                if (!(is.na(pre_mtch4)) & pre_mtch4 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3)))){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch4 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch4 + cnt])

                                                cnt = cnt - 1

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch4 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }
                                        while (can_be_num(cur_lst[pre_mtch4 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch4 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(is.null(pre_nb)) & !(is.null(post_nb)))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch4+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch4+cnt2):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx4 <- c(var_idx4, pre_mtch4)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="o", table=cur_lst)

                                }else{

                                        pre_mtch1 <- match(x="o", table=cur_lst[-var_idx1]) +length(var_idx1)

                                }

                                if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="e", table=cur_lst) 

                                }else{

                                        pre_mtch2 <- match(x="e", table=cur_lst[-var_idx2]) + length(var_idx2)

                                }

                                if (is.null(var_idx3)){

                                        pre_mtch3 <- match(x="!", table=cur_lst)

                                }else{

                                        pre_mtch3 <- match(x="!", table=cur_lst[-var_idx3]) + length(var_idx3)

                                }

                                if (is.null(var_idx4)){

                                        pre_mtch4 <- match(x="^", table=cur_lst)

                                }else{

                                        pre_mtch4 <- match(x="^", table=cur_lst[-var_idx4]) + length(var_idx4)

                                }

                        }
                        
                        pre_mtch1 <- match(x="*", table=cur_lst) 
                        
                        pre_mtch2 <- match(x="/", table=cur_lst) 

                        var_pres <- FALSE

                        if (!(is.null(var_idx1)) | !(is.null(var_idx2)) | !(is.null(var_idx3)) | !(is.null(var_idx4))){

                                var_pres <- TRUE

                        }

                        var_idx1 <- c()

                        var_idx2 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                                no_exc <- TRUE

                                no_exc2 <- TRUE

                                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch1 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch1 + cnt) > 2){

                                                        if (cur_lst[pre_mtch1 + cnt] == "-" & cur_lst[pre_mtch1 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch1 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch1 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch1 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch1+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch1+cnt2):length(cur_lst)])

                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                cur_lst[pre_mtch1+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                cur_lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }else if (pre_mtch1 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) > 0){

                                                        if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                cur_lst[pre_mtch1 + cnt] <- "+"

                                                        }

                                                }

                                        }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }else if (!(is.na(pre_mtch2))){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch2 + cnt) > 2){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-" & cur_lst[pre_mtch2 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch2 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch2 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch2 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch2+cnt2):length(cur_lst)])

                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch2+cnt] == "-"){

                                                                cur_lst[pre_mtch2+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch2+cnt] == "+"){

                                                                cur_lst[pre_mtch2+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch2+cnt))

                                                        }

                                                }else if (pre_mtch2 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) > 0){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                cur_lst[pre_mtch2 + cnt] <- "+"

                                                        }

                                                }

                                        }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="*", table=cur_lst)
                               
                                }else{

                                        pre_mtch1 <- match(x="*", table=cur_lst[-var_idx1]) + length(var_idx1)
                               
                                }

                               if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="/", table=cur_lst)
                               
                                }else{

                                        pre_mtch2 <- match(x="/", table=cur_lst[-var_idx2]) + length(var_idx2) 
                               
                                }

                        }

                        pre_mtch1 <- match(x="+", table=cur_lst) 
                        
                        pre_mtch2 <- match(x="-", table=cur_lst[3:length(cur_lst)]) + 2

                        if (length(var_idx1) > 0 | length(var_idx2) > 0){

                                var_pres <- TRUE

                        }

                        var_idx1 <- c()

                        var_idx2 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                                no_exc <- TRUE

                                no_exc2 <- TRUE

                                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                                pre_nb <- c()

                                                cnt = -1

                                                cnt2 = 1

                                                while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                        pre_nb <- c(pre_nb, cur_lst[pre_mtch1 + cnt])

                                                        cnt = cnt - 1

                                                }

                                                if (!(is.null(pre_nb))){

                                                        if ((pre_mtch1 + cnt) > 2){

                                                                if (cur_lst[pre_mtch1 + cnt] == "-" & cur_lst[pre_mtch1 + cnt - 1] == "+"){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        cnt = cnt - 1

                                                                }else if (cur_lst[pre_mtch1 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch1 + cnt - 1]))){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        no_exc <- FALSE

                                                                }else if (cur_lst[pre_mtch1 + cnt] == "_"){

                                                                        no_exc2 <- FALSE

                                                                }

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        cnt = cnt - 1

                                                        }

                                                }

                                                post_nb <- c()

                                                if (cur_lst[pre_mtch1 + cnt2] == "-"){

                                                        post_nb <- c(post_nb, "-")

                                                        cnt2 = cnt2 + 1

                                                }

                                                while (can_be_num(cur_lst[pre_mtch1 + cnt2])){

                                                        post_nb <- c(post_nb, cur_lst[pre_mtch1 + cnt2])

                                                        cnt2 = cnt2 + 1

                                                }

                                                if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                        cur_lst <- c(cur_lst[1:(pre_mtch1+cnt)], 
                                                        abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), 
                                                        cur_lst[(pre_mtch1 + cnt2):length(cur_lst)])

                                                        if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                                if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                        cur_lst[pre_mtch1+cnt] <- "+"

                                                                }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                        cur_lst[pre_mtch1+cnt] <- "-"

                                                                }else{

                                                                        cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                                }

                                                        }else if (pre_mtch1 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) > 0){

                                                                if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                        cur_lst[pre_mtch1 + cnt] <- "+"

                                                                }

                                                        }

                                                }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                                }

                                }else if (!(is.na(pre_mtch2))){
                                
                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch2 + cnt) > 2){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-" & cur_lst[pre_mtch2 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch2 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch2 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch2 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), cur_lst[(pre_mtch2 + cnt2):length(cur_lst)])
                                                
                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch2+cnt] == "-"){

                                                                cur_lst[pre_mtch2+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch2+cnt] == "+"){

                                                                cur_lst[pre_mtch2+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch2+cnt))

                                                        }

                                                }

                                        }else{

                                                var_idx2 <- c(var_idx2, pre_mtch2)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="+", table=cur_lst)
                               
                                }else{

                                        pre_mtch1 <- match(x="+", table=cur_lst[-var_idx1]) + length(var_idx1)
                               
                                }

                               if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="-", table=cur_lst[3:length(cur_lst)]) + 2
                               
                                }else{

                                        pre_mtch2 <- match(x="-", table=cur_lst[-c(1, 2, var_idx2)]) + length(var_idx2) + 2
                               
                                }

                        }
                        
                        if (length(var_idx1) > 0 | length(var_idx2) > 0){

                                var_pres <- TRUE

                        }

                        if (var_pres){

                                lst_par <- lst_par[-c(el, (el-1))] 

                                lst_pos <- lst_pos[-c(el, (el-1))] 

                                lst <- c(pre_lst, cur_lst, post_lst)

                                lst_pos[(el - 1):length(lst_pos)] <- lst_pos[(el - 1):length(lst_pos)] - (pre_lngthf - length(lst))

                                el = 2

                        }else{

                                nmrl <- TRUE

                                if (lst_pos[(el-1)] > 1){

                                        saved_pre_pos <- lst_pos[(el-1)] - 1

                                }else{

                                        nmrl <- FALSE

                                }

                                lst_par <- lst_par[-c(el, (el-1))] 

                                lst_pos <- lst_pos[-c(el, (el-1))] 

                                lst <- c(pre_lst, cur_lst[2:(length(cur_lst) - 1)], post_lst)                                

                                if (nmrl){

                                        if (lst[saved_pre_pos] == "-" & lst[(saved_pre_pos + 1)] == "-"){

                                                lst <- lst[-(saved_pre_pos + 1)]

                                                if (length(lst) > 3){

                                                        if (lst[saved_pre_pos - 1] %in% c("*", "/", "_")){

                                                                lst <- lst[-(saved_pre_pos)]

                                                        }else{

                                                                lst[saved_pre_pos] <- "+"

                                                        }

                                                }else{

                                                        lst <- lst[-1]

                                                }

                                        }

                                        if (lst[saved_pre_pos] == "+" & lst[(saved_pre_pos+1)] == "-"){

                                                lst <- lst[-saved_pre_pos]

                                        }

                                }

                                lst_pos[(el - 1):length(lst_pos)] <- lst_pos[(el - 1):length(lst_pos)] - (pre_lngthf - length(lst))

                                el = 2

                        }

                }else{

                        el = el + 1

                }

        }

        if (length(lst) > 2){

                if (lst[1] == "-" & lst[2] == "-"){

                        lst <- lst[-c(1, 2)]

                }else if (lst[1] == "+" & lst[2] == "-"){

                        lst <- lst[-1]

                }

                if (length(lst) > 3){

                        i = 2

                        while (i <= length(lst)){

                                if (lst[(i-1)] == "-" & lst[i] == "-"){

                                        lst <- lst[-i]

                                        lst[(i-1)] <- "+"

                                }else if (lst[(i-1)] == "+" & lst[i] == "-"){

                                        lst <- lst[-(i-1)]

                                }

                                i = i + 1

                        }

                }

        }

        pre_mtch1 <- match(x="o", table=lst)

        pre_mtch2 <- match(x="e", table=lst)

        pre_mtch3 <- match(x="!", table=lst)

        pre_mtch4 <- match(x="^", table=lst)

        var_idx1 <- c()

        var_idx2 <- c()

        var_idx3 <- c()

        var_idx4 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2)) | !(is.na(pre_mtch3)) | !(is.na(pre_mtch4))){

                cur_pre_lngth <- length(lst)

                if (!(is.na(pre_mtch1)) & pre_mtch1 > max(c(ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        base_ <- c()

                        cnt = 0

                        while (can_be_num(lst[pre_mtch1 + 3 + cnt])){

                                base_ <- c(base_, lst[pre_mtch1 + 3 + cnt])

                                cnt = cnt + 1

                        }

                        cur_nb <- c()

                        cnt = 4 + length(base_)

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch1 + cnt])

                                if ((pre_mtch1 + cnt) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt + 1

                                }

                        }

                        if (all(!(is.null(base_)), !(is.null(cur_nb)))){

                                if ((pre_mtch1-2) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch1-2)], 
                                                     log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse=""))), 
                                                     lst[(pre_mtch1+cnt):length(lst)])

                                }else if ((pre_mtch1-2) == 0 & !(no_btm)){

                                        lst <- log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch1-2)], 
                                                     log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse="")))) 


                                }else {

                                        lst <- c(log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse=""))), lst[(pre_mtch1+cnt):length(lst)])

                                }

                                if (log(x=as.numeric(paste(cur_nb, collapse="")), 
                                base=as.numeric(paste(base_, collapse=""))) < 0){

                                        if (lst[pre_mtch1+cnt] == "-"){

                                                lst[pre_mtch1+cnt] <- "+"

                                        }else if (lst[pre_mtch1+cnt] == "+"){

                                                lst[pre_mtch1+cnt] <- "-"

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                        }

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }

                if (!(is.na(pre_mtch2)) & pre_mtch2 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        cur_nb <- c()

                        cnt = 2

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt + 1

                                }

                        }

                        if (!(is.null(cur_nb))){

                                if ((pre_mtch2-1) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch2-1)], 
                                                     exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                     lst[(pre_mtch2+cnt):length(lst)])

                                }else if ((pre_mtch2-1) == 0 & !(no_btm)){

                                        lst <- exp(x=as.numeric(paste(cur_nb, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch2-1)], 
                                                     exp(x=as.numeric(paste(cur_nb, collapse="")))) 


                                }else {

                                        lst <- c(exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                     lst[(pre_mtch2+cnt):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx2 <- c(var_idx2, pre_mtch2)

                        }

                }

                if (!(is.na(pre_mtch3)) & pre_mtch3 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        cur_nb <- c()

                        cnt = -1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch3 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch3 + cnt])

                                if ((pre_mtch3 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(cur_nb))){

                                if ((pre_mtch3+cnt-1) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch3+cnt-1)], 
                                                     factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                     lst[(pre_mtch3+1):length(lst)])

                                }else if ((pre_mtch3+cnt-1) == 0 & !(no_btm)){

                                        lst <- factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch3+cnt-1)], 
                                                     factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse="")))) 


                                }else {

                                        lst <- c(factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                     lst[(pre_mtch3+1):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx3 <- c(var_idx3, pre_mtch3)

                        }

                }

                if (!(is.na(pre_mtch4)) & pre_mtch4 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3)))){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch4 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch4 + cnt])

                                if ((pre_mtch4 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch4 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }
                       
                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch4 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch4 + cnt2])

                                if ((pre_mtch4 + cnt2) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(is.null(pre_nb)), !(is.null(post_nb)))){

                                if ((pre_mtch4+cnt) > 1 & no_btm){

                                        lst <- c(lst[1:(pre_mtch4+cnt)], 
                                                     abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                     lst[(pre_mtch4+cnt2):length(lst)])

                                }else if ((pre_mtch4+cnt) == 1 & !(no_btm)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch4+cnt)], 
                                                     abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse="")))) 


                                }else {

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                     lst[(pre_mtch4+cnt2):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx4 <- c(var_idx4, pre_mtch4)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="o", table=lst)

                }else{

                        pre_mtch1 <- match(x="o", table=lst[-var_idx1]) + length(var_idx1)

                }

                if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="e", table=lst)

                }else{

                        pre_mtch2 <- match(x="e", table=lst[-var_idx2]) + length(var_idx2)


                }

                if (is.null(var_idx3)){

                        pre_mtch3 <- match(x="!", table=lst)

                }else{

                        pre_mtch3 <- match(x="!", table=lst[-var_idx3]) + length(var_idx3)

                }

                if (is.null(var_idx4)){

                        pre_mtch4 <- match(x="^", table=lst)

                }else{

                        pre_mtch4 <- match(x="^", table=lst[-var_idx4]) + length(var_idx4)

                }

        }

        pre_mtch1 <- match(x="*", table=lst) 
        
        pre_mtch2 <- match(x="/", table=lst) 

        var_idx1 <- c()

        var_idx2 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                no_exc <- TRUE

                no_exc2 <- TRUE
                
                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch1 + cnt])

                                if ((pre_mtch1 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch1 + cnt) > 1){

                                        if (lst[pre_mtch1 + cnt] == "-" & lst[pre_mtch1 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch1 + cnt] == "-" & !can_be_num(lst[pre_mtch1 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch1 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch1 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch1 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch1 + cnt2])

                                if (pre_mtch1 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch1+cnt) == 1 & (pre_mtch1+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch1+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                }else if ((pre_mtch1+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                }

                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch1+cnt) > 0){

                                                if (lst[pre_mtch1+cnt] == "-"){

                                                        lst[pre_mtch1+cnt] <- "+"

                                                }else if (lst[pre_mtch1+cnt] == "+"){

                                                        lst[pre_mtch1+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                        }

                                }else if (pre_mtch1 + cnt > 0 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) > 0){

                                        if (lst[pre_mtch1 + cnt] == "-"){

                                                lst[pre_mtch1 + cnt] <- "+"

                                        }

                                }

                        }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }else if (!(is.na(pre_mtch2))){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch2 + cnt) > 1){

                                        if (lst[pre_mtch2 + cnt] == "-" & lst[pre_mtch2 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch2 + cnt] == "-" & !can_be_num(lst[pre_mtch2 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch2 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch2 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch2 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch2 + cnt2])

                                if (pre_mtch2 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch2+cnt) == 1 & (pre_mtch2+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch2+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }else if ((pre_mtch2+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }

                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch2+cnt) > 0){

                                                if (lst[pre_mtch2+cnt] == "-"){

                                                        lst[pre_mtch2+cnt] <- "+"

                                                }else if (lst[pre_mtch2+cnt] == "+"){

                                                        lst[pre_mtch2+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                        }

                                }else if (pre_mtch2 + cnt > 0 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) > 0){

                                        if (lst[pre_mtch2 + cnt] == "-"){

                                                lst[pre_mtch2 + cnt] <- "+"

                                        }

                                }

                        }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="*", table=lst)
               
                }else{

                        pre_mtch1 <- match(x="*", table=lst[-var_idx1]) + length(var_idx1)
               
                }

               if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="/", table=lst)
               
                }else{

                        pre_mtch2 <- match(x="/", table=lst[-var_idx2]) + length(var_idx2) 
               
                }

        }

        pre_mtch1 <- match(x="+", table=lst) 
        
        pre_mtch2 <- match(x="-", table=lst[3:length(lst)]) + 2

        var_idx1 <- c()

        var_idx2 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                no_exc <- TRUE

                no_exc2 <- TRUE

                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                pre_nb <- c()

                                cnt = -1

                                cnt2 = 1

                                no_btm <- TRUE

                                while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                        pre_nb <- c(pre_nb, lst[pre_mtch1 + cnt])

                                        if ((pre_mtch1 + cnt) == 1){
        
                                                no_btm <- FALSE

                                        }else{

                                                cnt = cnt - 1

                                        }

                                }

                                if (!(is.null(pre_nb))){

                                        if ((pre_mtch1 + cnt) > 1){

                                                if (lst[pre_mtch1 + cnt] == "-" & lst[pre_mtch1 + cnt - 1] == "+"){

                                                        pre_nb <- c(pre_nb, "-")

                                                        cnt = cnt - 1

                                                }else if (lst[pre_mtch1 + cnt] == "-" & !can_be_num(lst[pre_mtch1 + cnt - 1])){

                                                        pre_nb <- c(pre_nb, "-")

                                                        no_exc <- FALSE

                                                }else if (lst[pre_mtch1 + cnt] == "_"){

                                                        no_exc2 <- FALSE

                                                }

                                        }else if (lst[pre_mtch1 + cnt] == "-"){

                                                        pre_nb <- c(pre_nb, "-")

                                                        cnt = cnt - 1

                                                        no_exc <- FALSE

                                        }

                                }

                                post_nb <- c()

                                if (lst[pre_mtch1 + cnt2] == "-"){

                                        post_nb <- c(post_nb, "-")

                                        cnt2 = cnt2 + 1

                                }

                                no_btm <- TRUE

                                while (can_be_num(lst[pre_mtch1 + cnt2]) & no_btm){

                                        post_nb <- c(post_nb, lst[pre_mtch1 + cnt2])

                                        if (pre_mtch1 + cnt2 == length(lst)){

                                                no_btm <- FALSE

                                        }else{

                                                cnt2 = cnt2 + 1

                                        }

                                }

                                if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                        if ((pre_mtch1+cnt) == 1 & (pre_mtch1+cnt2) == length(lst)){

                                                lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse="")))

                                        }else if ((pre_mtch1+cnt) == 1){

                                                lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                        }else if ((pre_mtch1+cnt2) == length(lst)){

                                                lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))))


                                        }else{

                                                lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                        }

                                        if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                if ((pre_mtch1+cnt) > 0){

                                                        if (lst[pre_mtch1+cnt] == "-"){

                                                                lst[pre_mtch1+cnt] <- "+"

                                                        }else if (lst[pre_mtch1+cnt] == "+"){

                                                                lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                }

                                        }else if (pre_mtch1 + cnt > 0 & as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse="")) > 0){

                                                if (lst[pre_mtch1 + cnt] == "-"){

                                                        lst[pre_mtch1 + cnt] <- "+"

                                                }

                                        }

                                }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                }

                }else if (!(is.na(pre_mtch2))){
                
                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch2 + cnt) > 1){

                                        if (lst[pre_mtch2 + cnt] == "-" & lst[pre_mtch2 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch2 + cnt] == "-" & !can_be_num(lst[pre_mtch2 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch2 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch2 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch2 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch2 + cnt2])

                                if (pre_mtch2 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch2+cnt) == 1 & (pre_mtch2+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch2+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }else if ((pre_mtch2+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }
                                
                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch2+cnt) > 0){

                                                if (lst[pre_mtch2+cnt] == "-"){

                                                        lst[pre_mtch2+cnt] <- "+"

                                                }else if (lst[pre_mtch2+cnt] == "+"){

                                                        lst[pre_mtch2+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                        }

                                }

                        }else{

                                var_idx2 <- c(var_idx2, pre_mtch2)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="+", table=lst)
               
                }else{

                        pre_mtch1 <- match(x="+", table=lst[-var_idx1]) + length(var_idx1)
               
                }

               if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="-", table=lst[3:length(lst)]) + 2
               
                }else{

                        pre_mtch2 <- match(x="-", table=lst[-c(1, 2, var_idx2)]) + length(var_idx2) + 2
               
                }

        }

        return(paste(lst, collapse=""))

}

#' calcall_var 
#' 
#' Does the same thing as calcall function but calculates the formula that have variables. The values of the variables have to be given in a list of vectors, see examples.
#'
#' @param inpt is the input formula, with the variables
#' @param var_name_v is the vector that contains the variables name in the order of apparition in the formula. If the variable appears multiple times in the formula, it has to be specified in this vector, see examples. 
#' @param var_val_l is the list containing the vectors containing the values of each variable, for each point you want to calculate. The vectors has to be given in the same order has the variable in var_name_v.
#' @examples
#'
#' print(calcall_var(inpt="(6+m*-(4-imp))+3/jp", var_name_v=c("m", "imp", "jp"), 
#'                   var_val_l=list(
#' 
#'                                  c(1:6), 
#' 
#'                                  c(3, 4, 2, 5, 6, 1),
#' 
#'                                  c(6:1))))
#'
#'  [1] "5.5"  "6.6"  "0.75" "11"   "17.5" "-9"
#'
#' print(calcall_var(inpt="(6+m*-(4-imp))+3/jp+jp", var_name_v=c("m", "imp", "jp", "jp"), 
#'                    var_val_l=list(
#'  
#'                                   c(1:6), 
#'  
#'                                   c(3, 4, 2, 5, 6, 1),
#'  
#'                                   c(6:1))))
#' 
#' [1] "11.5" "11.6" "4.75" "14" "19.5" "-8"
#'
#' @export

calcall_var <- function(inpt, var_name_v, var_val_l){

        calcall <- function(inpt){
        
        can_be_num <- function(x){

                regex_spe_detect <- function(inpt){

                        fillr <- function(inpt_v, ptrn_fill="\\.\\.\\.\\d"){
                          
                          ptrn <- grep(ptrn_fill, inpt_v)

                          while (length(ptrn) > 0){
                           
                            ptrn <- grep(ptrn_fill, inpt_v)

                            idx <- ptrn[1] 
                            
                            untl <- as.numeric(c(unlist(strsplit(inpt_v[idx], split="\\.")))[4]) - 1
                           
                            pre_val <- inpt_v[(idx - 1)]

                            inpt_v[idx] <- pre_val

                            if (untl > 0){
                            
                              for (i in 1:untl){
                                
                                inpt_v <- append(inpt_v, pre_val, idx)
                                
                              }
                              
                            }

                          ptrn <- grep(ptrn_fill, inpt_v)
                            
                          }
                          
                          return(inpt_v)
                          
                        }

                   inpt <- unlist(strsplit(x=inpt, split=""))

                   may_be_v <- c("[", "]", "{", "}", "-", "_", ".", "(", ")", "/", "%", "*", "^", "?", "$")

                   pre_idx <- unique(match(x=inpt, table=may_be_v))

                   pre_idx <- pre_idx[!(is.na(pre_idx))]

                   for (el in may_be_v[pre_idx]){

                           for (i in grep(pattern=paste("\\", el, sep=""), x=inpt)){

                                   inpt <- append(x=inpt, values="\\", after=(i-1))

                           }

                   }

                
                   return(paste(inpt, collapse=""))

            }
    
            if (typeof(x) == "double"){

                    return(TRUE)

            }else{

                vec_bool <- c()

                v_ref <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "0", ".")    

                v_wrk <- unlist(str_split(x, ""))

                alrd <- TRUE

                for (i in 1:length(v_wrk)){ 

                        if (v_wrk[i] == "." & alrd){ 

                                vec_bool <- append(vec_bool, 1) 

                                alrd <- FALSE

                        }else{

                                vec_bool <- append(vec_bool, sum(grepl(pattern=regex_spe_detect(v_wrk[i]), x=v_ref))) 

                        }

                }

                if (sum(vec_bool) == length(vec_bool)){

                        return(TRUE)

                }else{

                        return(FALSE)

                }

            }

        }

        lst <- unlist(strsplit(x=inpt, split=""))

        lst_par <- c()

        lst_par_calc <- c()

        lst_pos <- c()

        paires = 1

        pre_paires = 1

        pre_paires2 = 1

        if ((length(grep(x=lst, pattern="\\(")) * 2) > 0){

                for (i in 1:(length(grep(x=lst, pattern="\\(")) * 2)){ 

                        lst_par <- c(lst_par, 0)

                        lst_par_calc <- c(lst_par_calc, 0)

                        lst_pos <- c(lst_pos, 0)


                }

        }

        vec_ret <- c()

        par_ = 1

        lvl_par = 0

        for (el in 1:length(lst)){

           if (lst[el] == "("){

                   if (!(is.null(vec_ret))){

                           lst_par_calc[pre_paires2:pre_paires][-vec_ret] <- lst_par_calc[pre_paires2:pre_paires][-vec_ret] + 1

                   }else{

                           lst_par_calc[pre_paires2:pre_paires] <- lst_par_calc[pre_paires2:pre_paires] + 1

                   }

                   pre_paires = pre_paires + 1

                   pre_cls <- TRUE

                   lst_pos[par_] <- el

                   par_ = par_ + 1

                   lvl_par = lvl_par + 1

           }

           if (lst[el] == ")"){

                   lvl_par = lvl_par - 1

                   if (!(is.null(vec_ret))){

                        lst_par_calc[c(pre_paires2:pre_paires)][-vec_ret] <- lst_par_calc[pre_paires2:pre_paires][-vec_ret] - 1

                        pre_val <- lst_par_calc[pre_paires2:pre_paires][vec_ret]

                        lst_par_calc[pre_paires2:pre_paires][vec_ret] <- (-2)
                   
                   }else{

                        lst_par_calc[c(pre_paires2:pre_paires)] <- lst_par_calc[pre_paires2:pre_paires] - 1

                   }

                   if (!(is.null(vec_ret))){ 

                           pre_mtch <- match(x=c(0, -1), table=lst_par_calc[pre_paires2:pre_paires][-vec_ret])

                           lst_par_calc[pre_paires2:pre_paires][vec_ret] <- pre_val 

                   }else{

                           pre_mtch <- match(x=c(0, -1), table=lst_par_calc[pre_paires2:pre_paires])

                   }

                   cnt_par = 1

                   cnt2 = 0

                   if (!(is.null(vec_ret))){

                           vec_ret <- sort(vec_ret)

                           if (pre_mtch[1] >= min(vec_ret)){

                                cnt2 = 2

                                while (pre_mtch[1] > cnt_par & cnt2 <= length(vec_ret)){

                                        if ((vec_ret[cnt2] - vec_ret[(cnt2 - 1)]) > 1){

                                                cnt_par = cnt_par + (vec_ret[cnt2] - vec_ret[(cnt2 - 1)]) - 1

                                        }

                                        cnt2 = cnt2 + 1

                                }

                                if (pre_mtch[1] > cnt_par){

                                        cnt_par = length(vec_ret) / 2 + 1

                                }

                                cnt2 = cnt2 - 1

                           }

                   }

                   lst_par[pre_mtch[1] + (pre_paires2 - 1) + ifelse(cnt2 %% 2 == 0, cnt2, (cnt2 - 1))] <- paires 

                   lst_par[pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret)] <- paires 

                   if ((pre_mtch[1] + (pre_paires2 - 1)) == 1){

                        pre_paires2 = pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret) + 1

                        vec_ret <- c()

                        cnt_par = 0

                   } else if (lst_par_calc[(pre_mtch[1] + (pre_paires2 - 1) - 1)] == -1 & ifelse(is.null(vec_ret), TRUE, 
                                is.na(match(x=-1, table=lst_par_calc[pre_paires2:pre_paires][-vec_ret])))){

                        pre_paires2 = pre_mtch[2] + (pre_paires2 - 1) + length(vec_ret) + 1

                        vec_ret <- c()

                        cnt_par = 0

                   } else{

                        vec_ret <- c(vec_ret, (pre_mtch[1]) + ifelse(cnt2 %% 2 == 0, cnt2, (cnt2 - 1)), 
                                     (pre_mtch[2] + length(vec_ret)))

                   }

                   paires = paires + 1

                   pre_paires = pre_paires + 1

                   pre_cls <- FALSE

                   lst_pos[par_] <- el

                   par_ = par_ + 1

           }

        }

        el = 2

        vec_btwn <- c()

        while (length(lst_pos) > 0 & el <= length(lst_pos)){

                if (lst_par[el] == lst_par[(el - 1)]){

                        if (lst_pos[(el-1)] > 1){

                                pre_lst <- lst[1:(lst_pos[(el - 1)] - 1)]

                        }else{

                                pre_lst <- c()

                        }

                        cur_lst <- lst[lst_pos[(el - 1)]:lst_pos[el]]

                        if (lst_pos[el] < length(lst)){

                                post_lst <- lst[(lst_pos[el]+1):length(lst)] 

                        }else{

                                post_lst <- c()

                        }

                        if (cur_lst[2] == "-" & cur_lst[3] == "-"){

                                cur_lst <- cur_lst[-c(2, 3)]

                        }else if (cur_lst[2] == "+" & cur_lst[3] == "-"){

                                cur_lst <- cur_lst[-2]

                        }

                        i2 = 3

                        while (i2 < length(cur_lst)){

                                if (cur_lst[i2] == "-" & cur_lst[(i2 - 1)] == "-"){
                                        
                                        cur_lst <- cur_lst[-i2]

                                        cur_lst[(i2 - 1)] <- "+"

                                }

                                if (cur_lst[i2] == "-" & cur_lst[(i2 - 1)] == "+"){
                                        
                                        cur_lst <- cur_lst[-(i2 - 1)]

                                }

                                i2 = i2 + 1

                        }

                        pre_lngthf <- length(lst)

                        pre_mtch1 <- match(x="o", table=cur_lst)

                        pre_mtch2 <- match(x="e", table=cur_lst)

                        pre_mtch3 <- match(x="!", table=cur_lst)

                        pre_mtch4 <- match(x="^", table=cur_lst)

                        var_idx1 <- c()

                        var_idx2 <- c()

                        var_idx3 <- c()

                        var_idx4 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2)) | !(is.na(pre_mtch3)) | !(is.na(pre_mtch4))){

                                cur_pre_lngth <- length(cur_lst)

                                if (!(is.na(pre_mtch1)) & pre_mtch1 > max(c(ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        base_ <- c()

                                        cnt = 0

                                        while (can_be_num(cur_lst[pre_mtch1 + 3 + cnt])){

                                                base_ <- c(base_, cur_lst[pre_mtch1 + 3 + cnt])

                                                cnt = cnt + 1

                                        }

                                        cur_nb <- c()

                                        cnt = 4 + length(base_)

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch1 + cnt])

                                                cnt = cnt + 1

                                        }

                                        if (all(!(is.null(cur_nb)), !(is.null(base_)))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch1-2)], 
                                                             abs(log(x=as.numeric(paste(cur_nb, collapse="")), 
                                                base=as.numeric(paste(base_, collapse="")))), 
                                                cur_lst[(pre_mtch1 + cnt):length(cur_lst)])

                                                if (log(x=as.numeric(paste(cur_nb, collapse="")), 
                                                base=as.numeric(paste(base_, collapse=""))) < 0){

                                                        if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                cur_lst[pre_mtch1+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                cur_lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }

                                if (!(is.na(pre_mtch2)) & pre_mtch2 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        cur_nb <- c()

                                        cnt = 2

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt + 1

                                        }

                                        if (!(is.null(cur_nb))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2-1)], 
                                                             exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                             cur_lst[(pre_mtch2+cnt):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx2 <- c(var_idx2, pre_mtch2)

                                        }

                                }

                                if (!(is.na(pre_mtch3)) & pre_mtch3 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                                        cur_nb <- c()

                                        cnt = -1

                                        while (can_be_num(cur_lst[pre_mtch3 + cnt])){

                                                cur_nb <- c(cur_nb, cur_lst[pre_mtch3 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(cur_nb))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch3+cnt)], 
                                                factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                cur_lst[(pre_mtch3 + 1):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx3 <- c(var_idx3, pre_mtch3)

                                        }

                                }

                                if (!(is.na(pre_mtch4)) & pre_mtch4 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3)))){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch4 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch4 + cnt])

                                                cnt = cnt - 1

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch4 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }
                                        while (can_be_num(cur_lst[pre_mtch4 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch4 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(is.null(pre_nb)) & !(is.null(post_nb)))){

                                                cur_lst <- c(cur_lst[1:(pre_mtch4+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch4+cnt2):length(cur_lst)])

                                                if (!(is.null(var_idx1))){

                                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx2))){

                                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx3))){

                                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(cur_lst))
                                                
                                                }

                                                if (!(is.null(var_idx4))){

                                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(cur_lst))

                                                }

                                        }else{

                                                var_idx4 <- c(var_idx4, pre_mtch4)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="o", table=cur_lst)

                                }else{

                                        pre_mtch1 <- match(x="o", table=cur_lst[-var_idx1]) +length(var_idx1)

                                }

                                if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="e", table=cur_lst) 

                                }else{

                                        pre_mtch2 <- match(x="e", table=cur_lst[-var_idx2]) + length(var_idx2)

                                }

                                if (is.null(var_idx3)){

                                        pre_mtch3 <- match(x="!", table=cur_lst)

                                }else{

                                        pre_mtch3 <- match(x="!", table=cur_lst[-var_idx3]) + length(var_idx3)

                                }

                                if (is.null(var_idx4)){

                                        pre_mtch4 <- match(x="^", table=cur_lst)

                                }else{

                                        pre_mtch4 <- match(x="^", table=cur_lst[-var_idx4]) + length(var_idx4)

                                }

                        }
                        
                        pre_mtch1 <- match(x="*", table=cur_lst) 
                        
                        pre_mtch2 <- match(x="/", table=cur_lst) 

                        var_pres <- FALSE

                        if (!(is.null(var_idx1)) | !(is.null(var_idx2)) | !(is.null(var_idx3)) | !(is.null(var_idx4))){

                                var_pres <- TRUE

                        }

                        var_idx1 <- c()

                        var_idx2 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                                no_exc <- TRUE

                                no_exc2 <- TRUE

                                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch1 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch1 + cnt) > 2){

                                                        if (cur_lst[pre_mtch1 + cnt] == "-" & cur_lst[pre_mtch1 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch1 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch1 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch1 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch1 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch1+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch1+cnt2):length(cur_lst)])

                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                cur_lst[pre_mtch1+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                cur_lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }else if (pre_mtch1 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) > 0){

                                                        if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                cur_lst[pre_mtch1 + cnt] <- "+"

                                                        }

                                                }

                                        }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }else if (!(is.na(pre_mtch2))){

                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch2 + cnt) > 2){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-" & cur_lst[pre_mtch2 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch2 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch2 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch2 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), 
                                                cur_lst[(pre_mtch2+cnt2):length(cur_lst)])

                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch2+cnt] == "-"){

                                                                cur_lst[pre_mtch2+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch2+cnt] == "+"){

                                                                cur_lst[pre_mtch2+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch2+cnt))

                                                        }

                                                }else if (pre_mtch2 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) > 0){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                cur_lst[pre_mtch2 + cnt] <- "+"

                                                        }

                                                }

                                        }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="*", table=cur_lst)
                               
                                }else{

                                        pre_mtch1 <- match(x="*", table=cur_lst[-var_idx1]) + length(var_idx1)
                               
                                }

                               if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="/", table=cur_lst)
                               
                                }else{

                                        pre_mtch2 <- match(x="/", table=cur_lst[-var_idx2]) + length(var_idx2) 
                               
                                }

                        }

                        pre_mtch1 <- match(x="+", table=cur_lst) 
                        
                        pre_mtch2 <- match(x="-", table=cur_lst[3:length(cur_lst)]) + 2

                        if (length(var_idx1) > 0 | length(var_idx2) > 0){

                                var_pres <- TRUE

                        }

                        var_idx1 <- c()

                        var_idx2 <- c()

                        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                                no_exc <- TRUE

                                no_exc2 <- TRUE

                                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                                pre_nb <- c()

                                                cnt = -1

                                                cnt2 = 1

                                                while (can_be_num(cur_lst[pre_mtch1 + cnt])){

                                                        pre_nb <- c(pre_nb, cur_lst[pre_mtch1 + cnt])

                                                        cnt = cnt - 1

                                                }

                                                if (!(is.null(pre_nb))){

                                                        if ((pre_mtch1 + cnt) > 2){

                                                                if (cur_lst[pre_mtch1 + cnt] == "-" & cur_lst[pre_mtch1 + cnt - 1] == "+"){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        cnt = cnt - 1

                                                                }else if (cur_lst[pre_mtch1 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch1 + cnt - 1]))){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        no_exc <- FALSE

                                                                }else if (cur_lst[pre_mtch1 + cnt] == "_"){

                                                                        no_exc2 <- FALSE

                                                                }

                                                        }else if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                        pre_nb <- c(pre_nb, "-")

                                                                        cnt = cnt - 1

                                                        }

                                                }

                                                post_nb <- c()

                                                if (cur_lst[pre_mtch1 + cnt2] == "-"){

                                                        post_nb <- c(post_nb, "-")

                                                        cnt2 = cnt2 + 1

                                                }

                                                while (can_be_num(cur_lst[pre_mtch1 + cnt2])){

                                                        post_nb <- c(post_nb, cur_lst[pre_mtch1 + cnt2])

                                                        cnt2 = cnt2 + 1

                                                }

                                                if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                        cur_lst <- c(cur_lst[1:(pre_mtch1+cnt)], 
                                                        abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), 
                                                        cur_lst[(pre_mtch1 + cnt2):length(cur_lst)])

                                                        if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                                if (cur_lst[pre_mtch1+cnt] == "-"){

                                                                        cur_lst[pre_mtch1+cnt] <- "+"

                                                                }else if (cur_lst[pre_mtch1+cnt] == "+"){

                                                                        cur_lst[pre_mtch1+cnt] <- "-"

                                                                }else{

                                                                        cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch1+cnt))

                                                                }

                                                        }else if (pre_mtch1 + cnt > 1 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) > 0){

                                                                if (cur_lst[pre_mtch1 + cnt] == "-"){

                                                                        cur_lst[pre_mtch1 + cnt] <- "+"

                                                                }

                                                        }

                                                }else{

                                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                                }

                                }else if (!(is.na(pre_mtch2))){
                                
                                        pre_nb <- c()

                                        cnt = -1

                                        cnt2 = 1

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt])){

                                                pre_nb <- c(pre_nb, cur_lst[pre_mtch2 + cnt])

                                                cnt = cnt - 1

                                        }

                                        if (!(is.null(pre_nb))){

                                                if ((pre_mtch2 + cnt) > 2){

                                                        if (cur_lst[pre_mtch2 + cnt] == "-" & cur_lst[pre_mtch2 + cnt - 1] == "+"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "-" & !(can_be_num(cur_lst[pre_mtch2 + cnt - 1]))){

                                                                pre_nb <- c(pre_nb, "-")

                                                                no_exc <- FALSE

                                                        }else if (cur_lst[pre_mtch2 + cnt] == "_"){

                                                                no_exc2 <- FALSE

                                                        }

                                                }else if (cur_lst[pre_mtch2 + cnt] == "-"){

                                                                pre_nb <- c(pre_nb, "-")

                                                                cnt = cnt - 1

                                                }

                                        }

                                        post_nb <- c()

                                        if (cur_lst[pre_mtch2 + cnt2] == "-"){

                                                post_nb <- c(post_nb, "-")

                                                cnt2 = cnt2 + 1

                                        }

                                        while (can_be_num(cur_lst[pre_mtch2 + cnt2])){

                                                post_nb <- c(post_nb, cur_lst[pre_mtch2 + cnt2])

                                                cnt2 = cnt2 + 1

                                        }

                                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                                cur_lst <- c(cur_lst[1:(pre_mtch2+cnt)], 
                                                abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), cur_lst[(pre_mtch2 + cnt2):length(cur_lst)])
                                                
                                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                        if (cur_lst[pre_mtch2+cnt] == "-"){

                                                                cur_lst[pre_mtch2+cnt] <- "+"

                                                        }else if (cur_lst[pre_mtch2+cnt] == "+"){

                                                                cur_lst[pre_mtch2+cnt] <- "-"

                                                        }else{

                                                                cur_lst <- append(x=cur_lst, values="-", after=(pre_mtch2+cnt))

                                                        }

                                                }

                                        }else{

                                                var_idx2 <- c(var_idx2, pre_mtch2)

                                        }

                                }

                                if (is.null(var_idx1)){

                                        pre_mtch1 <- match(x="+", table=cur_lst)
                               
                                }else{

                                        pre_mtch1 <- match(x="+", table=cur_lst[-var_idx1]) + length(var_idx1)
                               
                                }

                               if (is.null(var_idx2)){

                                        pre_mtch2 <- match(x="-", table=cur_lst[3:length(cur_lst)]) + 2
                               
                                }else{

                                        pre_mtch2 <- match(x="-", table=cur_lst[-c(1, 2, var_idx2)]) + length(var_idx2) + 2
                               
                                }

                        }
                        
                        if (length(var_idx1) > 0 | length(var_idx2) > 0){

                                var_pres <- TRUE

                        }

                        if (var_pres){

                                lst_par <- lst_par[-c(el, (el-1))] 

                                lst_pos <- lst_pos[-c(el, (el-1))] 

                                lst <- c(pre_lst, cur_lst, post_lst)

                                lst_pos[(el - 1):length(lst_pos)] <- lst_pos[(el - 1):length(lst_pos)] - (pre_lngthf - length(lst))

                                el = 2

                        }else{

                                nmrl <- TRUE

                                if (lst_pos[(el-1)] > 1){

                                        saved_pre_pos <- lst_pos[(el-1)] - 1

                                }else{

                                        nmrl <- FALSE

                                }

                                lst_par <- lst_par[-c(el, (el-1))] 

                                lst_pos <- lst_pos[-c(el, (el-1))] 

                                lst <- c(pre_lst, cur_lst[2:(length(cur_lst) - 1)], post_lst)                                

                                if (nmrl){

                                        if (lst[saved_pre_pos] == "-" & lst[(saved_pre_pos + 1)] == "-"){

                                                lst <- lst[-(saved_pre_pos + 1)]

                                                if (length(lst) > 3){

                                                        if (lst[saved_pre_pos - 1] %in% c("*", "/", "_")){

                                                                lst <- lst[-(saved_pre_pos)]

                                                        }else{

                                                                lst[saved_pre_pos] <- "+"

                                                        }

                                                }else{

                                                        lst <- lst[-1]

                                                }

                                        }

                                        if (lst[saved_pre_pos] == "+" & lst[(saved_pre_pos+1)] == "-"){

                                                lst <- lst[-saved_pre_pos]

                                        }

                                }

                                lst_pos[(el - 1):length(lst_pos)] <- lst_pos[(el - 1):length(lst_pos)] - (pre_lngthf - length(lst))

                                el = 2

                        }

                }else{

                        el = el + 1

                }

        }

        if (length(lst) > 2){

                if (lst[1] == "-" & lst[2] == "-"){

                        lst <- lst[-c(1, 2)]

                }else if (lst[1] == "+" & lst[2] == "-"){

                        lst <- lst[-1]

                }

                if (length(lst) > 3){

                        i = 2

                        while (i <= length(lst)){

                                if (lst[(i-1)] == "-" & lst[i] == "-"){

                                        lst <- lst[-i]

                                        lst[(i-1)] <- "+"

                                }else if (lst[(i-1)] == "+" & lst[i] == "-"){

                                        lst <- lst[-(i-1)]

                                }

                                i = i + 1

                        }

                }

        }

        pre_mtch1 <- match(x="o", table=lst)

        pre_mtch2 <- match(x="e", table=lst)

        pre_mtch3 <- match(x="!", table=lst)

        pre_mtch4 <- match(x="^", table=lst)

        var_idx1 <- c()

        var_idx2 <- c()

        var_idx3 <- c()

        var_idx4 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2)) | !(is.na(pre_mtch3)) | !(is.na(pre_mtch4))){

                cur_pre_lngth <- length(lst)

                if (!(is.na(pre_mtch1)) & pre_mtch1 > max(c(ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        base_ <- c()

                        cnt = 0

                        while (can_be_num(lst[pre_mtch1 + 3 + cnt])){

                                base_ <- c(base_, lst[pre_mtch1 + 3 + cnt])

                                cnt = cnt + 1

                        }

                        cur_nb <- c()

                        cnt = 4 + length(base_)

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch1 + cnt])

                                if ((pre_mtch1 + cnt) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt + 1

                                }

                        }

                        if (all(!(is.null(base_)), !(is.null(cur_nb)))){

                                if ((pre_mtch1-2) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch1-2)], 
                                                     log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse=""))), 
                                                     lst[(pre_mtch1+cnt):length(lst)])

                                }else if ((pre_mtch1-2) == 0 & !(no_btm)){

                                        lst <- log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch1-2)], 
                                                     log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse="")))) 


                                }else {

                                        lst <- c(log(x=as.numeric(paste(cur_nb, collapse="")), base=as.numeric(paste(base_, collapse=""))), lst[(pre_mtch1+cnt):length(lst)])

                                }

                                if (log(x=as.numeric(paste(cur_nb, collapse="")), 
                                base=as.numeric(paste(base_, collapse=""))) < 0){

                                        if (lst[pre_mtch1+cnt] == "-"){

                                                lst[pre_mtch1+cnt] <- "+"

                                        }else if (lst[pre_mtch1+cnt] == "+"){

                                                lst[pre_mtch1+cnt] <- "-"

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                        }

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }

                if (!(is.na(pre_mtch2)) & pre_mtch2 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        cur_nb <- c()

                        cnt = 2

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt + 1

                                }

                        }

                        if (!(is.null(cur_nb))){

                                if ((pre_mtch2-1) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch2-1)], 
                                                     exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                     lst[(pre_mtch2+cnt):length(lst)])

                                }else if ((pre_mtch2-1) == 0 & !(no_btm)){

                                        lst <- exp(x=as.numeric(paste(cur_nb, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch2-1)], 
                                                     exp(x=as.numeric(paste(cur_nb, collapse="")))) 


                                }else {

                                        lst <- c(exp(x=as.numeric(paste(cur_nb, collapse=""))), 
                                                     lst[(pre_mtch2+cnt):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx2 <- c(var_idx2, pre_mtch2)

                        }

                }

                if (!(is.na(pre_mtch3)) & pre_mtch3 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch4), 0, pre_mtch4)))){

                        cur_nb <- c()

                        cnt = -1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch3 + cnt]) & no_btm){

                                cur_nb <- c(cur_nb, lst[pre_mtch3 + cnt])

                                if ((pre_mtch3 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(cur_nb))){

                                if ((pre_mtch3+cnt-1) > 0 & no_btm){

                                        lst <- c(lst[1:(pre_mtch3+cnt-1)], 
                                                     factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                     lst[(pre_mtch3+1):length(lst)])

                                }else if ((pre_mtch3+cnt-1) == 0 & !(no_btm)){

                                        lst <- factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch3+cnt-1)], 
                                                     factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse="")))) 


                                }else {

                                        lst <- c(factorial(x=as.numeric(paste(cur_nb[length(cur_nb):1], collapse=""))), 
                                                     lst[(pre_mtch3+1):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx3 <- c(var_idx3, pre_mtch3)

                        }

                }

                if (!(is.na(pre_mtch4)) & pre_mtch4 > max(c(ifelse(is.na(pre_mtch1), 0, pre_mtch1), 
                                                            ifelse(is.na(pre_mtch2), 0, pre_mtch2), 
                                                            ifelse(is.na(pre_mtch3), 0, pre_mtch3)))){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch4 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch4 + cnt])

                                if ((pre_mtch4 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch4 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }
                       
                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch4 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch4 + cnt2])

                                if ((pre_mtch4 + cnt2) == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(is.null(pre_nb)), !(is.null(post_nb)))){

                                if ((pre_mtch4+cnt) > 1 & no_btm){

                                        lst <- c(lst[1:(pre_mtch4+cnt)], 
                                                     abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                     lst[(pre_mtch4+cnt2):length(lst)])

                                }else if ((pre_mtch4+cnt) == 1 & !(no_btm)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse="")))

                                }else if (!(no_btm)){

                                        lst <- c(lst[1:(pre_mtch4+cnt)], 
                                                     abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse="")))) 


                                }else {

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) ** as.numeric(paste(post_nb, collapse=""))), 
                                                     lst[(pre_mtch4+cnt2):length(lst)])

                                }

                                if (!(is.null(var_idx1))){

                                        var_idx1 <- var_idx1 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx2))){

                                        var_idx2 <- var_idx2 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx3))){

                                        var_idx3 <- var_idx3 - (cur_pre_lngth - length(lst))
                                
                                }

                                if (!(is.null(var_idx4))){

                                        var_idx4 <- var_idx4 - (cur_pre_lngth - length(lst))

                                }

                        }else{

                                var_idx4 <- c(var_idx4, pre_mtch4)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="o", table=lst)

                }else{

                        pre_mtch1 <- match(x="o", table=lst[-var_idx1]) + length(var_idx1)

                }

                if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="e", table=lst)

                }else{

                        pre_mtch2 <- match(x="e", table=lst[-var_idx2]) + length(var_idx2)


                }

                if (is.null(var_idx3)){

                        pre_mtch3 <- match(x="!", table=lst)

                }else{

                        pre_mtch3 <- match(x="!", table=lst[-var_idx3]) + length(var_idx3)

                }

                if (is.null(var_idx4)){

                        pre_mtch4 <- match(x="^", table=lst)

                }else{

                        pre_mtch4 <- match(x="^", table=lst[-var_idx4]) + length(var_idx4)

                }

        }

        pre_mtch1 <- match(x="*", table=lst) 
        
        pre_mtch2 <- match(x="/", table=lst) 

        var_idx1 <- c()

        var_idx2 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                no_exc <- TRUE

                no_exc2 <- TRUE
                
                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch1 + cnt])

                                if ((pre_mtch1 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch1 + cnt) > 1){

                                        if (lst[pre_mtch1 + cnt] == "-" & lst[pre_mtch1 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch1 + cnt] == "-" & !can_be_num(lst[pre_mtch1 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch1 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch1 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch1 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch1 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch1 + cnt2])

                                if (pre_mtch1 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch1+cnt) == 1 & (pre_mtch1+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch1+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                }else if ((pre_mtch1+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                }

                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch1+cnt) > 0){

                                                if (lst[pre_mtch1+cnt] == "-"){

                                                        lst[pre_mtch1+cnt] <- "+"

                                                }else if (lst[pre_mtch1+cnt] == "+"){

                                                        lst[pre_mtch1+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                        }

                                }else if (pre_mtch1 + cnt > 0 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) * as.numeric(paste(post_nb, collapse=""))) > 0){

                                        if (lst[pre_mtch1 + cnt] == "-"){

                                                lst[pre_mtch1 + cnt] <- "+"

                                        }

                                }

                        }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }else if (!(is.na(pre_mtch2))){

                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch2 + cnt) > 1){

                                        if (lst[pre_mtch2 + cnt] == "-" & lst[pre_mtch2 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch2 + cnt] == "-" & !can_be_num(lst[pre_mtch2 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch2 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch2 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch2 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch2 + cnt2])

                                if (pre_mtch2 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch2+cnt) == 1 & (pre_mtch2+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch2+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }else if ((pre_mtch2+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }

                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch2+cnt) > 0){

                                                if (lst[pre_mtch2+cnt] == "-"){

                                                        lst[pre_mtch2+cnt] <- "+"

                                                }else if (lst[pre_mtch2+cnt] == "+"){

                                                        lst[pre_mtch2+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                        }

                                }else if (pre_mtch2 + cnt > 0 & (as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) / as.numeric(paste(post_nb, collapse=""))) > 0){

                                        if (lst[pre_mtch2 + cnt] == "-"){

                                                lst[pre_mtch2 + cnt] <- "+"

                                        }

                                }

                        }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="*", table=lst)
               
                }else{

                        pre_mtch1 <- match(x="*", table=lst[-var_idx1]) + length(var_idx1)
               
                }

               if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="/", table=lst)
               
                }else{

                        pre_mtch2 <- match(x="/", table=lst[-var_idx2]) + length(var_idx2) 
               
                }

        }

        pre_mtch1 <- match(x="+", table=lst) 
        
        pre_mtch2 <- match(x="-", table=lst[3:length(lst)]) + 2

        var_idx1 <- c()

        var_idx2 <- c()

        while (!(is.na(pre_mtch1)) | !(is.na(pre_mtch2))){

                no_exc <- TRUE

                no_exc2 <- TRUE

                if (!(is.na(pre_mtch1)) & pre_mtch1 < ifelse(is.na(pre_mtch2), length(lst), pre_mtch2)){

                                pre_nb <- c()

                                cnt = -1

                                cnt2 = 1

                                no_btm <- TRUE

                                while (can_be_num(lst[pre_mtch1 + cnt]) & no_btm){

                                        pre_nb <- c(pre_nb, lst[pre_mtch1 + cnt])

                                        if ((pre_mtch1 + cnt) == 1){
        
                                                no_btm <- FALSE

                                        }else{

                                                cnt = cnt - 1

                                        }

                                }

                                if (!(is.null(pre_nb))){

                                        if ((pre_mtch1 + cnt) > 1){

                                                if (lst[pre_mtch1 + cnt] == "-" & lst[pre_mtch1 + cnt - 1] == "+"){

                                                        pre_nb <- c(pre_nb, "-")

                                                        cnt = cnt - 1

                                                }else if (lst[pre_mtch1 + cnt] == "-" & !can_be_num(lst[pre_mtch1 + cnt - 1])){

                                                        pre_nb <- c(pre_nb, "-")

                                                        no_exc <- FALSE

                                                }else if (lst[pre_mtch1 + cnt] == "_"){

                                                        no_exc2 <- FALSE

                                                }

                                        }else if (lst[pre_mtch1 + cnt] == "-"){

                                                        pre_nb <- c(pre_nb, "-")

                                                        cnt = cnt - 1

                                                        no_exc <- FALSE

                                        }

                                }

                                post_nb <- c()

                                if (lst[pre_mtch1 + cnt2] == "-"){

                                        post_nb <- c(post_nb, "-")

                                        cnt2 = cnt2 + 1

                                }

                                no_btm <- TRUE

                                while (can_be_num(lst[pre_mtch1 + cnt2]) & no_btm){

                                        post_nb <- c(post_nb, lst[pre_mtch1 + cnt2])

                                        if (pre_mtch1 + cnt2 == length(lst)){

                                                no_btm <- FALSE

                                        }else{

                                                cnt2 = cnt2 + 1

                                        }

                                }

                                if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                        if ((pre_mtch1+cnt) == 1 & (pre_mtch1+cnt2) == length(lst)){

                                                lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse="")))

                                        }else if ((pre_mtch1+cnt) == 1){

                                                lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                        }else if ((pre_mtch1+cnt2) == length(lst)){

                                                lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))))


                                        }else{

                                                lst <- c(lst[1:(pre_mtch1+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch1 + cnt2):length(lst)])

                                        }

                                        if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                                if ((pre_mtch1+cnt) > 0){

                                                        if (lst[pre_mtch1+cnt] == "-"){

                                                                lst[pre_mtch1+cnt] <- "+"

                                                        }else if (lst[pre_mtch1+cnt] == "+"){

                                                                lst[pre_mtch1+cnt] <- "-"

                                                        }else{

                                                                lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                        }

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch1+cnt))

                                                }

                                        }else if (pre_mtch1 + cnt > 0 & as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) + as.numeric(paste(post_nb, collapse="")) > 0){

                                                if (lst[pre_mtch1 + cnt] == "-"){

                                                        lst[pre_mtch1 + cnt] <- "+"

                                                }

                                        }

                                }else{

                                        var_idx1 <- c(var_idx1, pre_mtch1)

                                }

                }else if (!(is.na(pre_mtch2))){
                
                        pre_nb <- c()

                        cnt = -1

                        cnt2 = 1

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt]) & no_btm){

                                pre_nb <- c(pre_nb, lst[pre_mtch2 + cnt])

                                if ((pre_mtch2 + cnt) == 1){

                                        no_btm <- FALSE

                                }else{

                                        cnt = cnt - 1

                                }

                        }

                        if (!(is.null(pre_nb))){

                                if ((pre_mtch2 + cnt) > 1){

                                        if (lst[pre_mtch2 + cnt] == "-" & lst[pre_mtch2 + cnt - 1] == "+"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                        }else if (lst[pre_mtch2 + cnt] == "-" & !can_be_num(lst[pre_mtch2 + cnt - 1])){

                                                pre_nb <- c(pre_nb, "-")

                                                no_exc <- FALSE

                                        }else if (lst[pre_mtch2 + cnt] == "_"){

                                                no_exc2 <- FALSE

                                        }

                                }else if (lst[pre_mtch2 + cnt] == "-"){

                                                pre_nb <- c(pre_nb, "-")

                                                cnt = cnt - 1

                                                no_exc <- FALSE

                                }

                        }

                        post_nb <- c()

                        if (lst[pre_mtch2 + cnt2] == "-"){

                                post_nb <- c(post_nb, "-")

                                cnt2 = cnt2 + 1

                        }

                        no_btm <- TRUE

                        while (can_be_num(lst[pre_mtch2 + cnt2]) & no_btm){

                                post_nb <- c(post_nb, lst[pre_mtch2 + cnt2])

                                if (pre_mtch2 + cnt2 == length(lst)){

                                        no_btm <- FALSE

                                }else{

                                        cnt2 = cnt2 + 1

                                }

                        }

                        if (all(!(c(is.null(pre_nb), is.null(post_nb)))) & no_exc2){

                                if ((pre_mtch2+cnt) == 1 & (pre_mtch2+cnt2) == length(lst)){

                                        lst <- abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse="")))

                                }else if ((pre_mtch2+cnt) == 1){

                                        lst <- c(abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }else if ((pre_mtch2+cnt2) == length(lst)){

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))))


                                }else{

                                        lst <- c(lst[1:(pre_mtch2+cnt)], abs(as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))), lst[(pre_mtch2 + cnt2):length(lst)])

                                }
                                
                                if ((as.numeric(paste(pre_nb[length(pre_nb):1], collapse="")) - as.numeric(paste(post_nb, collapse=""))) < 0 & no_exc){

                                        if ((pre_mtch2+cnt) > 0){

                                                if (lst[pre_mtch2+cnt] == "-"){

                                                        lst[pre_mtch2+cnt] <- "+"

                                                }else if (lst[pre_mtch2+cnt] == "+"){

                                                        lst[pre_mtch2+cnt] <- "-"

                                                }else{

                                                        lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                                }

                                        }else{

                                                lst <- append(x=lst, values="-", after=(pre_mtch2+cnt))

                                        }

                                }

                        }else{

                                var_idx2 <- c(var_idx2, pre_mtch2)

                        }

                }

                if (is.null(var_idx1)){

                        pre_mtch1 <- match(x="+", table=lst)
               
                }else{

                        pre_mtch1 <- match(x="+", table=lst[-var_idx1]) + length(var_idx1)
               
                }

               if (is.null(var_idx2)){

                        pre_mtch2 <- match(x="-", table=lst[3:length(lst)]) + 2
               
                }else{

                        pre_mtch2 <- match(x="-", table=lst[-c(1, 2, var_idx2)]) + length(var_idx2) + 2
               
                }

        }

        return(paste(lst, collapse=""))

        }
       
        rtn_v <- c()

        for (I in 1:length(unlist(var_val_l[1]))){

                cur_inpt <- inpt

                for (var in var_name_v){

                        cur_inpt <- sub(pattern=var, 

                                         replacement=unlist(var_val_l[match(x=var, table=var_name_v)])[I],

                                         x=cur_inpt)

                }

                rtn_v <- c(rtn_v, calcall(cur_inpt))

        }

        return(rtn_v)

}

#' common_tracks
#'
#' From a time serie, allow to get the most common route for each individual at a given depth (time - 1). Access the frequency value as an element from the output vector and the value itself (the path) as a name of its element, see examples.
#'
#' @param inpt_datf is the input time serie as a dataframe
#' @param col_target is the column name or number that refers to the value of each individual
#' @param id_col is the column name or number that refers to the individual (ids)
#' @param untl_last is the depth value
#' @examples
#'
#' datf_test <- data.frame("id" = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 5),
#'                          "city" = c("A", "C", "B", "B", "A", "C", "A", "C", "A", "C", "B", "A", "A", "E"))
#' 
#' print(individual_route(inpt_datf = datf_test, 
#'                        col_target = "city", 
#'                        id_col = "id",
#'                        untl_last = 2))
#' 
#' AC CA BA 
#'  2  1  2 
#' 
#' print(individual_route(inpt_datf = datf_test, 
#'                        col_target = "city", 
#'                        id_col = "id",
#'                        untl_last = 3))
#'
#' ACB  AC CAC  BA BAA 
#'   1   2   1   2   1 
#'
#' @export

individual_route <- function(inpt_datf, col_target, id_col, untl_last = 2){
  if (typeof(col_target) == "character"){
    col_target <- match(x = col_target, table = colnames(inpt_datf))
  }
  if (typeof(id_col) == "character"){
    id_col <- match(x = id_col, table = colnames(inpt_datf))
  }
  inpt_datf <- inpt_datf[, c(col_target, id_col)]
  mods_v <- unique(inpt_datf[, 1])
  rtn_v <- c()
  for (mod in mods_v){
    cur_ids <- c()
    for (i in unique(inpt_datf[, 2])){
      cur_vec <- inpt_datf[grep(pattern = i, x = inpt_datf[, 2]), 1]
      pre_mtch <- match(x = mod, table = cur_vec)
      if (!(is.na(pre_mtch))){
        if (pre_mtch == 1){
          cur_ids <- c(cur_ids, grep(pattern = i, x = inpt_datf[, 2]))
        }
      }
    }
    if (length(cur_ids) > 0){
      cur_datf <- inpt_datf[cur_ids,]
      un_individual <- unique(cur_datf[, 2])
      concat_v <- c(matrix(data = "", nrow = length(un_individual), ncol = 1))
      for (i in 1:length(un_individual)){
        cur_vec <- cur_datf[grep(pattern = un_individual[i], x = cur_datf[, 2]), 1]
        cnt = 1
        while (cnt <= untl_last & cnt <= length(cur_vec)){
          concat_v[i] <- paste0(concat_v[i], cur_vec[cnt])
          cnt = cnt + 1
        }
      }
      pre_lngth <- length(rtn_v) + 1
      for (i in unique(concat_v)){
        rtn_v <- c(rtn_v, sum(grepl(pattern = i, x = concat_v)))
      }
      names(rtn_v)[pre_lngth:(pre_lngth - 1 + length(unique(concat_v)))] <- unique(concat_v) 
      inpt_datf <- inpt_datf[-cur_ids,]
    }
  }
  return(rtn_v)
}

#' v_Rmach_fold
#'
#' Allow to create uniform sampling dataset for cross validation, 
#' train and test, see examples and variables
#'
#' @param inpt_datf is the input dataframe
#' @param train_prop is the training proportion 
#' @param n_fold is the number of distinc pair of training and test dataset that will be outputed
#'
#' @examples
#'
#' print(v_Rmach_fold(inpt_datf = iris[1:25,],
#'              train_prop = 0.7,
#'              n_fold = 4))
#' 
#' [[1]]
#' [[1]]$train
#'      Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 11            5.4         3.7          1.5         0.2  setosa           0
#' 17            5.4         3.9          1.3         0.4  setosa           0
#' 22            5.1         3.7          1.5         0.4  setosa           0
#' 19            5.7         3.8          1.7         0.3  setosa           0
#' 7             4.6         3.4          1.4         0.3  setosa           0
#' 6             5.4         3.9          1.7         0.4  setosa           0
#' 14            4.3         3.0          1.1         0.1  setosa           0
#' 7.1           4.6         3.4          1.4         0.3  setosa           0
#' 13            4.8         3.0          1.4         0.1  setosa           0
#' 17.1          5.4         3.9          1.3         0.4  setosa           0
#' 14.1          4.3         3.0          1.1         0.1  setosa           0
#' 23            4.6         3.6          1.0         0.2  setosa           0
#' 15            5.8         4.0          1.2         0.2  setosa           0
#' 1             5.1         3.5          1.4         0.2  setosa           0
#' 10            4.9         3.1          1.5         0.1  setosa           0
#' 14.2          4.3         3.0          1.1         0.1  setosa           0
#' 14.3          4.3         3.0          1.1         0.1  setosa           0
#' 3             4.7         3.2          1.3         0.2  setosa           0
#' 
#' [[1]]$test
#'    Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 6           5.4         3.9          1.7         0.4  setosa           1
#' 10          4.9         3.1          1.5         0.1  setosa           1
#' 22          5.1         3.7          1.5         0.4  setosa           1
#' 9           4.4         2.9          1.4         0.2  setosa           1
#' 21          5.4         3.4          1.7         0.2  setosa           1
#' 4           4.6         3.1          1.5         0.2  setosa           1
#' 3           4.7         3.2          1.3         0.2  setosa           1
#' 
#' 
#' [[2]]
#' [[2]]$train
#'      Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 21            5.4         3.4          1.7         0.2  setosa           0
#' 23            4.6         3.6          1.0         0.2  setosa           0
#' 12            4.8         3.4          1.6         0.2  setosa           0
#' 22            5.1         3.7          1.5         0.4  setosa           0
#' 3             4.7         3.2          1.3         0.2  setosa           0
#' 12.1          4.8         3.4          1.6         0.2  setosa           0
#' 15            5.8         4.0          1.2         0.2  setosa           0
#' 24            5.1         3.3          1.7         0.5  setosa           0
#' 12.2          4.8         3.4          1.6         0.2  setosa           0
#' 11            5.4         3.7          1.5         0.2  setosa           0
#' 15.1          5.8         4.0          1.2         0.2  setosa           0
#' 15.2          5.8         4.0          1.2         0.2  setosa           0
#' 6             5.4         3.9          1.7         0.4  setosa           0
#' 5             5.0         3.6          1.4         0.2  setosa           0
#' 7             4.6         3.4          1.4         0.3  setosa           0
#' 7.1           4.6         3.4          1.4         0.3  setosa           0
#' 4             4.6         3.1          1.5         0.2  setosa           0
#' 14            4.3         3.0          1.1         0.1  setosa           0
#' 
#' [[2]]$test
#'      Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 17            5.4         3.9          1.3         0.4  setosa           1
#' 15            5.8         4.0          1.2         0.2  setosa           1
#' 5             5.0         3.6          1.4         0.2  setosa           1
#' 5.1           5.0         3.6          1.4         0.2  setosa           1
#' 3             4.7         3.2          1.3         0.2  setosa           1
#' 23            4.6         3.6          1.0         0.2  setosa           1
#' 15.1          5.8         4.0          1.2         0.2  setosa           1
#' 
#' 
#' [[3]]
#' [[3]]$train
#'      Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 24            5.1         3.3          1.7         0.5  setosa           0
#' 9             4.4         2.9          1.4         0.2  setosa           0
#' 24.1          5.1         3.3          1.7         0.5  setosa           0
#' 20            5.1         3.8          1.5         0.3  setosa           0
#' 9.1           4.4         2.9          1.4         0.2  setosa           0
#' 18            5.1         3.5          1.4         0.3  setosa           0
#' 10            4.9         3.1          1.5         0.1  setosa           0
#' 18.1          5.1         3.5          1.4         0.3  setosa           0
#' 12            4.8         3.4          1.6         0.2  setosa           0
#' 5             5.0         3.6          1.4         0.2  setosa           0
#' 19            5.7         3.8          1.7         0.3  setosa           0
#' 2             4.9         3.0          1.4         0.2  setosa           0
#' 7             4.6         3.4          1.4         0.3  setosa           0
#' 23            4.6         3.6          1.0         0.2  setosa           0
#' 8             5.0         3.4          1.5         0.2  setosa           0
#' 17            5.4         3.9          1.3         0.4  setosa           0
#' 16            5.7         4.4          1.5         0.4  setosa           0
#' 2.1           4.9         3.0          1.4         0.2  setosa           0
#' 
#' [[3]]$test
#'     Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 24           5.1         3.3          1.7         0.5  setosa           1
#' 14           4.3         3.0          1.1         0.1  setosa           1
#' 8            5.0         3.4          1.5         0.2  setosa           1
#' 9            4.4         2.9          1.4         0.2  setosa           1
#' 5            5.0         3.6          1.4         0.2  setosa           1
#' 6            5.4         3.9          1.7         0.4  setosa           1
#' 9.1          4.4         2.9          1.4         0.2  setosa           1
#' 
#' 
#' [[4]]
#' [[4]]$train
#'      Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 22            5.1         3.7          1.5         0.4  setosa           0
#' 4             4.6         3.1          1.5         0.2  setosa           0
#' 1             5.1         3.5          1.4         0.2  setosa           0
#' 9             4.4         2.9          1.4         0.2  setosa           0
#' 4.1           4.6         3.1          1.5         0.2  setosa           0
#' 21            5.4         3.4          1.7         0.2  setosa           0
#' 14            4.3         3.0          1.1         0.1  setosa           0
#' 9.1           4.4         2.9          1.4         0.2  setosa           0
#' 3             4.7         3.2          1.3         0.2  setosa           0
#' 21.1          5.4         3.4          1.7         0.2  setosa           0
#' 20            5.1         3.8          1.5         0.3  setosa           0
#' 20.1          5.1         3.8          1.5         0.3  setosa           0
#' 23            4.6         3.6          1.0         0.2  setosa           0
#' 8             5.0         3.4          1.5         0.2  setosa           0
#' 9.2           4.4         2.9          1.4         0.2  setosa           0
#' 8.1           5.0         3.4          1.5         0.2  setosa           0
#' 15            5.8         4.0          1.2         0.2  setosa           0
#' 24            5.1         3.3          1.7         0.5  setosa           0
#' 
#' [[4]]$test
#'    Sepal.Length Sepal.Width Petal.Length Petal.Width Species test_status
#' 24          5.1         3.3          1.7         0.5  setosa           1
#' 23          4.6         3.6          1.0         0.2  setosa           1
#' 15          5.8         4.0          1.2         0.2  setosa           1
#' 4           4.6         3.1          1.5         0.2  setosa           1
#' 17          5.4         3.9          1.3         0.4  setosa           1
#' 3           4.7         3.2          1.3         0.2  setosa           1
#' 6           5.4         3.9          1.7         0.4  setosa           1
#'
#' @export

v_Rmach_fold <- function(inpt_datf, train_prop, n_fold){
  nb_train <- train_prop * nrow(inpt_datf)
  if (str_detect(pattern = "\\.", 
        string = nb_train)){
    nb_train <- round(x = nb_train, digits = 0)
  }
  if (nb_train == 1){
    return("Training number too high")
  }else if (nb_train == 0){
    return("Training number too low")
  }
  rtn_l <- list()
  for (I in 1:n_fold) {
    cur_datf <- cbind(inpt_datf[runif(n = nb_train, 
                min = 1, max = nrow(inpt_datf)),], 
                "test_status" = rep(x = 0, times = nb_train))
    cur_datf2 <- cbind(
                        inpt_datf[runif(n = (nrow(inpt_datf) - nb_train), 
                        min = 1, max = nrow(inpt_datf)),], 
                        "test_status" = rep(x = 1, times = (nrow(inpt_datf) - nb_train))
                  ) 
    rtn_l <- append(x = rtn_l, values = list(list("train" = cur_datf, "test" = cur_datf2)))
  }
  return(rtn_l)
}



