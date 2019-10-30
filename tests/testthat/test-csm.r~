library(testthat)
library(ape)
library(csm)

set.seed(98749)

test_that("birth_only_constant_crown",{

   # rough test of pure birth
   # exp growth with N(0) = 2
   p <- 5
   n <- 1000
   v <- mat.or.vec(n,1)
   cat(" performing simulations... \n")
   i <- 1
   while (i < n) { 
      s  <- CSMsim(p,p,0,0,0,crown_or_stem=0,pp=1,trait_step=0.01)

      valid = FALSE
      if ( s$run_status == "timeout") { 
         valid = TRUE
         v[i] <- nrow(s$csm_traits)
      }
      if ( s$run_status == "extinction") {
         if ( s$net_tip_count > 1  ) { 
            valid = TRUE 
         }
      }

      if (valid) {

        i <- i + 1
        cat(i, " ")
      }

      if (!valid) {
         cat(" ebc ")
      }

   }
   cat("  \n")
   # analytical expression
   a <- 2*exp(p)
   cat("expected ", a, "\n")
   cat("observed ", mean(v), "\n")
   expect_equal( mean(v) , a, tolerance=0.1, scale=a)

})

test_that("birth_death_constant_crown",{

   # rough test of pure birth
   # exp growth with N(0) = 1
   p1 <- 7
   p2 <- 2
   n <- 1000
   v <- mat.or.vec(n,1)
   cat(" performing simulations... \n")
   i <- 1
   while (i < n) { 
      s  <- CSMsim(p1,p1,0,p2,p2,crown_or_stem=0,pp=1,trait_step=0.01)

      valid = FALSE
      if ( s$run_status == "timeout") { 
         valid = TRUE
         v[i] <- nrow(s$csm_traits)
      }
      if ( s$run_status == "extinction") {
         if ( s$net_tip_count > 1  ) { 
            valid = TRUE 
         }
      }

      if (valid) {

        i <- i + 1
        cat(i, " ")
      }

      if (!valid) {
         cat(" ebc ", s$run_status, "\n")
      }

   }
   cat("  \n")
   a <- 2*exp(p1 - p2)
   cat("expected ", a, "\n")
   cat("observed ", mean(v), "\n")
   expect_equal( mean(v) , a, tolerance=0.1, scale=a)

})


test_that("test_randomness_parent_assignment",{

   # rough test of pure birth
   # exp growth with N(0) = 1
   n <- 2000
   v <- mat.or.vec(n,1)
   cat(" performing simulations... \n")
   i <- 1
   while (i < n) { 
      s  <- csm(3,3,0,0,0,0,0.1,0.1)

      valid = FALSE
      if ( s[[1]] == "timeout") { 
         valid = TRUE
      }
      if ( nrow(s[[2]]) < 5) {
         valid = FALSE 
      }

      if (valid) {
        v[i] <- s[[2]][5,3]
        i <- i + 1
        cat(i, " ")
      }

      if (!valid) {
         cat(" ex ")
      }

   }
   cat("  \n")
   a <- n/4
   cat("expected ", a, "\n")
   cat("observed ", length(which(v==0)), "\n")
   cat("observed ", length(which(v==1)), "\n")
   cat("observed ", length(which(v==2)), "\n")
   cat("observed ", length(which(v==3)), "\n")

   expect_equal( length(which(v==0)) , a, tolerance=0.1, scale=a)
   expect_equal( length(which(v==1)) , a, tolerance=0.1, scale=a)
   expect_equal( length(which(v==2)) , a, tolerance=0.1, scale=a)
   expect_equal( length(which(v==3)) , a, tolerance=0.1, scale=a)

})


test_that("test_choose_integer",{

   n <- 20000
   v <- mat.or.vec(n,1)
   vs <- mat.or.vec(n,1)
   m <- 12

   for (i in 1:n)  {
      v[i] <- choose_random_integer(m)  
      vs[i] <- sample(1:m)[1]      
   }

   a <- 1/m
   t1 <- length(which(v==1))/n
   t4 <- length(which(v==4))/n
   t12 <- length(which(v==12))/n

   expect_equal( t1 , a, tolerance=0.1, scale=a)
   expect_equal( t4 , a, tolerance=0.1, scale=a)
   expect_equal( t12 , a, tolerance=0.1, scale=a)

   sc <- sum(abs(table(v) - n/m))
   ss <- sum(abs(table(vs) - n/m))

   cat("  total departure from uniform \n")
   cat("  From", n, "rcpp util:", sc, "\n")
   cat("  From", n, "R sample :", ss, "\n")

   expect_equal( sc , ss, tolerance=1.5, scale=ss)

})



test_that("test_get_reaction",{

   n <- 10000
   v <- mat.or.vec(n,1)

   c1=1;c2=1;c3=1;c4=1;c5=1;

   for (i in 1:n)  {
      v[i] <- get_reaction(c1,c2,c3, c4,c5)  
   }

   a1 <- c1 / (c1+c2+c3+c4+c5)
   a3 <- c3 / (c1+c2+c3+c4+c5)
   a5 <- c5 / (c1+c2+c3+c4+c5)

   t1 <- length(which(v==1))/n
   t3 <- length(which(v==3))/n
   t5 <- length(which(v==5))/n

   expect_equal( t1 , a1, tolerance=0.1, scale=a1)
   expect_equal( t3 , a3, tolerance=0.1, scale=a3)
   expect_equal( t5 , a5, tolerance=0.1, scale=a5)


   c1=10;c2=0;c3=5;c4=1;c5=0;

   for (i in 1:n)  {
      v[i] <- get_reaction(c1,c2,c3, c4,c5)  
   }

   a1 <- c1 / (c1+c2+c3+c4+c5)
   a3 <- c3 / (c1+c2+c3+c4+c5)
   a5 <- c5 / (c1+c2+c3+c4+c5)

   t1 <- length(which(v==1))/n
   t3 <- length(which(v==3))/n
   t5 <- length(which(v==5))/n

   expect_equal( t1 , a1, tolerance=0.1, scale=a1)
   expect_equal( t3 , a3, tolerance=0.1, scale=a3)
   expect_equal( t5 , a5, tolerance=0.1, scale=a5)



})


