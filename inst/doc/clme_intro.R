### R code from vignette source 'clme_intro.Rnw'

###################################################
### code chunk number 1: clme_intro.Rnw:413-414
###################################################
options( prompt="R> ") 


###################################################
### code chunk number 2: clme_intro.Rnw:417-419
###################################################
library("CLME")
data("rat.blood") 


###################################################
### code chunk number 3: clme_intro.Rnw:442-447
###################################################
const <- list(order = "simple", node = 1, decreasing = FALSE)
hct1 <- clme(hct ~ time + temp + sex + (1|id), data = rat.blood, seed = 42, nsim = 1,
   constraints = const, levels = list(2, levels(rat.blood$time)))
summary(hct1)
plot( hct1, ci = TRUE, legendx = "bottomright", inset = 0.08)


###################################################
### code chunk number 4: clme_intro.Rnw:454-455
###################################################
plot( hct1, ci = TRUE, legendx = "bottomright", inset = 0.08)


###################################################
### code chunk number 5: clme_intro.Rnw:464-467
###################################################
hct2 <- clme(hct ~ time + temp + sex + (1|id), data = rat.blood, seed = 42, nsim = 1, 
   levels = list(2, levels(rat.blood$time)))
summary(hct2)


###################################################
### code chunk number 6: clme_intro.Rnw:472-473
###################################################
invisible( dev.off() )


###################################################
### code chunk number 7: clme_intro.Rnw:478-482
###################################################
hct3 <- clme(hct ~ time + temp + sex + (1|id), data = rat.blood, seed = 42, nsim = 1, 
   gfix = rat.blood$time, constraints = const, 
   levels = list(2, levels(rat.blood$time)))
summary( hct3 )


###################################################
### code chunk number 8: clme_intro.Rnw:489-495
###################################################
const <- list( order = "simple.tree" , node = 1 , decreasing = FALSE)
wbc <- clme(wbc ~ time + temp + sex + (1|id), data = rat.blood, seed = 42, nsim = 1,
   constraints = const, levels = list(2, levels(rat.blood$time)),
   tsf = w.stat )
summary(wbc)
plot(wbc, legend = "topleft", inset = 0.08)


###################################################
### code chunk number 9: clme_intro.Rnw:500-501
###################################################
plot(wbc, legend = "topleft", inset = 0.08)


###################################################
### code chunk number 10: clme_intro.Rnw:519-524
###################################################
data("fibroid")
race.age <- factor(paste0( fibroid$race, ".", fibroid$age ) , 
   levels = c("Black.Yng", "Black.Mid", "Black.Old", 
   "White.Yng", "White.Mid", "White.Old") )
fibroid$race.age <- race.age


###################################################
### code chunk number 11: clme_intro.Rnw:529-535
###################################################
initVol <- rep("small", nrow(fibroid))
idx1 <- (14000 <= fibroid$vol & fibroid$vol < 65000)
idx2 <- (65000 <= fibroid$vol)
initVol[idx1] <- "medium"
initVol[idx2] <- "large"
fibroid$initVol <- factor(initVol, levels = c("small", "medium", "large"))


###################################################
### code chunk number 12: clme_intro.Rnw:540-544
###################################################
const <- list()
const$A <- cbind( 2:6 , 1:5 )[-3, ]
const$B <- rbind( c(3, 1), c(6, 4) )
const


###################################################
### code chunk number 13: clme_intro.Rnw:553-563
###################################################
w.blk.wht <- function (theta, cov.theta, B, A, ...) {
   stats <- vector("numeric", length = nrow(B))
   ctd   <- diag(cov.theta)
   stats <- apply(B, 1, FUN = function(a, theta, cov, ctd) {
     std <- sqrt(ctd[a[1]] + ctd[a[2]] - 2 * cov.theta[a[1], a[2]])
     (theta[a[2]] - theta[a[1]])/std
   }, theta = theta, cov = cov.theta, ctd = ctd)
   names(stats) <- c("Black.Yng - Black.Old", "White.Yng - White.Old" )
   return(stats)
 }


###################################################
### code chunk number 14: clme_intro.Rnw:570-573
###################################################
fib <- clme(lfgr ~ race.age + initVol + (1|id), data = fibroid, seed = 42, nsim = 1,
   constraints = const, tsf = w.blk.wht, levels = list(10, levels(race.age)))
summary( fib )


###################################################
### code chunk number 15: clme_intro.Rnw:576-577
###################################################
invisible( dev.off() )


###################################################
### code chunk number 16: fibroidfig
###################################################
plot(x = 1, y = 0, col = 0, ylim = c(-6, 22), xlim = c(0.9, 3.1), 
   xlab = "", ylab = "Estimated Coefficient", xaxt = "n")
axis(side=1, at=1:3, 
   labels = c("Young (<35)", "Middle aged (35-44)", "Older (>45)") )
for( y1 in seq(-15, 25, 5) ){
   lines( x = c(0, 7), y = c(y1, y1), col = "grey", lty = 2 )
}
lty1 <- 1 + (fib$p.value.ind[1]<0.05)
lty2 <- 1 + (fib$p.value.ind[2]<0.05)
lty3 <- 1 + (fib$p.value.ind[3]<0.05)
lty4 <- 1 + (fib$p.value.ind[4]<0.05)
points(c(1, 2), fib$theta[1:2], col = 1, type = "l", lwd = 2 , lty = lty1)
points(c(2, 3), fib$theta[2:3], col = 1, type = "l", lwd = 2 , lty = lty2)
points(c(1, 2), fib$theta[4:5], col = 3, type = "l", lwd = 2 , lty = lty3)
points(c(2, 3), fib$theta[5:6], col = 3, type = "l", lwd = 2 , lty = lty4)
points(1:3, fib$theta[1:3], col = 1, cex = 1.5, pch = 21, bg = "white")
points(1:3, fib$theta[4:6], col = 3, cex = 1.5, pch = 24, bg = "white")
legend("bottom", lty = c(1, 1), pch = c(21, 24), col = c(1, 3), pt.bg = 0,
   pt.cex = 1.1, lwd = 2, inset = 0.03, cex = 0.9,
   legend = c("Blacks    ", "Whites"))


###################################################
### code chunk number 17: clme_intro.Rnw:607-608
###################################################
plot(x = 1, y = 0, col = 0, ylim = c(-6, 22), xlim = c(0.9, 3.1), 
   xlab = "", ylab = "Estimated Coefficient", xaxt = "n")
axis(side=1, at=1:3, 
   labels = c("Young (<35)", "Middle aged (35-44)", "Older (>45)") )
for( y1 in seq(-15, 25, 5) ){
   lines( x = c(0, 7), y = c(y1, y1), col = "grey", lty = 2 )
}
lty1 <- 1 + (fib$p.value.ind[1]<0.05)
lty2 <- 1 + (fib$p.value.ind[2]<0.05)
lty3 <- 1 + (fib$p.value.ind[3]<0.05)
lty4 <- 1 + (fib$p.value.ind[4]<0.05)
points(c(1, 2), fib$theta[1:2], col = 1, type = "l", lwd = 2 , lty = lty1)
points(c(2, 3), fib$theta[2:3], col = 1, type = "l", lwd = 2 , lty = lty2)
points(c(1, 2), fib$theta[4:5], col = 3, type = "l", lwd = 2 , lty = lty3)
points(c(2, 3), fib$theta[5:6], col = 3, type = "l", lwd = 2 , lty = lty4)
points(1:3, fib$theta[1:3], col = 1, cex = 1.5, pch = 21, bg = "white")
points(1:3, fib$theta[4:6], col = 3, cex = 1.5, pch = 24, bg = "white")
legend("bottom", lty = c(1, 1), pch = c(21, 24), col = c(1, 3), pt.bg = 0,
   pt.cex = 1.1, lwd = 2, inset = 0.03, cex = 0.9,
   legend = c("Blacks    ", "Whites"))


