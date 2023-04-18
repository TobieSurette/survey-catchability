library(gulf.data)
library(gulf.graphics)

# Load RV data:
x <- read.csv("results/tables/size-frequencies males RV.csv")
x <- x[x$cw <= 150, ]
years <- sort(unique(x$year))
xx <- matrix(NA, nrow = length(years), ncol = 150)
dimnames(xx) <- list(year = years, cw = 1:150)
for (i in 1:nrow(xx)){
   ix <- x$year == years[i]
   xx[i, as.character(x$cw[ix])] <- x$mean[ix]
}

# Read SCS data:
load("results/tables/size-frequencies males scs 2001-2022.rdata")
f <- apply(f, c(1,2), sum)
f <- f[,as.character(1:150)]
f <- f[rownames(xx), ]

# Heatmap plot SCS:
tiff(file = paste0("results/figures/size-frequencies/size-frequency male all 2001-2021 heatmap SCS - english.tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
zlim <- c(0, 400.0)
f[f >= zlim[2]] <- zlim[2]
image(as.numeric(rownames(f)), as.numeric(colnames(f)), f, 
      xlab = "", ylab = "", zlim = zlim, ylim = c(0, 140), yaxs = "i", 
      xaxt = "n", col = hcl.colors(150, "YlOrRd", rev = TRUE))
axis(1, at = seq(2000, 2022, by = 2))
hline(95, col = "red", lty = "dashed", lwd = 2)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.25, cex = 1.25)
box(col = "grey70")
dev.off()

# Heatmap plot SCS:
load("results/tables/size-frequencies males scs 2001-2022.rdata")
f <- apply(f, c(1,2), sum)
f <- f[,as.character(1:150)]
tiff(file = paste0("results/figures/size-frequencies/size-frequency male all 2001-2022 heatmap SCS - english.tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
zlim <- c(0, 400.0)
f[f >= zlim[2]] <- zlim[2]
image(as.numeric(rownames(f)), as.numeric(colnames(f)), f, 
      xlab = "", ylab = "", zlim = zlim, ylim = c(0, 140), yaxs = "i", 
      xaxt = "n", col = hcl.colors(150, "YlOrRd", rev = TRUE))
axis(1, at = seq(2000, 2022, by = 2))
hline(95, col = "red", lty = "dashed", lwd = 2)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.25, cex = 1.25)
box(col = "grey70")
dev.off()

# Heatmap plot RV:
tiff(file = paste0("results/figures/size-frequencies/size-frequency male all 2001-2021 heatmap RV - english.tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
zlim <- c(0, 1.25)
xx[xx >= zlim[2]] <- zlim[2]
image(as.numeric(rownames(xx)), as.numeric(colnames(xx)), xx, 
      xlab = "", ylab = "", zlim = zlim, ylim = c(0, 140), yaxs = "i", 
      xaxt = "n", col = hcl.colors(150, "YlOrRd", rev = TRUE))
axis(1, at = seq(2000, 2022, by = 2))
hline(95, col = "red", lty = "dashed", lwd = 2)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.25, cex = 1.25)
box(col = "grey70")
dev.off()

# Ratio plot:
load("results/tables/size-frequencies males scs 2001-2022.rdata")
f <- apply(f, c(1,2), sum)
f <- f[,as.character(1:150)]
f <- f[rownames(xx), ]
# Load RV data:
x <- read.csv("results/tables/size-frequencies males RV.csv")
x <- x[x$cw <= 150, ]
years <- sort(unique(x$year))
xx <- matrix(NA, nrow = length(years), ncol = 150)
dimnames(xx) <- list(year = years, cw = 1:150)
for (i in 1:nrow(xx)){
   ix <- x$year == years[i]
   xx[i, as.character(x$cw[ix])] <- x$mean[ix]
}
zlim <- c(-1.5, 10.5)
r <- log(f / xx)
r[r > zlim[2]] <- zlim[2]
image(as.numeric(rownames(xx)), as.numeric(colnames(xx)), r, 
      xlab = "", ylab = "", zlim = zlim, ylim = c(0, 140), yaxs = "i", 
      xaxt = "n", col = hcl.colors(150, "YlOrRd", rev = TRUE))
axis(1, at = seq(2000, 2022, by = 2))
hline(95, col = "red", lty = "dashed", lwd = 2)
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.25, cex = 1.25)
box(col = "grey70")

# Plot ratio by size:
smooth = TRUE
log = TRUE
legend = TRUE
file = paste0("results/figures/size-frequencies/SCS to RV ratios male ", ifelse(log, "log ", "") , ifelse(smooth, "smooth ", ""), "- english.tiff")
tiff(file = file, compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
if (log) fun <- base::log else fun <- function(x) x
plot(c(0, 140), fun(exp(c(-1.5, 8))), type = "n", xlab = "", ylab = "", xaxs = "i")
grid()
for (i in 1:length(years)){
   r <- log(f[i,] / xx[i,])
   y <- 1:150
   ix <- is.finite(r) & !is.na(r)
   r <- r[ix]
   y <- y[ix]
   r <- exp(r)
   if (!smooth) lines(y, fun(r), col = rainbow(length(years))[i])
   r <- log(r)
   m <- mgcv::gam(r ~ s(y))
   if (smooth) lines(y, fun(exp(predict(m))), lwd = 2, col = rainbow(length(years))[i])
}
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
if (log) mtext("log(ratio) (SCS/RV)", 2, 2.5, cex = 1.25) else mtext("Ratio (SCS/RV)", 2, 2.5, cex = 1.25)
pos <- "bottomright"
if (!log) pos <- "topleft"
legend(pos, legend = years, col = rainbow(length(years)), cex = 0.6, lwd = 2)
box(col = "grey60")
dev.off()

windows()
gbarplot(apply(xx, 2, sum), grid = TRUE)

gbarplot(apply(f, 2, sum), grid = TRUE)

gbarplot(log(apply(f, 2, sum)), grid = TRUE)

# GAM analysis:
r <- log(f / xx)
ix <- is.finite(r) & !is.na(r)
data <- data.frame(year = as.factor(repvec(as.numeric(rownames(r)), ncol = ncol(r))[ix]),
                   cw = repvec(as.numeric(colnames(r)), nrow = nrow(r))[ix],
                   ratio = r[ix])

model <- mgcv::gam(ratio ~ s(cw) + year, data = data)
res <- as.data.frame(predict(model, newdata = data.frame(year = sort(unique(data$year)), cw = 90), se.fit = TRUE))
res$year <- as.numeric(as.character(sort(unique(data$year))))

file = paste0("results/figures/size-frequencies/SCS to RV ratios male GAM - english.tiff")
tiff(file = file, compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
gbarplot(res$fit, res$year, xaxt = "n")
error.bar(res$year, lower = res$fit - 1.96 * res$se.fit, upper = res$fit + 1.96 * res$se.fit)
axis(1, at = seq(2000, 2022, by = 2))
mtext("Year", 1, 2.5, cex = 1.25)
mtext("SCS to RV ratio", 2, 2.5, cex = 1.25)
box(col = "grey60")
dev.off()

