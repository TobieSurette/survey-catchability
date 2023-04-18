library(gulf.data)
library(gulf.stats)
library(gulf.graphics)

years <- 2001:2022
sex <- 1

# Read data:
s <- read.scsset(year = years, survey = "regular", valid = 1) # Tow data.
b <- read.scsbio(year = years, survey = "regular", sex = sex)   # Biological data.
b$tow.id <- tow.id(b)

# Determine morphometric maturity:
for (i in 1:length(years)){
   print(years[i])
   ix <- year(b) == years[i]
   bb <- b[ix, ]
   z <- rep(NA, nrow(bb))
   z[which(bb$carapace.width < 30)] <- 0 # Crab smaller than 30mm are considered immature.
   theta <- fit.morphometry.scsbio(bb$carapace.width, bb$chela.height, z, sex = 1, discrete = years[i] < 1998) # Fit morphometric model.
   v <- morphometry.scsbio(bb$carapace.width, bb$chela.height, theta = theta, discrete = years[i] < 1998)$p_mature_posterior
   b$maturity[ix] <- v > 0.5
}
b$maturity <- b$maturity == "TRUE"

# Attach size-frequencies:
groups <- c("immature", "skip", "recruit", "residual")
res <- list()
for (i in 1:length(groups)) res[i] <- list(NULL)
names(res) <- groups

# Immatures:
tmp <- s
import(tmp, fill = 0) <- freq(b[which(!b$maturity & is.new.shell(b)),], by = c("date", "tow.id"))
res$immature <- tmp

# Skip:
tmp <- s
import(tmp, fill = 0) <- freq(b[which(!b$maturity & !is.new.shell(b)),], by = c("date", "tow.id"))
res$skip <- tmp

# Recruit:
tmp <- s
import(tmp, fill = 0) <- freq(b[which(b$maturity & is.new.shell(b)),], by = c("date", "tow.id"))
res$recruit <- tmp

# Residuals:
tmp <- s
import(tmp, fill = 0) <- freq(b[which(b$maturity & !is.new.shell(b)),], by = c("date", "tow.id"))
res$residual <- tmp

# Standardize variables:
for (i in 1:length(res)){
   fvars <- names(res[[i]])[gsub("[0-9]", "", names(res[[i]])) == ""]
   res[[i]][fvars] <- 10^6 * res[[i]][fvars] / repvec(res[[i]]$swept.area, ncol = length(fvars))
}

# Frequency variables:
fvars <- unique(unlist(lapply(res, names)))
fvars <- fvars[gsub("[0-9]", "", fvars) == ""]

# Buffer frequencies with zeroes:
for (i in 1:length(res)) res[[i]][setdiff(fvars, names(res[[i]]))] <- 0

# Calculate annual size-frequencies:
f <- array(NA, dim = c(length(years), length(fvars), length(groups)))
dimnames(f) <- list(year = years, size = fvars, group = groups)
for (i in 1:length(groups)){
   tmp <- aggregate(res[[i]][fvars], by = list(year = year(res[[i]])), mean)
   rownames(tmp) <- tmp$year
   f[,,i] <- as.matrix(tmp[, -1])
}

# Heat map:
clg()
if (sex == 1) file <- "males" else file <- "females"
cols <- hcl.colors(12, "YlOrRd", rev = TRUE)
for (i in 1:length(dimnames(f)$group)){
   file <- paste0("results/figures/size-frequencies/size-frequency ", sex(sex), " ", dimnames(f)$group[i], " ", min(year),"-", max(year), " heatmap - ", language, ".tiff")
   print(file)
   tiff(file = file, compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
   zlim <- c(0, 150)
   if (dimnames(f)$group[i] == "immature") zlim <- c(0, 400)
   ylim <- c(0, 130)
   breaks <- c(seq(zlim[1], zlim[2], len = length(cols)), max(f))
   image(as.numeric(dimnames(f)$year), as.numeric(dimnames(f)$size), f[,,i], zlim = zlim, ylim = xlim, 
         xlab = "", ylab = "", xaxt = "n",
         col = cols,
         breaks = breaks)
   hline(95, col = "red", lty = "dashed", lwd = 2)
   axis(1, at = seq(2000, 2025, by = 2))
   mtext("Year", 1, 2.75, cex = 1.25)
   mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
   box(col = "grey50")
   mtext(dimnames(f)$group[i], 3, 0.5, cex = 1.5)
   dev.off()
}

clg()
if (sex == 1) file <- "males" else file <- "females"
cols <- hcl.colors(12, "YlOrRd", rev = TRUE)
file <- paste0("results/figures/size-frequencies/size-frequency ", sex(sex), " ", "all", " ", min(year),"-", max(year), " heatmap - ", language, ".tiff")
print(file)
tiff(file = file, compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
zlim <- c(0, 400)
ylim <- c(0, 130)
breaks <- c(seq(zlim[1], zlim[2], len = length(cols)), max(f))
image(as.numeric(dimnames(f)$year), as.numeric(dimnames(f)$size), apply(f, c(1,2), sum), zlim = zlim, ylim = xlim, 
      xlab = "", ylab = "", xaxt = "n",
      col = cols,
      breaks = breaks)
hline(95, col = "red", lty = "dashed", lwd = 2)
axis(1, at = seq(2000, 2025, by = 2))
mtext("Year", 1, 2.75, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
box(col = "grey50")
#mtext(dimnames(f)$group[i], 3, 0.5, cex = 1.5)
dev.off()

save(f, file = "results/tables/size-frequencies males scs 2001-2022.rdata")

