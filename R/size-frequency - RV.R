library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)

sex <- 1
year <- 2001:2021
language <- language("en")

source(paste0(getwd(),"/R/size-frequencies/station.number.rvset.R"))

# Read set card:
set <- read.gulf.set(year, password = password)
set$year <- year(date(set))
set <- set[-which(set$experiment %in% c(3, 8)), ]
set <- set[which(set$stratum %in% 415:439), ]
set <- set[which(set$vessel.code != "C"), ]      # Remove Cartier samples.
set$station.number <- station.number.rvset(set)

# Read length card:
len <- read.gulf.len(year, species = 2526, password = password)
len$year <- year(date(len))
ix <- match(len[key(set)], set[key(set)])
len <- len[!is.na(ix), ]
ix <- ix[!is.na(ix)]
set$start.time <- (as.POSIXct(paste0(set$date, set$start.time)) - date(set)) / 3600
len$start.time <- set$start.time[ix]

# Split non-sexed equally into male and female
vars <- paste0("freq", 1:14)
len[len$sex == 0 , vars] <- 0.5 * len[len$sex == 0 , vars]
len$sex[len$sex == 0] <- sex
len <- len[len$sex == sex, ]

# Needler and Teleost night length-based correction:
ix <- which(((len$start.time < 7) | (len$start.time > 19)) & (tolower(len$vessel.code) %in% c("n","t")))
vars <- paste0("freq", 1:14)
if (length(ix) > 0){
   lens <- repvec(len$start.length, ncol = length(vars)) + repvec(0:(length(vars)-1), nrow = nrow(len))
   lens[lens < 45] <- 45
   lens[lens > 115] <- 115
   coef <- -1.857 + (0.017*lens[ix, ])
   len[ix, vars] <- len[ix, vars] / exp(coef)
}

# Frequencies:
f <- freq(len, by = key(set))
fvars <- setdiff(names(f), key(set))

# Set crab sizes less than 5 mm to zero:
f[, fvars[as.numeric(fvars) <= 5]] <- 0

# Attach length-frequencies to set card:
ix <- match(f[key(set)], set[key(set)])
f <- f[!is.na(ix), ]
ix <- ix[!is.na(ix)]
set[fvars] <- 0
set[ix, fvars] <- f[fvars]
set[fvars] <- 1.75 * set[fvars] / repvec(set$distance, ncol = length(fvars))

# Aggregate out comparative experiments:
set <- aggregate(set[fvars], by = set[c("year", key(set), "stratum")], mean)

# Get stratum area:
stratum.info <- as.data.frame(read.gulf.spatial("stratum"))
stratum.info <- stratum.info[which((stratum.info$label %in% unique(set$stratum)) & (stratum.info$region == "gulf")), ]
area <- stratum.info$t_units
names(area) <-  stratum.info$label
area <- area / sum(area)

# Calculate stratified mean:
res <- aggregate(set[fvars], by = set[c("year", "stratum")], mean, na.rm = TRUE)
tmp <- res[fvars]
tmp[is.na(tmp)] <- 0
res[fvars] <- tmp
res[fvars] <- res[fvars] * repvec(area[as.character(res$stratum)], ncol = length(fvars))
res <- aggregate(res[fvars], by = res["year"], sum)

# Calculate standard error about the mean:
res.sd <- aggregate(set[fvars], by = set[c("year", "stratum")], function(x) sd(x, na.rm = TRUE)/sum(!is.na(x)))
tmp <- res.sd[fvars]
tmp[is.na(tmp)] <- 0
res.sd[fvars] <- tmp
res.sd[fvars] <- res.sd[fvars] * repvec(area[as.character(res.sd$stratum)], ncol = length(fvars))
res.sd <- aggregate(res.sd[fvars], by = res.sd["year"], sum)

clg()
if (sex == 1){
   xlim = c(0, 130)
   ylim = c(0, 2.0)
   breaks <- c(seq(0, 2, len = 12), 10)
}else{
   xlim = c(0, 95)
   ylim = c(0, 2.5)
   breaks <- c(seq(0, 2, len = 12), 10)
}
m <- kronecker(matrix(1:6, ncol = 1), matrix(1, nrow = 6, ncol = 6))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(year)){
   plot(xlim, ylim, xaxt = "n",  xaxs = "i", yaxs = "i", type = "n");
   grid();
   gbarplot(res[i, fvars], width = 1, xaxt = "n",  xaxs = "i", add = TRUE, border = "grey40")
   error.bar(as.numeric(fvars), lower = as.numeric(res[i, fvars] - 1.96*res.sd[i, fvars]), upper = as.numeric(res[i, fvars] + 1.96*res.sd[i, fvars]))
   lines(c(95,95), par("usr")[3:4], col = "red", lty = "dashed", lwd = 1)
   box()
   if (i == floor(length(year)/2)) mtext("Crab per tow", 2, 2.5, at = 0, cex = 1.50)
   text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.8 * diff(par("usr")[3:4]), year[i])
}
axis(1)
mtext("Carapace width (mm)", 1, 3.0, cex = 1.50)

# Heat map:
if (sex == 1) file <- "males" else file <- "females"
tiff(file = paste0("results/figures/size-frequencies/size-frequency.", file, ".", min(year),"-", max(year), " RV heatmap - ", language, ".tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width =7.5)
rownames(z) <- year
colnames(z) <- fvars
image(year, as.numeric(fvars), z, zlim = ylim, ylim = xlim, 
      xlab = "", ylab = "", xaxt = "n",
      col = hcl.colors(12, "YlOrRd", rev = TRUE),
      breaks = breaks)
hline(95, col = "red", lty = "dashed", lwd = 2)
axis(1, at = seq(2000, 2025, by = 2))
mtext("Year", 1, 2.75, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
box(col = "grey50")
dev.off()

# Linearize for output:
tab <- NULL
for (i in 1:length(fvars)){
   tmp <- data.frame(res$year, as.numeric(fvars[i]), res[, fvars[i]], res.sd[, fvars[i]])
   names(tmp) <- c("year", "cw", "mean", "sd")
   tab <- rbind(tab, tmp)
}

# Output:
write.csv(tab, file = paste0("results/tables/size-frequencies ", file, " RV.csv"), row.names = FALSE)

