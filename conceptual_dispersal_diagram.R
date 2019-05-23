#### Create a conceptual diagram of dispersal for each time period ####
library(RColorBrewer)
library(plotly)
cols <- brewer.pal(12, "Paired")
col2rgb(cols) # paired color vector

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/PADEconnectivity/conceptual_dispersal_diagram.png", width=11, height=5, res=300, units="in")

par(
  mfrow = c(1,3),
  mar=c(1, 0.5, 1, 0.5), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=15, # point size, which is the font size
  xpd=NA
)

## Early ##
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,20), ylim = c(0,22))
text(0.2,22, 'A', cex = 1.2)

# Plot time period
text(10,21.5, 'Early: 1989-1993', cex = 0.8)

# Plot ingress sites
text(2.4, 19.4, bquote(underline("Ingress")), cex = 0.8)
text(2.4, 18.5, bquote(underline("site")), cex = 0.8)
text(2, 14,"NJ")
text(2, 12.5,"DE")
text(2, 10.5,"VA")
text(2, 7.3,"NC")

# Plot adult origin sites
text(17, 19.5, bquote(underline("Origin pop")), cex = 0.8)
text(17, 18,"A")
text(17, 16,"B")
text(17, 14,"C")
text(17, 12,"D")
text(17, 10,"E")
text(17, 8,"F")
text(17, 6,"G")
text(17, 4,"H")
text(17, 2,"I")
text(17, 0,"J")

arrows(16,8,3.6,7.3, length = 0.1, lwd = 1, col = rgb(253/255,191/255,111/255,0.4)) # E1 NC

arrows(16,10.1,3.6,7.5, length = 0.1, lwd = 1, col = rgb(31/255,120/255,180/255,0.4)) # E2 NC

arrows(16,9.9,3.6,7.4, length = 0.1, lwd = 1, col = rgb(227/255,26/255,28/255,0.4)) # E3 NC

arrows(16,14,3.6,14, length = 0.1, lwd = 1, col = rgb(251/255,154/255,153/255,0.4)) # E4 NJ

arrows(16,6.1,3.6,7.2, length = 0.1, lwd = 1, col = rgb(51/255,160/255,44/255,0.4)) # E5 NC

arrows(16,5.9,3.6,7.1, length = 0.1, lwd = 1, col = rgb(177/255,89/255,40/255,0.4)) # E6 NC

box()
legend('bottomleft', legend = c('1-5', '6-10', '11-15', '16-20', '21-25', '26-30'), col = 'black', lwd = c(1:6), cex = 0.8, title = '# of larvae', text.font = 1)

## Middle ##
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,20), ylim = c(0,22))
text(0.2,22, 'B', cex = 1.2)

# Plot time period
text(10,21, 'Middle: 1998-2002', cex = 0.8)

# Plot ingress sites
text(2.4, 19.4, bquote(underline("Ingress")), cex = 0.8)
text(2.4, 18.5, bquote(underline("site")), cex = 0.8)
text(2, 14,"NJ")
text(2, 12.5,"DE")
text(2, 10.5,"VA")
text(2, 7.3,"NC")

# Plot adult origin sites
text(17, 19.5, bquote(underline("Origin pop")), cex = 0.8)
text(17, 18,"A")
text(17, 16,"B")
text(17, 14,"C")
text(17, 12,"D")
text(17, 10,"E")
text(17, 8,"F")
text(17, 6,"G")
text(17, 4,"H")
text(17, 2,"I")
text(17, 0,"J")

# Draw arrows of where each pie chart ingressed to #
arrows(16,10,3.6,14.1, length = 0.1, lwd = 5, col = rgb(202/255,178/255,214/255,0.4)) # M1 NJ

arrows(16,9.9,3.6,13.8, length = 0.1, lwd = 1, col = rgb(255/255,127/255,0/255,0.4)) # M2 NJ
arrows(16,9.9,3.6,7.3, length = 0.1, lwd = 5, col = rgb(255/255,127/255,0/255,0.4)) # M2 NC

box()

## Late ##
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '', xlim = c(0,20), ylim = c(0,22))
text(0.2,22, 'C', cex = 1.2)

# Plot time period
text(10,21.5, 'Late: 2008-2012', cex = 0.8)

# Plot ingress sites
text(2.4, 19.4, bquote(underline("Ingress")), cex = 0.8)
text(2.4, 18.5, bquote(underline("site")), cex = 0.8)
text(2, 14,"NJ")
text(2, 12.5,"DE")
text(2, 10.5,"VA")
text(2, 7.3,"NC")

# Plot adult origin sites
text(17, 19.5, bquote(underline("Origin pop")), cex = 0.8)
text(17, 18,"A")
text(17, 16,"B")
text(17, 14,"C")
text(17, 12,"D")
text(17, 10,"E")
text(17, 8,"F")
text(17, 6,"G")
text(17, 4,"H")
text(17, 2,"I")
text(17, 0,"J")

# Draw arrows of where each pie chart ingressed to 
arrows(16,8.1,3.6,7.3, length = 0.1, lwd = 1, col = rgb(251/255,154/255,153/255,0.4)) # L3 NC
arrows(16,8.1,3.6,10.4, length = 0.1, lwd = 3, col = rgb(251/255,154/255,153/255,0.4)) # L3 VA
arrows(16,8.1,3.6,12.5, length = 0.1, lwd = 2, col = rgb(251/255,154/255,153/255,0.4)) # L3 DE
arrows(16,8.1,3.6,14.1, length = 0.1, lwd = 1, col = rgb(251/255,154/255,153/255,0.4)) # L3 NJ

arrows(16,7.9,3.6,12.2, length = 0.1, lwd = 6, col = rgb(166/255,206/255,227/255,0.4)) # L1 DE
arrows(16,7.9,3.6,7.1, length = 0.1, lwd = 1, col = rgb(166/255,206/255,227/255,0.4)) # L1 NC
arrows(16,7.9,3.6,13.9, length = 0.1, lwd = 1, col = rgb(166/255,206/255,227/255,0.4)) # L1 NJ

arrows(16,10,3.6,7.4, length = 0.1, lwd = 2, col = rgb(51/255,160/255,44/255,0.4)) # L2 NC
arrows(16,10,3.6,10.5, length = 0.1, lwd = 2, col = rgb(51/255,160/255,44/255,0.4)) # L2 VA
arrows(16,10,3.6,12.6, length = 0.1, lwd = 1, col = rgb(51/255,160/255,44/255,0.4)) # L2 DE

box()

dev.off()
