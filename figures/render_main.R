#!/usr/bin/env Rscript
# Figure coordinates from saved results and the exact frozen step response.
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))
source(file.path(here, "bj_style.R"))
load_json <- function(name) jsonlite::fromJSON(file.path(data_dir, paste0(name, ".json")))
joint <- load_json("classification")
resolution <- load_json("independent_gate_readout")
steps <- function(t, constants) {
  vapply(pmax(t, 0), function(x) {
    weights <- vapply(constants, function(z) prod(z / (z - constants[constants != z])), numeric(1))
    1 - sum(weights * exp(-x / constants))
  }, numeric(1))
}
pulse <- function(t, start, end, amplitude, constants) {
  amplitude * (steps(t - start, constants) - steps(t - end, constants))
}
t <- 0:300
q <- pulse(t, 20, 80, .3, 8)
a <- pulse(t, 20, 80, .3, c(8,25)) + pulse(t,140,200,.2,25)
x <- pulse(t, 20, 80, .3, c(8,25,20)) + pulse(t,140,200,.2,c(25,20))
f <- .32*q + .68*a
stopifnot(abs(max(abs(q-f))-joint$adversarial_twin$latent_gate_max_difference) < 1e-12)
long <- function(values, labels) data.frame(t=rep(t,length(values)),
   value=unlist(values), trace=factor(rep(labels,each=length(t)),levels=labels))
shaded <- function() list(
  annotate("rect", xmin=20,xmax=80,ymin=-Inf,ymax=Inf,fill=green_tint),
  annotate("rect", xmin=140,xmax=200,ymin=-Inf,ymax=Inf,fill=lilac_pale))
arms <- c(u="Gate input",v="Actin input",both="Both inputs")
channels <- c(x="Deformation",FRET="FRET",actin="Actin",joint="Joint")
bars <- do.call(rbind,lapply(names(arms),function(arm) do.call(rbind,lapply(names(channels),function(ch)
  data.frame(arm=arms[[arm]],channel=channels[[ch]],value=joint$classification[[arm]][[ch]]$accuracy)))))
bars$arm <- factor(bars$arm,levels=arms)
bars$channel <- factor(bars$channel,levels=channels)
p1 <- ggplot(bars,aes(arm,channel,fill=value)) +
  geom_tile(width=.96,height=.94,colour="white",linewidth=.5) +
  geom_text(aes(label=sprintf("%.2f%%",100*value)),family=figure_font,size=3.2,colour=ink) +
  scale_fill_gradient2(low=lilac_pale,mid="#F8FAFC",high=green_fill,
                       midpoint=.5,limits=c(0,1),guide="none") +
  scale_y_discrete(limits=rev(levels(bars$channel))) +
  labs(x=NULL,y=NULL,title="Accuracy by input and readout") + th() +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        plot.title=element_text(size=9.5,margin=margin(t=7,b=4)))
p2 <- ggplot(long(list(q,f),c("Model A","Model B")),aes(t,value,colour=trace,linetype=trace)) +
  shaded() + geom_line(linewidth=trace_width) +
  scale_colour_manual(values=c(forest,plum)) + scale_linetype_manual(values=c(1,2)) +
  labs(x="Time (assigned s)",y="Gate deviation") + th()
p3 <- ggplot(long(list(f,a,x),c("Mixed FRET","Actin","Deformation")),aes(t,value,colour=trace,linetype=trace)) +
  shaded() + geom_line(linewidth=trace_width) +
  scale_colour_manual(values=c(forest,plum,green_ink)) +
  scale_linetype_manual(values=comparison_lines) +
  labs(x="Time (assigned s)",y="Output deviation") + th() +
  theme(legend.key.width=grid::unit(6,"mm"),legend.key.spacing.x=grid::unit(2,"pt"),
        legend.text=element_text(size=8.2)) +
  guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1))
r1 <- data.frame(b=seq(0,1,.01))
p4 <- ggplot(r1,aes(b)) +
  geom_ribbon(aes(ymin=1-b,ymax=1,fill="Model A"),colour=forest,linewidth=.38) +
  geom_ribbon(aes(ymin=.32*(1-b),ymax=.32,fill="Model B"),colour=plum,linewidth=.38,alpha=.85) +
  geom_vline(xintercept=.68,linetype=2,colour=muted_ink,linewidth=reference_width) +
  scale_fill_manual(values=c("Model A"=green_fill,"Model B"=lilac_pale)) +
  scale_x_continuous(breaks=c(0,.2,.5,.68,1)) +
  labs(x="Upper bound on actin sensitivity",y="Coefficient") + th()
r2 <- subset(resolution$rows,tau==10)
r2$bound <- factor(sprintf("%.2f",r2$bound),levels=c("0.20","0.50","0.68"))
p5 <- ggplot(r2,aes(cadence,correct,colour=bound,shape=bound,linetype=bound)) +
  annotate("rect",xmin=-Inf,xmax=Inf,ymin=0,ymax=.05,fill=pink_tint) +
  geom_line(linewidth=trace_width) + geom_point(size=2.4) +
  scale_colour_manual(values=c(forest,plum,green_ink)) +
  scale_linetype_manual(values=comparison_lines) +
  scale_x_continuous(breaks=c(5,15,60)) +
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1.03),expand=expansion(mult=c(0,0))) +
  labs(x="Sampling interval (assigned s)",y="P(correct)") +
  th() + labs(colour="Actin bound",shape="Actin bound",linetype="Actin bound") +
  theme(legend.title=element_text(size=8.5,face="bold"))
order_counts <- read.csv(file.path(data_dir,"pulse_order_counts.csv"))
order_counts$order <- factor(order_counts$order, levels=rev(order_counts$order))
p6 <- ggplot(order_counts,aes(correct/total,order)) +
  geom_vline(xintercept=.5,colour=muted_ink,linetype=2,linewidth=reference_width) +
  geom_segment(aes(x=.5,xend=correct/total,yend=order),colour=green_ink,linewidth=.53) +
  geom_point(aes(fill=order,shape=order),size=3.25,stroke=mark_stroke,colour=ink) +
  geom_text(aes(label=paste0(correct,"/",total)),nudge_y=.22,family=figure_font,size=3.15,colour=ink) +
  scale_fill_manual(values=c("Actin then gate"=lilac_fill,"Gate then actin"=green_fill),guide="none") +
  scale_shape_manual(values=c(22,21),guide="none") +
  scale_x_continuous(breaks=seq(0,1,.25),labels=scales::label_percent(),
                     limits=c(0,1.10),expand=expansion(mult=c(.01,.01))) +
  labs(x="Joint accuracy",y=NULL,title="Fixed classifier under pulse reversal") + th() +
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
        plot.title=element_text(size=9.5,margin=margin(t=7,b=4)))

lag <- subset(resolution$rows, cadence==60 & bound<.68)
lag$bound <- factor(sprintf("%.2f",lag$bound))
p7 <- ggplot(lag,aes(tau,correct,colour=bound,shape=bound,linetype=bound)) +
  geom_line(linewidth=trace_width) + geom_point(size=2.4) +
  scale_colour_manual(values=c(forest,plum)) +
  scale_linetype_manual(values=comparison_lines[1:2]) +
  scale_x_continuous(breaks=c(0,10,30)) +
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1.03),expand=expansion(mult=c(0,0))) +
  labs(x="Reporter response time (assigned s)",y="P(correct) at 60 s",
       colour="Actin bound",shape="Actin bound",linetype="Actin bound") + th() +
  theme(legend.title=element_text(size=8.5,face="bold"))

res <- subset(resolution$rows, tau==10 & cadence==15)
outcomes <- rbind(data.frame(bound=res$bound, probability=res$correct, outcome="Correct"),
                  data.frame(bound=res$bound, probability=1-res$correct-res$wrong, outcome="Abstain"),
                  data.frame(bound=res$bound, probability=res$wrong, outcome="Wrong"))
outcomes$outcome <- factor(outcomes$outcome,levels=c("Wrong","Abstain","Correct"))
p8 <- ggplot(outcomes,aes(factor(bound),probability,fill=outcome)) +
  geom_col(width=.64,colour=ink,linewidth=.2) +
  scale_fill_manual(values=c("Correct"=green_fill,"Abstain"=pink_fill,"Wrong"=lilac_fill),
                    breaks=c("Correct","Abstain","Wrong")) +
  scale_y_continuous(labels=scales::label_percent(),limits=c(0,1.02),expand=expansion(mult=c(0,0))) +
  labs(x="Upper bound on actin sensitivity",y="Decision fraction") + th()

# Diagram marks are conceptual. Only the kinetics in the first schematic
# carry model values. All geometry and text are drawn as R vector primitives.
diagram <- function() {
  ggplot() +
    coord_cartesian(xlim=c(0,100),ylim=c(0,26),expand=FALSE) +
    theme_void(base_family=figure_font) + theme(plot.margin=margin(14,8,1,8))
}
field <- function(left,right,fill) {
  annotate("rect",xmin=left,xmax=right,ymin=1,ymax=25,fill=fill,
           colour=NA)
}
word <- function(x,y,label,size=3.1,colour=ink,bold=FALSE,hjust=.5) {
  annotate("text",x=x,y=y,label=label,family=figure_font,size=size,
           colour=colour,fontface=if (bold) "bold" else "plain",hjust=hjust,
           lineheight=1.04)
}
link <- function(x,y,xend,yend,colour=forest,arrowhead=TRUE,linetype=1) {
  annotate("segment",x=x,y=y,xend=xend,yend=yend,colour=colour,linetype=linetype,
           linewidth=.40,arrow=if (arrowhead) grid::arrow(length=grid::unit(1.45,"mm"),type="closed"))
}
node <- function(x,y,w,h,label,fill,edge=forest,size=3.05) {
  list(annotate("rect",xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2,
                fill=fill,colour=edge,linewidth=.31),
       word(x,y,label,size,edge,TRUE))
}

a1 <- diagram() +
  field(2,40,green_tint) + field(44,66,green_tint) + field(70,99,lilac_tint) +
  word(4,22.3,"Exchange the relaxation times",3.52,forest,TRUE,0) +
  word(46,22.3,"Shared observations",3.38,forest,TRUE,0) +
  word(72,22.3,"Separate classifier test",3.38,plum,TRUE,0) +
  word(c(16,31),17.2,c("Gate","Actin"),2.95,forest,TRUE) +
  word(5,c(12.6,6),c("A","B"),3.4,forest,TRUE) +
  node(16,12.6,11.2,4.4,"8 s", "white") +
  node(31,12.6,11.2,4.4,"25 s",green_pale) +
  node(16,6,11.2,4.4,"25 s",lilac_pale,plum) +
  node(31,6,11.2,4.4,"8 s","white") +
  link(21.9,c(12.6,6),25,c(12.6,6)) +
  link(40.1,13,43.6,13) +
  word(55,16.8,"Same inputs",2.9,muted_ink) +
  link(55,14.4,55,11.7) +
  word(55,8.8,"Mixed FRET\nActin · deformation",3.1,forest,TRUE) +
  word(72,c(13.5,7.2),c("Original","Reversed"),2.9,ink,FALSE,0) +
  node(87,13.5,7.3,4.3,"Gate",green_pale,size=2.78) +
  node(95,13.5,7.1,4.3,"Actin",lilac_pale,plum,2.75) +
  node(87,7.2,7.3,4.3,"Actin",lilac_pale,plum,2.75) +
  node(95,7.2,7.1,4.3,"Gate",green_pale,size=2.78) +
  word(84.5,17.5,"Weights held fixed",2.85,muted_ink)

# One schematic step response, with illustrative sample times rather than data.
reporter_response <- data.frame(x=42+seq(0,17,length.out=90))
reporter_response$y <- 6+8*(1-exp(-(reporter_response$x-42)/6))
sample_times <- data.frame(x=c(44,48,52,56,59))
sample_times$y <- 6+8*(1-exp(-(sample_times$x-42)/6))
a2 <- diagram() +
  field(2,33,green_tint) + field(36.5,63.5,green_tint) + field(67,99,lilac_tint) +
  word(4,22.3,"Calibrate reporter mixing",3.38,forest,TRUE,0) +
  word(38.2,22.3,"Filter and sample",3.38,forest,TRUE,0) +
  word(69,22.3,"Decide or abstain",3.38,plum,TRUE,0) +
  node(10,14.6,11.6,4.3,"Gate (q)",green_pale,size=2.78) +
  node(10,7.8,11.6,4.3,"Actin (a)",lilac_pale,plum,2.78) +
  link(16.1,14.6,24,11.6) +
  link(16.1,7.8,24,10.9,plum) +
  word(20,16.6,"1 − m",2.65,forest) +
  word(20,5.4,"m ≤ b",2.65,plum) +
  annotate("point",x=27.2,y=11.4,shape=21,fill="white",colour=forest,size=5.9,stroke=.45) +
  word(27.2,11.4,"+",3.4,forest,TRUE) +
  link(33,11.4,36.2,11.4) +
  link(41,6,60,6,muted_ink,FALSE) +
  link(42,6,42,16,muted_ink,FALSE) +
  geom_line(data=reporter_response,aes(x,y),colour=forest,linewidth=.56) +
  geom_segment(data=sample_times,aes(x=x,xend=x,y=6,yend=y),colour=green_ink,
               linewidth=.26,linetype=3) +
  geom_point(data=sample_times,aes(x,y),shape=21,fill=lilac_fill,colour=forest,
             size=1.9,stroke=.26) +
  word(50,17.5,"Response and acquisition",2.7,muted_ink) +
  link(63.5,11.4,66.7,11.4) +
  node(74,11.4,10,5,"B",lilac_pale,plum) +
  node(83,11.4,8,5,"Abstain",pink_tint,muted_ink,size=2.65) +
  node(92,11.4,10,5,"A",green_pale) +
  word(83,17.1,"Projected readout",2.9,muted_ink) +
  link(72,5.7,94,5.7,muted_ink)

# N-WASP stays upstream of the branched network. The dashed return describes
# the candidate feedback topology rather than a measured cellular pathway.
a3 <- diagram() +
  field(2,24,green_tint) + field(27,49,green_tint) +
  field(52,74,lilac_tint) + field(77,99,green_tint) +
  word(c(13,38,63,88),17.1,c("Local N-WASP","Accessible pool","Branched F-actin","Protrusion mechanics"),
       3.32,forest,TRUE) +
  annotate("segment",x=6,xend=20,y=7.5,yend=7.5,colour=green_ink,linewidth=.62) +
  annotate("point",x=c(8,11,16,19),y=c(10,12,10,12),shape=21,
           fill=green_fill,colour=forest,size=2.6,stroke=.33) +
  word(13,4.1,"Local amount, p",2.78,forest) +
  annotate("path",x=c(31,35,38,41,45),y=c(10,12.9,8.7,12.8,10.1),
           colour=forest,linewidth=.68) +
  annotate("point",x=38,y=8.7,shape=21,fill=lilac_fill,colour=plum,size=2.6,stroke=.35) +
  word(38,4.1,"Active pool, (p / V) × q",2.78,forest) +
  annotate("segment",x=c(55,61,63.5,66),xend=c(70,65,68.7,70.4),
           y=c(8,9.6,10.3,11),yend=c(12,13.6,13.9,14),
           colour=plum,linewidth=.76) +
  annotate("point",x=c(61,63.5,66),y=c(9.6,10.3,11),shape=21,
           fill=lilac_fill,colour=plum,size=1.8,stroke=.26) +
  word(63,4.1,"Arp2/3 · actin turnover",2.78,plum) +
  annotate("curve",x=81,xend=95,y=8,yend=8,curvature=-.52,
           colour=forest,linewidth=.67) +
  link(88,10.2,88,14.6,green_ink) +
  word(88,4.1,"Edge motion · load",2.78,forest) +
  link(c(24,49,74),10.6,c(26.7,51.7,76.7),10.6) +
  # Keep the label and return arrow on white space above the molecular fields.
  annotate("rect",xmin=2,xmax=99,ymin=19.55,ymax=25,fill="white",colour=NA) +
  annotate("curve",x=89,xend=39,y=19.7,yend=19.7,curvature=.045,linetype=2,
           colour=plum,linewidth=.43,
           arrow=grid::arrow(length=grid::unit(1.55,"mm"),type="closed")) +
  word(64,23.9,"Hypothesized mechanical feedback",3.0,plum,TRUE)

iso <- data.frame(p=seq(.25,1.4,length.out=201))
iso$q <- .25/iso$p
states <- data.frame(p=c(.5,1),q=c(.5,.25),
                     label=c("Low abundance","Low accessibility"))
p9 <- ggplot(iso,aes(p,q)) +
  geom_line(linewidth=trace_width,colour=forest) +
  geom_point(data=states,aes(fill=label,shape=label),colour=ink,stroke=mark_stroke,size=3) +
  scale_fill_manual(values=c("Low abundance"=green_fill,"Low accessibility"=lilac_fill)) +
  scale_shape_manual(values=c("Low abundance"=21,"Low accessibility"=22)) +
  annotate("text",x=.96,y=.7,label="Active pool = 0.25",family=figure_font,colour=ink,size=3.2) +
  scale_x_continuous(breaks=c(.25,.5,1,1.4)) +
  scale_y_continuous(breaks=seq(0,1,.25),limits=c(0,1.04)) +
  labs(x="Relative local N-WASP",y="Accessible fraction") + th() +
  theme(legend.key.width=grid::unit(4,"mm"))
optics <- load_json("optical_reference")
errors <- data.frame(stage=factor(c("Raw ratio","Optically corrected"),
                                levels=c("Raw ratio","Optically corrected")),
                    value=c(optics$median_absolute_relative_ratio_error_naive,
                            optics$median_absolute_relative_ratio_error_corrected))
p10 <- ggplot(errors,aes(value,stage)) +
  geom_segment(aes(x=0,xend=value,yend=stage),colour=muted_ink,linewidth=.50) +
  geom_point(aes(fill=stage),shape=21,colour=ink,size=3,stroke=mark_stroke) +
  geom_text(aes(label=scales::percent(value,accuracy=.01)),nudge_x=.009,hjust=0,colour=ink,
            family=figure_font,size=3.15) +
  scale_fill_manual(values=c(green_fill,lilac_fill),guide="none") +
  scale_x_continuous(breaks=seq(0,.25,.05),labels=scales::label_percent(),limits=c(0,.265),
                     expand=expansion(mult=0)) +
  labs(x="Median absolute relative optical ratio error",y=NULL) + th() +
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(14,8,8,13))

# Curve and marked points from the frozen equal-time feedback stability condition.
K <- data.frame(gain=seq(0,13,length.out=261))
stability <- function(gain) -1 + gain^(1/3)/2
examples <- data.frame(gain=c(3.6,10),rate=stability(c(3.6,10)), label=c("q* = 0.10","q* = 0.50"))
p11 <- ggplot(K,aes(gain,stability(gain))) +
  annotate("rect",xmin=8,xmax=13,ymin=0,ymax=.21,fill=lilac_tint) +
  geom_hline(yintercept=0,colour=muted_ink,linewidth=reference_width,linetype=2) +
  geom_vline(xintercept=8,colour=plum,linewidth=reference_width,linetype=3) +
  geom_line(colour=forest,linewidth=trace_width) +
  geom_point(data=examples,aes(y=rate,fill=label,shape=label),colour=ink,size=2.8,stroke=mark_stroke) +
  scale_fill_manual(values=c(green_fill,lilac_fill)) +
  scale_shape_manual(values=c(21,22)) +
  scale_x_continuous(breaks=c(0,4,8,12),expand=expansion(mult=c(.025,.02))) +
  coord_cartesian(ylim=c(-.57,.21)) + labs(x="Effective feedback gain, K",y="Leading real part") + th()

phase <- read.csv(file.path(data_dir, "S02_state.csv"))
chaos <- jsonlite::fromJSON(file.path(data_dir, "chaos_exponents.json"),
                           simplifyVector = FALSE)
exponents <- data.frame(setting=factor(c("r = 10","r = 28"),levels=c("r = 28","r = 10")),
                 exponent=c(chaos$cases$low_control$spectrum[[1]],chaos$cases$feedback$spectrum[[1]]))
phase_plot <- ggplot(phase,aes(q,m)) +
  geom_path(colour=forest,alpha=.55,linewidth=.31) +
  labs(x="Gate proxy",y="Feedback proxy",title="Autonomous state projection (r = 28)") +
  th()
exponent_plot <- ggplot(exponents,aes(exponent,setting,fill=setting)) +
  geom_segment(aes(x=0,xend=exponent,yend=setting),colour=muted_ink,linewidth=.48) +
  geom_vline(xintercept=0,colour=muted_ink,linewidth=reference_width,linetype=2) +
  geom_point(shape=21,size=3,colour=ink,stroke=mark_stroke) +
  geom_text(aes(label=sprintf("%.4f",exponent)),nudge_y=.18,family=figure_font,size=3.0,colour=ink) +
  scale_fill_manual(values=c("r = 10"=green_pale,"r = 28"=lilac_fill),guide="none") +
  coord_cartesian(xlim=c(-.8,1.18),ylim=c(.72,2.46)) +
  labs(x="Largest exponent / model time",y=NULL,title="Fixed point and chaotic benchmark") +
  th()
p12 <- aligned_row(list(phase_plot, exponent_plot))

figure1_data <- labelled_grid(list(p2, p3, p1, p6), ncol=2, labels=letters[2:5])
figure1 <- plot_grid(tag_panel(a1,"a"),figure1_data,ncol=1,rel_heights=c(.61,2))
figure2_data <- labelled_grid(list(p4, p5, p7, p8), ncol=2, labels=letters[2:5])
figure2 <- plot_grid(tag_panel(a2,"a"),figure2_data,tag_panel(p10,"f"),
                     ncol=1,rel_heights=c(.63,2,.56))
figure3_left <- align_plots(p9, phase_plot, align="v", axis="lr")
figure3_right <- align_plots(p11, exponent_plot, align="v", axis="lr")
figure3_data <- aligned_row(list(figure3_left[[1]], figure3_right[[1]]), c("b", "c"))
figure3_benchmark <- aligned_row(list(figure3_left[[2]], figure3_right[[2]]))
figure3 <- plot_grid(tag_panel(a3,"a"),figure3_data,tag_panel(figure3_benchmark,"d"),
                     ncol=1,rel_heights=c(.62,1,1))
export(figure1,"F01",height=162)
export(figure2,"F02",height=187)
export(figure3,"F03",height=167)

panels <- list(F01a=a1,F01b=p2,F01c=p3,F01d=p1,F01e=p6,
               F02a=a2,F02b=p4,F02c=p5,F02d=p7,F02e=p8,F02f=p10,
               F03a=a3,F03b=p9,F03c=p11,F03d=p12)
panel_dir <- file.path(output_dir, "panels")
dir.create(panel_dir,showWarnings=FALSE)
for (id in names(panels)) {
  schema <- id %in% c("F01a","F02a","F03a")
  wide <- schema || id %in% c("F02f","F03d")
  export(tag_panel(panels[[id]],substr(id,4,4)),id,
         width=if (wide) 180 else 90,height=if (schema) 38 else if (id=="F02f") 37 else 63,
         directory=panel_dir)
}
