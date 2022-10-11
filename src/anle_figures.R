
start_time = Sys.time()
cat(file=stderr(), 'Loading dependencies...'); flush.console()

options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(survival))
suppressMessages(library(MASS)); select=dplyr::select; summarize=dplyr::summarize
suppressMessages(library(magick))
suppressMessages(library(ordinal))

if(interactive()) {
  setwd('~/d/sci/src/anle138b/')
}

###############
#### FUNCTIONS & CONSTANTS
###############

cat(file=stderr(), 'done.\nLoading functions and constants...'); flush.console()


percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x < 0, '-', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot


vacuole_scores = tibble(value=c('VVV','VV','V','(V)','(V?)','((V))','-'),
                        score=c(7:1))
puncta_scores = tibble(value=c('PP','P','(P)','((P))','diffuse'),
                       score=c(5:1))


###############
#### OUTPUTS
###############


cat(file=stderr(), 'done.\nCreating output stream...'); flush.console()

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T



###############
#### Figure 1: RML
###############

cat(file=stderr(), 'done.\nCreating Figure 1...'); flush.console()

rml_meta = read_tsv('data/animals/rml_meta.tsv', col_types=cols())
rml_master = read_tsv('data/animals/rml_master.tsv', col_types=cols())
rml_flux = read_tsv('data/animals/rml_flux.tsv', col_types=cols())

rml_survfit = survfit(Surv(death_dpi, disease_status) ~ diet, data=rml_master)
rml_survdiff = survdiff(Surv(death_dpi, disease_status) ~ diet, data=rml_master)
rml_surv_p_value = 1-pchisq(rml_survdiff$chisq,df=1)
rml_survfit$diet = gsub('diet=','',names(rml_survfit$strata))
rml_survfit$color = rml_meta$color[match(rml_survfit$diet,rml_meta$diet)]


resx=300
png('display_items/figure-1.png',width=3.25*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(1:3, nrow=3, byrow=T)
layout(layout_matrix, heights=c(.6,1,1))


panel = 1



par(mar=c(0,1,2,0))
figure_1a = image_convert(image_read('data_working/anle138b_hires.png'),'png')
plot(as.raster(figure_1a))
mtext(LETTERS[panel], side=3, cex=2, adj = 0.05, line = -0.5)
panel = panel + 1

par(mar=c(3,4,2,1))
xlims = c(0,500)
ylims = c(0,1.05)
plot(rml_survfit, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, col=rml_survfit$color, lwd=2)
axis(side=1, at=0:5*100, labels=NA, cex.axis=0.8)
axis(side=1, at=0:5*100, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='dpi')
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, cex.axis=0.8, tck=-0.06)
mtext(side=2, line=2.5, text='survival')
par(xpd=T)
legend(x=325, y=1.05, rml_meta$diet, col=rml_meta$color, text.col=rml_meta$color, lwd=2, bty='n', cex=0.7)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

rml_master %>%
  group_by(diet) %>%
  summarize(.groups='keep',
            mean = mean(death_dpi),
            sd = sd(death_dpi)) %>%
  ungroup() %>%
  mutate(dpi_disp = paste0(formatC(mean,format='f',digits=0),'±',formatC(sd,format='f',digits=0))) -> rml_surv_smry

write(paste('RML survival anle138b vs. placebo: ',paste(rml_surv_smry$dpi_disp, collapse=' vs. '),
            ', an increase of ',percent(rml_surv_smry$mean[1]/rml_surv_smry$mean[2]-1),
            ', P = ',formatC(rml_surv_p_value,format='fg',digits=2),'\n',sep=''),text_stats_path,append=T)



rml_flux %>%
  mutate(approx_dpi = round(bli_dpi/7,0)*7) %>%
  inner_join(rml_master, by='animal') %>%
  inner_join(rml_meta, by='diet') %>%
  group_by(diet, color, approx_dpi) %>%
  summarize(.groups='keep',
            n=n(),
            mean = mean(flux),
            l95 = lower(flux),
            u95 = upper(flux)) %>%
  filter(n > 1) -> rml_flux_smry

xlims = c(0,500)
ylims = c(0, 1.2e7)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=0:5*100, labels=NA, cex.axis=0.8)
axis(side=1, at=0:5*100, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='dpi')
axis(side=2, at=0:3*5e6, labels=c('0','5M','10M','15M'), las=2, tck=-0.06, cex.axis=0.8)
axis(side=2, at=0:15*1e6, lwd=0, lwd.ticks=1, tck=-0.03, labels=NA)
mtext(side=2, line=3, text='flux (photons)')
for (coh in unique(rml_flux_smry$diet)) {
  subs = rml_flux_smry[rml_flux_smry$diet==coh,]
  points(x=subs$approx_dpi, y=subs$mean, col=subs$color, type='l', lwd=1.5)
  points(x=subs$approx_dpi, y=subs$mean, col=subs$color, pch=20, cex=0.6)
  polygon(x=c(subs$approx_dpi, rev(subs$approx_dpi)), y=c(subs$l95, rev(subs$u95)), border=NA, col=alpha(subs$color,ci_alpha))
}
rml_flux %>%
  inner_join(rml_master, by='animal') %>%
  inner_join(rml_meta, by=c('diet')) %>%
  arrange(animal, bli_dpi) -> rml_flux_indiv
for (this_animal in unique(rml_flux_indiv$animal)) {
  subs = subset(rml_flux_indiv, animal==this_animal)
  points(subs$bli_dpi, subs$flux, type='l', lwd=.25, col=subs$color)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)

rml_flux_indiv %>% group_by(animal, diet) %>% summarize(.groups='keep', maxflux = max(flux)) -> rml_max_bli
rml_max_bli %>% group_by(diet) %>% summarize(.groups='keep', mean=mean(maxflux)) -> rml_max_bli_smry
rml_bli_test_obj = ks.test(rml_max_bli$maxflux[rml_max_bli$diet=='placebo diet'], rml_max_bli$maxflux[rml_max_bli$diet=='anle138b'])

write(paste('RML max BLI flux anle138b vs. placebo: ',paste(rml_max_bli_smry$mean, collapse=' vs. '),
            ', a difference of ',percent(rml_max_bli_smry$mean[1]/rml_max_bli_smry$mean[2]-1),
            ', P = ',formatC(rml_bli_test_obj$p.value,format='fg',digits=2),'\n',sep=''),text_stats_path,append=T)


unnecessary_message = dev.off()

cat(file=stderr(), 'done.\nCreating Figure 2...'); flush.console()


###############
#### Figure 2: KI
###############


resx=300
png('display_items/figure-2.png', width=3.25*resx, height=6*resx, res=resx)


layout_matrix = matrix(c(1,1,1,
                         2,2,2,
                         3,4,5), nrow=3, byrow=T)
layout(layout_matrix, heights=c(1,1,2))
par(mar=c(4,4,3,2))
panel = 1

ki_meta = read_tsv('data/animals/ki_meta.tsv', col_types=cols())
ki_master = read_tsv('data/animals/ki_master.tsv', col_types=cols())
ki_flux = read_tsv('data/animals/ki_flux.tsv', col_types=cols())
ki_vacuoles = read_tsv('data/animals/ki_vacuoles.tsv', col_types=cols())
ki_puncta = read_tsv('data/animals/ki_puncta.tsv', col_types=cols())
region_meta = read_tsv('data/animals/region_meta.tsv', col_types=cols())

ki_diet_date = as.Date('2014-08-28')
range_diet_days = range(ki_diet_date - ki_master$dob)
mean_diet_days = mean(ki_diet_date - ki_master$dob)
write(paste('KI mice days from birth to drug diet: mean=',paste(round(as.numeric(mean_diet_days),1)),
            ', range= ',paste(range_diet_days, collapse=' to '),'\n',sep=''),text_stats_path,append=T)



# need to exclude any who died after 1/4 of planned harvests took place - otherwise makes it look as though a disproportionate number died
exclude_after_day = quantile(ki_master$days[ki_master$death_comments=='planned harvest'], .25)

ki_master$exclude_from_curves = ki_master$days > exclude_after_day & ki_master$death_status==1

ki_survfit  = survfit(Surv(days, death_status) ~ gt, data=subset(ki_master, tx=='placebo diet' & !exclude_from_curves))
ki_survdiff = survdiff(Surv(days, death_status) ~ gt, data=subset(ki_master, tx=='placebo diet' & !exclude_from_curves))
ki_surv_p_value = 1-pchisq(ki_survdiff$chisq,df=2)
ki_survfit$gt = gsub('gt=','',names(ki_survfit$strata))
ki_survfit$color = ki_meta$color[match(ki_survfit$gt,ki_meta$gt)]


xlims = c(0,800)
ylims = c(0,1.05)
plot(ki_survfit, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, col=ki_survfit$color, lwd=2)
axis(side=1, at=0:8*100, labels=NA, cex.axis=0.8)
axis(side=1, at=0:8*100, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days')
axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, cex.axis=0.8, tck=-0.06)
mtext(side=2, line=2.5, text='survival',cex=0.8)
par(xpd=T)
legend('bottomleft', ki_meta$gt_full[ki_meta$tx=='placebo diet'], col=ki_meta$color[ki_meta$tx=='placebo diet'], text.col=ki_meta$color[ki_meta$tx=='placebo diet'], lwd=2, bty='n', cex=0.7)
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

write(paste('KI survival WT vs. FFI vs. CJD: P = ',formatC(ki_surv_p_value,format='fg',digits=2),'\n',sep=''),text_stats_path,append=T)



ki_flux %>%
  mutate(approx_dpi = round(bli_dpi/30.44,0)*30.44) %>%
  inner_join(ki_master, by='animal') %>%
  inner_join(ki_meta, by=c('tx','gt')) %>%
  group_by(gt, tx, color, approx_dpi) %>%
  summarize(.groups='keep',
            n=n(),
            mean = mean(flux),
            l95 = lower(flux),
            u95 = upper(flux)) %>%
  filter(n > 1) -> ki_flux_smry

ki_flux_smry %>%
  filter(tx=='placebo diet') -> ki_ctrl_flux_smry

xlims = c(0,800)
ylims = c(0, 1.2e7)
plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
axis(side=1, at=0:8*100, labels=NA, cex.axis=0.8)
axis(side=1, at=0:8*100, lwd=0, line=-0.5)
mtext(side=1, line=1.5, text='days')
axis(side=2, at=0:3*5e6, labels=c('0','5M','10M','15M'), las=2, tck=-0.06, cex.axis=0.8)
axis(side=2, at=0:15*1e6, lwd=0, lwd.ticks=1, tck=-0.03, labels=NA)
# axis(side=2, at=0:3*5e5, labels=c('0','500K','1M','1.5M'), las=2, tck=-0.06, cex.axis=0.8)
# axis(side=2, at=0:15*1e5, lwd=0, lwd.ticks=1, tck=-0.03, labels=NA)
mtext(side=2, line=2.5, text='flux (photons)',cex=0.8)
for (coh in unique(ki_ctrl_flux_smry$gt)) {
  subs = ki_ctrl_flux_smry[ki_ctrl_flux_smry$gt==coh,]
  points(x=subs$approx_dpi, y=subs$mean, col=subs$color, type='l', lwd=1.5)
  points(x=subs$approx_dpi, y=subs$mean, col=subs$color, pch=20, cex=0.6)
  polygon(x=c(subs$approx_dpi, rev(subs$approx_dpi)), y=c(subs$l95, rev(subs$u95)), border=NA, col=alpha(subs$color,ci_alpha))
}
ki_flux %>%
  inner_join(ki_master, by='animal') %>%
  inner_join(ki_meta, by=c('gt','tx')) %>%
  arrange(animal, bli_dpi) -> ki_flux_indiv
for (this_animal in unique(ki_master$animal[ki_master$tx=='placebo diet'])) {
  subs = subset(ki_flux_indiv, animal==this_animal)
  points(subs$bli_dpi, subs$flux, type='l', lwd=.25, col=subs$color)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1


gt_offset = 0.14

ki_vacuoles %>%
  inner_join(ki_master, by='animal') %>%
  inner_join(ki_meta, by=c('tx','gt')) %>%
  inner_join(region_meta, by='region') %>%
  mutate(gt_fct = factor(gsub('WT','_WT',gt), ordered=F)) %>% # make WT the reference group in models
  mutate(tx_fct = factor(gsub('placebo diet','_placebo diet',tx), ordered=F)) %>% # make placebo the reference group in models
  mutate(reg_fct = factor(region, ordered=F)) %>%
  mutate(animal_fct = factor(animal, ordered=F)) %>%
  mutate(vacuole_fct = factor(vacuole_score, ordered=T)) %>% # make an ordinal variable for models
  select(gt, tx, animal, region, short, long, x, color, vacuole_score, vacuole_fct, reg_fct, gt_fct, tx_fct, animal_fct) -> all_vacuoles

all_vacuoles %>%
  filter(tx=='placebo diet') %>%
  mutate(x_offset=case_when(gt == 'WT' ~ -gt_offset,
                            gt == 'FFI' ~ 0,
                            gt == 'CJD' ~ gt_offset)) %>%
  mutate(y_final = max(x) - x + 1) -> ctrl_vacuoles

ctrl_vacuoles %>%
  group_by(gt, color, region, short, long, y_final, x_offset, vacuole_score) %>%
  summarize(.groups='keep', n=n()) -> ctrl_vac_histo

ctrl_vacuoles %>%
  group_by(gt, color, region, short, long, y_final) %>%
  summarize(.groups='keep',
            n = n(),
            med = median(vacuole_score),
            liqr = quantile(vacuole_score, .25),
            uiqr = quantile(vacuole_score, .75),
            mean = mean(vacuole_score),
            l95 = lower(vacuole_score),
            u95 = upper(vacuole_score)) %>%
  ungroup() -> ctrl_vac_smry

ctrl_vac_smry$p = as.numeric(NA)
for (i in 1:nrow(ctrl_vac_smry)) {
  if (ctrl_vac_smry$gt[i]=='WT') {
    ctrl_vac_smry$p[i] = 1
  } else {
    control_values = ctrl_vacuoles$vacuole_score[ctrl_vacuoles$region == ctrl_vac_smry$region[i] & ctrl_vacuoles$gt=='WT']
    test_values    = ctrl_vacuoles$vacuole_score[ctrl_vacuoles$region == ctrl_vac_smry$region[i] & ctrl_vacuoles$gt==ctrl_vac_smry$gt[i]]
    test_obj = suppressWarnings(wilcox.test(control_values, test_values))
    ctrl_vac_smry$p[i] = replace_na(test_obj$p.value, 1)
  }
}
ctrl_vac_smry$pfdr = as.numeric(NA)
ctrl_vac_smry$pfdr[ctrl_vac_smry$gt != 'WT'] = p.adjust(ctrl_vac_smry$p[ctrl_vac_smry$gt != 'WT'], method='fdr')
ctrl_vac_smry$psymb = case_when( ctrl_vac_smry$pfdr < 0.001 ~ '***',
                                 ctrl_vac_smry$pfdr < 0.01 ~ '**',
                                 ctrl_vac_smry$pfdr < 0.05 ~ '*',
                                 TRUE ~ '')

ctrl_vac_smry %>%
  group_by(region, short, long) %>%
  summarize(.groups='keep', y=mean(y_final)) %>% ungroup() -> y_meta

# legend panel
ylims = range(y_meta$y) + c(-0.1, 1.0)
xlims = c(0,1)
par(mar=c(4,0,3,0.25))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4, line=0, adj=1, at=y_meta$y+0.5, text=gsub('_','\n',y_meta$long), las=2, cex=0.6)

par(mar=c(4,0,3,1))
ylims = range(y_meta$y) + c(-0.1, 1.0)
xlims = range(all_vacuoles$vacuole_score) + c(-0.6, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
xbar = gt_offset/2
ybar = 0.85
rect(xleft=ctrl_vac_histo$vacuole_score-xbar+ctrl_vac_histo$x_offset, xright=ctrl_vac_histo$vacuole_score+xbar+ctrl_vac_histo$x_offset,
     ybottom=ctrl_vac_histo$y_final, ytop=ctrl_vac_histo$y_final + ybar*ctrl_vac_histo$n/8,
     col=ctrl_vac_histo$color, border=NA)
abline(h=unique(ctrl_vac_histo$y_final), lwd=.125)
par(xpd=T)
text(x=max(xlims)*.97,y=max(ylims)*1.01,labels='N mice',adj=0,pos=3,srt=270,font=3,cex=0.6)
arrow_y = 0.7
arrows(x0=2,x1=0.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
text(x=c(1.25,6.75),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(1,1),cex=0.8)
#text(x=c(2,6),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(4,2),cex=0.7)
arrows(x0=6,x1=7.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
par(xpd=F)
mtext(side=1, line=2, text='vacuoles', cex=0.8)
psymb_offset = 0.01
mtext(side=4, at=ctrl_vac_smry$y_final+psymb_offset, line=0.25, text=ctrl_vac_smry$psymb, col=ctrl_vac_smry$color, las=2)
for (i in 1:nrow(y_meta)) {
  axis(side=4, at=y_meta$y[i] + c(0:8)*ybar/8, tck=0.01, labels=NA, lwd=0.25)
  segments(x0=unique(all_vacuoles$vacuole_score)+0.5,x1=unique(all_vacuoles$vacuole_score)+0.5,y0=y_meta$y[i],y1=y_meta$y[i]+ybar/8*.3,lwd=0.125)
}
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1

vac_gt_model = suppressWarnings(clmm2(vacuole_fct ~ reg_fct + gt_fct, random=animal_fct, data=subset(all_vacuoles, tx=='placebo diet'), Hess=T))
vac_gt_model_cjd_est = suppressWarnings(summary(vac_gt_model)$coefficients['gt_fctCJD','Estimate'])
vac_gt_model_cjd_p = suppressWarnings(summary(vac_gt_model)$coefficients['gt_fctCJD','Pr(>|z|)'])
vac_gt_model_ffi_est = suppressWarnings(summary(vac_gt_model)$coefficients['gt_fctFFI','Estimate'])
vac_gt_model_ffi_p = suppressWarnings(summary(vac_gt_model)$coefficients['gt_fctFFI','Pr(>|z|)'])

write(paste('KI vacuole CLMM2 model: ',
            'WT vs. CJD P = ',formatC(vac_gt_model_cjd_p,format='e',digits=1), ', location coefficient = ',formatC(vac_gt_model_cjd_est,format='e',digits=1),
            ', WT vs. FFI P = ',formatC(vac_gt_model_ffi_p,format='e',digits=1), ', location coefficient = ',formatC(vac_gt_model_ffi_est,format='e',digits=1),
            '\n',sep=''),text_stats_path,append=T)


ki_puncta %>%
  inner_join(ki_master, by='animal') %>%
  inner_join(ki_meta, by=c('tx','gt')) %>%
  inner_join(region_meta, by='region') %>%
  mutate(gt_fct = factor(gsub('WT','_WT',gt), ordered=F)) %>% # make WT the reference group in models
  mutate(tx_fct = factor(gsub('placebo diet','_placebo diet',tx), ordered=F)) %>% # make placebo the reference group in models
  mutate(reg_fct = factor(region, ordered=F)) %>%
  mutate(animal_fct = factor(animal, ordered=F)) %>%
  mutate(puncta_fct = factor(puncta_score, ordered=T)) %>% # make an ordinal variable for models
  select(gt, tx, animal, region, x, color, puncta_score, puncta_fct, gt_fct, tx_fct, reg_fct, animal_fct) -> all_puncta

all_puncta %>%
  filter(tx=='placebo diet') %>%
  mutate(x_offset=case_when(gt == 'WT' ~ -gt_offset,
                            gt == 'FFI' ~ 0,
                            gt == 'CJD' ~ gt_offset)) %>%
  mutate(y_final = max(x) - x + 1) -> ctrl_puncta

ctrl_puncta %>%
  group_by(gt, color, region, y_final, x_offset, puncta_score) %>%
  summarize(.groups='keep', n=n()) -> ctrl_punc_histo

ctrl_puncta %>%
  group_by(gt, color, region, y_final) %>%
  summarize(.groups='keep',
            n = n(),
            med = median(puncta_score),
            liqr = quantile(puncta_score, .25),
            uiqr = quantile(puncta_score, .75),
            mean = mean(puncta_score),
            l95 = lower(puncta_score),
            u95 = upper(puncta_score)) %>%
  ungroup() -> ctrl_punc_smry

ctrl_punc_smry$p = as.numeric(NA)
for (i in 1:nrow(ctrl_punc_smry)) {
  if (ctrl_punc_smry$gt[i]=='WT') {
    ctrl_punc_smry$p[i] = 1
  } else {
    control_values = ctrl_puncta$puncta_score[ctrl_puncta$region == ctrl_punc_smry$region[i] & ctrl_puncta$gt=='WT']
    test_values    = ctrl_puncta$puncta_score[ctrl_puncta$region == ctrl_punc_smry$region[i] & ctrl_puncta$gt==ctrl_punc_smry$gt[i]]
    test_obj = suppressWarnings(wilcox.test(control_values, test_values))
    ctrl_punc_smry$p[i] = replace_na(test_obj$p.value, 1)
  }
}
ctrl_punc_smry$pfdr = as.numeric(NA)
ctrl_punc_smry$pfdr[ctrl_punc_smry$gt != 'WT'] = p.adjust(ctrl_punc_smry$p[ctrl_punc_smry$gt != 'WT'], method='fdr')
ctrl_punc_smry$psymb = case_when(ctrl_punc_smry$pfdr < 0.001 ~ '***',
                                 ctrl_punc_smry$pfdr < 0.01 ~ '**',
                                 ctrl_punc_smry$pfdr < 0.05 ~ '*',
                                 TRUE ~ '')

ctrl_punc_smry %>%
  group_by(region) %>%
  summarize(.groups='keep', y=mean(y_final)) -> y_meta

ctrl_punc_histo$x_offset = case_when(ctrl_punc_histo$gt=='WT' ~ -0.14,
                                     ctrl_punc_histo$gt=='FFI' ~ 0,
                                     ctrl_punc_histo$gt=='CJD' ~ 0.14)
ctrl_punc_histo$y_final = max(region_meta$x[match(ctrl_punc_histo$region, region_meta$region)]) - region_meta$x[match(ctrl_punc_histo$region, region_meta$region)] + 1

ylims = range(y_meta$y) + c(-0.1, 1.0)
xlims = range(all_puncta$puncta_score) + c(-0.5, 0.5)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
#mtext(side=2, line=0.25, at=y_meta$y, text=gsub('_',' ',y_meta$region), las=2)
xbar = 0.07
ybar = 0.85
rect(xleft=ctrl_punc_histo$puncta_score-xbar + ctrl_punc_histo$x_offset, xright=ctrl_punc_histo$puncta_score+xbar+ ctrl_punc_histo$x_offset,
     ybottom=ctrl_punc_histo$y_final, ytop=ctrl_punc_histo$y_final + ybar*ctrl_punc_histo$n/8,
     col=ctrl_punc_histo$color, border=NA)
abline(h=unique(ctrl_punc_histo$y_final), lwd=.125)
for (i in 1:nrow(y_meta)) {
  axis(side=4, at=y_meta$y[i] + c(0:8)*ybar/8, tck=0.01, labels=NA, lwd=0.25)
  segments(x0=unique(all_puncta$puncta_score)+0.5,x1=unique(all_puncta$puncta_score)+0.5,y0=y_meta$y[i],y1=y_meta$y[i]+ybar/8*.3,lwd=0.125)
}
#mtext(side=4, at=5, line=0.5, text='N mice', font=3, cex=0.6)
par(xpd=T)
text(x=max(xlims)*.97,y=max(ylims)*1.01,labels='N mice',adj=0,pos=3,srt=270,font=3,cex=0.6)
arrow_y = 0.7
arrows(x0=2.5,x1=0.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
text(x=c(1.5,4.5),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(1,1),cex=0.8)
#text(x=c(2,6),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(4,2),cex=0.7)
arrows(x0=3.5,x1=5.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
par(xpd=F)
mtext(side=1, line=2, text='PrP puncta', cex=0.8)
mtext(side=4, at=ctrl_punc_smry$y_final+psymb_offset, line=0.25, text=ctrl_punc_smry$psymb, col=ctrl_punc_smry$color, las=2)
mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
panel = panel + 1



punc_gt_model = suppressWarnings(clmm2(puncta_fct ~ reg_fct + gt_fct, random=animal_fct, data=subset(all_puncta, tx=='placebo diet'), Hess=T))
punc_gt_model_cjd_p = suppressWarnings(summary(punc_gt_model)$coefficients['gt_fctCJD','Pr(>|z|)'])
punc_gt_model_ffi_p = suppressWarnings(summary(punc_gt_model)$coefficients['gt_fctFFI','Pr(>|z|)'])

write(paste('KI puncta CLMM2 model: ',
            'WT vs. CJD P = ',formatC(punc_gt_model_cjd_p,format='e',digits=1),
            ', WT vs. FFI P = ',formatC(punc_gt_model_ffi_p,format='e',digits=1),
            '\n',sep=''),text_stats_path,append=T)


unnecessary_message = dev.off()

cat(file=stderr(), 'done.\nCreating Figure 4...'); flush.console()

###############
#### Figure 3 will be manually assembled histo
###############


###############
#### Figure 4: anle138b in KI
###############


resx=300
png('display_items/figure-4.png', width=6.5*resx, height=6*resx, res=resx)

layout_matrix = matrix(c(1,1,1,2,2,2,3,3,3,
                         4,4,4,5,5,5,6,6,6,
                         7,8,9,10,11,12,13,14,15), nrow=3, byrow=T)
layout(layout_matrix, heights=c(1,1,2))
par(mar=c(3,4,3,2))
panel = 1

tx_meta = tibble(tx=c('anle138b','placebo diet'),color=c('#37BC61','#A9A9A9'),lty=c(1,3))

for (this_gt in c('WT','FFI','CJD')) {
  ki_survfit = survfit(Surv(days, death_status) ~ tx, data=subset(ki_master, gt==this_gt & !exclude_from_curves))
  ki_survdiff = survdiff(Surv(days, death_status) ~ tx, data=subset(ki_master, gt==this_gt & !exclude_from_curves))
  
  ki_surv_p_value = 1-pchisq(ki_survdiff$chisq,df=1)
  ki_survfit$tx = gsub('tx=','',names(ki_survfit$strata))
  ki_survfit$lty = tx_meta$lty[match(ki_survfit$tx,tx_meta$tx)]
  ki_survfit$color = ki_meta$color[ki_meta$gt==this_gt][1] #tx_meta$color[match(ki_survfit$tx,tx_meta$tx)]
  
  xlims = c(0,800)
  ylims = c(0,1.05)
  plot(ki_survfit, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, col=ki_survfit$color, lty=ki_survfit$lty, lwd=2)
  par(new=T)
  plot(ki_survfit, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F, col=ki_survfit$color, lwd=.25)
  axis(side=1, at=0:8*100, labels=NA, cex.axis=0.8)
  axis(side=1, at=0:8*100, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='days')
  axis(side=2, at=0:4/4, labels=percent(0:4/4), las=2, cex.axis=0.8, tck=-0.06)
  mtext(side=2, line=2.5, text='survival',cex=0.8)
  par(xpd=T)
  #legend('bottomleft', tx_meta$tx, col=tx_meta$color, text.col=tx_meta$color, lwd=2, bty='n', cex=0.7)
  legend('bottomleft', tx_meta$tx, col=ki_survfit$color, text.col=ki_survfit$color, lty=ki_survfit$lty, lwd=2, bty='n', cex=0.7)
  legend('bottomleft', c('',''), col=ki_survfit$color, text.col=ki_survfit$color, lwd=.25, bty='n', cex=0.7)
  par(xpd=F)
  mtext(side=3, line=0, text=ki_meta$gt_full[ki_meta$gt==this_gt][1], col=ki_meta$color[ki_meta$gt==this_gt][1])
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  write(paste('KI survival anle138b vs. placebo in ',this_gt,': P = ',formatC(ki_surv_p_value,format='fg',digits=2),'\n',sep=''),text_stats_path,append=T)

}

ki_survdiff_overall = survdiff(Surv(days, death_status) ~ tx + gt, data=subset(ki_master, !exclude_from_curves))
ki_surv_overall_p_value = 1-pchisq(ki_survdiff_overall$chisq,df=5)
write(paste('KI survival anle138b vs. placebo overall: P = ',formatC(ki_surv_overall_p_value,format='fg',digits=2),'\n',sep=''),text_stats_path,append=T)

for (this_gt in c('WT','FFI','CJD')) {
  
  ki_flux %>%
    mutate(approx_dpi = round(bli_dpi/30.44,0)*30.44) %>%
    inner_join(ki_master, by='animal') %>%
    inner_join(ki_meta, by=c('tx','gt')) %>%
    filter(gt == this_gt) %>%
    group_by(gt, tx, color, approx_dpi) %>%
    summarize(.groups='keep',
              n=n(),
              mean = mean(flux),
              l95 = lower(flux),
              u95 = upper(flux)) %>%
    filter(n > 1) -> ki_flux_smry
  
  ki_flux_smry$color = tx_meta$color[match(ki_flux_smry$tx, tx_meta$tx)]
  
  xlims = c(0,800)
  ylims = c(0, 1.5e6)
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', axes=F, ann=F)
  axis(side=1, at=0:8*100, labels=NA, cex.axis=0.8)
  axis(side=1, at=0:8*100, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='days')
  axis(side=2, at=0:3*5e5, labels=NA, tck=-0.06)
  axis(side=2, at=0:3*5e5, labels=c('0','500K','1M','1.5M'), las=2, lwd=0, line=-0.5, cex.axis=0.8)
  axis(side=2, at=0:15*1e5, lwd=0, lwd.ticks=1, tck=-0.03, labels=NA)
  mtext(side=2, line=2.5, text='flux (photons)',cex=0.8)
  for (coh in unique(ki_flux_smry$tx)) {
    subs = ki_flux_smry[ki_flux_smry$tx==coh,]
    points(x=subs$approx_dpi, y=subs$mean, col=subs$color, type='l', lwd=1.5)
    points(x=subs$approx_dpi, y=subs$mean, col=subs$color, pch=20, cex=0.6)
    polygon(x=c(subs$approx_dpi, rev(subs$approx_dpi)), y=c(subs$l95, rev(subs$u95)), border=NA, col=alpha(subs$color,ci_alpha))
  }
  for (this_animal in unique(ki_master$animal[ki_master$gt==this_gt])) {
    subs = subset(ki_flux_indiv, animal==this_animal)
    subs$color = tx_meta$color[match(subs$tx, tx_meta$tx)]
    points(subs$bli_dpi, subs$flux, type='l', lwd=.25, col=subs$color)
  }
  
  # par(xpd=T)
  # legend(x=600,y=max(ylims)*1.2,tx_meta$tx,col=tx_meta$color,text.col=tx_meta$color,lwd=2,bty='n',cex=0.7)
  # par(xpd=F)
  par(xpd=T)
  legend(x=0,y=max(ylims)*1.1,tx_meta$tx,col=tx_meta$color,text.col=tx_meta$color,lwd=2,bg='white',box.lwd=0,cex=0.7)
  par(xpd=F)
  mtext(side=3, line=0, text=ki_meta$gt_full[ki_meta$gt==this_gt][1], col=ki_meta$color[ki_meta$gt==this_gt][1])
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  panel = panel + 1
}

for (this_gt in unique(ki_meta$gt)) {
  
  tx_offset = 0.14
  
  all_vacuoles %>%
    filter(gt==this_gt) %>%
    mutate(x_offset=case_when(tx=='placebo diet' ~ -tx_offset/2,
                              tx=='anle138b' ~ tx_offset/2)) %>%
    mutate(y_final = max(x) - x + 1) %>%
    mutate(color = case_when(tx=='placebo diet' ~ tx_meta$color[tx_meta$tx=='placebo diet'],
                             tx=='anle138b' ~ tx_meta$color[tx_meta$tx=='anle138b'])) -> this_vacuoles
  
  this_vacuoles %>%
    group_by(tx, color, region, short, long, y_final, x_offset, vacuole_score) %>%
    summarize(.groups='keep', n=n()) -> this_vac_histo
  
  this_vacuoles %>%
    group_by(tx, color, region, short, long, y_final) %>%
    summarize(.groups='keep',
              n = n(),
              med = median(vacuole_score),
              liqr = quantile(vacuole_score, .25),
              uiqr = quantile(vacuole_score, .75),
              mean = mean(vacuole_score),
              l95 = lower(vacuole_score),
              u95 = upper(vacuole_score)) %>%
    ungroup() -> this_vac_smry
  
  this_vac_smry$p = as.numeric(NA)
  for (i in 1:nrow(this_vac_smry)) {
    if (this_vac_smry$tx[i]=='placebo diet') {
      this_vac_smry$p[i] = 1
    } else {
      control_values = this_vacuoles$vacuole_score[this_vacuoles$region == this_vac_smry$region[i] & this_vacuoles$tx=='placebo diet']
      test_values    = this_vacuoles$vacuole_score[this_vacuoles$region == this_vac_smry$region[i] & this_vacuoles$tx==this_vac_smry$tx[i]]
      test_obj = suppressWarnings(wilcox.test(control_values, test_values))
      this_vac_smry$p[i] = replace_na(test_obj$p.value, 1)
    }
  }
  this_vac_smry$pfdr = as.numeric(NA)
  this_vac_smry$pfdr[this_vac_smry$tx != 'placebo diet'] = p.adjust(this_vac_smry$p[this_vac_smry$tx != 'placebo diet'], method='fdr')
  this_vac_smry$psymb = case_when( this_vac_smry$pfdr < 0.001 ~ '***',
                                   this_vac_smry$pfdr < 0.01 ~ '**',
                                   this_vac_smry$pfdr < 0.05 ~ '*',
                                   TRUE ~ '')
  
  this_vac_smry %>%
    group_by(region, short, long) %>%
    summarize(.groups='keep', y=mean(y_final)) %>% ungroup() -> y_meta
  
  # legend panel
  ylims = range(y_meta$y) + c(-0.1, 1.0)
  xlims = c(0,1)
  par(mar=c(4,0,3,0.25))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  mtext(side=4, line=0, adj=1, at=y_meta$y+0.5, text=gsub('_','\n',y_meta$long), las=2, cex=0.6)
  par(xpd=T)
  legend(x=0,y=max(ylims)*1.1,tx_meta$tx,col=tx_meta$color,text.col=tx_meta$color,lwd=2,bg='white',box.lwd=0,cex=0.7)
  par(xpd=F)
  
  par(mar=c(4,0,3,1))
  ylims = range(y_meta$y) + c(-0.1, 1.0)
  xlims = range(all_vacuoles$vacuole_score) + c(-0.65, 0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  xbar = gt_offset/2
  ybar = 0.85
  rect(xleft=this_vac_histo$vacuole_score-xbar+this_vac_histo$x_offset, xright=this_vac_histo$vacuole_score+xbar+this_vac_histo$x_offset,
       ybottom=this_vac_histo$y_final, ytop=this_vac_histo$y_final + ybar*this_vac_histo$n/8,
       col=this_vac_histo$color, border=NA)
  abline(h=unique(this_vac_histo$y_final), lwd=.125)
  par(xpd=T)
  text(x=max(xlims)*.97,y=max(ylims)*1.01,labels='N mice',adj=0,pos=3,srt=270,font=3,cex=0.6)
  arrow_y = 0.7
  arrows(x0=2.5,x1=0.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
  text(x=c(1.5,6.5),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(1,1),cex=0.8)
  #text(x=c(2,6),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(4,2),cex=0.7)
  arrows(x0=5.5,x1=7.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
  
  mtext(side=1, line=2, text='vacuoles', cex=0.8)
  psymb_offset = 0.01
  mtext(side=4, at=this_vac_smry$y_final+psymb_offset, line=0.25, text=this_vac_smry$psymb, col=this_vac_smry$color, las=2)
  for (i in 1:nrow(y_meta)) {
    axis(side=4, at=y_meta$y[i] + c(0:8)*ybar/8, tck=0.01, labels=NA, lwd=0.25)
    segments(x0=unique(all_vacuoles$vacuole_score)+0.5,x1=unique(all_vacuoles$vacuole_score)+0.5,y0=y_meta$y[i],y1=y_meta$y[i]+ybar/8*.3,lwd=0.125)
  }
  
  mtext(side=3, line=1.5, at=max(xlims)*.75, adj=0, text=ki_meta$gt_full[ki_meta$gt==this_gt][1], col=ki_meta$color[ki_meta$gt==this_gt][1])
  par(xpd=F)
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.1, line = 0.5)
  panel = panel + 1
  
  
  all_puncta %>%
    filter(gt==this_gt) %>%
    mutate(x_offset=case_when(tx=='placebo diet' ~ -tx_offset/2,
                              tx=='anle138b' ~ tx_offset/2)) %>%
    mutate(y_final = max(x) - x + 1) %>%
    mutate(color = case_when(tx=='placebo diet' ~ tx_meta$color[tx_meta$tx=='placebo diet'],
                             tx=='anle138b' ~ tx_meta$color[tx_meta$tx=='anle138b'])) -> this_puncta
  
  this_puncta %>%
    group_by(tx, color, region, y_final, x_offset, puncta_score) %>%
    summarize(.groups='keep', n=n()) -> this_punc_histo
  
  this_puncta %>%
    filter(region %in% region_meta$region) %>%
    group_by(tx, color, region, y_final) %>%
    summarize(.groups='keep',
              n = n(),
              med = median(puncta_score),
              liqr = quantile(puncta_score, .25),
              uiqr = quantile(puncta_score, .75),
              mean = mean(puncta_score),
              l95 = lower(puncta_score),
              u95 = upper(puncta_score)) %>%
    ungroup() -> this_punc_smry
  
  this_punc_smry$p = as.numeric(NA)
  for (i in 1:nrow(this_punc_smry)) {
    if (this_punc_smry$tx[i]=='placebo diet') {
      this_punc_smry$p[i] = 1
    } else {
      control_values = this_puncta$puncta_score[this_puncta$region == this_punc_smry$region[i] & this_puncta$tx=='placebo diet']
      test_values    = this_puncta$puncta_score[this_puncta$region == this_punc_smry$region[i] & this_puncta$tx==this_punc_smry$tx[i]]
      test_obj = suppressWarnings(wilcox.test(control_values, test_values))
      this_punc_smry$p[i] = replace_na(test_obj$p.value, 1)
    }
  }
  this_punc_smry$pfdr = as.numeric(NA)
  this_punc_smry$pfdr[this_punc_smry$tx != 'placebo diet'] = p.adjust(this_punc_smry$p[this_punc_smry$tx != 'placebo diet'], method='fdr')
  this_punc_smry$psymb = case_when(this_punc_smry$pfdr < 0.001 ~ '***',
                                   this_punc_smry$pfdr < 0.01 ~ '**',
                                   this_punc_smry$pfdr < 0.05 ~ '*',
                                   TRUE ~ '')
  
  this_punc_smry %>%
    group_by(region) %>%
    summarize(.groups='keep', y=mean(y_final)) -> y_meta
  
  ylims = range(y_meta$y) + c(-0.1, 1.0)
  xlims = range(all_puncta$puncta_score) + c(-0.6, 0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
  #mtext(side=2, line=0.25, at=y_meta$y, text=gsub('_',' ',y_meta$region), las=2)
  xbar = 0.07
  ybar = 0.85
  rect(xleft=this_punc_histo$puncta_score-xbar + this_punc_histo$x_offset, xright=this_punc_histo$puncta_score+xbar+ this_punc_histo$x_offset,
       ybottom=this_punc_histo$y_final, ytop=this_punc_histo$y_final + ybar*this_punc_histo$n/8,
       col=this_punc_histo$color, border=NA)
  abline(h=unique(this_punc_histo$y_final), lwd=.125)
  par(xpd=T)
  text(x=max(xlims)*.97,y=max(ylims)*1.01,labels='N mice',adj=0,pos=3,srt=270,font=3,cex=0.6)
  arrow_y = 0.7
  arrows(x0=2,x1=0.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
  text(x=c(1.25,4.75),y=c(arrow_y,arrow_y),labels=c('none','severe'),pos=c(1,1),cex=0.8)
  arrows(x0=4,x1=5.5,y0=arrow_y,y1=arrow_y,code=2,angle=30,length=0.03)
  par(xpd=F)
  mtext(side=1, line=2, text='PrP puncta', cex=0.8)
  for (i in 1:nrow(y_meta)) {
    axis(side=4, at=y_meta$y[i] + c(0:8)*ybar/8, tck=0.01, labels=NA, lwd=0.25)
    segments(x0=unique(all_puncta$puncta_score)+0.5,x1=unique(all_puncta$puncta_score)+0.5,y0=y_meta$y[i],y1=y_meta$y[i]+ybar/8*.3,lwd=0.125)
  }
  psymb_offset = 0.01
  mtext(side=4, at=this_punc_smry$y_final+psymb_offset, line=0.25, text=this_punc_smry$psymb, col=this_punc_smry$color, las=2)
  
}


tx_punc_model = suppressWarnings(clmm2(puncta_fct ~ reg_fct + tx_fct + gt_fct, random=animal_fct, data=subset(all_puncta, gt %in% c('FFI','CJD')), Hess=T))
tx_punc_est = suppressWarnings(summary(tx_punc_model)$coefficients['tx_fctanle138b','Estimate'])
tx_punc_p = suppressWarnings(summary(tx_punc_model)$coefficients['tx_fctanle138b','Pr(>|z|)'])

tx_vac_model = suppressWarnings(clmm2(vacuole_fct ~ reg_fct + tx_fct + gt_fct, random=animal_fct, data=subset(all_vacuoles, gt %in% c('FFI','CJD')), Hess=T))
tx_vac_est = suppressWarnings(summary(tx_vac_model)$coefficients['tx_fctanle138b','Estimate'])
tx_vac_p = suppressWarnings(summary(tx_vac_model)$coefficients['tx_fctanle138b','Pr(>|z|)'])

write(paste('KI CLMM2 model for FFI & CJD, anle138b vs. placebo: vacuoles P = ',formatC(tx_vac_p,format='e',digits=1),', location coefficient = ',formatC(tx_vac_est,format='e',digits=1),
            '; puncta P = ',formatC(tx_punc_p,format='e',digits=1),', location coefficient = ',formatC(tx_punc_est,format='e',digits=1),
            '\n',sep=''),text_stats_path,append=T)


unnecessary_message = dev.off()

cat(file=stderr(), paste0('done.\nAll tasks completed in ',round(as.numeric(Sys.time() - start_time),1),' seconds.\n')); flush.console()






