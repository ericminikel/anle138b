options(stringsAsFactors=F)
library(tidyverse)
library(janitor)
library(openxlsx)
library(stringr)
library(survival)
setwd('~/d/sci/src/anle138b/')

#### RML STUDY

read.xlsx('received/promising/RMLinocGFAPlucDec2014.xlsx', sheet='Master') %>% 
  clean_names() %>%
  filter(harvest_date=='age to endpoint') %>%
  mutate(animal=as.integer(number),
         sex=toupper(substr(sex,1,1)),
         inoculation_date = as.Date(inoc_date, origin='1900-01-01')-2,
         diet = tolower(diet)) %>%
  select(animal, sex, inoculation_date, diet) -> rml_master

read_tsv('data/rml_surv.tsv')  -> rml_surv

rml_master %>%
  select(animal, inoculation_date) %>%
  inner_join(rml_surv,by='animal') %>%
  mutate(dod = as.character(as.Date(dod, format='%m/%d/%y'))) %>%
  mutate(dpi = as.numeric(as.Date(dod) - as.Date(inoculation_date)),
         status=1) %>%
  mutate(inoculum='RML') %>%
  select(animal, inoculum, sex, diet, death_dpi=dpi, disease_status=status, inoculation_date, death_date=dod, death_notes=cause) -> survdata

write_tsv(survdata, 'data/animals/rml_master.tsv')

all_files = list.files('received/drive-download-20220506T210132Z-001/', full.name=T)
potential_quant_files = all_files[grepl('(quantification [0-9]{6}\\.xlsx)|([0-9]{6} quan\\.xlsx)',all_files) & !grepl('~',all_files)]
quant_data = tibble(date=character(0), animal=integer(0), flux=numeric(0), adjusted=logical(0))
for (qfile in potential_quant_files) {
  date_text = str_extract(gsub('.*\\/','',qfile),'[0-9]{6}')
  if (date_text=='051105') { date_text = '051115' } # fix one typo
  date = as.Date(date_text, format='%m%d%y')
  mod = read.xlsx(qfile, sheet=1, colNames=F) # Sheet1 in the "quantification" files or ori in the "quan" files
  total_flux_col = grep('Total Flux',mod)
  final_flux_col = grep('final flux',mod)
  if (length(total_flux_col) > 0) {
    # if more than 1 option for data col, take the one with the most N of animal sessions
    n_matches = numeric(length(total_flux_col))
    for (j in 1:length(total_flux_col)) {
      n_matches[j] = length(grep('Total Flux',mod[,total_flux_col[j]]))
    }
    best_j = which.max(n_matches)
    flux_col = total_flux_col[best_j]
    flux_type = 'total'
  } else if (length(final_flux_col) > 0) {
    flux_col = final_flux_col
    flux_type = 'final'
  } else {
    next
  }
  rows_with_data = which(!is.na(suppressWarnings(as.numeric(mod[,flux_col]))))
  for (row in rows_with_data) {
    animal = as.integer(str_extract(mod[row,'X1'],'[0-9]{5}'))
    flux = as.numeric(mod[row,flux_col])
    quant_data = rbind(quant_data, cbind(date=as.character(date), animal=animal, flux=flux, adjusted=(flux_type=='final')))
  }
}
quant_data = quant_data[!is.na(quant_data$animal) & quant_data$animal %in% survdata$animal,]

survdata$animal = as.character(survdata$animal)
quant_data %>%
  inner_join(survdata, by='animal') %>%
  mutate(bli_dpi = as.numeric(as.Date(date) - as.Date(inoculation_date))) %>%
  select(animal, bli_dpi, flux) %>%
  arrange(bli_dpi, animal) -> rml_flux

write_tsv(rml_flux, 'data/animals/rml_flux.tsv')

rml_meta = tibble(diet=c('anle138b','placebo diet'),
                  color=c('#0001CD','#A9A9A9'))


write_tsv(rml_meta, 'data/animals/rml_meta.tsv')


### KI STUDY

meta = tibble(gt=rep(c('WT','FFI','CJD'),2), 
              tx=rep(c('placebo','anle138b'),each=3),
              color=alpha(rep(c('#000000','#D95F02','#7570B3'),2), rep(c(1,.35),each=3)))

write_tsv(meta, 'data/animals/ki_meta.tsv')

master = read_tsv('data_working/master.tsv')
master$days = as.numeric(as.Date(master$dod, '%m/%d/%y') - as.Date(master$dob, '%m/%d/%y'))

master %>%
  mutate(animal = as.character(id),
         gt = gsub(' .*','',group),
         tx = gsub('.+ ','',group),
         dob = as.Date(dob, format='%m/%d/%y'),
         dod = as.Date(dod, format='%m/%d/%y'),
         death_status = ifelse(death_cause=='planned harvest',0,1)) %>%
  mutate(tx = gsub('Anle','anle138b',tx)) %>%
  mutate(tx = gsub('placebo','placebo diet',tx)) %>%
  rename(death_comments=death_cause) %>%
  select(animal, gt, tx, group, dob, dod, days, death_status, death_comments) -> ki_master

write_tsv(ki_master, 'data/animals/ki_master.tsv')

flux_dates_raw = read.xlsx('received/promising/quantification flux_fixed.xlsx', sheet='final flux', startRow=1)[1,]
flux_raw = read.xlsx('received/promising/quantification flux_fixed.xlsx', sheet='final flux', startRow=3) %>% clean_names()

flux_dates = as.Date(str_extract(unlist(flux_dates_raw[1,11:33]),'[0-9]+\\/[0-9]+\\/[0-9]+'), format='%m/%d/%y')
flux_raw$dob = as.Date(flux_raw$dob, origin = '1900-01-01') - 2 # for some reason, subtract 2 days match up with how Excel renders it by default - 2014-06-27 and 2014-07-01 are hte first two unique values
flux_raw$diet_date = as.Date(flux_raw$diet_date, origin = '1900-01-01') - 2 # same
table(flux_raw[,c('prp','treatment')]) # 6 cohorts of 10 each
colnames(flux_raw)[11:33] = as.character(flux_dates)

flux_raw$gt = gsub('3F4 ','',gsub(' HOZ','',flux_raw$prp))
flux_raw$tx = gsub(' ','',tolower(flux_raw$treatment))

flux_raw %>%
  select(animal=number, sex,gt, tx, dob, diet_date, 11:33) %>%
  pivot_longer(cols=`2014-09-07`:`2016-02-26`) %>%
  rename(flux_date=name, flux_value=value) %>%
  filter(!is.na(flux_value) & flux_value > 1e4) %>%
  mutate(days = as.numeric(as.Date(flux_date) - as.Date(dob))) -> flux

all_files = list.files('received/drive-download-20220506T210132Z-001/', full.name=T)
potential_quant_files = all_files[grepl('051016|042016|042916',all_files) & !grepl('~',all_files)]
quant_data = tibble(date=character(0), animal=integer(0), flux=numeric(0), adjusted=logical(0))
for (qfile in potential_quant_files) {
  date_text = str_extract(gsub('.*\\/','',qfile),'[0-9]{6}')
  date = as.Date(date_text, format='%m%d%y')
  potential_data_sheets = getSheetNames(qfile)
  potential_data_sheets = potential_data_sheets[!(potential_data_sheets %in% c('graph','together'))]
  for (potential_data_sheet in potential_data_sheets) {
    xlsx_raw = read.xlsx(qfile, sheet=potential_data_sheet, colNames=F)
    if (grepl('#',xlsx_raw[2,2])) {
      animal = as.integer(str_extract(xlsx_raw[2,2],'[0-9]{5}'))
    } else if (grepl('#',xlsx_raw[1,1])) {
      animal = as.integer(str_extract(xlsx_raw[1,1],'[0-9]{5}'))
    } else if (grepl('#',xlsx_raw[3,2])) {
      animal = as.integer(str_extract(xlsx_raw[3,2],'[0-9]{5}'))
    } else {
      animal = as.integer(NA)
    }
    total_flux_col = grep('Total Flux',xlsx_raw)
    if (length(total_flux_col) > 0) {
      fluxes = xlsx_raw[3:20,total_flux_col[1]]
      max_flux = max(as.numeric(fluxes), na.rm=T)
      quant_data = rbind(quant_data,cbind(date=as.character(date),animal,flux=max_flux,adjusted=F))
    }
    # second_animal_startrow = grep('^2nd',xlsx_raw[,1])
    # if (length(second_animal_startrow) > 0) {
    #   # turns out it is just re-imaging same animal, no new data
    # }
  }
}
quant_data %>% filter(!is.na(animal)) %>% arrange(date, animal)

quant_data %>%
  mutate(animal = as.character(animal)) %>%
  inner_join(ki_master, by='animal') %>%
  mutate(bli_dpi = as.Date(date) - as.Date(dob)) %>%
  select(animal, bli_dpi, flux) -> flux_from_quant_data

flux %>%
  select(animal, bli_dpi=days, flux=flux_value) -> flux_from_flux_sheet

# mean(as.numeric(flux_from_quant_data$flux), na.rm=T) # 14449111
# mean(flux_from_flux_sheet$flux) # 504341
# so turns out something was *very* different about the quant_data readings. we'll exclude those then.

rbind(flux_from_flux_sheet) %>%
  arrange(bli_dpi, animal) -> flux_all

# check no dups - great!
# flux_all %>%
#   group_by(animal, bli_dpi) %>%
#   summarize(.groups='keep', n=n()) %>%
#   ungroup() %>%
#   arrange(desc(n))

write_tsv(flux_all, 'data/animals/ki_flux.tsv')

### HISTO

histo_header1 = read.xlsx('histo/FFICJDanalysis.xlsx',sheet='Sheet1', startRow=2, rows = 2) %>% clean_names()
histo = read.xlsx('histo/FFICJDanalysis.xlsx',sheet='Sheet1', startRow=3) %>% clean_names()
colnames(histo) = c('group','animal',paste(rep(colnames(histo_header1),each=2),colnames(histo)[3:ncol(histo)],sep='_'))
colnames(histo) = gsub('_[0-9]*$','',colnames(histo))

last_group = histo$group[1]
for (i in 1:nrow(histo)) {
  if (is.na(histo$group[i])) histo$group[i] = last_group
  else last_group = histo$group[i]
}
histo = histo[!is.na(suppressWarnings(as.integer(histo$animal))),]


histo %>%
  mutate(across(everything(), as.character)) %>%
  pivot_longer(cols=-c(group,animal)) %>%
  mutate(stain = ifelse(grepl('h_e',name),'h_e',ifelse(grepl('x3f4',name),'x3f4','')),
         region = gsub('_h_e|_x3f4','',name)) -> histo_long

# vacuole = factor(c('VVV','VV','V','(V)','(V?)','((V))','-'),ordered=T)
# puncta = factor(c('PPP','PP','P','(P)','((P))','diffuse'),ordered=T)

vacuole_scores = tibble(value=c('VVV','VV','V','(V)','(V?)','((V))','-'),
                        score=c(7:1))
puncta_scores = tibble(value=c('PP','P','(P)','((P))','diffuse'),
                       score=c(5:1))

histo_long %>%
  mutate(gt = gsub(' .*','',group),
         tx = gsub('.+ ','',group)) %>%
  mutate(value = case_when(grepl('diffuse',value) ~ 'diffuse',
                           value %in% c('V, patch') ~ 'V',
                           value %in% c('((V?))') ~ '((V))',
                           value %in% c('a few V?','(few V?)') ~ '(V?)',
                           value %in% c('(P), threads','(puncta)','3F4 positive puncta') ~ '(P)',
                           TRUE ~ value)) %>%
  mutate(valid = value %in% c(vacuole_scores$value, puncta_scores$value)) %>%
  select(gt, tx, animal, region, stain, value, valid) %>%
  filter(region != 'striatum') -> histo_data # only 1 score for striatum, impossible to do anything with it



histo_data %>%
  filter(stain=='h_e') %>%
  inner_join(vacuole_scores, by=c('value')) %>%
  select(animal, region, vacuole_score=score) -> vacuoles

histo_data %>%
  filter(stain=='x3f4') %>%
  inner_join(puncta_scores, by=c('value')) %>%
  select(animal, region, puncta_score=score) -> puncta

region_meta = tibble(region=unique(histo_data$region),
                     x=1:length(unique(histo_data$region)))


write_tsv(vacuoles, 'data/animals/ki_vacuoles.tsv')
write_tsv(puncta, 'data/animals/ki_puncta.tsv')
write_tsv(region_meta, 'data/animals/region_meta.tsv')


gt_meta = tibble(region=c('WT','CJD','FFI'),
                     x=1:length(unique(histo_data$gt)))

tx_meta = tibble(region=unique(histo_data$tx),
                 x=1:length(unique(histo_data$tx)))

summary(lm(vacuole_score ~ gt + tx + region, data=vacuoles))
summary(lm(puncta_score ~ gt + tx + region, data=puncta))

polr_model = polr(as.factor(vacuole_score) ~ gt + tx + region, data=vacuoles, Hess=T)
summary(polr_model)



vacuoles %>%
  group_by(gt, tx, region) %>%
  summarize(.groups='keep',
            median = median(vacuole_score),
            iqrl = quantile(vacuole_score, .25),
            iqru = quantile(vacuole_score, .75))
  

