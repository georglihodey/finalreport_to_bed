##### Вообще основной код по трансформации c прекрасного канала https://www.youtube.com/@GenomicsBootCamp/videos
######## Illumina Final Reports to PLINK files https://www.youtube.com/watch?v=5_pNby7l2dE&t=1s
######## How to solve the SNP data merge error - TOP and FORWARD allele coding https://www.youtube.com/watch?v=RhHI0lwXDNI

library(data.table)
library(stringr)

# ==============================================================================
# КОНФИГУРАЦИЯ И ПУТИ
# ==============================================================================
DATA_PATH <- "Путь до ваших FinalReport'ов"
WORK_DIR  <- "РАБОЧАЯ ДИРЕКТОРИЯ"
setwd(WORK_DIR)

dir.create("./fin/rs_bed", recursive = TRUE, showWarnings = FALSE)
dir.create("./fin/comb", recursive = TRUE, showWarnings = FALSE)
dir.create("./fin/transversion", recursive = TRUE, showWarnings = FALSE)
dir.create("./arch_a/raw", recursive = TRUE, showWarnings = FALSE)

#### Предлагаемые шаблоны
con_plink_start <- "~/plink/plink --cow --nonfounders --allow-no-sex "
con_plink_end   <- " --make-bed --out "

# ==============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ==============================================================================
comp_map <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C", 
              "0" = "0", "I" = "D", "D" = "I", "PLUS" = "MINUS", "MINUS" = "PLUS")

vector_nucleotide_complement2 <- function(vec) {
  sapply(strsplit(vec, "/"), function(x) unname(comp_map[x]))
}

vector_nucleotide_complement3 <- function(vec) {
  sapply(strsplit(vec, "/"), function(x) paste(comp_map[x], collapse = "/"))
}

vector_nucleotide_complement4 <- function(vec) {
  sapply(strsplit(vec, "/"), function(x) paste(comp_map[x][c(2, 1)], collapse = "/"))
}

### Распознавание патерна CHR-POS-Allele1-Allele2-rsSNP
CHR_POS_Allele1_Allele2_rsSNP_dt <- function(bim) {
  chr_map <- c(as.character(0:31), "33", "X" = "30", "Y" = "31", "MT" = "33")
  names(chr_map)[1:33] <- c(as.character(0:31), "33")
  
  # Парсинг CHR-POS-...
  parts_dash <- tstrsplit(bim$V2, "-", keep = 1:2)
  chr_dash <- unname(chr_map[parts_dash[[1]]])
  pos_dash <- suppressWarnings(as.numeric(parts_dash[[2]]))
  valid_dash <- !is.na(chr_dash) & !is.na(pos_dash)
  if (any(valid_dash)) {
    bim[valid_dash, `:=`(V1 = as.numeric(chr_dash[valid_dash]), V4 = pos_dash[valid_dash])]
  }
  
  # Парсинг CHR:POS...
  cleaned_v2 <- gsub("^addit_|_ilmndup\\d+$", "", bim$V2)
  if (any(grepl(":", cleaned_v2))) {
    parts_colon <- tstrsplit(cleaned_v2, ":", keep = 1:2)
    chr_colon <- unname(chr_map[parts_colon[[1]]])
    pos_colon <- suppressWarnings(as.numeric(parts_colon[[2]]))
    valid_colon <- !is.na(chr_colon) & !is.na(pos_colon)
    if (any(valid_colon)) {
      bim[valid_colon, `:=`(V1 = as.numeric(chr_colon[valid_colon]), V4 = pos_colon[valid_colon])]
    }
  }
  return(bim)
}

forward_To_top_issue_dt <- function(bim, SNPchimp) {
  sc_sub <- SNPchimp[bim, on = .(SNP_name = V2), nomatch = NULL]
  sc_sub <- sc_sub[Alleles_A_B_FORWARD != Alleles_A_B_TOP]
  if (nrow(sc_sub) == 0) return(bim)
  
  sc_sub <- unique(sc_sub, by = "SNP_name")
  
  bim_sub <- bim[V2 %in% sc_sub$SNP_name]
  bim_sub[, comp_V6 := comp_map[V6]]
  bim_non_comp <- bim_sub[V5 != comp_V6]
  
  sc_non_comp <- sc_sub[SNP_name %in% bim_non_comp$V2]
  merged_nc <- bim_non_comp[sc_non_comp, on = .(V2 = SNP_name)]
  
  cond_forward <- merged_nc[(V2 == Al_AF | V5 == Al_BF | V6 == Al_AF | V6 == Al_BF), .N]
  cond_top     <- merged_nc[(V5 == Al_AT | V5 == Al_BT | V6 == Al_AT | V6 == Al_BT), .N]
  
  if (cond_forward < cond_top) {
    # Коррекция мономорфных/нулевых
    target_snps <- bim[V2 %in% sc_sub$SNP_name & V5 == "0", V2]
    if (length(target_snps) > 0) {
      sc_targets <- sc_sub[SNP_name %in% target_snps]
      bim[sc_targets, on = .(V2 = SNP_name), V5 := fifelse(V6 == i.Al_AT, i.Al_BT, fifelse(V6 == i.Al_BT, i.Al_AT, V5))]
    }
  } else {
    # Смена ориентации аллелей
    bim[sc_sub, on = .(V2 = SNP_name), `:=`(
      V6 = fifelse(V6 == i.Al_AF, i.Al_AT, fifelse(V6 == i.Al_BF, i.Al_BT, V6)),
      V5 = fifelse(V6 == i.Al_AF, i.Al_BT, fifelse(V6 == i.Al_BF, i.Al_AT, V5))
    )]
  }
  bim[, comp_V6 := NULL]
  return(bim)
}

### Отлавливаем и исправляем потенциальные ошибки
### в ходе работы финального цикла
PLINK_error_try_catch <- function(log, bed_file, bed_merged, dir_i, dir_o, sufx_f, file_o) {
  code_p <- '1'
  
  # Несовпадение позиций или хромосом
  has_mult_pos <- grepl('Multiple positions seen for variant', log)
  has_mult_chr <- grepl('Multiple chromosomes seen for variant', log)
  
  if (has_mult_pos || has_mult_chr) {
    pattern <- if (has_mult_pos) "Multiple positions seen for variant" else "Multiple chromosomes seen for variant"
    warning_lines <- grep(pattern, unlist(strsplit(log, "\n")), value = TRUE)
    snp_warnings <- unique(gsub("'", "", unlist(regmatches(warning_lines, gregexpr("'([^']+)'", warning_lines)))))
    
    snp_in_man  <- intersect(snp_warnings, manifest$Name)
    snp_not_man <- setdiff(snp_warnings, manifest$Name)
    
    for (bpath in c(bed_file, bed_merged)) {
      bim_dt <- fread(paste0(bpath, '.bim'), header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
      if (length(snp_in_man) > 0) {
        bim_dt[manifest, on = .(V2 = Name), `:=`(V1 = i.Chr, V4 = i.MapInfo)]
        fwrite(bim_dt, paste0(bpath, '.bim'), sep = '\t', col.names = FALSE)
      }
    }
    
    if (length(snp_not_man) > 0) {
      bim_f <- fread(paste0(bed_file, '.bim'), header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
      bim_m <- fread(paste0(bed_merged, '.bim'), header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
      
      # Отбираем валидные координаты из merged для обновления file
      m_valid <- bim_m[V2 %in% snp_not_man & V1 != "0" & V4 != 0, .(V2, new_V1 = V1, new_V4 = V4)]
      if (nrow(m_valid) > 0) {
        bim_f[m_valid, on = .(V2), `:=`(
          V1 = fifelse(V1 == "0" | V4 == 0, i.new_V1, V1),
          V4 = fifelse(V1 == "0" | V4 == 0, i.new_V4, V4)
        )]
        fwrite(bim_f, paste0(bed_file, '.bim'), sep = '\t', col.names = FALSE)
      }
      
      # Отбираем валидные координаты из file для обновления merged
      f_valid <- bim_f[V2 %in% snp_not_man & V1 != "0" & V4 != 0, .(V2, new_V1 = V1, new_V4 = V4)]
      if (nrow(f_valid) > 0) {
        bim_m[f_valid, on = .(V2), `:=`(
          V1 = fifelse(V1 == "0" | V4 == 0, i.new_V1, V1),
          V4 = fifelse(V1 == "0" | V4 == 0, i.new_V4, V4)
        )]
        fwrite(bim_m, paste0(bed_merged, '.bim'), sep = '\t', col.names = FALSE)
      }
    }
    code_p <- paste0(code_p, '_2')
  }
  
  # 3+ аллеля
  if (grepl("3\\+ alleles present", log)) {
    missnp_path <- paste0(dir_o, file_o, "-merge.missnp")
    if (grepl("_flip", log)) {
      if (file.exists(missnp_path)) {
        tri_all <- fread(missnp_path, header = FALSE)[[1]]
        tri_in_man <- intersect(tri_all, manifest$Name)
        tri_not_man <- setdiff(tri_all, manifest$Name)
        
        for (bpath in c(bed_file, bed_merged)) {
          bim_dt <- fread(paste0(bpath, '.bim'), header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
          sub_m <- manifest[Name %in% tri_in_man]
          
          # Проверка на совпадение аллелей
          bim_dt[sub_m, on = .(V2 = Name), is_match := {
            comp_snp <- vector_nucleotide_complement3(i.SNP)
            (paste0(V5, "/", V6) == i.SNP | paste0(V6, "/", V5) == i.SNP |
                paste0(V5, "/", V6) == comp_snp | paste0(V6, "/", V5) == comp_snp)
          }]
          bim_dt[V2 %in% tri_in_man & is_match == FALSE, V2 := paste0(V2, '_dup')]
          bim_dt[, is_match := NULL]
          fwrite(bim_dt, paste0(bpath, '.bim'), sep = '\t', col.names = FALSE)
        }
        
        if (length(tri_not_man) > 0) {
          fwrite(data.table(tri_not_man), missnp_path, col.names = FALSE)
          system(paste0(con_plink_start, " --bfile ", bed_file, " --exclude ", missnp_path, con_plink_end, bed_file))
        }
        code_p <- paste0(code_p, '_3')
      }
    } else {
      system(paste0(con_plink_start, " --bfile ", bed_file, " --flip ", missnp_path, con_plink_end, bed_file, sufx_f))
      code_p <- paste0(code_p, '_4')
    }
  }
  return(code_p)
}

### Определение типа чипа для исправления проблемы FORWARD и TOP
standardize_chip_strand_dt <- function(bim, snpchimp_dt) {
  # 1. Точное определение модели чипа по максимальному пересечению SNP_name
  chip_counts <- snpchimp_dt[SNP_name %in% bim$V2, .N, by = chip_name]
  if (nrow(chip_counts) == 0) {
    return(list(chip = "UNKNOWN", flip_list = character()))
  }
  
  detected_chip <- chip_counts[order(-N)][1, chip_name]
  sub_chimp <- snpchimp_dt[chip_name == detected_chip]
  
  # 2. Находим маркеры с расхождением спецификаций цепей
  diff_snps <- sub_chimp[Alleles_A_B_FORWARD != Alleles_A_B_TOP]
  merged_diag <- merge(bim, diff_snps, by.x = "V2", by.y = "SNP_name")
  
  if (nrow(merged_diag) == 0) {
    return(list(chip = detected_chip, flip_list = character()))
  }
  
  # Сравниваем наблюдаемые аллели в bim с эталонами FORWARD и TOP
  merged_diag[, `:=`(
    obs1 = paste0(V5, "/", V6),
    obs2 = paste0(V6, "/", V5)
  )]
  
  cnt_forward <- merged_diag[obs1 == Alleles_A_B_FORWARD | obs2 == Alleles_A_B_FORWARD, .N]
  cnt_top     <- merged_diag[obs1 == Alleles_A_B_TOP     | obs2 == Alleles_A_B_TOP, .N]
  
  # 3. Если большинство маркеров соответствует FORWARD, формируем список флипа для приведения к TOP
  snps_to_flip <- character()
  if (cnt_forward > cnt_top * 2) {
    snps_to_flip <- diff_snps[SNP_name %in% bim$V2, SNP_name]
  }
  
  return(list(chip = detected_chip, flip_list = snps_to_flip))
}

# ==============================================================================
# ЗАГРУЗКА БАЗ ДАННЫХ И МАНИФЕСТОВ
# ==============================================================================
manifest <- fread('~/Documents/cow/Manifests/unite_manifest_umd3.txt', sep = ',')
manifest <- manifest[, .(Name, SNP, Chr = as.character(Chr), MapInfo = as.numeric(MapInfo))]
manifest[Chr == "X", Chr := "30"]
manifest[Chr == "Y", Chr := "31"]
manifest[Chr == "MT", Chr := "33"]
manifest[, Chr := as.numeric(Chr)]
setkey(manifest, Name)

synonims_snp <- fread("~/Documents/cow/Manifests/synonymus_snp3.txt", sep = ',')

SNPchimp <- fread("~/Documents/cow/Manifests/SNPchimp/SNPchimp.tsv", sep = '\t')
SNPchimp[, `:=`(
  Al_AF = tstrsplit(Alleles_A_B_FORWARD, '/', keep = 1)[[1]],
  Al_BF = tstrsplit(Alleles_A_B_FORWARD, '/', keep = 2)[[1]],
  Al_AT = tstrsplit(Alleles_A_B_TOP, '/', keep = 1)[[1]],
  Al_BT = tstrsplit(Alleles_A_B_TOP, '/', keep = 2)[[1]]
)]
setkey(SNPchimp, SNP_name)


# ==============================================================================
# ШАГ 1: КОНВЕРТАЦИЯ FINALREPORT -> PLINK BINARY
# ==============================================================================
fr <- gsub('\\.txt$', '', list.files(DATA_PATH, pattern = '\\.txt$'))

for (FinalReport_File in fr) {
  file_path <- paste0(DATA_PATH, FinalReport_File, '.txt')
  
  # Читаем файл; если присутствует секция [Data], начинаем парсинг с неё
  FinalReport <- tryCatch(
    fread(file_path, sep = "\t", skip = "[Data]"),
    error = function(e) fread(file_path, sep = "\t")
  )
  
  # Нормализация имен колонок: удаляем пробелы, дефисы, точки, подчеркивания и переводим в нижний регистр
  clean_names <- tolower(gsub("[ _\\.-]", "", names(FinalReport)))
  setnames(FinalReport, clean_names)
  
  # 1. Приведение ключевых идентификаторов
  if ("snpname" %in% names(FinalReport)) {
    # Уже в порядке
  } else if ("snp" %in% names(FinalReport)) {
    setnames(FinalReport, "snp", "snpname")
  }
  
  if ("sampleid" %in% names(FinalReport)) {
    # Уже в порядке
  } else if ("sample" %in% names(FinalReport)) {
    setnames(FinalReport, "sample", "sampleid")
  }
  
  # 2. Выбор аллелей по приоритету: Top -> Forward -> AB
  if (all(c("allele1top", "allele2top") %in% names(FinalReport))) {
    setnames(FinalReport, c("allele1top", "allele2top"), c("Allele1", "Allele2"))
  } else if (all(c("allele1forward", "allele2forward") %in% names(FinalReport))) {
    setnames(FinalReport, c("allele1forward", "allele2forward"), c("Allele1", "Allele2"))
  } else if (all(c("allele1ab", "allele2ab") %in% names(FinalReport))) {
    setnames(FinalReport, c("allele1ab", "allele2ab"), c("Allele1", "Allele2"))
  } else {
    stop(paste("Не найдены подходящие пары аллелей в файле:", FinalReport_File))
  }
  
  # 3. Приведение координат Chr и Position
  if ("mapinfo" %in% names(FinalReport) && !"position" %in% names(FinalReport)) {
    setnames(FinalReport, "mapinfo", "position")
  }
  
  if (!all(c("chr", "position") %in% names(FinalReport))) {
    FinalReport[manifest, on = .(snpname = Name), `:=`(chr = i.Chr, position = i.MapInfo)]
  }
  
  FinalReport[is.na(chr), chr := 0]
  FinalReport[is.na(position), position := 0]
  
  # 4. Фильтрация по GC Score (если колонка присутствует)
  gc_col <- intersect(c("gcscore", "gc"), names(FinalReport))
  if (length(gc_col) > 0) {
    FinalReport[, (gc_col[1]) := suppressWarnings(as.numeric(get(gc_col[1])))]
    if (sum(!is.na(FinalReport[[gc_col[1]]])) > 0 && max(FinalReport[[gc_col[1]]], na.rm = TRUE) > 0) {
      FinalReport <- FinalReport[get(gc_col[1]) > 0.30]
    }
  }
  
  # Замена пропусков на 0
  FinalReport[Allele1 %in% c('-', '', 'NA') | is.na(Allele1), Allele1 := '0']
  FinalReport[Allele2 %in% c('-', '', 'NA') | is.na(Allele2), Allele2 := '0']
  
  # 5. Экспорт в PLINK LGEN/MAP/FAM
  # FAM
  fam_dt <- unique(FinalReport[, .(sampleid)])
  fam_dt[, `:=`(FID = FinalReport_File, IID = sampleid, sire = 0, dam = 0, sex = 0, phenotype = -9)]
  fwrite(fam_dt[, .(FID, IID, sire, dam, sex, phenotype)], 
         paste0('./', FinalReport_File, ".fam"), sep = '\t', col.names = FALSE)
  
  # LGEN
  lgen_dt <- FinalReport[, .(FID = FinalReport_File, sampleid, snpname, Allele1, Allele2)]
  fwrite(lgen_dt, paste0('./', FinalReport_File, ".lgen"), sep = '\t', col.names = FALSE)
  
  # MAP
  map_dt <- unique(FinalReport, by = "snpname")[, .(chr, snpname, morgan = 0, position)]
  fwrite(map_dt, paste0('./', FinalReport_File, ".map"), sep = '\t', col.names = FALSE)
  
  # Сборка бинарного набора
  system(paste0(con_plink_start, " --lfile ./", FinalReport_File, 
                " --make-bed --out ./arch_a/raw/", FinalReport_File))
  file.remove(paste0('./', FinalReport_File, c(".lgen", ".map", ".fam")))
}

# Копирование исходников в fin/
for (bed in fr) {
  file.copy(paste0("./", bed, ".bed"), paste0("./fin/", bed, "_fin.bed"), overwrite = TRUE)
  file.copy(paste0("./", bed, ".fam"), paste0("./fin/", bed, "_fin.fam"), overwrite = TRUE)
  file.copy(paste0("./", bed, ".bim"), paste0("./fin/", bed, "_fin.bim"), overwrite = TRUE)
}

# ==============================================================================
# ШАГ 2: ВОССТАНОВЛЕНИЕ КООРДИНАТ ЧЕРЕЗ BCFTOOLS
# ==============================================================================

for (bed in fr) {
  bim_path <- paste0('./fin/', bed, '_fin.bim')
  bim <- fread(bim_path, header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
  
  # Извлекаем rsID строго регулярным выражением
  bim[, rsID := str_extract(V2, "rs\\d+")]
  
  # Находим маркеры, у которых хромосома или позиция не определены (0)
  rs_zero <- na.omit(unique(bim[V1 == "0" | V4 == 0, rsID]))
  
  if (length(rs_zero) > 0) {
    rs_query_file <- paste0("~/Documents/cow/rs_database/rsIDs_", bed, ".txt")
    vcf_out_file  <- paste0("~/Documents/cow/rs_database/rsID_subset_", bed, ".vcf")
    
    fwrite(data.table(rs_zero), rs_query_file, col.names = FALSE)
    
    # Извлечение координат из эталонного VCF через bcftools
    system(paste0(
      "~/miniconda3/envs/samtools/bin/bcftools view -i 'ID=@", rs_query_file, "' ",
      "~/Documents/cow/rs_database/9913_GCA_000003055.5_current_ids.vcf.gz -o ", vcf_out_file
    ))
    
    vcf <- tryCatch(
      fread(vcf_out_file, skip = 45, select = 1:3, col.names = c('CHROM', 'POS', 'ID'),
            colClasses = list(character = c(1, 3), numeric = 2)),
      error = function(e) NULL
    )
    
    if (!is.null(vcf) && nrow(vcf) > 0) {
      vcf <- unique(vcf, by = "ID")
      # Обновляем V1 и V4 только там, где они были нулевыми
      bim[vcf, on = .(rsID = ID), `:=`(
        V1 = fifelse(V1 == "0", as.character(i.CHROM), V1),
        V4 = fifelse(V4 == 0, i.POS, V4)
      )]
    }
    
    unlink(c(rs_query_file, vcf_out_file))
  }
  
  bim[, rsID := NULL]
  
  # Сохраняем с разделителем таб
  fwrite(bim, bim_path, sep = '\t', col.names = FALSE)
}

# ==============================================================================
# ШАГ 3: ГАРМОНИЗАЦИЯ МАРКЕРОВ, ПОЗИЦИЙ И СИНОНИМОВ
# ==============================================================================
all_synonyms <- unique(na.omit(unlist(synonims_snp)))

for (bed in fr) {
  bim <- fread(paste0('./fin/', bed, '_fin.bim'), header = FALSE, colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
  bim <- CHR_POS_Allele1_Allele2_rsSNP_dt(bim)
  strand_check <- standardize_chip_strand_dt(bim, SNPchimp)
  message(sprintf("Датасет %s: определен чип %s. Маркеров к переводу цепи: %d", 
                  bed, strand_check$chip, length(strand_check$flip_list)))
  
  # Сопоставление безымянных/несоответствующих маркеров по ChrPos
  bim[, ChrPos := paste0(V1, '_', V4)]
  target_snps <- bim[V1 != "0" & V4 != 0 & !V2 %in% manifest$Name & !V2 %in% all_synonyms]
  
  if (nrow(target_snps) > 0) {
    m_valid <- unique(manifest[Chr != 0 & MapInfo != 0, .(Name, SNP, ChrPos = paste0(Chr, '_', MapInfo))], by = "ChrPos")
    matched <- m_valid[target_snps, on = .(ChrPos), nomatch = NULL]
    
    if (nrow(matched) > 0) {
      matched[, `:=`(
        al1 = tstrsplit(SNP, '/', keep = 1)[[1]],
        al2 = tstrsplit(SNP, '/', keep = 2)[[1]]
      )]
      matched[, comp_al1 := comp_map[al1]]
      matched[, comp_al2 := comp_map[al2]]
      
      # Векторизованный перенос имен
      valid_match <- matched[
        (V5 == al1 | V5 == al2 | V5 == comp_al1 | V5 == comp_al2 |
           V6 == al1 | V6 == al2 | V6 == comp_al1 | V6 == comp_al2)
      ]
      bim[valid_match, on = .(V2), V2 := i.Name]
    }
  }
  
  # Синонимы
  for (col in names(synonims_snp)[-1]) {
    match_dt <- synonims_snp[!is.na(get(col)), .(Target = Name_1, Current = get(col))]
    bim[match_dt, on = .(V2 = Current), V2 := i.Target]
  }
  
  # Исправление нулевых аллелей (мономорфные маркеры)
  zero_al <- bim[V5 == "0" & V2 %in% manifest$Name]
  if (nrow(zero_al) > 0) {
    m_sub <- manifest[zero_al$V2]
    m_sub[, `:=`(
      al1 = tstrsplit(SNP, '/', keep = 1)[[1]],
      al2 = tstrsplit(SNP, '/', keep = 2)[[1]]
    )]
    m_sub[, comp_al1 := comp_map[al1]]
    m_sub[, comp_al2 := comp_map[al2]]
    
    bim[m_sub, on = .(V2 = Name), `:=`(
      V6 = fifelse(V6 == i.al1 | V6 == i.comp_al1, i.al1, fifelse(V6 == i.al2 | V6 == i.comp_al2, i.al2, V6)),
      V5 = fifelse(V6 == i.al1 | V6 == i.comp_al1, i.al2, fifelse(V6 == i.al2 | V6 == i.comp_al2, i.al1, V5))
    )]
    
    # Сохраняем маркеры, оставшиеся неразрешенными и трансверсии
    unresolved <- bim[V5 == "0" | V2 %in% m_sub[comp_al1 == al2,Name]]
    if (nrow(unresolved) > 0) {
      fwrite(unresolved, paste0('./fin/transversion/', bed, '_fin_rs.txt'), sep = '\t', col.names = FALSE)
    }
  }
  
  # Исправление Chr / Pos из манифеста
  bim[manifest, on = .(V2 = Name), `:=`(
    V1 = fifelse(V4 == 0 | ChrPos != paste0(i.Chr, '_', i.MapInfo), as.character(i.Chr), V1),
    V4 = fifelse(V4 == 0 | ChrPos != paste0(i.Chr, '_', i.MapInfo), i.MapInfo, V4)
  )]
  
  # Forward to Top
  bim[, ChrPos := NULL]
  bim <- forward_To_top_issue_dt(bim, SNPchimp)
  # --------------------------------------------------------------------------
  # Обработка дубликатов ID маркеров перед сборкой PLINK
  # --------------------------------------------------------------------------
  dup_ids <- bim$V2[duplicated(bim$V2)]
  
  if (length(dup_ids) > 0) {
    # Если остались маркеры с одинаковым именем, но разными координатами/аллелями,
    # переименовываем повторные вхождения, добавляя суффикс _dup
    if (anyDuplicated(bim$V2)) {
      bim[, dup_idx := seq_len(.N), by = V2]
      bim[dup_idx > 1, V2 := paste0(V2, "_dup", dup_idx - 1)]
      bim[, dup_idx := NULL]
    }
  }
  
  # Проверка ориентации цепей по SNPchimp
  strand_check <- standardize_chip_strand_dt(bim, SNPchimp)
  message(sprintf("Датасет %s: определен чип %s. Маркеров к переводу цепи: %d", 
                  bed, strand_check$chip, length(strand_check$flip_list)))
  
  
  # Финальная запись набора файлов
  out_prefix <- paste0('./fin/rs_bed/', bed, '_fin_rs')
  fwrite(bim, paste0(out_prefix, '.bim'), sep = '\t', col.names = FALSE)
  
  file.copy(paste0("./fin/", bed, "_fin.bed"), paste0(out_prefix, ".bed"), overwrite = TRUE)
  file.copy(paste0("./fin/", bed, "_fin.fam"), paste0(out_prefix, ".fam"), overwrite = TRUE)
  
  if (length(strand_check$flip_list) > 0) {
    flip_file <- paste0('./fin/rs_bed/', bed, '_chimp_flip.txt')
    fwrite(data.table(strand_check$flip_list), flip_file, col.names = FALSE)
    
    system(paste0(con_plink_start, " --bfile ", out_prefix, 
                  " --flip ", flip_file, con_plink_end, out_prefix))
    system(paste0("rm ./fin/rs_bed/*~ "))
    
    unlink(flip_file)
  } else {
    system(paste0(con_plink_start, " --bfile ", out_prefix, con_plink_end, out_prefix))
    system(paste0("rm ./fin/rs_bed/*~ "))
  }
}

# ==============================================================================
# ШАГ 4: РАЗРЕШЕНИЕ КОМПЛЕМЕНТАРНЫХ ЗАМЕН (A/T И G/C) ЧЕРЕЗ MAF И ЛОКАЛЬНЫЙ LD
# ==============================================================================

#### Данная часть самая спорная.
#### В случае объединения разных датосетов мы рекомендуем удалять палиндромные ОНП


# 1. Поиск всех палиндромных маркеров в манифесте
palindromic_snps <- manifest[SNP %in% c("A/T", "T/A", "C/G", "G/C"), Name]

# 2. Определение эталонного датасета (с максимальным количеством ОНП)
sample_counts <- sapply(fr, function(x) {
  fam_file <- paste0("./fin/rs_bed/", x, "_fin_rs.bim")
  if (file.exists(fam_file)) length(readLines(fam_file)) else 0
})

ref_dataset <- names(which.max(sample_counts))
ref_prefix  <- paste0("./fin/rs_bed/", ref_dataset, "_fin_rs")

# Расчет частот эталона
system(paste0(con_plink_start, " --bfile ", ref_prefix, " --freq --out ", ref_prefix))
ref_frq <- fread(paste0(ref_prefix, ".frq"))
ref_bim <- fread(paste0(ref_prefix, ".bim"), header = FALSE, 
                 colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
setnames(ref_bim, c("CHR", "SNP", "CM", "POS", "A1", "A2"))

# Сохраняем список асимметричных маркеров эталона с известной ориентацией
ref_asym_snps <- manifest[!SNP %in% c("A/T", "T/A", "C/G", "G/C") & Chr != 0 & MapInfo != 0, Name]
ref_asym_snps <- intersect(ref_asym_snps, ref_bim$SNP)

ambiguous_log_list <- list()

MAF_threshold <- 0.1

# 3. Обработка датасетов
for (ds in fr) {
  if (ds == ref_dataset) next
  
  target_prefix <- paste0("./fin/rs_bed/", ds, "_fin_rs")
  system(paste0(con_plink_start, " --bfile ", target_prefix, " --freq --out ", target_prefix))
  
  tgt_frq <- fread(paste0(target_prefix, ".frq"))
  tgt_bim <- fread(paste0(target_prefix, ".bim"), header = FALSE,
                   colClasses = list(character = c(1, 2, 5, 6), numeric = c(3, 4)))
  setnames(tgt_bim, c("CHR", "SNP", "CM", "POS", "A1", "A2"))
  
  common_palindromes <- intersect(tgt_bim$SNP, palindromic_snps)
  common_palindromes <- intersect(common_palindromes, ref_bim$SNP)
  if (length(common_palindromes) == 0) next
  
  # Сравнение частот аллелей
  comp_table <- merge(
    tgt_frq[SNP %in% common_palindromes, .(SNP, A1_tgt = A1, A2_tgt = A2, MAF_tgt = MAF)],
    ref_frq[SNP %in% common_palindromes, .(SNP, A1_ref = A1, A2_ref = A2, MAF_ref = MAF)],
    by = "SNP"
  )
  # Нормализация частот к единому аллелю (A1_ref или A2_ref, если A1_ref == 0)
  comp_table[, `:=`(
    # Если в референсе A1 == 0 (мономорфен), базовым берем A2_ref
    base_ref_allele = fifelse(A1_ref == "0", A2_ref, A1_ref),
    freq_ref_base   = fifelse(A1_ref == "0", 1 - MAF_ref, MAF_ref)
  )]
  
  # Вычисляем частоту того же базового аллеля в тестовом датасете
  comp_table[, freq_tgt_base := fifelse(
    A1_tgt == base_ref_allele, MAF_tgt,
    fifelse(A2_tgt == base_ref_allele, 1 - MAF_tgt, NA_real_)
  )]
  
  # Если аллель вообще не найден в целевом файле (например, в обоих A1 == 0)
  comp_table[is.na(freq_tgt_base) & A1_tgt == "0" & A2_tgt == base_ref_allele, freq_tgt_base := 1.0]
  comp_table[is.na(freq_tgt_base) & A1_tgt == "0" & A2_tgt != base_ref_allele, freq_tgt_base := 0.0]
  comp_table[, freq_diff := abs(freq_tgt_base - freq_ref_base)]
  
  to_flip   <- character()
  to_remove <- character()
  
  # Вместо вектора to_remove создаем таблицу для сбора сомнительных сайтов текущего датасета
  to_log_dt <- data.table(
    DataSet = character(),
    SNP = character(),
    freq_ref_base = numeric(),
    freq_tgt_base = numeric(),
    freq_diff = numeric(),
    Reason = character()
  )
  
  # Зона: Неоднозначные маркеры, требующие LD-валидации
  ambiguous_snps <- comp_table[freq_diff > MAF_threshold, SNP]
  
  if (length(ambiguous_snps) > 0) {
    ld_snps_to_flip <- character()
    
    for (focal_snp in ambiguous_snps) {
      focal_info <- ref_bim[SNP == focal_snp]
      chr_target <- focal_info$CHR
      pos_target <- focal_info$POS
      cur_row    <- comp_table[SNP == focal_snp]
      cur_diff   <- cur_row$freq_diff
      
      flanking_snps <- ref_bim[
        CHR == chr_target & 
          abs(POS - pos_target) <= 500000 & 
          SNP != focal_snp & 
          SNP %in% ref_asym_snps & 
          SNP %in% tgt_bim$SNP, 
        SNP
      ]
      
      if (length(flanking_snps) < 3) {
        if (cur_diff > MAF_threshold) {
          to_log_dt <- rbind(to_log_dt, data.table(
            DataSet = ds, SNP = focal_snp, 
            freq_ref_base = cur_row$freq_ref_base, 
            freq_tgt_base = cur_row$freq_tgt_base, 
            freq_diff = cur_diff, 
            Reason = "Low_flanking_markers"
          ))
        }
        next
      }
      
      # Выгрузка генотипов окна и расчет корреляций
      snps_subset <- c(focal_snp, flanking_snps)
      fwrite(data.table(snps_subset), "./fin/tmp_snplist.txt", col.names = FALSE)
      
      system(paste0(con_plink_start, " --silent --bfile ", ref_prefix, 
                    " --extract ./fin/tmp_snplist.txt --recode A --out ./fin/tmp_ref"))
      system(paste0(con_plink_start, " --silent --bfile ", target_prefix, 
                    " --extract ./fin/tmp_snplist.txt --recode A --out ./fin/tmp_tgt"))
      
      raw_ref <- fread("./fin/tmp_ref.raw")
      raw_tgt <- fread("./fin/tmp_tgt.raw")
      
      cols_ref <- names(raw_ref)[-(1:6)]
      cols_tgt <- names(raw_tgt)[-(1:6)]
      
      focal_col_ref <- grep(paste0("^", focal_snp, "_"), cols_ref, value = TRUE)
      focal_col_tgt <- grep(paste0("^", focal_snp, "_"), cols_tgt, value = TRUE)
      
      if (length(focal_col_ref) == 0 || length(focal_col_tgt) == 0) {
        to_log_dt <- rbind(to_log_dt, data.table(
          DataSet = ds, SNP = focal_snp, 
          freq_ref_base = cur_row$freq_ref_base, 
          freq_tgt_base = cur_row$freq_tgt_base, 
          freq_diff = cur_diff, 
          Reason = "Focal_SNP_missing_in_recode"
        ))
        next
      }
      
      focal_allele_ref <- sub(paste0("^", focal_snp, "_"), "", focal_col_ref)
      focal_allele_tgt <- sub(paste0("^", focal_snp, "_"), "", focal_col_tgt)
      
      vec_focal_ref <- raw_ref[[focal_col_ref]]
      vec_focal_tgt <- raw_tgt[[focal_col_tgt]]
      
      if (focal_allele_tgt != focal_allele_ref) {
        vec_focal_tgt <- 2 - vec_focal_tgt
      }
      
      flank_cols_ref <- cols_ref[cols_ref != focal_col_ref]
      flank_cols_tgt <- cols_tgt[cols_tgt != focal_col_tgt]
      common_flanking_exact <- intersect(flank_cols_ref, flank_cols_tgt)
      
      if (length(common_flanking_exact) < 2) {
        snps_ref_clean <- gsub("_.*", "", flank_cols_ref)
        snps_tgt_clean <- gsub("_.*", "", flank_cols_tgt)
        shared_snps <- intersect(snps_ref_clean, snps_tgt_clean)
        
        if (length(shared_snps) < 2) {
          to_log_dt <- rbind(to_log_dt, data.table(
            DataSet = ds, SNP = focal_snp, 
            freq_ref_base = cur_row$freq_ref_base, 
            freq_tgt_base = cur_row$freq_tgt_base, 
            freq_diff = cur_diff, 
            Reason = "No_shared_flanking_markers"
          ))
          next
        }
        
        for (s in shared_snps) {
          c_ref <- flank_cols_ref[snps_ref_clean == s]
          c_tgt <- flank_cols_tgt[snps_tgt_clean == s]
          al_ref <- sub(".*_", "", c_ref)
          al_tgt <- sub(".*_", "", c_tgt)
          if (al_ref != al_tgt) {
            raw_tgt[, (c_tgt) := 2 - get(c_tgt)]
            setnames(raw_tgt, c_tgt, c_ref)
          }
        }
        common_flanking_exact <- flank_cols_ref[snps_ref_clean %in% shared_snps]
      }
      
      mat_flank_ref <- as.matrix(raw_ref[, ..common_flanking_exact])
      mat_flank_tgt <- as.matrix(raw_tgt[, ..common_flanking_exact])
      
      r_ref <- cor(vec_focal_ref, mat_flank_ref, use = "pairwise.complete.obs")
      r_tgt <- cor(vec_focal_tgt, mat_flank_tgt, use = "pairwise.complete.obs")
      
      strong_ld <- which(!is.na(r_ref) & !is.na(r_tgt) & abs(r_ref) >= MAF_threshold)
      
      if (length(strong_ld) < 2) {
        if (cur_diff > MAF_threshold) {
          to_log_dt <- rbind(to_log_dt, data.table(
            DataSet = ds, SNP = focal_snp, 
            freq_ref_base = cur_row$freq_ref_base, 
            freq_tgt_base = cur_row$freq_tgt_base, 
            freq_diff = cur_diff, 
            Reason = "Low_LD_correlation"
          ))
        }
        next
      }
      
      r_prod <- r_ref[strong_ld] * r_tgt[strong_ld]
      prop_negative <- mean(r_prod < 0)
      ld_concordance <- sum(r_prod)
      
      if (prop_negative >= 0.75 && ld_concordance < -0.25) {
        ld_snps_to_flip <- c(ld_snps_to_flip, focal_snp)
      } else if (prop_negative <= 0.25 && ld_concordance > 0.10) {
        # Надежное совпадение ориентации цепей
      } else {
        if (cur_diff > MAF_threshold) {
          to_log_dt <- rbind(to_log_dt, data.table(
            DataSet = ds, SNP = focal_snp, 
            freq_ref_base = cur_row$freq_ref_base, 
            freq_tgt_base = cur_row$freq_tgt_base, 
            freq_diff = cur_diff, 
            Reason = "Ambiguous_LD_signal"
          ))
        }
      }
    }
    to_flip <- c(to_flip, ld_snps_to_flip)
  }
  
  # Применяем переворот для однозначно определенных маркеров
  if (length(to_flip) > 0) {
    fwrite(data.table(to_flip), "./fin/tmp_flip.txt", col.names = FALSE)
    system(paste0(con_plink_start, " --bfile ", target_prefix, 
                  " --flip ./fin/tmp_flip.txt ", con_plink_end, target_prefix))
  }
  
  # Сохраняем логи текущего датасета в общий накопитель
  if (nrow(to_log_dt) > 0) {
    ambiguous_log_list[[ds]] <- to_log_dt
  }
  
  # Очистка временных файлов шага
  unlink(c("./fin/tmp_snplist.txt", "./fin/tmp_ref.*", "./fin/tmp_tgt.*", "./fin/tmp_flip.txt"))
}

#### Записываем ОНП на которые стоит обратить внимание
write.table(rbindlist(ambiguous_log_list), paste0(dir_o, 'AttentionSNP.tsv'),
            col.names = T, row.names = F, sep = '\t', quote = F)



# ==============================================================================
# ШАГ 5: ИТЕРАТИВНЫЙ MERGE ЧЕРЕЗ PLINK
# ==============================================================================
dir_i  <- './fin/rs_bed/'
dir_o  <- './fin/comb/'
sufx_i <- '_fin_rs'
sufx_f <- '_flip'
file_o <- 'unite'

if (length(fr) >= 2) {
  for (i in seq_along(fr)) {
    bed <- fr[i]
    if (i == 1) next
    
    if (i == 2) {
      bed_file   <- paste0(dir_i, fr[1], sufx_i)
      bed_merged <- paste0(dir_i, fr[2], sufx_i)
    } else {
      bed_merged <- paste0(dir_i, file_o)
      bed_file   <- paste0(dir_i, bed, sufx_i)
    }
    
    system(paste0(con_plink_start, " --bfile ", bed_file, " --bmerge ", bed_merged, " ", con_plink_end, dir_o, file_o))
    log <- paste(readLines(paste0(dir_o, file_o, ".log")), collapse = "\n")
    
    code_p <- PLINK_error_try_catch(log, bed_file, bed_merged, dir_i, dir_o, sufx_f, file_o)
    
    if (grepl('_', code_p)) {
      if (grepl('4', code_p)) {
        bed_file_flipped <- paste0(bed_file, sufx_f)
        system(paste0(con_plink_start, " --bfile ", bed_file_flipped, " --bmerge ", bed_merged, " ", con_plink_end, dir_o, file_o))
        log <- paste(readLines(paste0(dir_o, file_o, ".log")), collapse = "\n")
        code_p <- PLINK_error_try_catch(log, bed_file_flipped, bed_merged, dir_i, dir_o, sufx_f, file_o)
        
        if (grepl('3', code_p)) {
          system(paste0(con_plink_start, " --bfile ", bed_file_flipped, " --bmerge ", bed_merged, " ", con_plink_end, dir_o, file_o))
        }
      } else if (grepl('2', code_p)) {
        system(paste0(con_plink_start, " --bfile ", bed_file, " --bmerge ", bed_merged, " ", con_plink_end, dir_o, file_o))
      }
    }
    
    file.copy(paste0(dir_o, file_o, c(".bim", ".bed", ".fam")), 
              paste0(dir_i, file_o, c(".bim", ".bed", ".fam")), 
              overwrite = TRUE)
  }
}
