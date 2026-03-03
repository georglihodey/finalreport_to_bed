### Распознавание патерна CHR-POS-Allele1-Allele2-rsSNP
CHR_POS_Allele1_Allele2_rsSNP <- function(bim) {
  # Преобразование каждого элемента вектора
  chr_to_chr <- setNames(c(0:31,33, 30,31,33),
                         c(c(0:31,33),c("X", "Y", "MT")))
  
  result <- as.data.frame.matrix(str_split(bim$V2, "-",simplify = T))
  result$V1 <- sapply(result$V1, function(x) (chr_to_chr[x]))
  result$V2 <- as.numeric(result$V2)
  rroows <- (c(!is.na(result$V1))&c(!is.na(result$V2)))
  if (any(rroows)){
    bim$V1[rroows] <- result$V1[rroows]
    bim$V4[rroows] <- result$V2[rroows]
  }
  
  result <- as.data.frame.matrix(str_split(gsub("^addit_|_ilmndup\\d+$", "", bim$V2), ":",simplify = T))
  if(ncol(result)>1){
    result$V1 <- sapply(result$V1, function(x) (chr_to_chr[x]))
    result$V2 <- as.numeric(result$V2)
    rroows <- (c(!is.na(result$V1))&c(!is.na(result$V2)))
    if (any(rroows)){
      bim$V1[rroows] <- result$V1[rroows]
      bim$V4[rroows] <- result$V2[rroows]
    }
  }
  return(bim)
  
}

##### Я хз вообще зачем так много "одинаковых" функций, просто так получилось
vector_nucleotide_complement2 <- function(vec) {
  # Таблица соответствия нуклеотидов
  complement_table <- setNames(c("T", "A", "G", "C", '0', "I", "D","PLUS","MINUS"),
                               c("A", "T", "C", "G", '0', "D", "I","MINUS","PLUS"))
  
  # Преобразование каждого элемента вектора
  result <- sapply(strsplit(vec, "/"), function(x) complement_table[x])
  
  return(result)
}
vector_nucleotide_complement3 <- function(vec) {
  # Таблица соответствия нуклеотидов
  complement_table <- setNames(c("T", "A", "G", "C", '0', "I", "D","PLUS","MINUS"),
                               c("A", "T", "C", "G", '0', "D", "I","MINUS","PLUS"))
  
  # Преобразование каждого элемента вектора
  result <- sapply(strsplit(vec, "/"), function(x) paste(complement_table[x], collapse="/"))
  
  return(result)
}
vector_nucleotide_complement4 <- function(vec) {
  # Таблица соответствия нуклеотидов
  complement_table <- setNames(c("T", "A", "G", "C", '0', "I", "D","PLUS","MINUS"),
                               c("A", "T", "C", "G", '0', "D", "I","MINUS","PLUS"))
  
  # Преобразование каждого элемента вектора
  result <- sapply(strsplit(vec, "/"), function(x) paste(complement_table[x][c(2,1)], collapse="/"))
  
  return(result)
}
forward_To_top_issue <- function(bim, SNPchimp) {
  SNPchimp_tmp <- SNPchimp[SNPchimp$SNP_name%in%bim$V2,]
  SNPchimp_tmp <- SNPchimp_tmp[SNPchimp_tmp$Alleles_A_B_FORWARD!=SNPchimp_tmp$Alleles_A_B_TOP,]
  SNPchimp_tmp <- SNPchimp_tmp[order(SNPchimp_tmp$SNP_name),]
  SNPchimp_tmp <- SNPchimp_tmp[!duplicated(SNPchimp_tmp$SNP_name),]
  bim_tmp <- bim[bim$V2%in%SNPchimp_tmp$SNP_name,]
  bim_tmp <- bim_tmp[order(bim_tmp$V2),]
  bim_tmp_non_comp <- bim_tmp[bim_tmp$V5!=vector_nucleotide_complement3(bim_tmp$V6),]
  SNPchimp_non_comp <- SNPchimp_tmp[SNPchimp_tmp$SNP_name%in%bim_tmp_non_comp$V2,]
  
  SNPchimp_tmp$Al_AF <- str_split(SNPchimp_tmp$Alleles_A_B_FORWARD,'\\/', simplify = T)[,1]
  SNPchimp_tmp$Al_BF <- str_split(SNPchimp_tmp$Alleles_A_B_FORWARD,'\\/', simplify = T)[,2]
  SNPchimp_tmp$Al_AT <- str_split(SNPchimp_tmp$Alleles_A_B_TOP,'\\/', simplify = T)[,1]
  SNPchimp_tmp$Al_BT <- str_split(SNPchimp_tmp$Alleles_A_B_TOP,'\\/', simplify = T)[,2]
  
  if(nrow(bim_tmp_non_comp[(bim_tmp_non_comp$V2==SNPchimp_non_comp$Al_AF|bim_tmp_non_comp$V5==SNPchimp_non_comp$Al_BF|
                            bim_tmp_non_comp$V6==SNPchimp_non_comp$Al_AF|bim_tmp_non_comp$V6==SNPchimp_non_comp$Al_BF),])<
     nrow(bim_tmp_non_comp[(bim_tmp_non_comp$V5==SNPchimp_non_comp$Al_AT|bim_tmp_non_comp$V5==SNPchimp_non_comp$Al_BT|
                            bim_tmp_non_comp$V6==SNPchimp_non_comp$Al_AT|bim_tmp_non_comp$V6==SNPchimp_non_comp$Al_BT),])){
    bim_tmp <- bim_tmp[bim_tmp$V5=='0',]
    for (snp in bim_tmp$V2){
      SNPchimp_tmp_tmp <- SNPchimp_tmp[SNPchimp_tmp$SNP_name==snp,]
      bim_Nucl <- bim_tmp[bim_tmp$V2==snp,c(5,6)]
      if (bim_Nucl[2]==SNPchimp_tmp_tmp$Al_AT){
        bim$V5[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_BT
      } else if (bim_Nucl[2]==SNPchimp_tmp_tmp$Al_BT){
        bim$V5[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_AT
      }else{
        next
      }
      return(bim)
    }
  }else {
    for (snp in SNPchimp_tmp$SNP_name){
      SNPchimp_tmp_tmp <- SNPchimp_tmp[SNPchimp_tmp$SNP_name==snp,]
      bim_tmp[bim_tmp$V2==snp,]
      bim_Nucl <- bim_tmp[bim_tmp$V2==snp,c(5,6)]
      if (bim_Nucl[2]==SNPchimp_tmp_tmp$Al_AF){
        bim$V6[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_AT
        bim$V5[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_BT
      } else if (bim_Nucl[2]==SNPchimp_tmp_tmp$Al_BF){
        bim$V6[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_BT
        bim$V5[bim$V2==snp] <- SNPchimp_tmp_tmp$Al_AT
      }else{
        next
      }
    }
    return(bim)
  }
}

##### Вообще основной код по трансформации c прекрасного канала https://www.youtube.com/@GenomicsBootCamp/videos

######## Illumina Final Reports to PLINK files https://www.youtube.com/watch?v=5_pNby7l2dE&t=1s
######## How to solve the SNP data merge error - TOP and FORWARD allele coding https://www.youtube.com/watch?v=RhHI0lwXDNI

########################
# Transform final report to PLINK files
#########################
# Note: This is *one* of the possible solutions
#       Any other way could be used that creates the same file structure
#########################

# Clear workspace"~/Documents/cow/Manifests/synonims_snp.txt""~/Documents/cow/Manifests/synonims_snp.txt""~/Documents/cow/Manifests/synonims_snp.txt"
#rm(list = ls())
# Set working directory
#load packages
library(tidyverse)

DATA_PATH <- "Путь до ваших FinalReport'ов"

###Задаём рабочую директорию
setwd("РАБОЧАЯ ДИРЕКТОРИЯ")

###Загружаем манифест
manifest <- read.table('~/ПУТЬ ДО ФАЙЛА/unite_manifest_umd3.txt', header = T,
                             sep=',')[c('Name','SNP','Chr', 'MapInfo')]

manifest$Chr[manifest$Chr=="X"] <- 30
manifest$Chr[manifest$Chr=="Y"] <- 31
manifest$Chr[manifest$Chr=="MT"] <- 33
manifest$Chr <- as.numeric(manifest$Chr)
###Загружаем таблицу синонимов
synonims_snp <- read.table("~/ПУТЬ ДО ФАЙЛА/synonymus_snp3.txt", header = T,
                           sep=',')

###Загружаем SNPchimp файл с Forward-Top трансформацией
SNPchimp <- read.table("~/ПУТЬ ДО ФАЙЛА/SNPchimp.tsv", header=T, sep = '\t')

con_plink_start <- "~/plink/plink --cow --nonfounders --allow-no-sex "
con_plink_end <- " --make-bed --out  "


# read in final report file
#
# base R - puts . in column names, instead of spaces


###Читаем полный список файналрепортов
fr <- list.files(DATA_PATH)[grep('.txt',list.files(DATA_PATH))]
fr <- gsub('.txt','',fr)
###Задаём нужные колонки
nes_col <- c('SNPName', 'SampleID', 'Allele1Forward', 'Allele2Forward', 'Allele1Top', 'Allele2Top',  'GCScore', 'Chr', 'Position')

###Цикл по переработки файналрепортов в BED
###Делаем это через трансформацию FinalReport>LGEN>BED
for (FinalReport_File in fr[c(9:12)]){
  FinalReport <- read.csv(paste0(DATA_PATH, FinalReport_File,'.txt'), sep = "\t", header = T)
  names(FinalReport) <- gsub('\\.','',names(FinalReport))
  if(length(grep('Position', names(FinalReport)))==0){
    manifest_tmp <- manifest[manifest$Name%in%FinalReport$SNPName,]
    FinalReport <- left_join(FinalReport, manifest_tmp[c('Name','Chr', 'MapInfo')], join_by("SNPName"=="Name"))
    FinalReport$Chr[is.na(FinalReport$Chr)] <- 0
    FinalReport$MapInfo[is.na(FinalReport$MapInfo)] <- 0
    names(FinalReport) <- gsub('MapInfo','Position',names(FinalReport))
  }
  
  nes_col1 <- nes_col[nes_col%in%names(FinalReport)]
  FinalReport <- FinalReport[,nes_col1]
  
  if ("Allele1Top"%in%nes_col1){
    appendix <- 'Allele1Top'
    names(FinalReport) <- gsub('Top','',names(FinalReport))
    
  } else if ("Allele1Forward"%in%nes_col1){
    appendix <- 'Allele1Forward'
    names(FinalReport) <- gsub('Forward','',names(FinalReport))
    
  }else if ("Allele1AB"%in%nes_col1){
    appendix <- 'Allele1AB'
    names(FinalReport) <- gsub('AB','',names(FinalReport))
  }
  
  if(sum(FinalReport$GCScore, na.rm = T)!=0){
    FinalReport  <- subset(FinalReport, FinalReport$GCScore>0.30)
  }
  
  FinalReport$Allele1[FinalReport$Allele1=='-']=0
  FinalReport$Allele2[FinalReport$Allele2=='-']=0
  
  FinalReport %>%
    distinct(`SampleID`) %>%
    mutate(FID = FinalReport_File, sire = 0, dam = 0, sex = 0, phenotype = -9) %>%
    relocate(`SampleID`, .after = FID) %>%
    write.table(paste0('./', FinalReport_File,".fam"), col.names = F,row.names = F,
                quote = F, sep ='\t')
  # Lgen file
  FinalReport %>%
    mutate(FID = FinalReport_File) %>%
    dplyr::select(`FID`, `SampleID`, `SNPName`, `Allele1`, `Allele2`) %>%
    write.table(paste0('./', FinalReport_File,".lgen"), col.names = F,row.names = F,
                quote = F, sep ='\t')  # Map file
  FinalReport %>%
    distinct(`SNPName`, .keep_all = TRUE) %>%
    mutate(morgan = 0) %>%
    dplyr::select(`Chr`, `SNPName`, morgan, `Position`) %>%
    write.table(paste0('./', FinalReport_File,".map"), col.names = F,row.names = F,
                quote = F, sep ='\t')
  
  system(paste0(con_plink_start, " --lfile ./",
                FinalReport_File,con_plink_end, "  ./arch/", 
                FinalReport_File, ""))
  system(paste0("rm ./",FinalReport_File, ".lgen"))
  system(paste0("rm ./",FinalReport_File, ".map"))
  system(paste0("rm ./",FinalReport_File, ".fam"))
}

rm(FinalReport, FinalReport_File)

###Тут можно подгрузить уже ранее сделанные BED файлы
#fr <- list.files("./arch/")[grep('.bed',list.files("./arch/"))]
#fr <- gsub('.bed','',fr)


###Здесь и далее будет сделаны папки в которых будут храниться промежуточные файлы 
###(созданные промежуточные файлы не удаляться в ходе выполнения скрипта!)
system('mkdir fin')
for (bed in fr){
  system(paste0("cp ./",bed,".bed ./fin/",bed,"_fin.bed"))
  system(paste0("cp ./",bed,".fam ./fin/",bed,"_fin.fam"))
  system(paste0("cp ./",bed,".bim ./fin/",bed,"_fin.bim"))
}

####Если в именах маркера есть rs и нет позиций востанавлниваем позиции
for (bed in fr){
  bim <- read.table(paste0('./fin/',bed,'_fin.bim'))
  bim$rsID <- (str_extract(bim$V2, "rs\\d+"))
  rsID <- na.omit(unique(unique(bim$rsID[bim$V1==0])))
  if(is_empty(rsID)){next}
  write.table(rsID, "./rsIDs.txt", col.names = F, row.names = F, quote = F)
  ###Нужен доступ к bcftools
  system("bcftools view -i 'ID=@~/ПУТЬ ДО ФАЙЛА/rsIDs.txt' ~/ПУТЬ ДО ФАЙЛА/9913_GCA_000003055.5_current_ids.vcf.gz -o ~/ПУТЬ ДО ФАЙЛА/rsID_subset.vcf")
  transversion_tmp <- tryCatch(
    read.table('~/Documents/cow/rs_database/rsID_subset.vcf', skip = 45)[c(1:3)],
    error = function(e) NULL
  )
  if (is.null(transversion_tmp)) next
  vcf <- read.table('~/ПУТЬ ДО ФАЙЛА/rsID_subset.vcf', skip = 45)[c(1:3)]
  names(vcf) <- c('CHROM','POS','ID')
  for (rs in vcf$ID){
    bim$V1[bim$rsID==rs] <- vcf$CHROM[vcf$ID==rs][1]
    bim$V4[bim$rsID==rs] <- vcf$POS[vcf$ID==rs][1]
  }
  write.table(bim[-7],paste0('./fin/',bed,'_fin.bim'), col.names = F, row.names = F, quote = F, sep = ' ')
}

#### Приводим (по крайней мере пытаемся) к единообразию все маркеры из разных датасетов
for (bed in fr){
  bim <- read.table(paste0('./fin/',bed,'_fin.bim'))
  bim <- CHR_POS_Allele1_Allele2_rsSNP(bim)

  #Проблема с не теми названиями маркеров
  bim_zero <- bim
  bim_zero <- bim_zero[c(bim_zero$V1!=0&bim_zero$V4!=0),]
  bim_zero$ChrPos <- paste0(bim_zero$V1,'_',bim_zero$V4)
  bim_zero <- bim_zero[!bim_zero$V2%in%manifest$Name,]
  bim_zero <- bim_zero[!bim_zero$V2%in%na.omit(as.vector(as.matrix.data.frame(synonims_snp))),]
  
  manifest_zero <- manifest[c(manifest$MapInfo!=0&manifest$Chr!=0),]
  manifest_zero$ChrPos <- paste0(manifest_zero$Chr,'_',manifest_zero$MapInfo)
  manifest_zero <- manifest_zero[!duplicated(manifest_zero$ChrPos),]
  manifest_zero <- manifest_zero[!manifest_zero$Name%in%bim$V2,]

  bim_zero <- na.omit(left_join(bim_zero,manifest_zero, join_by('ChrPos')))
  bim_zero <- bim_zero[bim_zero$V2!=bim_zero$Name,]
  bim_zero <- na.omit(bim_zero[!duplicated(bim_zero$MapInfo),])

  if(nrow(bim_zero)!=0){
    for (snp in bim_zero$Name){
      bim_zero_tmp <- bim_zero[bim_zero$Name==snp,]
      if ((bim_zero_tmp$V5==0&bim_zero_tmp$V6==0)|
          !(any(c(bim_zero_tmp$V5,bim_zero_tmp$V6)%in%str_split(bim_zero_tmp$SNP,'\\/',simplify = T))|
          any(c(bim_zero_tmp$V5,bim_zero_tmp$V6)%in%vector_nucleotide_complement3(str_split(bim_zero_tmp$SNP,'\\/',simplify = T))))){
        next
      }
      if ((bim_zero_tmp$V5==0|bim_zero_tmp$V6==0)){
          if ((bim_zero_tmp$V5=='I'|bim_zero_tmp$V6=='I')|(bim_zero_tmp$V5=='D'|bim_zero_tmp$V6=='D')&
              (bim_zero_tmp$SNP=='I/D'|bim_zero_tmp$SNP=='D/I')){
            bim$V2[bim$V2==bim_zero_tmp$V2] <- bim_zero_tmp$Name
        } else if (str_split(bim_zero_tmp$SNP,'\\/',simplify = T)[,1]!=
                   vector_nucleotide_complement3(str_split(bim_zero_tmp$SNP,'\\/',simplify = T)[,2])){
          snp_nucl <- str_split(bim_zero_tmp$SNP,'\\/',simplify = T)
          if (any(snp_nucl%in%bim_zero_tmp$V6)){
            bim$V5[bim$V2==bim_zero_tmp$V2] <- snp_nucl[!snp_nucl%in%bim_zero_tmp$V6]
            bim$V6[bim$V2==bim_zero_tmp$V2] <- snp_nucl[snp_nucl%in%bim_zero_tmp$V6]
            bim$V2[bim$V2==bim_zero_tmp$V2] <- bim_zero_tmp$Name
          } else {
            bim$V5[bim$V2==bim_zero_tmp$V2] <- snp_nucl[!snp_nucl%in%vector_nucleotide_complement3(bim_zero_tmp$V6)]
            bim$V6[bim$V2==bim_zero_tmp$V2] <- snp_nucl[snp_nucl%in%vector_nucleotide_complement3(bim_zero_tmp$V6)]
            bim$V2[bim$V2==bim_zero_tmp$V2] <- bim_zero_tmp$Name
          }
        }
      }
      else{
        if ( all(c(bim_zero_tmp$V5,bim_zero_tmp$V6)%in%str_split(bim_zero_tmp$SNP,'\\/',simplify = T))|
             all(c(bim_zero_tmp$V5,bim_zero_tmp$V6)%in%vector_nucleotide_complement3(str_split(bim_zero_tmp$SNP,'\\/',simplify = T))) ){
          bim$V2[bim$V2==bim_zero_tmp$V2] <- bim_zero_tmp$Name
        } else{
          bim$V2[bim$V2==bim_zero_tmp$V2] <- paste0(bim_zero_tmp$Name,'_dup')
        }
      }
    }
  }  
  
  #Проблема с нулевыми аллелями
  bim_zero <- bim[bim$V5==0,]
  if (nrow(bim_zero)!=0){
    manifest_zero <- manifest[manifest$Name%in%bim_zero$V2,]
    if (nrow(manifest_zero)!=0){
      manifest_zero <- manifest_zero[
        vector_nucleotide_complement3(str_split(manifest_zero$SNP,'\\/',simplify = T)[,1]) !=
          (str_split(manifest_zero$SNP,'\\/',simplify = T)[,2]),
      ]
      for (snp in manifest_zero$Name){
        manifest_zero_tmp <- manifest_zero[manifest_zero$Name==snp,]
        
        bim[bim$V2==snp,]
        
        if(bim$V6[bim$V2==snp]==str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,1]|
           bim$V6[bim$V2==snp]==vector_nucleotide_complement3(str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,1])){
          bim$V6[bim$V2==snp] <- str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,1]
          bim$V5[bim$V2==snp] <- str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,2]
        } else if ( bim$V6[bim$V2==snp]==str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,2]|
                    bim$V6[bim$V2==snp]==vector_nucleotide_complement3(str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,2]) ) {
          bim$V6[bim$V2==snp] <- str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,2]
          bim$V5[bim$V2==snp] <- str_split(manifest_zero_tmp$SNP,'\\/',simplify = T)[,1]
        }
      }
      bim_zero <- bim[bim$V5==0,]
      if (nrow(bim_zero)!=0){
        write.table(bim_zero, paste0('./fin/transversion/',bed,'_fin_rs.txt'), col.names = F, row.names = F, quote = F, sep = '\t')
      }
      rm(bim_zero, manifest_zero_tmp,manifest_zero, manifest_tmp)
    }
  }
  #Проблема с нулевыми хромосомами и позициями
  bim_zero <- bim[bim$V4==0,]
  if (nrow(bim_zero)!=0){
    manifest_zero <- manifest[manifest$Name%in%bim_zero$V2,]
    manifest_zero <- manifest_zero[manifest_zero$Chr!=0,]
    for (snp in manifest_zero$Name){
      manifest_zero_tmp <- manifest_zero[manifest_zero$Name==snp,]
      bim$V1[bim$V2==snp] <- manifest_zero_tmp$Chr
      bim$V4[bim$V2==snp] <- manifest_zero_tmp$MapInfo
    }
    rm(bim_zero, manifest_zero_tmp,manifest_zero)
    
  }
  #Проблема с не теми хромосомами и позициями
  bim_zero <- bim
  bim_zero$ChrPos <- paste0(bim_zero$V1,'_',bim_zero$V4)
  manifest_zero <- manifest[manifest$Name%in%bim_zero$V2,]
  manifest_zero <- manifest_zero[manifest_zero$MapInfo!=0,]
  manifest_zero$ChrPos <- paste0(manifest_zero$Chr,'_',manifest_zero$MapInfo)
  manifest_zero <- left_join(manifest_zero,bim_zero, join_by('Name'=='V2'))
  manifest_zero <- manifest_zero[manifest_zero$ChrPos.x!=manifest_zero$ChrPos.y,]
  if(nrow(manifest_zero)!=0){
    for (snp in manifest_zero$Name){
      manifest_zero_tmp <- manifest_zero[manifest_zero$Name==snp,]
      bim$V1[bim$V2==snp] <- manifest_zero_tmp$Chr
      bim$V4[bim$V2==snp] <- manifest_zero_tmp$MapInfo
    }
  }
  
  # Проблема синонимов
  synonims_snp_tmp <- synonims_snp[synonims_snp$Name_2%in%bim$V2,]
  for (nr in 3:ncol(synonims_snp)){
    synonims_snp_tmp <- rbind(synonims_snp_tmp, synonims_snp[synonims_snp[[nr]]%in%bim$V2,])
  }
  rm(nr)
  synonims_snp_tmp <- unique(synonims_snp_tmp)
  if (length(synonims_snp_tmp$Name_1%in%bim$V2)!=0){
    synonims_snp_tmp <- synonims_snp_tmp[!synonims_snp_tmp$Name_1%in%bim$V2,]
  }
  
  if (nrow(synonims_snp_tmp)!=0){
    for (nr in 1:nrow(synonims_snp_tmp)){
      synonims_snp_tmp1 <- synonims_snp_tmp[nr,]
      snp_tmp1 <- synonims_snp_tmp1[which(synonims_snp_tmp1%in%bim$V2)][[1]]
      bim$V2[bim$V2==snp_tmp1] <- synonims_snp_tmp1[[1]]
    }
    rm(synonims_snp_tmp, synonims_snp_tmp1,snp_tmp1,nr)
  }
  
  #Исправляем Forward на Top
  bim <- forward_To_top_issue(bim, SNPchimp)
  
                 write.table(bim, paste0('./fin/rs_bed/',bed,'_fin_rs.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
  system(paste0("cp ./fin/",bed,"_fin.bed ./fin/rs_bed/",bed,"_fin_rs.bed"))
  system(paste0("cp ./fin/",bed,"_fin.fam ./fin/rs_bed/",bed,"_fin_rs.fam"))
  
  system(paste0(con_plink_start, " --bfile ./fin/rs_bed/",
                bed, "_fin_rs ",con_plink_end ," ./fin/rs_bed/", 
                bed, "_fin_rs"))
  

}

rm(bim, bim_zero,manifest_zero, manifest_zero_tmp,manifestOpa)


###### Существует проблема с комплементарными заменами A/T и G/C здесь попытка её разрешить
#### Это чатсь требует внимания и лучше глазками побегать, особенно если у вас в датасете 1 образец
transversion <- data.frame()
for (FinalReport_File in fr){
  transversion_tmp <- read.table(paste0('./fin/transversion/',FinalReport_File, '_fin_rs.txt'))
  transversion_tmp$DataSet <- FinalReport_File
  transversion <- rbind(transversion, transversion_tmp)
}
rm(transversion_tmp)
attention_snp <- c()
for (snp in unique(transversion$V2)){
  transversion_tmp <- transversion[transversion$V2==snp,]
  if (length(table(transversion_tmp$V6))==1){
    next
  } else {
    attention_snp <- c(attention_snp,snp)
  }
}
transversion <- transversion[transversion$V2%in%attention_snp,]

manifest$SNP2 <- vector_nucleotide_complement4(manifest$SNP)
transversion_snp <- manifest$Name[manifest$SNP2==manifest$SNP]
transversion <- transversion[transversion$V2%in%transversion_snp,]

freq <- data.frame()
for (FinalReport_File in fr){
  system(paste0(con_plink_start," --bfile ./fin/rs_bed/",
                FinalReport_File, "_fin_rs --freq --out ./fin/rs_bed/", 
                FinalReport_File, "_fin_rs"))
  freq_tmp <- read.table(paste0("./fin/rs_bed/",FinalReport_File, "_fin_rs.frq"), header = T)
  freq_tmp$Data <- FinalReport_File
  freq <- rbind(freq, freq_tmp)
}

freq <- freq[freq$SNP%in%unique(c(transversion$V2)),]
freq <- freq[!is.na(freq$MAF),]
#freq <- freq[freq$NCHROBS>2,]
snp <- "1_41756715"
for (snp in unique(freq$SNP)){
  freq_tmp <- freq[freq$SNP==snp,]
  if (length(table(freq_tmp$A2))>1) {
    ssaz <- freq_tmp$A2!=names(which.max(table(freq_tmp$A2)))
    problem <- freq_tmp$Data[!freq_tmp$A2%in%str_split(manifest$SNP[manifest$Name==snp],'\\/',simplify = T)]
    if (is_empty(problem)) {problem <- ''}
    
    for (FinalReport_File in freq_tmp$Data[ssaz]) {
      if (FinalReport_File%in%problem){
        bim <- read.table(paste0('./fin/rs_bed/',FinalReport_File,'_fin_rs.bim'))
        bim$V2[bim$V2==snp] <- paste0(snp,'_dup')
        write.table(bim, paste0('./fin/rs_bed/',FinalReport_File,'_fin_rs.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
      }else {
        bim <- read.table(paste0('./fin/rs_bed/',FinalReport_File,'_fin_rs.bim'))
        if (''%in%problem){nuclist <- freq_tmp$A2}else{nuclist <- freq_tmp$A2[!freq_tmp$Data%in%problem]}
        max_nuc <- names(which.max(table(nuclist)))
        bim$V6[bim$V2==snp] <- max_nuc
        nuclist <- nuclist[nuclist!=max_nuc]
        if (is_empty(nuclist)){next}
        bim$V5[bim$V2==snp] <- names(which.min(table(nuclist)))
        write.table(bim, paste0('./fin/rs_bed/',FinalReport_File,'_fin_rs.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
      }
    }
  }
}


rm(bim, freq, freq_tmp,transversion, transversion_snp,transversion_tmp)

dir_i <- './fin/rs_bed/'
dir_o <- './fin/comb/'
sufx_i <- '_fin_rs'
sufx_f <- '_flip'
sufx_t <- '_tmp'
file_o <- 'ural'

for (bed in fr){
  if (bed==fr[2]){
    
    bed_file <- paste0(dir_i, fr[1],sufx_i)
    bed_merged <- paste0(dir_i, fr[2],sufx_i)
    
    system(paste0(
      con_plink_start, " --bfile ",bed_file,
      
      " --bmerge ",bed_merged," ",
      
      con_plink_end, dir_o,file_o
    )
    )
    log <- read_file(paste0(dir_o,file_o,".log"))
    
    code_p <- PLINK_error_try_catch(log, bed_file, bed_merged, dir_i, dir_o, sufx_f,file_o)
    if (!is_empty(grep('\\_',code_p)==1)){
      if (!is_empty(grep('4',code_p))){
        bed_file <- paste0(dir_i,fr[1],sufx_i, sufx_f)
        system(
          paste0(
            con_plink_start, " --bfile ",bed_file,
            " --bmerge ",bed_merged," ",
            con_plink_end, " ", dir_o,file_o
          )
        )
        log <- read_file(paste0(dir_o,file_o,".log"))
        code_p <- PLINK_error_try_catch(log, bed_file, bed_merged, dir_i, dir_o, sufx_f,file_o)
        if (!is_empty(grep('3',code_p))){
          system(
            paste0(
              con_plink_start, " --bfile ",bed_file,
              " --bmerge ",bed_merged," ",
              con_plink_end, " ", dir_o,file_o
            )
          )
        }
        
      } else if (!is_empty(grep('2',code_p))&is_empty(grep('4',code_p))&is_empty(grep('3',code_p))){
        system(
          paste0(
            con_plink_start, " --bfile ",bed_file,
            
            " --bmerge ",bed_merged," ",
            
            con_plink_end, dir_o,file_o
          )
        )
      } 
    }
    
    system(paste0("cp ", dir_o, file_o,".bim ",dir_i,file_o,".bim"))
    system(paste0("cp ", dir_o, file_o,".bed ",dir_i,file_o,".bed"))
    system(paste0("cp ", dir_o, file_o,".fam ",dir_i,file_o,".fam"))
    
  }else if (bed==fr[1]){
    next
  }else{
    
    bed_merged <- paste0(dir_i, file_o)
    bed_file <- paste0(dir_i, bed,sufx_i)
    
    
    system(
      paste0(
        con_plink_start," --bfile ", bed_file, 
        " --bmerge ", bed_merged," ", 
        con_plink_end, " ",dir_o,file_o
      )
    )
    
    log <- read_file(paste0(dir_o,file_o,".log"))
    
    code_p <- PLINK_error_try_catch(log, bed_file, bed_merged, dir_i, dir_o, sufx_f,file_o)

    if (!is_empty(grep('\\_',code_p)==1)){
      if (!is_empty(grep('4',code_p))){
        bed_file <- paste0(dir_i,bed,sufx_i, sufx_f)
        system(
          paste0(
            con_plink_start, " --bfile ",bed_file,
            " --bmerge ",bed_merged," ",
            con_plink_end, " ", dir_o,file_o
          )
        )
        log <- read_file(paste0(dir_o,file_o,".log"))
        code_p <- PLINK_error_try_catch(log, bed_file, bed_merged, dir_i, dir_o, sufx_f,file_o)
        if (!is_empty(grep('3',code_p))){
          system(
            paste0(
              con_plink_start, " --bfile ",bed_file,
              " --bmerge ",bed_merged," ",
              con_plink_end, " ", dir_o,file_o
            )
          )
        }
        
      } else if (!is_empty(grep('2',code_p))&is_empty(grep('4',code_p))&is_empty(grep('3',code_p))){
        system(
          paste0(
            con_plink_start, " --bfile ",bed_file,
            
            " --bmerge ",bed_merged," ",
            
            con_plink_end, dir_o,file_o
          )
        )
      } 
    }
    
    system(paste0("cp ", dir_o, file_o,".bim ",dir_i,file_o,".bim"))
    system(paste0("cp ", dir_o, file_o,".bed ",dir_i,file_o,".bed"))
    system(paste0("cp ", dir_o, file_o,".fam ",dir_i,file_o,".fam"))
  }
}

PLINK_error_try_catch <- function(log, bed_file, bed_merged, dir_i, dir_o, sufx_f,file_o) {
  code_p <- '1'
  if (length(grep('Multiple positions seen for variant', log))!=0){
    warning_lines <- unlist(strsplit(log, "\n"))
    warning_lines <- warning_lines[grepl("Multiple positions seen for variant", warning_lines)]
    snp_from_warnings <- regmatches(warning_lines, gregexpr("'([^']+)'", warning_lines))
    snp_from_warnings <- unlist(snp_from_warnings)
    snp_from_warnings <- gsub("'", "", snp_from_warnings)
    snp_from_warnings_not <- snp_from_warnings[!snp_from_warnings%in%manifest$Name]
    
    snp_from_warnings <- snp_from_warnings[snp_from_warnings%in%manifest$Name]
    for (bed_cor in c(bed_file, bed_merged)){
      bim_cor <- read.table(paste0(bed_cor,'.bim'))
      for (snp_warnings in snp_from_warnings){
        bim_cor$V1[bim_cor$V2==snp_warnings] <- manifest$Chr[manifest$Name==snp_warnings]
        bim_cor$V4[bim_cor$V2==snp_warnings] <- manifest$MapInfo[manifest$Name==snp_warnings]
        
      }
      write.table(bim_cor, paste0(bed_cor,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    }
    bim_cor_file <- read.table(paste0(bed_file,'.bim'))
    bim_cor_merged <- read.table(paste0(bed_merged,'.bim'))
    
    for (snp_warnings in snp_from_warnings_not){
      if (bim_cor_file$V1[bim_cor_file$V2==snp_warnings]=='0'|bim_cor_file$V4[bim_cor_file$V2==snp_warnings]=='0'){
        bim_cor_file$V1[bim_cor_file$V2==snp_warnings]==bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings]
        bim_cor_file$V4[bim_cor_file$V2==snp_warnings]==bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings]
      } else if (bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings]=='0'|bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings]=='0'){
        bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings]==bim_cor_file$V1[bim_cor_file$V2==snp_warnings]
        bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings]==bim_cor_file$V4[bim_cor_file$V2==snp_warnings]
      }
    }
    write.table(bim_cor_file, paste0(bed_cor,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    write.table(bim_cor_merged, paste0(bed_merged,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    
    code_p <- paste0(code_p,'_2')
  }
  if (length(grep('Multiple chromosomes seen for variant', log))!=0){
    warning_lines <- unlist(strsplit(log, "\n"))
    warning_lines <- warning_lines[grepl("Multiple chromosomes seen for variant", warning_lines)]
    snp_from_warnings <- regmatches(warning_lines, gregexpr("'([^']+)'", warning_lines))
    snp_from_warnings <- unlist(snp_from_warnings)
    snp_from_warnings <- gsub("'", "", snp_from_warnings)
    snp_from_warnings_not <- snp_from_warnings[!snp_from_warnings%in%manifest$Name]
    
    snp_from_warnings <- snp_from_warnings[snp_from_warnings%in%manifest$Name]

    for (bed_cor in c(bed_file, bed_merged)){
      bim_cor <- read.table(paste0(bed_cor,'.bim'))
      
      
      for (snp_warnings in snp_from_warnings){
        bim_cor$V1[bim_cor$V2==snp_warnings] <- manifest$Chr[manifest$Name==snp_warnings]
        bim_cor$V4[bim_cor$V2==snp_warnings] <- manifest$MapInfo[manifest$Name==snp_warnings]
      }
      write.table(bim_cor, paste0(bed_cor,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    }
    bim_cor_file <- read.table(paste0(bed_file,'.bim'))
    bim_cor_merged <- read.table(paste0(bed_merged,'.bim'))
    
    for (snp_warnings in snp_from_warnings_not){
      if (bim_cor_file$V1[bim_cor_file$V2==snp_warnings]=='0'|bim_cor_file$V4[bim_cor_file$V2==snp_warnings]=='0'){
        bim_cor_file$V1[bim_cor_file$V2==snp_warnings] <- bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings]
        bim_cor_file$V4[bim_cor_file$V2==snp_warnings] <- bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings]
      } else if (bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings]=='0'|bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings]=='0'){
        bim_cor_merged$V1[bim_cor_merged$V2==snp_warnings] <- bim_cor_file$V1[bim_cor_file$V2==snp_warnings]
        bim_cor_merged$V4[bim_cor_merged$V2==snp_warnings] <- bim_cor_file$V4[bim_cor_file$V2==snp_warnings]
      }
    }
    write.table(bim_cor_file, paste0(bed_cor,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    write.table(bim_cor_merged, paste0(bed_merged,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
    
    code_p <- paste0(code_p,'_2')
  }
  
  if (length(grep("3\\+ alleles present", log))!=0){
    if (length(grep("_flip", log))!=0){
      TRI_alleles <- read.table(paste0(dir_o,file_o,"-merge.missnp") )[[1]]
      TRI_alleles_not <- TRI_alleles[!TRI_alleles%in%manifest$Name]
      TRI_alleles <- TRI_alleles[TRI_alleles%in%manifest$Name]
      bed_cor <- bed_merged
      for (bed_cor in c(bed_file, bed_merged)){
        bim_cor <- read.table(paste0(bed_cor,'.bim'))
        TRI_allele <- 'BTB-00385685'
        for (TRI_allele in TRI_alleles){
          bim_cor_tmp <- bim_cor[bim_cor$V2==TRI_allele,]
          manifest_tmp <- manifest[manifest$Name==TRI_allele,]
          if (all(c(bim_cor_tmp$V5,bim_cor_tmp$V6)%in%vector_nucleotide_complement2(manifest_tmp$SNP))|
              all(c(bim_cor_tmp$V5,bim_cor_tmp$V6)%in%str_split(manifest_tmp$SNP, '\\/',simplify = T ))){
            next
          } else {
            bim_cor$V2[bim_cor$V2==TRI_allele] <- paste0(TRI_allele,'_dup')
          }
        }
        write.table(bim_cor, paste0(bed_cor,'.bim'), col.names = F, row.names = F, quote = F, sep = '\t')
      }
      
      if (!is_empty(TRI_alleles_not)){
        write.table(paste0(dir_o,file_o,"-merge.missnp"), col.names = F, row.names = F, quote = F)
        
        system(paste0(
          con_plink_start," --bfile ", bed_file,
          " --exclude ",dir_o,file_o, "-merge.missnp ", 
          con_plink_end, " ", bed_file
        )
        )
        
      }
      code_p <- paste0(code_p,'_3')
      
      
    }else{
      system(
        paste0(
          con_plink_start, " --bfile ",bed_file, 
          
          " --flip ",dir_o,file_o,"-merge.missnp  ", 
          
          con_plink_end, " ",bed_file,sufx_f
        )
      )
      code_p <- paste0(code_p,'_4')
    }
  }
  return(code_p)
}
