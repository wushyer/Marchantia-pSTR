######build reference for MP######
perl -ne 'if (/^>(\S+)/) { close OUT; open OUT, ">$1.fasta" } print OUT' /MpTak_v7.1.fa
#run trf
for chrom in $(seq 1 8) U V
do
echo MP/raw_fasta/chr$chrom.fasta MP/trf_results 5
done | xargs -L 1 -P 30 ./scripts/run_TRF.sh
#filter over 6pn str
for chrom in $(seq 1 8) U V
do
echo scripts/fix_trf_output.py MP/trf_results/chr$chrom.fasta MP/fixed_trf_results/chr$chrom.fasta
done | xargs -L 1 -P 40 python
#3. reformt str

files=""
for chrom in $(seq 1 8) U V
do
files="$files,MP/fixed_trf_results/chr$chrom.fasta"
done
files=`echo $files | sed "s/,//"`
python scripts/trf_parser.py $files > filtered_repeats.MP.bed
bedtools sort -i filtered_repeats.MP.bed > filtered_repeats.MP.sorted.bed
#merge overlapping STRs into single entries and filter repeats that fail merging
python scripts/analyze_overlaps.py filtered_repeats.MP.sorted.bed pass.MP fail.MP
#then remove any entries within 10bp of a failed merge region:
bedtools window -w 10 -a pass.MP -b fail.MP -v > pass.MP.r2
#polish
bedtools merge -i pass.MP.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.MP.r3
bedtools merge -i pass.MP.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g" >> pass.MP.r3
#final
cat pass.MP.r3 | bedtools sort | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "MP_STR_"NR, $5}' > MP.hipstr_reference.bed
rm -r MP
rm fail.MP filtered_repeats.MP.bed filtered_repeats.MP.sorted.bed
rm pass.MP pass.MP.r2 pass.MP.r3

######run HipSTR#####

module load gcc/8.3.0
export LD_LIBRARY_PATH=/software/2020/software/gcccore/8.3.0/lib64:$LD_LIBRARY_PATH

/scratch-cbe/users/shuangyang.wu/str/WI-Ce-STRs/scripts/HipSTR_nf/HipSTR-0.7/HipSTR --bams bgld1.TAK1.markdup.bam,Bonn2.TAK1.markdup.bam,Bonn3.TAK1.markdup.bam,Bonn6.TAK1.markdup.bam,Cam1.TAK1.markdup.bam,Cam2.TAK1.markdup.bam,Canada1.TAK1.markdup.bam,Croatia.TAK1.markdup.bam,Czech1.TAK1.markdup.bam,Denmark1.TAK1.markdup.bam,Ehrenhausen3.TAK1.markdup.bam,Ehrenhausen6.TAK1.markdup.bam,France1.TAK1.markdup.bam,Germany1.TAK1.markdup.bam,Grafenegg1.TAK1.markdup.bam,Grafenegg2.TAK1.markdup.bam,Ireland10.TAK1.markdup.bam,Ireland12.TAK1.markdup.bam,Ireland13.TAK1.markdup.bam,Ireland18.TAK1.markdup.bam,Ireland19.TAK1.markdup.bam,Ireland1.TAK1.markdup.bam,Ireland21.TAK1.markdup.bam,Ireland22.TAK1.markdup.bam,Ireland2.TAK1.markdup.bam,Ireland36.TAK1.markdup.bam,Ireland37.TAK1.markdup.bam,Ireland38.TAK1.markdup.bam,Ireland3.TAK1.markdup.bam,Ireland5.TAK1.markdup.bam,Ireland6.TAK1.markdup.bam,Ireland7.TAK1.markdup.bam,Japan10.TAK1.markdup.bam,Japan12.TAK1.markdup.bam,Japan13.TAK1.markdup.bam,Japan16.TAK1.markdup.bam,Japan17.TAK1.markdup.bam,Japan18.TAK1.markdup.bam,Japan19.TAK1.markdup.bam,Japan1.TAK1.markdup.bam,Japan20.TAK1.markdup.bam,Japan21.TAK1.markdup.bam,Japan22.TAK1.markdup.bam,Japan26.TAK1.markdup.bam,Japan27.TAK1.markdup.bam,Japan28.TAK1.markdup.bam,Japan29.TAK1.markdup.bam,Japan2.TAK1.markdup.bam,Japan30.TAK1.markdup.bam,Japan4.TAK1.markdup.bam,Japan5.TAK1.markdup.bam,Japan6.TAK1.markdup.bam,Japan7.TAK1.markdup.bam,Japan9.TAK1.markdup.bam,MaG1.TAK1.markdup.bam,MiBa1.TAK1.markdup.bam,MiBa2.TAK1.markdup.bam,Norwich1.TAK1.markdup.bam,Norwich2.TAK1.markdup.bam,Norwich4.TAK1.markdup.bam,Oxford12.TAK1.markdup.bam,Oxford15.TAK1.markdup.bam,Oxford4.TAK1.markdup.bam,Oxford6.TAK1.markdup.bam,Schubert13.TAK1.markdup.bam,Schubert2.TAK1.markdup.bam,Schubert4.TAK1.markdup.bam,Schubert5.TAK1.markdup.bam,Schubert6.TAK1.markdup.bam,Schubert8.TAK1.markdup.bam,Schubert9.TAK1.markdup.bam,Sopron1.TAK1.markdup.bam,Tak1F7.TAK1.markdup.bam,Tak2F7.TAK1.markdup.bam,US1.TAK1.markdup.bam,Zurich14.TAK1.markdup.bam,Zurich6.TAK1.markdup.bam,Zurich9.TAK1.markdup.bam --fasta shuangyang.wu/ref/TakV7/TakV7.fa --regions /scratch-cbe/users/shuangyang.wu/str/HipSTR-references-1/MP.hipstr_reference.update.bed --chrom chr6 --haploid-chrs chr6 --str-vcf str_calls_6.vcf.gz

bcftools concat chr*.vcf.gz | bcftools sort -Oz -o STR_all_raw.vcf.gz

python2 /scratch-cbe/users/shuangyang.wu/str/WI-Ce-STRs/scripts/HipSTR_nf/HipSTR-0.7/scripts/filter_haploid_vcf.py --vcf STR_all_raw.vcf.gz --min-call-qual 0.9 --max-call-flank-indel 0.15 --max-call-stutter 0.15 > STR_filtered.vcf.gz

bcftools view STR_filtered.vcf.gz -Oz -o STR_all_filtered.vcf.gz

bcftools view STR_all_filtered.vcf.gz | bcftools filter -i "F_MISSING<0.1" -Oz -o STR_all_filtered_Fmiss01.vcf.gz

#####Table S1 general information of str and pstr#####

less -S filtered_repeats.MP.sorted.bed|awk  '{print $1"_"$2"_"$3"_"$5"\t"$0}' > format2
rm MP.hipstr_reference.bed
less -S MP.hipstr_reference.update.bed |awk  '{print $1"_"$2"_"$3"_"$NF"\t"$0}' > format1
perl shuangyang.wu/script/merge_files.pl -r format1 -n 1 -i format2 -c 1 -o out
less -S out |cut -f 2-8,17 >STR.all.info

RESULT:STR type annotation: (tools3) [shuangyang.wu@clip-login-1 HipSTR-references]$ nohup python 3new-hipstrANNOTATE.py  --gff shuangyang.wu/ref/TakV7/MpTak_v7.1.gff --te_bed ../1-hipstr/TE.bed --compartment ../1-hipstr/chrom_regions.tsv -t 10  STR.all.info STR.all.annotation.update.info  &
######str gene gff annotation#####
/scratch-cbe/users/shuangyang.wu/str/1-hipstr
module load cufflinks/2.2.1-foss-2018b
gffread MpTak_v7.1.gff -T -o my.gtf
module load kent_tools/20190507-linux.x86_64
gtfToGenePred  -genePredExt my.gtf genome_refGene.txt
shuangyang.wu/software/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile MpTak_v7.1.fa genome_refGene.txt --out genome_refGeneMrna.fa
shuangyang.wu/software/annovar/convert2annovar.pl -format vcf4old STR_all_filtered_Fmiss01.vcf.gz  >vcf.annovar.inputi
shuangyang.wu/software/annovar/annotate_variation.pl --geneanno --neargene 2000 --dbtype refGene -out tes-1 --buildver genome vcf.annovar.inputi shuangyang.wu/ref/TakV7/genome/

  
  
  bcftools query -f '%CHROM:%POS\t%ID\t%REF\t%ALT\t%INFO/BPDIFFS\t[%GT\t]\n' STR_all_filtered_Fmiss01.vcf.gz \
| awk '
BEGIN {
    OFS="\t";
    print "Variant_POS","ID","REF","ALT","BPDIFFS","allele_count","REF_count","ALT_count"
}
{
    pos=$1; id=$2; ref_allele=$3; alt_allele=$4; bpdiffs=$5;
    ref=0; alt=0; total=0;

    for(i=6; i<=NF; i++) {
        split($i, gt, "[|/]");
        for(j in gt) {
            if(gt[j] != "." && gt[j] != "") {
                total++;
                if(gt[j] == 0) ref++;
                else alt++;
            }
        }
    }
    print pos, id, ref_allele, alt_allele, bpdiffs, total, ref, alt
}' > hipstr_summary.tsv
#update this into a weigh script to do statistics
RESULT:python 1-hipstr.statistics-1.py #hipstr_summary_final.tsv

RESULT: perl shuangyang.wu/script/merge_files.pl -r ../HipSTR-references/STR.all.annotation.update.info -n 6 -i ../1-hipstr/hipstr_summary_final.tsv -c 2 -i ../1-hipstr/STR.pstr.info -c 1 -o out
less -S out |sort -k1,1 -k2,2n > out1

Update: 1-hipstr.statistics-2.py:less -S hipstr_summary_all_filtered.tsv |cut -f1,9-12 > expansionOrcontraction.tsv
python str_fraction.py STR_all_filtered.vcf.gz > STR_fraction.tsv &
  
#TABLE S1 20250806-STR_pSTR.XLS

####TABLE S1 in desktop STR/

############# Figure  1    ###############
#          pSTR distribution             #
##########################################
(tools3) [shuangyang.wu@clip-login-0 Table]$ less -S 20250813STR_pSTR.tsv|awk '{print $19}'|sort|uniq -c
14825 No
1 polymorphic
47315 Yes
(tools3) [shuangyang.wu@clip-login-0 Table]$ less -S 20250813STR_pSTR.tsv|awk '{print $12}'|sort|uniq -c
33055 No
1 polymorphicFM
29085 Yes


library(tidyverse)
library(ggplot2)
 
theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

ancestry.colours <- c("gold2","plum4", "darkorange1",   "lightskyblue2", 
                      "springgreen4", "lightpink2",  "deepskyblue4", 
                      "yellow3",  "yellow4",  
                      'black','red2', 'cornflowerblue', 'magenta', 'darkolivegreen4', 
                      'indianred1', 'tan4', 'darkblue', 'yellowgreen', "tan1",
                      'darkgray', 'wheat4', '#DDAD4B', 'chartreuse','seagreen1',
                      'moccasin', 'mediumvioletred', 'cadetblue1',"darkolivegreen1" ,"#7CE3D8",
                      "gainsboro","#E69F00","#009E73", "#F0E442", "sienna4", "#0072B2", 
                      "mediumpurple4","#D55E00", "burlywood3","gray51","#CC79A7","gray19", "firebrick") 

period_size_color <- c("1"="#b2182b","2"="#ef8a62","3"="#fddbc7",
                       "4"= "#e6f598","5" = "#99d594", "6" = "#3288bd" )

fisher_enrich_func <- function(pt_deg){
  sis <- pt_deg$sig_n  
  silp <- pt_deg$n - sis  
  fis <- pt_deg$sig_total - sis   
  filp <- pt_deg$total - pt_deg$sig_total - silp  
  ftp <- fisher.test(matrix(c(sis,silp,fis,filp), 2, 2), alternative='greater')
  return(ftp$p.value)
}

table_s1 <- data.table::fread("~/Desktop/STR//20250813STR_pSTR.tsv") 
table_s1_poly <- table_s1 %>% 
  dplyr::filter(polymorphicFM=="Yes")

str_dist_polym <- table_s1_poly %>% 
  dplyr::filter(!Chr=="chrU")  %>%
  dplyr::filter(!Chr=="chrV")

#domain_count_polym <- str_dist_polym  %>% 
  #dplyr::group_by(Chr,domain,domain_start ,domain_end ) %>% 
  #dplyr::count() %>%
  #dplyr::mutate(npm=n*1e6/(domain_end+1-domain_start),
                #Pos=(domain_end-domain_start)/2+domain_start) 
#manual revise chr3
domain_count_polym<-read.table("~/Desktop/STR/domain.tsv",header = T)
fig_1a <-  ggplot() + 
  geom_histogram(data=str_dist_polym, aes(x=start/1e6), bins = 50,fill="gray69") +
  geom_line(data=domain_count_polym, aes(x=Pos/1e6,y=npm/2),color="red",size=0.5 ) +
  geom_point(data=domain_count_polym, aes(x=Pos/1e6,y=npm/2),fill="red",shape=25,size=2 ) +
  scale_y_continuous( name = "Number of\npSTRs", sec.axis = sec_axis( trans=~.*2, name="Number of\npSTRs / Mb") )+
  facet_grid(.~Chr,scales = "free", space="free") +
  theme_cust +
  xlab("Genomic Position (Mb)") 

fig_1a <- ggplot() + 
  geom_histogram(
    data = str_dist_polym,
    aes(x = start / 1e6),
    bins = 50,
    fill = "gray69"
  ) +
  geom_line(
    data = domain_count_polym,
    aes(x = Pos / 1e6, y = npm),
    color = "#1B9E77",   # 柔和绿，更优雅
    size = 0.8,
    lineend = "round"
  ) +
  
  # 美化后的趋势点
  geom_point(
    data = domain_count_polym,
    aes(x = Pos / 1e6, y = npm),
    shape = 21,           # 实心圆点
    size = 2,
    fill = "#1B9E77",     # 同色填充
    color = "black",      # 外边框微微突出
    stroke = 0.3,
    alpha = 0.8
  ) +
  
  # 双坐标轴设置（保持原样）
  scale_y_continuous(
    name = "Number of\npSTRs",
    sec.axis = sec_axis(~ . * 2, name = "Number of\npSTRs / Mb")
  ) +
  
  # 染色体分面
  facet_grid(. ~ Chr, scales = "free", space = "free") +
  
  # 自定义主题
  theme_cust +
  xlab("Genomic Position (Mb)") +
  theme(
    strip.text = element_text(size = 13),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

polystr_motif_count  <- table_s1_poly %>% 
  dplyr::mutate( motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(motif_geno,motif_length ) %>% 
  dplyr::count() %>% 
  dplyr::ungroup() %>% 
  dplyr::top_n(10, n)  %>% 
  dplyr::arrange(  n )

polystr_motif_count$motif <- factor(polystr_motif_count$motif_geno, levels = polystr_motif_count$motif_geno )


fig_1b <- ggplot(polystr_motif_count,aes(x=motif,y=n ,fill=factor(motif_length))) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=period_size_color) +
  theme_cust +
  coord_flip() +
  theme(legend.position = "none") +
  labs(y="Number of sites",
       x="STR motif") +
  scale_y_continuous(breaks=c(0, 250,500,750,1000,1250,1500),limits = c(0,1800) )

###### fig_1c ######

PolySTR_region <- table_s1_poly %>% 
  dplyr::group_by(gfeature,motif_length) %>% 
  dplyr::add_count(name = "STRbyPS_REGION")%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::add_count(name = "STRby_REGION") %>% 
  dplyr::mutate(ps2region=100*STRbyPS_REGION/STRby_REGION) %>% 
  dplyr::distinct(gfeature,motif_length,STRbyPS_REGION,STRby_REGION,ps2region) %>% 
  dplyr::mutate(STRby_REGIONs=ifelse(motif_length==6,STRby_REGION,NA))


PolySTR_region$gfeatures<- factor(PolySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","TE","intergenic"))

#plot

fig_1c <- ggplot(PolySTR_region,aes(x=gfeatures,y=ps2region,fill=factor(motif_length))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  coord_flip() +
  theme_cust +
  xlab("Genomic features")+
  ylab("Percent of pSTRs (%)")  +
  labs(fill="Motif\nlength")+
  scale_fill_manual(values=period_size_color) +
  geom_text(aes(label=STRby_REGIONs),y=115,size = 10*5/14)+
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,125) ) +
  theme(legend.position = "left")+ 
  guides(fill = guide_legend(nrow = 6)) 

###### fig_1d ######

enrich_PolySTR_region_stats <- PolySTR_region %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(sig_n=STRbyPS_REGION,
                sig_total=STRby_REGION) %>% 
  dplyr::group_by(motif_length) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>% 
  dplyr::group_by(motif_length,gfeature)  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.1D")%>% 
  dplyr::mutate(group_factor="gfeature, motif_length",
                group_factor_catogory=paste0(gfeature,  ", ", motif_length),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrich_PolySTR_region <- enrich_PolySTR_region_stats %>% 
  dplyr::filter(fisherp<0.05) %>% 
  dplyr::mutate(logp=-log10(fisherp))  

enrich_PolySTR_region$gfeatures<- factor(enrich_PolySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","TE","intergenic"))


fig_1d <-ggplot(enrich_PolySTR_region,
                aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F)+
  theme_cust+
  theme( legend.position = "none")+
  scale_color_manual(values=period_size_color) +
  scale_x_continuous(breaks = c(0,100,200 ),labels = c("0","150", "200" )  )   +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 

###### fig_1e ######

enrichMotif_polySTR_region_stats <-  table_s1_poly %>% 
  dplyr::mutate( motif_geno=motif_geno_fwd ) %>% 
  dplyr::group_by(gfeature, motif_geno ) %>% 
  dplyr::count(name = "sig_n") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(sig_total=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by( motif_geno) %>% 
  dplyr::mutate(n=sum(sig_n)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(total=sum(sig_n)) %>%
  dplyr::group_by( gfeature,motif_geno )  %>% 
  dplyr::do(data.frame(x=fisher_enrich_func(.)))  %>%
  dplyr::rename(fisherp = x) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(fisherp_adj=p.adjust(fisherp,method="bonferroni")) %>% 
  dplyr::mutate(figure="FIG.1E")%>% 
  dplyr::mutate(group_factor="gfeature, motif_geno",
                group_factor_catogory=paste0(gfeature,  ", ", motif_geno),
                method="one-sided Fisher's Exact test",
                padjustment="BF")  

enrichMotif_polySTR_region <- enrichMotif_polySTR_region_stats %>% 
  dplyr::filter(fisherp_adj<0.05)   %>%  
  dplyr::mutate(logp=-log10(fisherp_adj)) %>% 
  dplyr::mutate(motif_length=nchar(motif_geno))%>% 
  dplyr::arrange(  motif_length )  %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::top_n(3, logp) 


enrichMotif_polySTR_region$gfeatures<- factor(enrichMotif_polySTR_region$gfeature,levels = c("promoter","enhancer","5'UTR","CDS","intron","3'UTR","TE","intergenic"))
fig_1e <- ggplot(enrichMotif_polySTR_region,
                 aes(y=gfeatures,x=logp,color=factor(motif_length)))+
  ggbeeswarm::geom_beeswarm(cex=2,groupOnX = F) +
  ggrepel::geom_text_repel(aes(label = motif_geno),  nudge_x = 3,segment.linetype=6,
                           max.overlaps=Inf,size=10*5/14,
                           box.padding = 0.5) +
  theme_cust+
  theme( legend.position = "none") +
  scale_color_manual(values=period_size_color) +
  ylab("Genomic\nfeatures")+
  xlab(expression(-log[10](italic(p)))) 
#for review
library(dplyr)
library(ggplot2)

## reviewer-requested plot:
## proportion of each motif-length category distributed across genomic features

library(dplyr)
library(tidyr)
library(ggplot2)

feature_order <- c(
  "intergenic",
  "TE",
  "3'UTR",
  "intron",
  "CDS",
  "5'UTR",
  "promoter"
)

feature_cols <- c(
  "intergenic" = "#E78AC3",
  "TE"         = "#8DA0CB",
  "3'UTR"      = "#1F78B4",
  "intron"     = "#66C2A5",
  "CDS"        = "#33A02C",
  "5'UTR"      = "#E6AB02",
  "promoter"   = "#E31A1C"
)

motif_length_order <- c(1, 2, 3, 4, 5, 6)

PolySTR_motif_feature <- table_s1_poly %>%
  filter(!Chr %in% c("chrU", "chrV")) %>%
  filter(!is.na(gfeature), !is.na(motif_length)) %>%
  mutate(
    motif_length = as.integer(motif_length),
    gfeature = as.character(gfeature),
    gfeature = case_when(
      gfeature %in% c("5UTR", "5_UTR") ~ "5'UTR",
      gfeature %in% c("3UTR", "3_UTR") ~ "3'UTR",
      TRUE ~ gfeature
    )
  ) %>%
  group_by(motif_length, gfeature) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(motif_length) %>%
  mutate(
    total_motif_length = sum(n),
    percent = 100 * n / total_motif_length
  ) %>%
  ungroup() %>%
  complete(
    motif_length = motif_length_order,
    gfeature = feature_order,
    fill = list(n = 0, percent = 0)
  ) %>%
  mutate(
    motif_length = factor(motif_length, levels = motif_length_order),
    gfeature = factor(gfeature, levels = feature_order)
  )

# 检查每个 motif length 是否都等于 100
PolySTR_motif_feature %>%
  group_by(motif_length) %>%
  summarise(total_percent = sum(percent), .groups = "drop") %>%
  print()

fig_review_motif_feature <- ggplot(
  PolySTR_motif_feature,
  aes(x = motif_length, y = percent, fill = gfeature)
) +
  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.25,
    position = "stack"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = feature_cols,
    breaks = feature_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = "Motif length",
    y = "Proportion of pSTRs (%)",
    fill = "Genomic\nfeature"
  ) +
  theme_cust +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

fig_review_motif_feature




###### fig 1 #####

fig_1bc <-  cowplot::plot_grid(fig_1b, fig_1c,
                               labels = c('', 'C'), 
                               rel_widths =  c(1.3,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)

fig_1de <-  cowplot::plot_grid(fig_1d,fig_1e,
                               labels = c('', 'E'), 
                               rel_widths =  c(1.2 ,2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "t",
                               nrow = 1)


fig_1 <-  cowplot::plot_grid(fig_1a, fig_1bc,fig_1de,  
                             labels = c('A', 'B',"D" ), 
                             rel_heights =  c(1,1.5,1 ),
                             label_size = 12, 
                             label_fontfamily="Helvetica",
                             axis = "lr",
                             nrow = 3)

ggsave(fig_1, filename = paste( "~/Desktop/STR/Fig_1_final.pdf",sep = ""), units = "mm",height = 160, width = 170)



#for review PCA 19 climate factor

library(tidyverse)
library(data.table)
library(ggplot2)
library(ggrepel)

# =========================
# 1. Read climate data
# =========================

pheno <- fread("~/Desktop/STR/pheno.csv") %>%
  as.data.frame()

# If the first column is unnamed, rename it as accession
colnames(pheno)[1] <- "accession"

bio_vars <- paste0("bio", 1:19)

clim_df <- pheno %>%
  dplyr::select(accession, all_of(bio_vars)) %>%
  dplyr::distinct()

# =========================
# 2. Climate PCA
# =========================

clim_pca <- prcomp(
  clim_df %>% dplyr::select(all_of(bio_vars)),
  center = TRUE,
  scale. = TRUE
)

# Variance explained
pca_var <- data.frame(
  PC = paste0("PC", seq_along(clim_pca$sdev)),
  variance_explained = clim_pca$sdev^2 / sum(clim_pca$sdev^2),
  cumulative_variance = cumsum(clim_pca$sdev^2 / sum(clim_pca$sdev^2))
)

pc1_pct <- round(pca_var$variance_explained[1] * 100, 1)
pc2_pct <- round(pca_var$variance_explained[2] * 100, 1)

write.table(
  pca_var,
  "Table_SXX_climate_PCA_variance_explained.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# 3. PCA scores
# =========================

pca_scores <- as.data.frame(clim_pca$x) %>%
  dplyr::mutate(accession = clim_df$accession)

# Define population group
# Japan: Japan*, Tak*, MiBa*, SopronMiBa
# Others are colored as Europe here, including samples clustering with Europe
pca_scores <- pca_scores %>%
  dplyr::mutate(
    population = dplyr::case_when(
      grepl("^Japan|^Tak", accession) ~ "Japan",
      TRUE ~ "Europe"
    )
  )

write.table(
  pca_scores,
  "Table_SXX_climate_PCA_scores.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# =========================
# 4. Plot 1: PCA scores
# =========================

fig_climate_pca_scores <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, label = accession, color = population)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = 0.3,
    color = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.3,
    color = "grey60"
  ) +
  geom_point(size = 2.2, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 2.6,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    segment.color = "grey70",
    segment.size = 0.2
  ) +
  scale_color_manual(
    values = c(
      "Europe" = "#2C7BB6",
      "Japan" = "#D7191C"
    )
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)"),
    color = NULL
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

fig_climate_pca_scores



# =========================
# 5. PCA loadings
# =========================

pca_loadings <- as.data.frame(clim_pca$rotation) %>%
  dplyr::mutate(variable = rownames(.)) %>%
  dplyr::select(variable, everything())

write.table(
  pca_loadings,
  "Table_SXX_climate_PCA_loadings.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Top GF variables
top_bio <- c("bio5", "bio4", "bio8", "bio10", "bio7")

pca_loadings_plot <- pca_loadings %>%
  dplyr::mutate(
    top_GF_variable = ifelse(
      variable %in% top_bio,
      "Top GF variable",
      "Other variable"
    )
  )

# =========================
# 6. Plot 2: PCA loadings
# =========================

fig_climate_pca_loading <- ggplot(
  pca_loadings_plot,
  aes(x = PC1, y = PC2, label = variable, color = top_GF_variable)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = 0.3,
    color = "grey60"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    size = 0.3,
    color = "grey60"
  ) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey70",
    segment.size = 0.25
  ) +
  scale_color_manual(
    values = c(
      "Other variable" = "#F8766D",
      "Top GF variable" = "#00BFC4"
    )
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = paste0("PC1 loading (", pc1_pct, "%)"),
    y = paste0("PC2 loading (", pc2_pct, "%)"),
    color = NULL
  ) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

fig_climate_pca_loading

#review cofound effects

/data2/wsst/SW_HS/str-fulldata/1-hipstr/pSTR-length

Rscript run_candidate_pSTR_PC_controlled_model.R

cat run_candidate_pSTR_PC_controlled_model.fig.r
library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
setwd("~/Desktop/STR/")

library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(patchwork)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

res <- fread("Table_SXX_candidate_pSTR_climate_association_PC1_PC2_controlled.tsv")

res <- res %>%
  mutate(
    short_label = paste0(CHROM, ":", POS)
  )

pA <- ggplot(res, aes(x = beta_climate_only, y = beta_PCcontrolled)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = significant_Bonferroni_005), shape = 21, color = "black", size = 3.2) +
  ggrepel::geom_text_repel(
    aes(label = short_label),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey70",
    segment.linewidth = 0.25
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#D55E00")) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  labs(
    x = "Climate effect before PC correction",
    y = "Climate effect after PC1/PC2 correction",
    fill = "Bonferroni\nsignificant\nafter PC control"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

pA

# ============================================================
# Panel B: significance before vs after PC correction
# ============================================================

plot_df <- res %>%
  mutate(
    locus_label = paste0(CHROM, ":", POS, " (", climate, ")"),
    locus_label = forcats::fct_reorder(locus_label, logp_PCcontrolled)
  )

pB <- ggplot(plot_df) +
  geom_segment(
    aes(
      y = locus_label,
      yend = locus_label,
      x = logp_climate_only,
      xend = logp_PCcontrolled
    ),
    color = "grey70",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = bonf_cutoff,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.4
  ) +
  geom_point(
    aes(x = logp_climate_only, y = locus_label),
    color = "#4C78A8",
    size = 3
  ) +
  geom_point(
    aes(x = logp_PCcontrolled, y = locus_label, fill = significant_Bonferroni_005),
    shape = 21,
    color = "black",
    size = 3.2
  ) +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#D55E00")) +
  theme_classic(base_size = 12) +
  labs(
    x = expression(-log[10](P)),
    y = NULL,
    fill = "Bonferroni\nsignificant\nafter PC control"
  ) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

# ============================================================
# Combine and save
# ============================================================

fig <- pA + pB +
  plot_layout(widths = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 15, face = "bold")
  )

fig

ggsave(
  "Fig_SXX_PCcontrolled_candidate_pSTRs.pdf",
  fig,
  width = 12.5,
  height = 5.2
)

ggsave(
  "Fig_SXX_PCcontrolled_candidate_pSTRs.png",
  fig,
  width = 12.5,
  height = 5.2,
  dpi = 300
)



############# Figure  2    ###############
#         expansion contraciton          #
##########################################

###### fig_2a ######  
#polymorphicFM : FM0.1 variant
#pSTR: homo allelic

library(ggplot2)
library(dplyr)
table_s1 <- data.table::fread("~/Desktop/STR/20250813STR_pSTR.tsv") 
table_s1_poly <- table_s1 %>% 
  dplyr::filter(MultiAllelic=="Yes")


bp_diff_BPD <- table_s1_poly %>% 
  dplyr::select(ref_STR,Chr,start,BPDIFFS) %>% 
  dplyr::arrange(Chr,start) %>% 
  splitstackshape::cSplit("BPDIFFS",",", direction = "long",sep = ",")  
bp_diff_BPD$BPDIFFS <- as.numeric(as.character(bp_diff_BPD$BPDIFFS))
fig_2a <- ggplot(bp_diff_BPD,aes(BPDIFFS))   + 
  geom_histogram(color="black", fill="white",bins =50,size=0.2) +
  theme_cust +
  ylab("Number of\nalleles") +
  xlab("Base-pair difference")

expansion_contractionS <- table_s1_poly %>% 
  dplyr::select(ref_STR ,Expansion,Contraction,motif_length) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length) %>% 
  dplyr::filter(!score==0) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length) %>% 
  dplyr::mutate(lg_sc=log10(abs(score))) %>% 
  dplyr::mutate(mean_lg_sc=mean(lg_sc),
                mean_sc=mean(score),
                median_sc=median(score)) %>% 
  dplyr::mutate(id=row_number(),
                mean_lg_sc=ifelse(id==1,round(mean_lg_sc,digits = 2),NA))


fig_2b <- ggplot(expansion_contractionS ,aes(score,fill=diff)) +
  geom_histogram(color="black", 
                 bins =100,size=0.2) +
  theme_cust +
  xlab("Contraction / Expansion score")+
  ylab("Number of\nSTRs")   +xlim(-1,2)+
  scale_fill_manual(values=c("#E7B800", "#00AFBB")) +
  theme( legend.position = "none") 

library(ggsci)
fig_2bNew<-ggplot(expansion_contractionS ,aes(score,fill=diff)) +
  geom_histogram(color="black", 
                 bins =100,size=0.2) +
  theme_cust +
  xlab("Contraction / Expansion score")+
  ylab("Number of\nSTRs")   +xlim(-1,1)+
  scale_fill_jco() +
  theme( legend.position = "none")
###### fig_2c ###### 

motif_length.labs <- c("Mono-STR", "Di-STR","Tri-STR","Tetra-STR","Penta-STR","Hexa-STR")
names(motif_length.labs) <- c("1", "2", "3", "4", "5", "6")

fig_2c <- ggplot(expansion_contractionS,aes(x=diff,y=lg_sc ,color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=mean_lg_sc),size=0.5,color="red") +
  theme_cust +
  ylab("Contraction /\nExpansion\ntransformed scores") +
  ggpubr::stat_compare_means(  label = "p.signif", method = "wilcox.test", 
                               label.x=1.35,
                               symnum.args = list(cutpoints = c(0, 0.00001, 0.0001, 0.001,  1), 
                                                  symbols = c("****","***","**",  "ns")) ,
                               label.y = 0.2) +
  facet_grid(.~motif_length,scales="free", 
             labeller = labeller(motif_length = motif_length.labs)) +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        legend.position = "none")  + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 0.5))

fig_2cNew <- ggplot(expansion_contractionS,aes(x=diff,y=lg_sc ,color=diff))+
  geom_quasirandom()+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=mean_lg_sc),size=0.5,color="red") +
  theme_cust +
  ylab("Contraction /\nExpansion\ntransformed scores") +
  ggpubr::stat_compare_means(  label = "p.signif", method = "wilcox.test", 
                               label.x=1.35,
                               symnum.args = list(cutpoints = c(0, 0.00001, 0.0001, 0.001,  1), 
                                                  symbols = c("****","***","**",  "ns")) ,
                               label.y = 0.2) +
  facet_grid(.~motif_length,scales="free", 
             labeller = labeller(motif_length = motif_length.labs)) +
  scale_color_jco() +
  theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
        legend.position = "none")  + 
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 0.5))





###### fig_2d ###### 

expand_frac <- table_s1_poly %>% 
  dplyr::select(ref_STR ,Contraction_fraction,Expansion_fraction,motif_length) %>% 
  tidyr::gather(diff,frac,-ref_STR,-motif_length) %>% 
  dplyr::mutate(frac=frac/100) %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(diff,motif_length)%>% 
  dplyr::mutate(mean_frac=mean(frac),
                median_frac=median(frac)) %>% 
  dplyr::mutate(id=row_number(),
                mean_frac=ifelse(id==1,round(mean_frac,digits = 2),NA)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate( diff=sub("(.*)(_frac)","\\1",diff))



fig_2d <- ggplot(expand_frac,aes(x=diff,y=frac ,color=diff))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=abs(mean_frac)),size=0.5,color="red") +
  theme_cust +
  ylab("Allele frequency") +
  xlab("Variation to median allele")+
  ggpubr::stat_compare_means(label = "p.signif",
                             label.y = 0.015 ,  label.x=1.35,
                             symnum.args = list(cutpoints = c(0,  0.001,  1), 
                                                symbols = c("****",   "ns")),
                             method = "wilcox.test" ) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_manual(values=c("#E7B800", "#00AFBB")) +
  theme(axis.text.x = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.02)) 

fig_2dNew <- ggplot(expand_frac,aes(x=diff,y=frac ,color=diff))+
  geom_violin()+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=diff,y=abs(mean_frac)),size=0.5,color="red") +
  theme_cust +
  ylab("Allele frequency") +
  xlab("Variation to median allele")+
  ggpubr::stat_compare_means(label = "p.signif",
                             label.y = 0.013 ,  label.x=1.35,
                             symnum.args = list(cutpoints = c(0,  0.001,  1), 
                                                symbols = c("****",   "ns")),
                             method = "wilcox.test" ) +
  facet_grid(.~motif_length,scales="free") +
  scale_color_jco() +
  theme(axis.text.x = element_blank(), 
        strip.text = element_blank(),
        legend.position = "none") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.015))




###### fig_2e ###### 

expansion_contractionS_di <- table_s1_poly %>% 
  dplyr::filter(!grepl("/",motif_geno)) %>% 
  dplyr::mutate(motif_geno=motif_geno_fwd ) %>%
  dplyr::select(ref_STR ,Expansion,Contraction,motif_length,motif_geno) %>% 
  tidyr::gather(diff,score,-ref_STR,-motif_length,-motif_geno) %>% 
  dplyr::filter(motif_length==2) %>% 
  dplyr::mutate(diff=ifelse(score==0,"substitution",diff)) %>% 
  dplyr::group_by(motif_geno,diff,motif_length) %>% 
  dplyr::count()  %>% 
  dplyr::group_by(motif_geno) %>% 
  dplyr::mutate(count_motif=sum(n)) %>% 
  dplyr::mutate(frac=100*n/count_motif) %>% 
  dplyr::mutate(count_motif2=ifelse(diff=="substitution",count_motif,NA),
                heading=ifelse(count_motif2==618,"The total number of di-STRs with each motif", NA)) %>% 
  dplyr::mutate(diff=sub("(.*)(_score)","\\1",diff)) 


fig_2e <- ggplot(expansion_contractionS_di,aes(x=motif_geno,y=frac,fill=factor(diff))) + 
  geom_bar(stat='identity', position = position_stack(reverse = TRUE)) +
  theme_cust +
  scale_fill_jco() +
  xlab("Motif")+
  ylab("Percent of\nALT alleles (%)")  +
  labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=105,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,120) )+
  geom_text(  aes(label=heading) ,y=117,size = 10*5/14  ) 


fig_2eNew<-ggplot(expansion_contractionS_di,aes(x=motif_geno,y=frac,fill=factor(diff))) +
  geom_bar(stat = "identity", aes(fill = factor(diff)),
           width = 0.5, color = "white") +
  labs(fill  = NULL,
       x = "Motif",
       y="Percent of\nALT alleles (%)") +labs(fill="Mutations") +
  geom_text(aes(label=count_motif2),y=105,size = 10*5/14) +
  scale_y_continuous(breaks=c(0, 50, 100),limits = c(0,120) )+
  geom_text(  aes(label=heading) ,y=117,size = 10*5/14  ) +
  scale_fill_jco() +
  theme_cust +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(color = "black"),
    plot.title = element_text(size = 12),
    plot.caption = element_text(hjust = 0)
  )



###### fig 2 #####

fig2ab <-  cowplot::plot_grid(fig_2a, fig_2bNew,
                              labels = c('', 'B'), 
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "tb",
                              nrow = 1)

fig2cd <-  cowplot::plot_grid(fig_2cNew,fig_2dNew,
                              labels = c('', 'D'), 
                              rel_heights =  c(1.1,1 ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              axis = "lr",
                              align = "v",
                              nrow = 2)

fig2 <-  cowplot::plot_grid(fig2ab, fig2cd,fig_2eNew,
                            labels = c('A', 'C' , "E"), 
                            rel_heights =  c(1,2,1.2),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            axis = "lr",
                            nrow = 3)


ggsave(fig2, filename = paste( "~/Desktop/STR//20250813_Fig2.pdf",sep = ""), units = "mm",height = 180, width = 170)



############# Figure  3    ###############
#             popgenet                    #
##########################################

###### fig_3a ######

table_s1 <- data.table::fread("~/Desktop/STR/20250813STR_pSTR.tsv") 
table_s1_poly <- table_s1 %>% 
  dplyr::filter(polymorphicFM=="Yes")


bp_diff_BPD <- table_s1_poly %>% 
  dplyr::select(ref_STR,Chr,start,BPDIFFS) %>% 
  dplyr::arrange(Chr,start) %>% 
  splitstackshape::cSplit("BPDIFFS",",", direction = "long",sep = ",")  
bp_diff_BPD$BPDIFFS <- as.numeric(as.character(bp_diff_BPD$BPDIFFS))



num_allele <- bp_diff_BPD %>% # fig_2a
  dplyr::group_by(ref_STR) %>% 
  dplyr::count() %>% 
  dplyr::mutate(n=n+1)

fig_3a <- ggplot(num_allele,aes(n))   + 
  geom_histogram(color="black", fill="white",bins=25,size=0.3) +
  theme_cust +
  xlab("Number of alleles\nper pSTR") +
  ylab("Number of pSTRs") + theme( axis.title.y = ggplot2::element_text(size=10,  color = "black", hjust = 0.2 )) 
#scale_y_continuous(breaks=c(0, 3000, 1500), limits=c(0, 3000))+ 
#scale_x_continuous(breaks=c(5,10,21) )  +
###### fig_3b ######
python 4-maf.py STR_all_filtered_Fmiss01.vcf All.major_allele.tsv
less -S All.major_allele.tsv > TableS_MAF.tsv
majorAF_Expected <- data.table::fread("~/Desktop/STR/TableS_MAF.tsv") 

fig_3b <- ggplot()   + 
  ggplot2::stat_density(data=majorAF_Expected,aes(x=Major_Allele_Freq ), geom="line",position = "identity",
                        size=0.8 ) +
  theme_cust +
  theme(legend.position = c(.4,.6) ,
        legend.spacing  = unit(0.1, 'cm'),
        legend.margin=margin(1,1,1,1)) + 
  ylab("Density") +
  xlab("Major\nallele frequency") +
  scale_x_continuous(breaks=c(0, 0.5, 1),limits = c(0,1) )


###PCA of STR
# 1. Load packages
install.packages("vcfR")
library(vcfR)

# 2. Read HipSTR VCF
vcf <- read.vcfR("~/Desktop/STR/STR_all_filtered_Fmiss01.vcf")

# 3. Extract GT field (haploid)
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# 4. Get REF/ALT sequences
ref <- getFIX(vcf)[, "REF"]
alt <- getFIX(vcf)[, "ALT"]

# 5. Map allele index to length
allele_lengths <- lapply(seq_along(ref), function(i) {
  alts <- unlist(strsplit(alt[i], ","))
  lengths <- c(nchar(ref[i]), nchar(alts))
  names(lengths) <- 0:(length(lengths)-1)
  lengths
})

# 6. Convert haploid GT to numeric matrix
geno_num <- matrix(NA, nrow = ncol(gt), ncol = nrow(gt),
                   dimnames = list(colnames(gt), rownames(gt)))

for (i in seq_len(nrow(gt))) {       # loop over loci
  for (j in seq_len(ncol(gt))) {     # loop over samples
    g <- gt[i, j]
    if (is.na(g) || g == ".") next
    # Haploid: one allele index only
    lengths <- allele_lengths[[i]][g]
    geno_num[j, i] <- lengths
  }
}

# 7. Run PCA (remove loci with any NA)
geno_complete <- geno_num[, colSums(is.na(geno_num)) == 0]

var_per_locus <- apply(geno_complete, 2, var, na.rm = TRUE)
geno_pca <- geno_complete[, var_per_locus > 0]

# Now run PCA
pca <- prcomp(geno_pca, center = TRUE, scale. = TRUE)
# 8. Inspect and plot
summary(pca)  # variance explained
plot(pca$x[,1], pca$x[,2],
     xlab = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"),
     pch = 19, col = "blue")
text(pca$x[,1], pca$x[,2], labels = rownames(pca$x), pos = 3)

pc1 <- pca$x[,1]
pc2 <- pca$x[,2]

"", ""

# 定义颜色规则
colors <- ifelse(pc1 < -0, "#00AFBB",
                 ifelse(pc1 > 0, "#E7B800", "grey80"))

# 绘图
plot(pc1, pc2,
     col = colors,
     pch = 19,
     xlab = paste0("PC1 (", round(100 * summary(pca)$importance[2,1], 1), "%)"),
     ylab = paste0("PC2 (", round(100 * summary(pca)$importance[2,2], 1), "%)"))


library(ggplot2)

# Create data frame for ggplot
df_pca <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pc1,
  PC2 = pc2
)

# Define colors (same logic as before)
df_pca$Color <- ifelse(df_pca$PC1 < 0, "#00AFBB",
                       ifelse(df_pca$PC1 > 0, "#E7B800", "grey80"))

# Make ggplot
fig_3c<-ggplot(df_pca, aes(x = PC1, y = PC2, color = Color)) +
  geom_point(size = 3) +
  scale_color_identity() +
  labs(
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
  ) +
  theme_cust 

table=read_tsv("~/Desktop/STR/STR.PCA")
country_colours <-
  c("Austria" = "#00A087",
    "Canada" = "#3C5488",
    "Croatia" = "#4DBBD5",
    "Czech" = "#9DAAC4",
    "Denmark" = "#E64B35",
    "UK" = "#8c2d04",
    "France"= "#762a83",
    "Germany"="#ec7014",
    "Hungary"="#fee090",
    "Ireland"="#c51b7d",
    "Japan"="black",
    "Switzerland"="#f4a582",
    "US" = "#5e4fa2"
  )
fig_3c<-ggplot(table,aes(x=PC1,y=PC2,color=Country),alpha=1)+geom_point(size=3, alpha=1)+theme_cust+scale_colour_manual(values = country_colours)+xlab("PCA1 7.1%")+ylab("PCA2 2.6%")

pi_data <- read_tsv("~/Desktop/STR//Europe.Asia.simplify.genetic.new.diversity")
library(tidyr)

pi_data_long <- pi_data %>%
  pivot_longer(
    cols = c(Europe_He, Japan_He),  # columns to reshape
    names_to = "pop",        # new column for names
    values_to = "Diversity"                 # new column for values
  )



fig_3d<-pi_data_long %>%
  ggplot(aes(POS * 20000, Diversity, 
             color = factor(pop, levels = c("Europe_He", "Japan_He"), labels = c("Europe", "Japan")),
             fill  = factor(pop, levels = c("Europe_He", "Japan_He"), labels = c("Europe", "Japan")))) +
  geom_smooth(span = 0.4, se = FALSE, method = "loess") +
  facet_grid(. ~ CHROM, scales = "free") +
  scale_color_manual(values = c("#FDD876", "#87c9c3")) +
  scale_fill_manual(values = c("#FDD876", "#87c9c3")) +
  theme_bw()+
  theme(
    legend.box = NULL,
    strip.background = element_blank(),
    legend.background = element_rect(colour = NA),
    legend.key = element_blank(),
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),axis.title.y = element_text(size = 8)
  ) +
  labs(x = "Genomic Position (Mb)", y = "pSTR Nei Diversity", color = "Population", fill = "Population")

ggplot(pi_data_long, aes(CHROM,Diversity,color=pop),guide="none")+geom_boxplot()+scale_color_aaas()+facet_grid(.~CHROM,scales = "free")+stat_compare_means()+theme_bw()+labs(y = "Diversity")

fig3abc <-  cowplot::plot_grid(fig_3a, fig_3b,fig_3c,
                               labels = c('', 'B','C'), 
                               rel_widths =  c(0.9,0.9 , 1.2),
                               label_size = 12, 
                               label_fontfamily="Helvetica",
                               axis = "tb",
                               align = "h",
                               nrow = 1)


fig3 <-  cowplot::plot_grid(fig3abc, fig_3d,   
                            labels = c('A' ,'D' ), 
                            rel_heights =  c(1.2,1.5),
                            label_size = 12, 
                            label_fontfamily="Helvetica",
                            align = "v",
                            nrow = 2)

ggsave(fig3, filename = paste( "~/Desktop/STR/20250813_Fig3.pdf",sep = ""), units = "mm",height = 180, width = 170)

ggplot(pi_data_long, aes(CHROM,Diversity,color=pop),guide="none")+geom_boxplot()+scale_color_aaas()+facet_grid(.~CHROM,scales = "free")+stat_compare_means()+theme_bw()+labs(y = "Diversity")

#Add new in fig3
table_s1 <- data.table::fread("~/Desktop/STR/20250813STR_pSTR.tsv") 
table_s1_poly <- table_s1 %>% dplyr::filter(polymorphicFM=="Yes")
pstr_fea <- table_s1_poly %>% dplyr::select(ref_STR,gfeature)
majorAF_Expected <- data.table::fread("~/Desktop/STR/TableS_MAF.Diversity") 
majorAF_Expected_fea <- majorAF_Expected %>% 
  dplyr::left_join(pstr_fea)%>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(id=row_number(),
                mean_het=round(mean(hets),digits = 3),
                mean_het2=ifelse(id==1,round(mean(hets),digits = 2),NA)) %>% 
  dplyr::ungroup()

majorAF_Expected_fea$gfeatures<- factor(majorAF_Expected_fea$gfeature,levels = c("promoter","5'UTR","CDS","intron","3'UTR","TE","intergenic"))

fig_S6a <- ggplot(majorAF_Expected_fea,aes(x=gfeatures,y=hets))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2)+ 
  geom_point(aes(x=gfeatures,y=mean_het2),size=0.5,color="red")+
  theme_cust +
  xlab("Genomic features")+
  ylab("Genetic diversity")+
  ylim(0,1) +
  theme(axis.text.x =  element_blank(),
        axis.title.x= element_blank() ) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 0.9,  p.adjust.method = "bonferroni", 
                             ref.group = "CDS")




# repeat_var <- table_s1_poly %>% 
#   dplyr::select(ref_STR, Chr, start, BPDIFFS, motif_length, gfeature) %>% 
#   dplyr::arrange(Chr, start) %>% 
#   splitstackshape::cSplit("BPDIFFS", ",", direction = "long", sep = ",") %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(BPDIFFS = as.numeric(BPDIFFS),   # <-- convert to numeric
#                 repeatN_var = abs(BPDIFFS) / motif_length) %>% 
#   dplyr::group_by(gfeature, ref_STR) %>% 
#   dplyr::mutate(mean_repeatN_var = mean(repeatN_var)) %>% 
#   dplyr::distinct(gfeature, ref_STR, mean_repeatN_var) %>% 
#   dplyr::group_by(gfeature) %>% 
#   dplyr::mutate(id = row_number(),
#                 gmean_repeatN_var = ifelse(id == 1, round(mean(mean_repeatN_var), digits = 2), NA)) %>% 
#   dplyr::ungroup()
# 
# 
# 
# repeat_var$gfeatures<- factor(repeat_var$gfeature,levels = c("promoter","5'UTR","CDS","intron","3'UTR","TE","intergenic"))
# 
# #plot
# 
# fig_S6b <- ggplot(repeat_var,aes(x=gfeatures,y=mean_repeatN_var))+
#   geom_violin(width=1.1,size=0.2)+
#   geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
#   geom_point(aes(x=gfeatures,y=gmean_repeatN_var),size=0.5,color="red") +
#   theme_cust +
#   xlab("Genomic features")+
#   ylab("Mean repeat\nnumber variance") +
#   theme(axis.text.x =  element_blank(),
#         axis.title.x= element_blank() ) +
#   ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 14,  
#                              p.adjust.method = "bonferroni", 
#                              symnum.args = list(cutpoints = c(0, 0.0001,0.05,  1), 
#                                                 symbols = c("****","*",  "ns")),
#                              ref.group = "CDS")+
#   ylim(0,15)




###### fig_S6c  ######


poly_gc <- table_s1_poly  %>% 
  dplyr::mutate(countA=stringr::str_count(motif_geno, "A"),
                countG=stringr::str_count(motif_geno, "G"),
                countC=stringr::str_count(motif_geno, "C"),
                countT=stringr::str_count(motif_geno, "T"),
                contentGC=100*(countG+countC)/(countG+countC+countA+countT)) %>% 
  dplyr::group_by(gfeature) %>% 
  dplyr::mutate(id=row_number(),
                gmean_gc=ifelse(id==1,round(mean(contentGC),digits = 2),NA)) %>% 
  dplyr::ungroup() 

poly_gc$gfeatures<- factor(poly_gc$gfeature,levels = c("promoter","5'UTR","CDS","intron","3'UTR","TE","intergenic"))

#plot

fig_S6c <- ggplot(poly_gc,aes(x=gfeatures,y=contentGC))+
  geom_violin(width=1.1,size=0.2)+
  geom_boxplot(width=0.08, color="black",fill="grey", alpha=0.2,outlier.shape = NA,size=0.2) + 
  geom_point(aes(x=gfeatures,y=gmean_gc),size=0.5,color="red") +
  theme_cust +
  xlab("Genomic features")+
  ylab("GC content (%)")+
  theme(axis.text.x =  element_text(size=10,   color = "black", angle = 45, hjust = 1)) +
  ggpubr::stat_compare_means(label = "p.signif", method = "wilcox.test",label.y = 101,
                             p.adjust.method = "bonferroni", 
                             symnum.args = list(cutpoints = c(0, 0.00001,0.0001,0.05,  1), 
                                                symbols = c("****","***","**",  "ns")),
                             ref.group = "CDS")+
  ylim(0,105)


fig_S6 <-  cowplot::plot_grid(fig_S6a,fig_S6c,  
                              labels = c('A' ,'B'), 
                              rel_heights =  c(1,1.5 ),
                              label_size = 12, 
                              label_fontfamily="Helvetica",
                              align = "v",
                              nrow = 2)

ggsave(fig_S6, filename = paste( "~/Desktop/STR//20250183_Fig3-add.pdf",sep = ""), units = "mm",height = 170, width = 170)

######fig 4#####

setwd("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM")
library(gradientForest)
library(MaizePal)

gfData <- read.table("1884.GF",header =T,sep="\t",row.names = "Pops")  # 加载环境数据、等位基因频率
candidate <- gfData[,grep("chr",names(gfData))] #提取包含candSNPs的列，等位基因频率
present <- gfData[,c(1,2,grep("bio",names(gfData)))]  # 提取第一列、第二列和包含bio的列，即坐标和生物气候数据
bioclimatic <- paste("bio",1:19,sep = "")  # 生成向量(bio_1, ..., bio_19)
maxLevel <- log2(0.368*nrow(candidate)/2)# 固定公式 log2（0.368*居群数/2）
gf_candidate <- gradientForest(cbind(present[,bioclimatic], candidate),  predictor.vars=colnames(present[,bioclimatic]),response.vars=colnames(candidate), ntree=500,maxLevel=maxLevel, trace=T, corr.threshold=0.50)  #决策树默认500，相关性阈值0.5/0.7
#plot(gf_candidate, plot.type = "Overall.Importance", col=c(rep("grey",15),MaizePal::maize_pal("HighlandMAGIC", 4) ),las=2,cex.names=0.8)
bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)]  

most_cand <- names(bio_cand[1]) 

barplot(bio_cand,las=2,cex.names=0.8,col=c(MaizePal::maize_pal("JimmyRed",5),rep("grey",14)),ylab="Weighted importance")


temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Max Temperature of Warmest Month (BIO5)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

#points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg=rev(id_wcol),cex=1.5) 

#points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg=id_ccol,cex=1.5) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 





most_cand <- names(bio_cand[2]) 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Temperature Seasonality (BIO4)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 


most_cand <- names(bio_cand[3]) 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Mean Temperature of Wettest Quarter (BIO8)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5)


most_cand <- names(bio_cand[4]) 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Mean Temperature of Warmest Quarter (BIO10)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5)


most_cand <- names(bio_cand[5]) 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Temperature Annual Range (BIO7)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5)




###QQ plot
theme_cust <- theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

df <- read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio5.pvalue.gwas.result", header = TRUE, sep = "\t")

# 构建 QQ plot 的数据
n <- nrow(df)
expected <- -log10(ppoints(n))  # 期望分位数
observed <- -log10(sort(df$P))  # 观测分位数排序

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# 绘图
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(shape = 19, size = 0.4, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(
    x = "Expected quantile BIO5 pvalue ",
    y = "Observed BIO5 pvalue"
  ) +
  theme_cust ->BIO5_QQ

df <- read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio4.pvalue.gwas.result", header = TRUE, sep = "\t")

# 构建 QQ plot 的数据
n <- nrow(df)
expected <- -log10(ppoints(n))  # 期望分位数
observed <- -log10(sort(df$P))  # 观测分位数排序

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# 绘图
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(shape = 19, size = 0.4, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(
    x = "Expected quantile BIO4 pvalue ",
    y = "Observed BIO4 pvalue"
  ) +
  theme_cust ->BIO4_QQ

df <- read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio7.pvalue.gwas.result", header = TRUE, sep = "\t")

# 构建 QQ plot 的数据
n <- nrow(df)
expected <- -log10(ppoints(n))  # 期望分位数
observed <- -log10(sort(df$P))  # 观测分位数排序

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# 绘图
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(shape = 19, size = 0.4, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(
    x = "Expected quantile BIO7 pvalue ",
    y = "Observed BIO7 pvalue"
  ) +
  theme_cust ->BIO7_QQ


df <- read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio8.pvalue.gwas.result", header = TRUE, sep = "\t")

# 构建 QQ plot 的数据
n <- nrow(df)
expected <- -log10(ppoints(n))  # 期望分位数
observed <- -log10(sort(df$P))  # 观测分位数排序

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# 绘图
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(shape = 19, size = 0.4, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(
    x = "Expected quantile BIO8 pvalue ",
    y = "Observed BIO8 pvalue"
  ) +
  theme_cust ->BIO8_QQ


df <- read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio10.pvalue.gwas.result", header = TRUE, sep = "\t")

# 构建 QQ plot 的数据
n <- nrow(df)
expected <- -log10(ppoints(n))  # 期望分位数
observed <- -log10(sort(df$P))  # 观测分位数排序

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# 绘图
ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(shape = 19, size = 0.4, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  labs(
    x = "Expected quantile BIO10 pvalue ",
    y = "Observed BIO10 pvalue"
  ) +
  theme_cust ->BIO10_QQ



bio5gwas<-read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio5.pvalue.gwas.result",header = T)
bio5gwas$sig <- ifelse(-log10(bio5gwas$P) > 5, "sig", "ns")

bio5Plot<-ggplot(bio5gwas, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_cust+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )+ggtitle("Max Temperature of Warmest Month (BIO5)")

bio7gwas<-read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio7.pvalue.gwas.result",header = T)
bio7gwas$sig <- ifelse(-log10(bio7gwas$P) > 5, "sig", "ns")

bio7Plot<-ggplot(bio7gwas, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_cust+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )+ggtitle("Temperature Annual Range (BIO7)")

bio4gwas<-read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio4.pvalue.gwas.result",header = T)
bio4gwas$sig <- ifelse(-log10(bio4gwas$P) > 5, "sig", "ns")

bio4Plot<-ggplot(bio4gwas, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_cust+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )+ggtitle("Temperature Seasonality (BIO4)")

bio8gwas<-read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio8.pvalue.gwas.result",header = T)
bio8gwas$sig <- ifelse(-log10(bio8gwas$P) > 5, "sig", "ns")

bio8Plot<-ggplot(bio8gwas, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_cust+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )+ggtitle("Mean Temperature of Wettest Quarter (BIO8)")

bio10gwas<-read.table("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/bio10.pvalue.gwas.result",header = T)
bio10gwas$sig <- ifelse(-log10(bio10gwas$P) > 5, "sig", "ns")

bio10Plot<-ggplot(bio10gwas, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_cust+
  theme(
    axis.text.x = element_blank(),
    legend.position = "none",
    panel.spacing = unit(0.1, "line")
  )+ggtitle("Mean Temperature of Warmest Quarter (BIO10)")

bio5Plot / bio4Plot / bio7Plot / bio8Plot / bio10Plot

bio_vars <- c("bio4", "bio5", "bio7", "bio8", "bio10")

# 循环处理每个变量
for (bio in bio_vars) {
  # 读取 GWAS 文件
  gwas_file <- paste0("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/", bio, ".pvalue.gwas.result")
  gwas_data <- read.table(gwas_file, header = TRUE)
  
  # 添加显著性标记
  gwas_data$sig <- ifelse(-log10(gwas_data$P) > 5, "sig", "ns")
  
  # 筛选显著位点
  sig_sites <- subset(gwas_data, sig == "sig")
  
  # 输出显著位点到文本文件
  out_file <- paste0("/scratch-cbe/users/shuangyang.wu/str/Table/ppopp/LFMM/", bio, "_sig_sites.txt")
  write.table(sig_sites[, c("CHR", "BP", "P")], 
              file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message("Finished ", bio)
}

library(ggplot2)
library(gggenes)

### 输入文件
gff_file <- "MpBHLH42.gff.2"          # gff文件
variant_file <- "var.id"     # 变异文件（至少两列：CHR, POS）
target_gene <- "Mp5g09710"         # 目标基因ID

### 读入 GFF
gff <- read.table(gff_file, sep="\t", header=FALSE, quote="")
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# 只保留目标基因相关的 feature
gff_gene <- gff[grep(target_gene, gff$attributes), ]

# 保留 CDS / UTR
conserve <- c("CDS","five_prime_UTR","three_prime_UTR")
gff_gene <- subset(gff_gene, type %in% conserve)

# 提取 gene name
gene_name <- sub(".*Name=([^;]+).*", "\\1", gff_gene$attributes[1])

# 整理数据
gene_df <- data.frame(
  molecule = gff_gene$seqid,
  type = gff_gene$type,
  start = gff_gene$start/1000,   # 转成 Kb
  end = gff_gene$end/1000,
  strand = gff_gene$strand,
  gene = gene_name
)

# strand 转成方向
gene_df$direction <- ifelse(gene_df$strand == "+", 1, -1)

### 读入 variants
variants <- read.table(variant_file, header=TRUE)  # 必须包含 CHR, POS
variants_sub <- subset(variants, CHR == unique(gene_df$molecule) &
                         POS >= min(gff_gene$start) &
                         POS <= max(gff_gene$end))
variants_sub$POS <- variants_sub$POS/1000  # 转 Kb

### 配色
type_color <- c("CDS"="#8dd3c7", "five_prime_UTR"="#ffffb3", "three_prime_UTR"="#bebada")

### 绘图
p <- ggplot(gene_df, aes(xmin=start, xmax=end, y=gene, fill=type, forward=direction)) +
  geom_gene_arrow() +
  scale_fill_manual(values=type_color[unique(gene_df$type)]) +
  geom_vline(data=variants_sub, aes(xintercept=POS), color="red", linetype="dashed", size=0.6) +
  xlab(paste0(unique(gene_df$molecule)," (Kb)")) + ylab("") +
  theme_genes() +
  theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

print(p)

#str lenght gwas

library(vcfR)
vcf <- read.vcfR("/scratch-cbe/users/shuangyang.wu/str-1/Table/STR_all_filtered_Fmiss01.vcf")
geno <- extract.gt(vcf) # matrix of genotypes
# Convert alleles to repeat lengths
# Assuming format: 10/12 etc. (repeat counts)
geno_numeric <- apply(geno, c(1,2), function(x){
  if(is.na(x)) return(NA)
  alleles <- strsplit(x, "/")[[1]]
  mean(as.numeric(alleles)) # average length per individual
})
write.csv(geno_numeric, "/scratch-cbe/users/shuangyang.wu/str-1/Table/pSTR_genotypes_numeric.csv", row.names = TRUE)


###qqplot for gemma
setwd("/scratch-cbe/users/shuangyang.wu/str-1/Table/gemma-strLength/output")

library(ggplot2)
library(readr)

# --- PARAMETERS ---
results_dir <- "./"
files <- list.files(results_dir, pattern = "\\.assoc.txt$", full.names = TRUE)

# --- FUNCTION: make QQ plot ---
make_qqplot <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Use Wald test p-values (you can switch to p_lrt or p_score if preferred)
  df <- df %>% filter(!is.na(p_wald))
  
  # Expected -log10(p) under uniform(0,1)
  n <- nrow(df)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(df$p_wald))
  
  qq_df <- data.frame(expected, observed)
  
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
      title = paste0("QQ plot: ", gsub("\\.assoc.*", "", basename(file))),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_minimal(base_size = 14)
  
  # Save PDF per phenotype
  outname <- paste0(gsub("\\.assoc.*", "", basename(file)), "_qqplot.pdf")
  ggsave(outname, plot = p, width = 5, height = 5)
  
  return(p)
}

# --- LOOP over files and generate QQ plots ---
plots <- lapply(files, make_qqplot)
library(ggplot2)
library(readr)

# --- PARAMETERS ---
results_dir <- "./"
files <- list.files(results_dir, pattern = "\\.assoc.txt$", full.names = TRUE)

# --- FUNCTION: make QQ plot ---
make_qqplot <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  # Use Wald test p-values (you can switch to p_lrt or p_score if preferred)
  df <- df %>% filter(!is.na(p_wald))
  
  # Expected -log10(p) under uniform(0,1)
  n <- nrow(df)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(df$p_wald))
  
  qq_df <- data.frame(expected, observed)
  
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(
      title = paste0("QQ plot: ", gsub("\\.assoc.*", "", basename(file))),
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) +
    theme_minimal(base_size = 14)
  
  # Save PDF per phenotype
  outname <- paste0(gsub("\\.assoc.*", "", basename(file)), "_qqplot.pdf")
  ggsave(outname, plot = p, width = 5, height = 5)
  
  return(p)
}

# --- LOOP over files and generate QQ plots ---
plots <- lapply(files, make_qqplot)

####gwas result loop
# Load packages
# Load packages
library(ggplot2)
library(dplyr)

# --- Load GEMMA result ---
file <- "pheno_5_results.assoc.txt"
df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# --- Check column names ---
colnames(df)
# Should include: "chr", "rs", "ps", "p_wald"

# --- Rename for standard GWAS plotting ---
df <- df %>%
  dplyr::rename(
    CHR = chr,
    BP  = ps,
    RS  = rs,
    P   = p_wald
  ) %>%
  dplyr::filter(!is.na(P))  # Remove NA p-values

# --- QQ Plot ---
n <- nrow(df)
expected <- -log10(ppoints(n))
observed <- -log10(sort(df$P))
qq_df <- data.frame(expected, observed)

qqplot <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = "Expected -log10(P)", y = "Observed -log10(P)",
       title = "QQ plot: pheno_5") +
  theme_minimal(base_size = 14)

ggsave("pheno_5_qqplot.pdf", qqplot, width = 5, height = 5)

# --- Manhattan Plot ---
df$sig <- ifelse(-log10(df$P) > 5, "sig", "ns")  # significance threshold 1e-5

manhattan_plot <- ggplot(df, aes(x = BP, y = -log10(P), color = sig)) +
  geom_point(size = 1) +
  facet_grid(. ~ CHR, scales = "free") +
  scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
  ylab(expression(-log[10](italic(P)))) +
  xlab("Genomic position (Mb)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.1, "line")) +
  ggtitle("Manhattan Plot: pheno_5")

ggsave("pheno_5_manhattan.pdf", manhattan_plot, width = 10, height = 5)

# --- Extract Significant Loci ---
sig_loci <- df %>% filter(sig == "sig") %>% select(RS, CHR, BP, P)
write.table(sig_loci, "pheno_5_significant_loci.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- Done ---
cat("QQ plot, Manhattan plot, and significant loci exported.\n")



library(ggplot2)
library(dplyr)

# --- PARAMETERS ---
results_dir <- "./"   # folder containing GEMMA .assoc.txt files
files <- list.files(results_dir, pattern = "\\.assoc\\.txt$", full.names = TRUE)

# --- LOOP over all files ---
for(file in files) {
  
  # Extract phenotype name from file
  pheno_name <- gsub("\\.assoc\\.txt$", "", basename(file))
  
  # --- Load GEMMA result ---
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # --- Rename for standard GWAS plotting ---
  df <- df %>%
    dplyr::rename(
      CHR = chr,
      BP  = ps,
      RS  = rs,
      P   = p_wald
    ) %>%
    dplyr::filter(!is.na(P))
  
  # --- QQ Plot ---
  n <- nrow(df)
  expected <- -log10(ppoints(n))
  observed <- -log10(sort(df$P))
  qq_df <- data.frame(expected, observed)
  
  qqplot <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(x = "Expected -log10(P)", y = "Observed -log10(P)",
         title = paste0("QQ plot: ", pheno_name)) +
    theme_minimal(base_size = 14)
  
  ggsave(paste0(pheno_name, "_qqplot.pdf"), qqplot, width = 5, height = 5)
  
  # --- Manhattan Plot ---
  df$sig <- ifelse(-log10(df$P) > 5, "sig", "ns")  # significance threshold 1e-5
  
  manhattan_plot <- ggplot(df, aes(x = BP, y = -log10(P), color = sig)) +
    geom_point(size = 1) +
    facet_grid(. ~ CHR, scales = "free") +
    scale_color_manual(values = c("ns" = "gray69", "sig" = "red")) +
    ylab(expression(-log[10](italic(P)))) +
    xlab("Genomic position (Mb)") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1, "line")) +
    ggtitle(paste0("Manhattan Plot: ", pheno_name))
  
  ggsave(paste0(pheno_name, "_manhattan.pdf"), manhattan_plot, width = 10, height = 5)
  
  # --- Extract Significant Loci ---
  sig_loci <- df %>% filter(sig == "sig") %>% select(RS, CHR, BP, P)
  write.table(sig_loci, paste0(pheno_name, "_significant_loci.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Processed: ", pheno_name, "\n")
}

cat("All files processed. QQ plots, Manhattan plots, and significant loci exported.\n")


#CMplot
cm_str=read.table("~/Desktop/STR/CM.data2",header = T)
SNPs <-  cm_str[cm_str$Bio5 < 1e-5 |
                  cm_str$Bio10 < 1e-5 |
                  cm_str$Bio4 < 1e-5 |cm_str$Bio7 < 1e-5|cm_str$Bio8 < 1e-5,1]
CMplot(cm_str,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
       threshold=1e-5,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,file="jpg",file.name="new-circos-3",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)


CMplot(cm_str,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),bin.size=1e6,outward=FALSE,file="jpg",file.name="new-circos",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

CMplot(cm_str,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),bin.size=1e6,outward=FALSE,file="jpg",file.name="new-circos-1",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=FALSE)


CMplot(cm_str,type="p",plot.type="c",r=2,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),highlight.cex=0.5,signal.line=1,signal.col=c("red","green"),chr.den.col=c("#d8b365","#f5f5f5","#5ab4ac"),bin.size=1e6,outward=FALSE,file="pdf",file.name="new-circos-1",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=FALSE)


####LD
python maf.py hipstr.vcf.gz bi_allelic_maxmaf_2.vcf.gz
python vcf_refalt_to_one_base.py -i bi_allelic_maxmaf_2.vcf.gz -o bi_allelic_maxmaf_4.vcf
change hipstr result to single base
python vcf_refalt_to_one_base.py -i bi_allelic_maxmaf_2.vcf.gz -o bi_allelic_maxmaf_4.vcf

#chunk bed
bedtools merge -i str.bed > new.str.bed
#get targt vcf then calculate the LD
# chunk hipstr.bed run using child bed
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.xaa.191.bed
bedtools intersect -a ../change.fmt.vcf -b xaa.191.bed -header >xaa.191.bed.snp.vcf
shuangyang.wu/software/PopLDdecay/PopLDdecay -InVCF xaa.191.bed.snp.vcf -OutStat xaa.191.bed.stat.gz -MaxDist 50 -Miss 0.1 -OutType 3


cat  out.xaa.*|awk '$5>0'|sed 's/:/\t/'|cut -f1-5 > SNP-STR.LD &

less -S SNP-STR.LD|awk '{print $0"\t"$1":"$3}' > SNP-STR.LD.1

perl shuangyang.wu/script/merge_files.pl -r SNP-STR.LD.1 -n 6 -i ../../hipstr.id -c 1 -o xxxx
less -S xxxx |awk '!($7>0)'|cut -f1-5 > SNP-STR.LD.final


/scratch-cbe/users/shuangyang.wu/str-1/LD/poplddecay/chunk/reform
less -S SNP-STR.LD.final  ##final result

import pandas as pd
import numpy as np

# 读取数据（假设空格或制表符分隔）
df = pd.read_csv('SNP-STR.LD.final', sep='\s+')

# 显示前几行确认数据读取正确
print("原始数据前几行:")
print(df.head())

# 创建距离区间（1bp间隔）
min_dist = df['Dist'].min()
max_dist = df['Dist'].max()
bins = np.arange(min_dist, max_dist + 2)  # 创建从最小距离到最大距离+1的区间

# 为每个SNP对分配距离区间
df['Dist_Bin'] = pd.cut(df['Dist'], bins=bins, include_lowest=True, right=False)

# 使用区间的左边界作为分组依据
df['Dist_Group'] = df['Dist_Bin'].apply(lambda x: x.left).astype(int)

# 分组聚合计算
grouped = df.groupby('Dist_Group').agg({
  'r^2': ['mean', 'sum', 'count']
}).reset_index()

# 平整化列名
grouped.columns = ['Dist', 'Mean_r^2', 'Sum_r^2', 'NumberPairs']

# 添加D'相关列（全部设为0）
grouped["Mean_D'"] = 0.0
grouped["Sum_D'"] = 0.0

# 重新排列列顺序以匹配PopLDdecay格式
final_df = grouped[['Dist', 'Mean_r^2', "Mean_D'", 'Sum_r^2', "Sum_D'", 'NumberPairs']]

# 格式化数值显示
final_df['Mean_r^2'] = final_df['Mean_r^2'].round(4)
final_df['Sum_r^2'] = final_df['Sum_r^2'].round(4)

# 显示结果
print("\nPopLDdecay格式数据前几行:")
print(final_df.head(10))

# 保存结果
final_df.to_csv('poplddecay_format_output_SNP_STR.txt', sep='\t', index=False)
print("\n结果已保存到 'poplddecay_format_output_SNP_STR.txt'")

# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
setwd("/scratch-cbe/users/shuangyang.wu/str-1/LD/poplddecay/chunk/reform")
# 读取两类LD数据（请替换为您的实际文件路径）
snp_snp_data <- read_tsv("snp-snp.stat") %>%
  mutate(Type = "SNP-SNP")

str_snp_data <- read_tsv("poplddecay_format_output_SNP_STR.txt") %>%
  mutate(Type = "STR-SNP")

# 合并数据
combined_data <- bind_rows(snp_snp_data, str_snp_data)

# 查看数据结构
head(combined_data)

# 创建自定义颜色方案
type_colors <- c("SNP-SNP" = "#E41A1C", "STR-SNP" = "#377EB8")

# 基础LD衰减图（散点+平滑曲线）

# 如果数据点过多，可以使用分箱平均值
binned_data <- combined_data %>%
  mutate(Dist_bin = cut(Dist, breaks = seq(0, max(Dist), by = 1000), include.lowest = TRUE)) %>%
  group_by(Dist_bin, Type) %>%
  summarise(
    Mean_r2 = mean(Mean_r2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(Mid_dist = as.numeric(sub("[\\[\\(](.*),.*", "\\1", Dist_bin)) + 500)

# 使用分箱数据绘图
ggplot(binned_data, aes(x = Mid_dist, y = Mean_r2, color = Type, group = Type)) +
  geom_line(size = 1.2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = type_colors) +
  labs(
    x = "Distance (bp)",
    y = expression(paste("r"^2)),
    color = "LD Type") +
  theme_cust



####

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# 设置工作目录
setwd("/scratch-cbe/users/shuangyang.wu/str-1/LD/poplddecay/chunk/reform")

# 读取两类原始LD数据
snp_snp_data <- read_tsv("snp_snp_raw.txt", col_names = c("chr", "Site1", "Site2", "r2", "Dist")) %>%
  mutate(Type = "SNP-SNP")

str_snp_data <- read_tsv("str_snp_raw.txt", col_names = c("chr", "Site1", "Site2", "r2", "Dist")) %>%
  mutate(Type = "STR-SNP")

# 合并数据
combined_data <- bind_rows(snp_snp_data, str_snp_data)

# 查看数据结构
head(combined_data)

# 创建自定义颜色方案
type_colors <- c("SNP-SNP" = "#E41A1C", "STR-SNP" = "#377EB8")

# 设置距离分箱大小（根据数据范围调整）
bin_size <- 1000  # 1kb分箱

# 创建距离分箱
combined_data <- combined_data %>%
  mutate(
    Dist_bin = cut(Dist, 
                   breaks = seq(0, max(Dist, na.rm = TRUE) + bin_size, by = bin_size),
                   include.lowest = TRUE),
    Mid_dist = as.numeric(sub("[\\[\\(](.*),.*", "\\1", Dist_bin)) + bin_size/2
  )

# 计算每个分箱和类型的均值
mean_data <- combined_data %>%
  group_by(Mid_dist, Type) %>%
  summarise(
    Mean_r2 = mean(r2, na.rm = TRUE),
    .groups = "drop"
  )

# 绘制箱线图并连接均值点
p <- ggplot(combined_data, aes(x = factor(Mid_dist), y = r2, fill = Type)) +
  # 绘制箱线图
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.8) +
  # 添加均值点
  geom_point(data = mean_data, 
             aes(x = factor(Mid_dist), y = Mean_r2, group = Type), 
             color = "black", size = 2, shape = 18) +
  # 连接均值点
  geom_line(data = mean_data, 
            aes(x = factor(Mid_dist), y = Mean_r2, group = Type, color = Type),
            size = 1.2) +
  # 设置颜色
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  # 设置坐标轴标签和标题
  labs(
    x = "Distance (bp)",
    y = expression(paste("r"^2)),
    title = "LD Decay with Distance",
    subtitle = paste("Binned by", bin_size, "bp intervals"),
    fill = "LD Type",
    color = "LD Type"
  ) +
  # 设置主题
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

# 显示图形
print(p)

# 保存图形
ggsave("ld_decay_boxplot_from_raw.png", width = 12, height = 7, dpi = 300)

# 可选：如果x轴标签过于密集，可以只显示部分标签
# 创建一个函数来选择要显示的标签
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

# 应用选择性标签
p + scale_x_discrete(breaks = every_nth(5))

# 保存选择性标签的图形
ggsave("ld_decay_boxplot_selective_labels.png", width = 12, height = 7, dpi = 300)

# 可选：使用分面图更清晰地比较两类LD
p_faceted <- ggplot(combined_data, aes(x = factor(Mid_dist), y = r2, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.8) +
  geom_point(data = mean_data, 
             aes(x = factor(Mid_dist), y = Mean_r2), 
             color = "black", size = 2, shape = 18) +
  geom_line(data = mean_data, 
            aes(x = factor(Mid_dist), y = Mean_r2, group = Type),
            color = "darkblue", size = 1.2) +
  facet_wrap(~ Type, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = type_colors) +
  labs(
    x = "Distance (bp)",
    y = expression(paste("r"^2)),
    title = "LD Decay with Distance by Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  scale_x_discrete(breaks = every_nth(5))

print(p_faceted)
ggsave("ld_decay_boxplot_faceted.png", width = 10, height = 10, dpi = 300)




####update length gwas ###

cm_str=read.table("~/Desktop/STR_Paper/revision/CM.data2.length",header = T)

SNPs <-  cm_str[cm_str$Bio5 < 1e-5 |
                  cm_str$Bio10 < 1e-5 |
                  cm_str$Bio4 < 1e-5 |cm_str$Bio7 < 1e-5|cm_str$Bio8 < 1e-5,1]
CMplot(cm_str,type="p",plot.type="m",LOG10=TRUE,highlight=SNPs,highlight.type="l",
       threshold=1e-5,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),signal.cex=1.2, signal.col="red", highlight.col="grey",highlight.cex=0.7,file="jpg",file.name="new-circos-3",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)


CMplot(cm_str,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),bin.size=1e6,outward=FALSE,file="jpg",file.name="new-circos",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

CMplot(cm_str,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),bin.size=1e6,outward=FALSE,file="jpg",file.name="new-circos-1",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=FALSE)


CMplot(cm_str,type="p",plot.type="c",r=2,col=c("grey30","grey60"),
       threshold=c(1e-5),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red","blue"),highlight.cex=0.5,signal.line=1,signal.col=c("red","green"),chr.den.col=c("#d8b365","#f5f5f5","#5ab4ac"),bin.size=1e6,outward=FALSE,file="pdf",file.name="new-circos-1",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10,cir.axis.grid=FALSE)


CMplot(
  cm_str,
  type = "p",
  plot.type = "c",
  r = 2,
  col = c("grey30", "grey60"),
  threshold = c(1e-5),
  cir.chr.h = 1.5,
  amplify = TRUE,
  threshold.lty = c(1, 2),
  threshold.col = c("red", "blue"),
  highlight.cex = 0.5,
  signal.line = 1,
  signal.col = c("red", "green"),
  chr.den.col = NULL,     # 🔑 remove density plot
  outward = FALSE,
  cir.axis.grid = FALSE,
  file.output = FALSE,   # show in RStudio
  verbose = TRUE
)


library(ggplot2)
library(gggenes)

### 输入文件
gff_file <- "Mp8g08100.gff"          # gff文件
variant_file <- "var.updat.id"     # 变异文件（至少两列：CHR, POS）
target_gene <- "Mp8g08100"         # 目标基因ID

### 读入 GFF
gff <- read.table(gff_file, sep="\t", header=FALSE, quote="")
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# 只保留目标基因相关的 feature
gff_gene <- gff[grep(target_gene, gff$attributes), ]

# 保留 CDS / UTR
conserve <- c("CDS","five_prime_UTR","three_prime_UTR")
gff_gene <- subset(gff_gene, type %in% conserve)

# 提取 gene name
gene_name <- sub(".*Name=([^;]+).*", "\\1", gff_gene$attributes[1])

# 整理数据
gene_df <- data.frame(
  molecule = gff_gene$seqid,
  type = gff_gene$type,
  start = gff_gene$start/1000,   # 转成 Kb
  end = gff_gene$end/1000,
  strand = gff_gene$strand,
  gene = gene_name
)

# strand 转成方向
gene_df$direction <- ifelse(gene_df$strand == "+", 1, -1)

### 读入 variants
variants <- read.table(variant_file, header=TRUE)  # 必须包含 CHR, POS
variants_sub <- subset(variants, CHR == unique(gene_df$molecule) &
                         POS >= min(gff_gene$start) &
                         POS <= max(gff_gene$end))
variants_sub$POS <- variants_sub$POS/1000  # 转 Kb

### 配色
type_color <- c("CDS"="#8dd3c7", "five_prime_UTR"="#ffffb3", "three_prime_UTR"="#bebada")

### 绘图
p <- ggplot(gene_df, aes(xmin=start, xmax=end, y=gene, fill=type, forward=direction)) +
  geom_gene_arrow() +
  scale_fill_manual(values=type_color[unique(gene_df$type)]) +
  geom_vline(data=variants_sub, aes(xintercept=POS), color="red", linetype="dashed", size=0.6) +
  xlab(paste0(unique(gene_df$molecule)," (Kb)")) + ylab("") +
  theme_genes() +
  theme_bw() + 
  theme(plot.title = ggplot2::element_text(size=10,  color = "black"),
        legend.text = ggplot2::element_text(size=10,  color = "black"),
        legend.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.title =  ggplot2::element_text(size=10,  color = "black"),
        axis.text =  ggplot2::element_text(size=10,  color = "black"),
        strip.text = ggplot2::element_text(size=10, vjust = 1,  color = "black"),
        strip.background = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(),
        text = ggplot2::element_text(family="Helvetica"))

print(p)



