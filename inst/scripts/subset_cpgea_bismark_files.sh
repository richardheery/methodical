wget https://download.cncb.ac.cn/gsa-human/HRA000099/processed/T100_WGBS_mC_Identification.stat.gz
wget https://download.cncb.ac.cn/gsa-human/HRA000099/processed/N100_WGBS_mC_Identification.stat.gz
wget https://download.cncb.ac.cn/gsa-human/HRA000099/processed/N103_WGBS_mC_Identification.stat.gz
wget https://download.cncb.ac.cn/gsa-human/HRA000099/processed/T103_WGBS_mC_Identification.stat.gz

zcat N100_WGBS_mC_Identification.stat.gz | awk 'BEGIN{OFS="\t"} $1 == 11 && $2 >= 67578812 && $2 <= 67588812 {$1="chr"$1; print}' | gzip > N100_WGBS_mC_Identification_GSTP1.CX_report.txt.gz
zcat T100_WGBS_mC_Identification.stat.gz | awk 'BEGIN{OFS="\t"} $1 == 11 && $2 >= 67578812 && $2 <= 67588812 {$1="chr"$1; print}' | gzip > T100_WGBS_mC_Identification_GSTP1.CX_report.txt.gz
zcat N103_WGBS_mC_Identification.stat.gz | awk 'BEGIN{OFS="\t"} $1 == 11 && $2 >= 67578812 && $2 <= 67588812 {$1="chr"$1; print}' | gzip > N103_WGBS_mC_Identification_GSTP1.CX_report.txt.gz
zcat T103_WGBS_mC_Identification.stat.gz | awk 'BEGIN{OFS="\t"} $1 == 11 && $2 >= 67578812 && $2 <= 67588812 {$1="chr"$1; print}' | gzip > T103_WGBS_mC_Identification_GSTP1.CX_report.txt.gz

