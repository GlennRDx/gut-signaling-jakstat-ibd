# DESeq2 DGE datasets
dsq_eec = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_EEC.csv")
dsq_eep = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_EEC_Progenitor.csv")
dsq_ent = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Enterocyte.csv")
dsq_enp = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Enterocyte_Progenitor.csv")
dsq_gob = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Goblet.csv")
dsq_gop = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Goblet_Progenitor.csv")
dsq_isc = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_ISC.csv")
dsq_pan = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Paneth.csv")
dsq_pap = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Paneth_Progenitor.csv")
dsq_tuf = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Tuft.csv")
dsq_tup = read.csv("/home/glennrdx/Documents/Research_Project/scRNAseq-MSc-Analysis/3. upstream_analysis/crypt/differential_expression_deseq2/diff_exp_CD_vs_HFD_Tuft_Progenitor.csv")

dsq_list = list(dsq_isc = dsq_isc,
               dsq_enp = dsq_enp,
               dsq_ent = dsq_ent,
               dsq_gob = dsq_gob,
               dsq_tuf = dsq_tuf,
               dsq_eec = dsq_eec)