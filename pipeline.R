################## 一. WGS —> Genotype ####
#我们的AD全基因组测序数据来自https://www.synapse.org/#!Synapse:syn11724002
#使用plink处理基因型数据，使格式符合后续分析要求。
#plink的安装和使用参考https://www.cog-genomics.org/plink2/
plink --vcf genotype.data.vcf --allow-extra-chr --make-bed --out genotype.data


################## 二. Genotype —> Polygenic risk score, PRS ####

##Disease-related GWAS data(GWAS) + Genotype data(.bim) →→→ Disease risk scores for each individual


##1. QC of Base Data#####

#Quality Control (QC) of GWAS data参考 https://choishingwan.github.io/PRS-Tutorial/
#QCed GWAS data文件可用于下游分析


##2. PRScs分析#####

#PRS for each individual参考 https://github.com/getian107/PRScs

python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE 
--out_dir=OUTPUT_DIR [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR 
                      --chrom=CHROM --beta_std=BETA_STD --seed=SEED]

##所有染色体PRScs计算完成后，合并所有染色体的PRScs结果
cd OUTPUT_DIR
cat *.txt > PRScs.result
plink --bfile genotype.data --score OUTPUT_DIR/PRScs.result 2 4 6 header sum --out OUTPUT_DIR/PRScs.final.result




################## 三. Genotype —> gene expression ####
#使用MetaXcan（https://github.com/hakyimlab/PrediXcan）的Predict.py包
#（https://github.com/hakyimlab/MetaXcan/wiki/Individual-level-PrediXcan:-introduction,-tutorials-and-manual），
#结合利用GTEx预测的snp与gene表达量关系文件（https://zenodo.org/record/3842289#.YrvrM7FBVYA和https://zenodo.org/record/3859065#.YrSgRrFBVYA）（TWAS的前半部分）

python Predict.py --model_db_path PATH_TO_MODEL_DB --vcf_genotypes genotype.data.vcf --vcf_mode genotyped 
--prediction_output prediction_OUTPUT_DIR/predict.txt
--prediction_summary_output prediction_OUTPUT_DIR/summary.txt
--verbosity 9
--throw



#####in R
setwd(work.dir)
library(GSVA)
library(psych)
library(HHG)


################## 四. gene expression —> GSVA ####

##1. 载入基因表达矩阵#####
geneexpression <- read.table("prediction_OUTPUT_DIR/predict.txt",header = TRUE, row.names = 1)

##2. 载入每个通路对应的基因集#####
load("geneSets.Rdata")

##3. 计算#####
gsva_es <- gsva(geneexpression, geneSets, method='gsva', kcdf='Gaussian', abs.ranking=TRUE, mx.diff=TRUE, min.sz=2)





################## 五. Association between disease risk (PRS) and pathways (GSVA) ####

##1. 载入PRS score#####
PRS_result <- read.table("OUTPUT_DIR/PRScs.final.result.profile",header = T,row.names = NULL)
PRS_result1 <- data.frame(SampleID = PRS_result$SampleID,PRScore=PRS_result$SCORESUM)

##2. 载入GSVA score#####
gsva_es_result <- as.data.frame(t(gsva_es)) #样本在行，通路在列
gsva_es_result$SampleID <- rownames(gsva_es_result)

##3. 合并矩阵以计算相关性#####
cor_data <- merge(PRS_result1,gsva_es_result,by="SampleID")

##4. 去除离群值#####
cor_data1 <- apply(cor_data[,-1], 2, function(x){winsor(x)})
cor_data2 <- data.frame(SampleID=cor_data[,1],cor_data1)

##5. 置换计算HHG的经验p值#####
ALLpath_permutation_p = list()

N_Large = SAMPLE.NUMBER #the number of samples
NullTable_for_N_Large_MXL_tables = Fast.independence.test.nulltable(N_Large, variant = 'ADP-EQP',nr.atoms = 40, nr.perm=1000)
#nr.atoms: Brill (2016) suggests a minimum of 40 atoms, with an increase of up to 60 for alternatives which are more difficult to detect (on the expense of computational complexity
#如果通路数量多，推荐使用并行计算，参考R包doParallel
##并行计算
library(doParallel)
# 设置核心数
num_cores <- 108

# 开始并行计算
registerDoParallel(num_cores)

results <- foreach(i = 3:ncol(cor_data2), .combine = rbind) %dopar% {
  permutation_HHG_result = data.frame(matrix(ncol = 5))
  colnames(permutation_HHG_result) = c("pathway","HHG.test.statistic","HHG.pvalue","permutation.cycles.number","permutation.empirical.p")
  permutation_HHG_result$pathway[1] = colnames(cor_data2)[i]
  
  #通路原始HHG p值
  X_Large = cor_data2[,2]
  Y_Large = cor_data2[,i]
  
  set.seed(1234) #每次的运行结果都不同,设置set.seed解决
  ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables) #HHG
  permutation_HHG_result$HHG.test.statistic[1] = ADP_EQP_ML_Result$MinP
  permutation_HHG_result$HHG.pvalue[1] = ADP_EQP_ML_Result$MinP.pvalue
  #MinP: The test statistic when the combining type is "MinP".
  #MinP.pvalue: The p-value when the combining type is "MinP"
  HHG_res_p = ADP_EQP_ML_Result$MinP.pvalue  #原本的p值
  
  
  #置换的p值结果
  result.table = data.frame(matrix(ncol = 3))
  colnames(result.table) = c("pathway","permutation.test.statistic","permutation.p")
  #置换p值permutation.p
  for (k in 1:10000){
    if(k<=1000|nrow(subset(result.table,permutation.p<HHG_res_p))<15){  #每个通路置换1千次以上1万次以下，且小于原pvalue的个数不低于15个
      #随机打乱，以毫秒时间为随机种子。
      set.seed(as.numeric(gsub("[^0-9]", "", format(Sys.time(), "%H:%M:%OS3"))))
      X_Large = sample(cor_data2[,2])
      Y_Large = sample(cor_data2[,i])
      
      set.seed(1234)
      ADP_EQP_ML_Result = Fast.independence.test(X_Large,Y_Large, NullTable_for_N_Large_MXL_tables) #HHG
      
      hhg_permutation_res = data.frame(colnames(cor_data2)[i],ADP_EQP_ML_Result$MinP,ADP_EQP_ML_Result$MinP.pvalue)
      colnames(hhg_permutation_res) = c("pathway","permutation.test.statistic","permutation.p")
      
      result.table = rbind(result.table,hhg_permutation_res)
    }
  }
  result.table = result.table[-1,]
  ALLpath_permutation_p[[colnames(cor_data2)[i]]] = result.table[,3]
  
  ##经验p值empirical.p
  permutation_HHG_result$permutation.cycles.number[1] = nrow(result.table)
  smallerp = subset(result.table,permutation.p<HHG_res_p) #小于原本p值的置换p值的个数
  permutation_HHG_result$permutation.empirical.p[1] = nrow(smallerp)/nrow(result.table)  #除以总的置换p值个数。得到经验p值
  
  return(list(ALLpath_permutation_p, permutation_HHG_result))
}
# 结束并行计算
stopImplicitCluster()

#合并结果
ALLpath_permutation_p <- results[, 1]
for(i in 1:length(ALLpath_permutation_p)){
  names(ALLpath_permutation_p)[[i]] <- names(ALLpath_permutation_p[[i]])
}

permutation_HHG_result1 <- results[, 2]
permutation_HHG_result1 <- do.call(rbind,permutation_HHG_result1)


##计算置换的经验q值
permutation_HHG_result1 = permutation_HHG_result1[order(permutation_HHG_result1$permutation.empirical.p),]
permutation_HHG_result1$permutation.empirical.q = p.adjust(permutation_HHG_result1$permutation.empirical.p,method = "BH")


##6. 计算每个通路对应的p值阈值#####
#找出经验q值最接近0.05时，每个通路对应的经验p值。
permutation_HHG_result2 = subset(permutation_HHG_result1,permutation.empirical.q!=0) #有些通路置换10000次以后置换p值都大于原本的p值
empirical_p = permutation_HHG_result2[which.min(abs(permutation_HHG_result2$permutation.empirical.q - 0.05)),"permutation.empirical.p"]
empirical_p = as.numeric(empirical_p)

#然后用这个经验p值去找置换p值，如果是经验p值是0.09，就是找1000个置换p值里面的第90个（从小往大排序)。
#找出的这个置换p值就是基因1的p阈值。
permutation_HHG_result1$threshold.p = ""
for (m in 1:nrow(permutation_HHG_result1)) {
  cycles_number = as.numeric(permutation_HHG_result1$permutation.cycles.number[m])
  order_p = ceiling(empirical_p*cycles_number) #置换p值的第..个
  #置换p值排序
  order_permutation_p = ALLpath_permutation_p[[permutation_HHG_result1$pathway[m]]][[1]]
  order_permutation_p1 = order_permutation_p[order(order_permutation_p)]
  permutation_HHG_result1$threshold.p[m] = order_permutation_p1[order_p]  ##第..个对应的置换p值
}
sig_res = subset(permutation_HHG_result1,HHG.pvalue<=threshold.p)






################## 六. Structural equation modeling (SEM) analysis ####

#使用通路的GSVA score的相关性矩阵，利用SPSS和lisrel构建结构方程模型



