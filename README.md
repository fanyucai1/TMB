##   TMB：

    TMB的定义如下：TMB was defined as the number of somatic, coding, base substitution, and indel mutations per megabase of genome examined. 可参考文献[1]
    
    panel size的问题对于TMB的计算很重要，不能太小，最小0.5吧，也有建议0.8以上，不过总的结论是不能太小,越大越接近外显子的结果.
    
    panel的设计初衷有时候不仅仅是call变异，也有SV、基因融合等等，这些区域在计算TMB的时候不应该考虑在内,对于变异位点，建议将同义突变也考虑进去。
    
    在测序深度上一般的panel是500X，为了较准确的发现低频突变，最低也不要低于200X。
    
    在call变异过程中VAF FoundationOne CDx and Oncomine assays (https://www.accessdata.fda.gov/cdrh_docs/pdf17/P170019B.pdf)与MSK-IMPACT都是大于等于5%，但是热点突变是2%[2-4]

    关于热点突变（文献里也指驱动突变）建议还是要去掉的，貌似不是去掉整个驱动基因而是驱动突变或者热点突变,关于热点突变有的文献粗略的定义cosmic数据库中CNT>=50的位点，这些位点你可以通过下载cosmic数据库的VCF获得[5]

    无论是正常样本还是肿瘤样本，germline位点一般VAF都集中在50%与100%附近，一般情况下在肿瘤样本中46%的germline位点的VAF位于40%-60%之间，35%的germline位点的VAF大于95%，约有19%的后选位点不在这个范围内（依据AMP/ASCO/CAP guidelines）[6]
    
    另外变异位点也是考虑体细胞突变，如何区分体细胞和遗传性突变，要么你就是配对样本，如果你是tumor only那么可以参考Foundation Medicine’s FoundationOne panel 开发的SGZ（https://github.com/jsunfmi/SGZ），不过是需要一定的样本作为训练集合的[5,7]
 
    TMB有明确证据证明在 non-small cell lung carcinoma (NSCLC)是biomarker，但并不是在所有癌种中都适用,泛癌种的TMB表现可以参考[9-10]

##  参考文献

1.  Chalmers Z R, Connelly C F, Fabrizio D, et al. Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden[J]. Genome medicine, 2017, 9(1): 34.

2.  Buchhalter I, Rempel E, Endris V, et al. Size matters: Dissecting key parameters for panel‐based tumor mutational burden analysis[J]. International journal of cancer, 2019, 144(4): 848-858.

3.  Büttner R, Longshore J W, López-Ríos F, et al. Implementing TMB measurement in clinical practice: considerations on assay requirements[J]. ESMO open, 2019, 4(1): e000442.

4.  Zehir A, Benayed R, Shah R H, et al. Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients[J]. Nature medicine, 2017, 23(6): 703.

5.  Wood D E, White J R, Georgiadis A, et al. A machine learning approach for somatic mutation discovery[J]. Science translational medicine, 2018, 10(457): eaar7939.

6.  Montgomery N D, Selitsky S R, Patel N M, et al. Identification of Germline Variants in Tumor Genomic Sequencing Analysis[J]. The Journal of molecular diagnostics: JMD, 2018, 20(1): 123-125.

7.  Sun J X, He Y, Sanford E, et al. A computational approach to distinguish somatic vs. germline origin of genomic alterations from deep sequencing of cancer specimens without a matched normal[J]. PLoS computational biology, 2018, 14(2): e1005965.

8.  Xu Z, Dai J, Wang D, et al. Assessment of tumor mutation burden calculation from gene panel sequencing data[J]. OncoTargets and therapy, 2019, 12: 3401.

9.  Fancello L, Gandini S, Pelicci P G, et al. Tumor mutational burden quantification from targeted gene panels: major advancements and challenges[J]. Journal for immunotherapy of cancer, 2019, 7(1): 183.

10.  Samstein R M, Lee C H, Shoushtari A N, et al. Tumor mutational load predicts survival after immunotherapy across multiple cancer types[J]. Nature genetics, 2019, 51(2): 202-206.


##  TMB.py

    依据文献[1,8]也就是FoundationOne CDx (F1CDx)计算TMB，主要是要考虑同义突变，防止样本噪音

### 识别癌症基因（oncogenes）与抑癌基因（TSG，suppressor genes）

    根据文献[1,8]stop-gain mutations in tumor suppressor gene 应该去掉，以下两个数据库有基因的分类

    Oncokb database:    https://oncokb.org/cancerGenes

    cosmic database:    https://cancer.sanger.ac.uk/census#cl_search

### 识别热点突变

    这里的热点可以定义在cosmic数据库中已经被认定为somatic的突变，在cosmic数据库中CosmicMutantExport.tsv包含了这样的信息

### 识别dbsnp已知的germline和common 位点

    下载最新的dbsnp
    ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/
    这里的germline位点（SAO=1）和common位点（COMMON=1）,这里common SNP is one that has at least one 1000Genomes population with a MAF >= 1% and for which 2 or more founders contribute to that minor allele frequency

### 人群频率
    
    依据文献[1]也就是FoundationOne CDx (F1CDx)计算TMB的方法是去掉了two or more counts in the ExAC database were not counted，这里我们选取的人群频率数据库阈值VAF大于1%的过滤掉
    需要说明的是有些人群频率数据库还包含一些人群子库，建议在过滤的时候只要有1个子库大于1%就过滤掉该位点
    
### SGZ[7]

    这是FoundationOne CDx (F1CDx)针对自己的panel设计的针对单Tumor only来识别somatic和germline,由于没有公开数据的原因我们很难使用该软件。
    但是有一点是清楚的，该软件主要还是针对VAF频率高于50%的位点是否是somatic位点，在相关文献和资料中提到，在肿瘤纯度大于50%的样本中是有可能存在的，这样的样本一般呈现的是TMB high。
    所以我们这里建议如果在你的最终结果中去掉频率大于50%的位点，如果TMB仍大于20，则建议保留这些位点，如果小于20则建议去掉这些位点。特别说明本脚本计算的是TMB对应的变异位点数，需要最后再除以你的code区域（或者coding区域覆盖度大于50X的区域）



