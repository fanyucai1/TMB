import re
import sys
import subprocess

cosmic_genelist="Cancer_Gene_Census.tsv"##download cosmic https://cancer.sanger.ac.uk/census#cl_search
OncoKB_Cancer_Gene_List="cancerGeneList.txt"##download https://oncokb.org/cancerGenes
cosmic_anno="/data/Database/COSMIC/release_v88/CosmicMutantExport.tsv"##download cosmic
cosmic_vcf="/data/Database/COSMIC/release_v88/CosmicCodingMuts.vcf"##download cosmic
annovar="/software/docker_tumor_base/Resource/Annovar/"
dbsnp_germline="/data/Database/hg19/dbsnp/germline.vcf"#(SAO=1)
common_snp="/data/Database/hg19/dbsnp/dbsnp.common.vcf"#(COMMON=1)
database = ['1000g2015aug_all', 'ExAC_ALL', 'esp6500siv2_all','genome_AF','exome_AF','exome_AF_popmax','genome_AF_popmax',
            'ExAC_AFR','ExAC_AMR','ExAC_EAS','ExAC_FIN','ExAC_NFE','ExAC_OTH','ExAC_SAS']
def run(vcf):
    germline = {}
    #######################################format vcf identify the VAF
    infile = open(vcf, "r")
    outfile = open("tmp.vcf", "w")
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for line in infile:
        line = line.strip()
        if not line.startswith("#"):
            array = line.split("\t")
            tmp = array[-2].split(":")
            info = array[-1].split(":")
            a, c =array[4].split(","), []
            for k in range(len(tmp)):
                if tmp[k] == "AF" or tmp[k] == "VAF":
                    c = info[k].split(",")  # AF
                else:
                    pass
            if len(a) == 1:
                outfile.write("%s\t%s\t%s\t%s\t%s\t.\t.\tVAF=%s"
                              % (array[0], array[1], array[2], array[3], array[4], c[0]))
                outfile.write("\n")
            else:
                for i in range(len(a)):
                    ALT = a[i]
                    outfile.write("%s\t%s\t%s\t%s\t%s\t.\t.\tVAF=%s"
                                  % (array[0], array[1], array[2], array[3], ALT, c[i]))
                    outfile.write("\n")
    infile.close()
    outfile.close()
    #############################step1:Known germline alterations in dbSNP(SAO=1) and common snp
    infile=open(dbsnp_germline,"r")
    for line in infile:
        if not line.startswith("#"):
            line = line.strip()
            array = line.split("\t")
            if re.search(r',',array[4]):
                newsite=array[4].split(",")
                for i in newsite:
                    tmp =array[0] + "_" + array[1] + "_" + array[3] + "_" +i
                    germline[tmp]=1
            else:
                tmp = array[0] + "_" + array[1] + "_" + array[3] + "_" + array[4]
                germline[tmp] = 1
    infile.close()
    infile = open(common_snp, "r")
    for line in infile:
        if not line.startswith("#"):
            line = line.strip()
            array = line.split("\t")
            if re.search(r',',array[4]):
                newsite=array[4].split(",")
                for i in newsite:
                    tmp =array[0] + "_" + array[1] + "_" + array[3] + "_" +i
                    germline[tmp]=1
            else:
                tmp = array[0] + "_" + array[1] + "_" + array[3] + "_" + array[4]
                germline[tmp] = 1
    infile.close()
    ###########################step2:identify the tumor suppressor genes
    infile=open(cosmic_genelist,"r")
    TSG,num,name={},0,{}
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        num+=1
        if num==1:
            for i in range(len(array)):
                name[array[i]]=i
        else:
            if array[name['Role in Cancer']]=="TSG":
                TSG[array[0]]=1
    infile.close()
    infile=open(OncoKB_Cancer_Gene_List,"r")
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        if array[5]=="Yes":
            TSG[array[0]]=1
    infile.close()
    ########################step3:identify the somatic site in cosmic
    cosmicID,num,id,status,somatic={},0,0,0,{}
    infile=open(cosmic_anno,"r")
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        num+=1
        if num==1:
            for i in range(len(array)):
                if array[i]=="Mutation ID":
                    id=i
                if array[i]=="Mutation somatic status":
                    status=i
        else:
            if array[status]=="Confirmed somatic variant":
                cosmicID[array[id]]=1
    infile.close()
    infile=open(cosmic_vcf,"r")
    for line in infile:
        if not line.startswith("#"):
            line = line.strip()
            array = line.split("\t")
            if array[2] in cosmicID:
                tmp="chr"+array[0]+"_"+array[1]+"_"+array[3]+"_"+array[4]
                somatic[tmp]=1
    infile.close()
    ###########################step4:filter somatic site in cosmic and known germline site in dbsnp
    infile=open(vcf,"r")
    outfile=open("tmp2.vcf","w")
    for line in infile:
        line = line.strip()
        array = line.split("\t")
        if line.startswith("#"):
            outfile.write("%s\n"%(line))
        else:
            tmp = array[0] + "_" + array[1] + "_" + array[3] + "_" + array[4]
            if not tmp in somatic:
                if not tmp in germline:
                    outfile.write("%s\n" % (line))
    outfile.close()
    ###########################step5:annotate vcf using annovar and delete stop-gain mutations in tumor suppressor genes
    par = " -protocol refGene,cytoBand,exac03,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eas,gnomad211_exome,gnomad211_genome,cosmic88_coding,clinvar_20190305,intervar_20180118"
    par += " -operation g,r,f,f,f,f,f,f,f,f,f "
    par += " -nastring . -polish "
    subprocess.check_call("perl %s/table_annovar.pl tmp2.vcf %s/humandb -buildver hg19 -out TMB -remove %s -vcfinput " % (annovar, annovar, par), shell=True)
    ###########################step6:counts TMB
    TMB,TMB_50,num,AF=0,0,0,[]
    infile=open("TMB.hg19_multianno.txt","r")
    for line in infile:
        line=line.strip()
        array=line.split("\t")
        num+=1
        if num==1:
            for i in range(len(array)):
                if array[i] in database:
                    AF.append(i)
        else:
            result="True"
            for i in AF:
                if array[i]!="." and float(array[i])>=0.01:
                    result = "false"
            if array[6] in TSG and array[8]=="stopgain":
                result="false"
            if not re.search(r'exonic', array[5]):
                result = "false"
            pattern=re.compile(r'VAF=([0-9.]+)')
            var=pattern.findall(line)
            if result=="True":
                TMB+=1
                if float(var[0])<0.5:
                    TMB_50+=1
    infile.close()
    print("TMB is %s and TMB_50(vaf<0.5) is %s"%(TMB,TMB_50))

if __name__=="__main__":
    if len(sys.argv)!=2:
        print("Usage:python3 %s input.vcf"%(sys.argv[0]))
        print("Email:fanyucai1@126.com")
    else:
        vcf=sys.argv[1]
        run(vcf)
