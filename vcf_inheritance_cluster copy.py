import allel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.options.mode.chained_assignment = None
import matplotlib.patches as mpatches
from sklearn.cluster import DBSCAN
import os
import argparse
import logging
from datetime import datetime
import csv

NOTE={"SVTYPE":"Type of structural variant",\
      "END":"End position of the variant described in this record",\
      "caller":"caller",\
      "mapper":"mapper",\
      "center":"center of SV",\
      "signal":"half of length",\
      "RE":"Number of read support this record",\
      "SVLEN":"Difference in length between REF and ALT alleles",
      "GT":"Genotype"}
CHROMOSOME=[str(x) for x in range(1,30)]+["X","MT"]
MAPPER = ["minimap2","ngmlr","DEL","INS","DUP","INV"]
LABELS={
    "cuteSV":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "sniffles":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "svim":["CHROM","POS","ID","REF","variants/ALT","END","SVTYPE","SVLEN","SUPPORT","samples","GT"],
    "nanoSV":["CHROM","POS","ID","REF","variants/ALT","END","SVTYPE","SVLEN","DV","samples","GT"],
    "overlap":["CHROM","POS","ID","REF","ALT_1","END","SVTYPE","SVLEN","GT","length"],
    "noisy":["CHROM","POS","ID","REF","ALT_1","END","SVTYPE","SVLEN","GT","length"]
}

ALTS={
    "cuteSV": {"INS":"INS","DEL":"DEL","DUP":"DUP","INV":"UNO","BND":"BND"},
    "sniffles": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","INVDUP":"DUP","TRA":"UNO"},
    "svim": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","DUP:TANDEM":"DUP","DUP:INT":"DUP","BND":"BND"},
    "nanoSV": {"DEL":"DEL","DUP":"DUP","BND":"BND","INS":"INS"},
    "overlap":{"DEL":"DEL","DUP":"DUP","INV":"INV","INS":"INS"},
    "noisy":{"DEL":"DEL","DUP":"DUP","INV":"INV","INS":"INS"}
}

class VcfGeno:

    def __init__(self,mapper,caller,vcf):
        self.mapper = mapper
        self.caller = caller
        self.vcf = vcf
        self.subvcf = vcf_get_loc(self.vcf)

    def get_subvcf(self):
        ref = ALTS[self.caller].keys()
        for sv in ref:
            self.subvcf[sv] = vcf_get_loc(self.vcf[self.vcf["SVTYPE"]==sv])

def vcf_get_loc(vcf):
    try:
        vcf["signal"] = abs(vcf["SVLEN"])//2
    except:
        vcf["signal"] = vcf["length"] //2
    vcf["center"] = vcf["POS"] + vcf["signal"]
    return vcf

class Caster:

    def __init__(self,combo,svs,mvcf,ovcf):
        self.combo = combo
        self.svs = svs
        self.vcf = {}
        self.match_record = None
        for chr in CHROMOSOME:
            msub = mvcf[mvcf["CHROM"]==chr]
            osub = ovcf[ovcf["CHROM"]==chr]
            for sv in self.svs:
                key = "{}_{}".format(chr,sv)
                v1 = msub[msub["SVTYPE"] == sv]
                v2 = osub[osub["SVTYPE"] == sv]
                self.vcf[key] = pd.concat([v1,v2],axis=0)
        self.overlap = {}
        self.precise = pd.DataFrame()
        self.precise50 = pd.DataFrame()
        self.valided = pd.DataFrame()

    def cast(self,size=0,eps=10):
        all_true_shoot = 0
        all_false_shoot = 0
        all_miss_shoot = 0
        all_true_call = 0
        all_false_call = 0
        output = [] #pd.DataFrame(columns=["key","true","false","all","precise"])
        match_records = []
        checked_groups = []
        valided = []
        for chr in CHROMOSOME:
            for sv in self.svs:
                true_shoot = 0
                false_shoot = 0
                miss_shoot = 0
                true_call = 0
                false_call = 0
                key = "{}_{}".format(chr,sv)
                pts = get_pts(self.vcf[key])
                if len(pts) != 0:
                    groups = DBSCAN(eps=eps,min_samples=2).fit_predict(pts)
                    df = pd.DataFrame({"group":groups})
                    #print(df)

                    self.vcf[key]["group"] = df.loc[:,"group"].values
                    #print(self.vcf[key]["GT"])
                    mask = self.vcf[key]["GT"].apply(lambda x: x != [-1,-1])
                    self.vcf[key]["match"] = "False"
                    self.vcf[key].loc[self.vcf[key]["group"] != -1,"match"] = "True"
                    match_record = self.vcf[key][["CHROM","POS","signal","match","label","caller","mapper","SVTYPE"]]
                    match_records.append(match_record)
                    self.overlap[key] = self.vcf[key][mask]


                    mother = self.overlap[key][
                        (self.overlap[key]["label"] == "mother") & (self.overlap[key]["signal"] >= size)]
                    off = self.overlap[key][
                        (self.overlap[key]["label"] == "offspring") & (self.overlap[key]["signal"] >= size)]
                    gtmask = mother["GT"].apply(lambda x: x == [1,1])
                    mref = mother[gtmask]
                    gtmask = off["GT"].apply(lambda x: x == [1, 1])
                    oref = off[gtmask]
                    if len(mref) == 0 or len(oref) == 0:
                        continue
                    mgroups = {x for x in mref["group"].tolist() if x != -1}
                    ogroups = {x for x in oref["group"].tolist()if x != -1}
                    print("mother: {}; child: {}".format(len(mgroups),len(ogroups)))
                    for group in mgroups:
                        gts=[]
                        for x in mref[mref["group"] == group]["GT"].tolist():
                            if x not in gts:
                                gts.append(x)
                        suboff = off[off["group"] == group]["GT"].tolist()
                        #print("the size : off: {}, mother: {}".format(len(suboff),len(mother[mother["group"] == group]["GT"].tolist())))
                        if [1, 1] in gts and [0,0] not in gts:
                            false_call += suboff.count([0,0])
                            true_call += suboff.count([1,1])
                            if len(suboff) == 0:
                                miss_shoot += 1
                            elif [0, 0] in suboff:
                                false_shoot += 1
                            else:
                                valided.append(off[off["group"] == group])
                                valided.append(mref[mref["group"] == group])
                                true_shoot += 1

                    if true_shoot != 0 or false_shoot != 0 or miss_shoot != 0:
                        tmp_all = true_shoot + false_shoot + miss_shoot
                        tmp_precise = true_shoot / (true_shoot + false_shoot + miss_shoot)

                        print(
                            "currently the precise for chr{} in mother is {}, {} record missed, true call: {}, false call: {}.".format(
                                key, tmp_precise, miss_shoot, true_call, false_call))

                    for group in ogroups:
                        gts = []
                        for x in oref[oref["group"] == group]["GT"].tolist():
                            if x not in gts:
                                gts.append(x)
                        submother = mother[mother["group"] == group]["GT"].tolist()
                        if [1, 1] in gts and [0,0] not in gts:
                            false_call += submother.count([0,0])
                            true_call += submother.count([1,1])
                            if len(submother) == 0:
                                miss_shoot += 1
                            elif [0, 0] in submother:
                                false_shoot += 1
                            else:
                                valided.append(oref[oref["group"] == group])
                                valided.append(mother[mother["group"] == group])
                                true_shoot += 1

                if true_shoot != 0 or false_shoot != 0 or miss_shoot != 0:
                    tmp_all = true_shoot+false_shoot+miss_shoot
                    tmp_precise = true_shoot/(true_shoot+false_shoot+miss_shoot)
                    all_false_shoot += false_shoot
                    all_true_shoot += true_shoot
                    all_miss_shoot += miss_shoot
                    all_true_call += true_call
                    all_false_call += false_call

                    output.append([key,true_shoot,false_shoot,tmp_all,tmp_precise,miss_shoot,true_call,false_call])
                    print("currently the precise for chr{} is {}, {} record missed, true call: {}, false call: {}.".format(key,tmp_precise,miss_shoot,true_call,false_call))
        self.valided = pd.concat(valided,axis=0)
        all_gt = self.valided["GT"]
        all_gt = ["/".join([str(x[0]),str(x[1])]) for x in all_gt]
        self.valided["GT"] = all_gt
        self.valided.drop_duplicates(inplace=True)
        self.match_record = pd.concat(match_records,axis=0)
        presice = 0
        if all_false_call != 0 or all_true_call != 0 or all_miss_shoot != 0:
            presice = all_true_shoot/(all_true_shoot+all_false_shoot+all_miss_shoot)
        all = all_true_shoot+all_false_shoot+all_miss_shoot
        for sv in self.svs:
            all_sv = [sv,0,0,0,0,0,0,0]
            for record in output:
                if sv in record[0]:
                    all_sv[1] += record[1]
                    all_sv[2] += record[2]
                    all_sv[5] += record[5]
                    all_sv[6] += record[6]
                    all_sv[7] += record[7]
            all_sv[3] = all_sv[1] + all_sv[2] + all_sv[4]
            if all_sv[3] != 0:
                all_sv[4] = all_sv[1]/all_sv[3]
            output.append(all_sv)
        output.append(["all",all_true_shoot,all_false_shoot,all,presice,all_miss_shoot,all_true_call,all_false_call])
        print(len(output))
        print("get pricise: ",presice)
        self.precise=pd.DataFrame(np.array(output),columns=["key","true","false","all","precise","missed","true_call","false_call"])

    def output(self,locat,sample,eps=10,size=0):
        out = []
        for chr in CHROMOSOME:
            for sv in self.svs:
                key = "{}_{}".format(chr,sv)
                if key in self.overlap.keys():
                    out.append(self.overlap[key])
        final = pd.concat(out,axis=0)
        allel.write_vcf("{}/{}/vcf/{}_valided.vcf".format(locat,sample,self.combo),self.valided,description=NOTE)
        allel.write_vcf("{}/{}/vcf/{}_overlap_geno.vcf".format(locat,sample,self.combo), final,description=NOTE)
        self.precise.to_csv("{}/{}/genotype/{}_precise.txt".format(locat,sample,self.combo), sep="\t")
        #allel.write_vcf("./{}/overlap_geno_50.vcf".format(locat), final)
        #self.precise50.to_csv("{}/genotype/{}_{}_precise_eps{}_50.txt".format(locat,sample,self.combo,eps),sep="\t")
        self.match_record.to_csv("{}/{}/match_record/{}_match_record.txt".format(locat,sample,self.combo),sep="\t")
        return final


def get_pts(vcf):

    pts = vcf[["center","signal"]].values
    return pts



def main():
    parser = argparse.ArgumentParser(description="test")
    req_grp = parser.add_argument_group(title='Required')
    #req_grp.add_argument('-i', '--input', type=str, help="Input VCF folder.", required=True)
    req_grp.add_argument('-1','--mother',type=str,help="Input the path of mother.",required=True)
    req_grp.add_argument('-2','--offspring',type=str,help="Input the path of offspring.",required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    req_grp.add_argument('-s', '--threshold', type=int, help="filter threshold", default=5)
    req_grp.add_argument('-n', '--sample',type=str,help="input the name of group",default="sample")
    req_grp.add_argument('-e', '--eps', type=int, help="clustering eps", default=10)
    req_grp.add_argument('-m', '--min_length', type=int, help="min sv length", default=0)
    #    parser.add_argument('-l', '--log', type=str, help="Log file name")
    mother_names = []
    offspring_names = []
    args = parser.parse_args()
    eps = int(args.eps)
    size = args.min_length
    if size == 0:
        name_size = ""
    else:
        name_size = "MinLen_"+str(size)
    if args.output[0] == "/":
        locat = "/"+args.output.strip('/')
    else:
        locat = args.output.strip('/')
    sample = "{}_eps{}{}".format(args.sample,str(eps),str(name_size))
    #vcf_path = args.input.strip('/') + '/'
    mother_path = args.mother
    offspring_path = args.offspring

    vcfs_mother = [x for x in os.listdir(mother_path) if '.vcf' in x]
    for vcf_file in vcfs_mother:
        print("Now read {}".format(vcf_file))
        full_path = mother_path + vcf_file
        basename = vcf_file.split(".")[0]
        caller = ""
        mapper = "basic"
        for i in LABELS.keys():
            if i in basename:
                caller = i
        for m in MAPPER:
            if m in basename:
                mapper = m
        print("mapper: {}, caller: {}".format(mapper, caller))
        mother_names.append((full_path, caller, mapper))

    vcfs_offspring = [x for x in os.listdir(offspring_path) if '.vcf' in x]
    for vcf_file in vcfs_offspring:
        print("Now read {}".format(vcf_file))
        full_path = offspring_path + vcf_file
        basename = vcf_file.split(".")[0]
        caller = ""
        mapper = "basic"
        for i in LABELS.keys():
            if i in basename:
                caller = i
        for m in MAPPER:
            if m in basename:
                mapper = m
        print("mapper: {}, caller: {}".format(mapper, caller))
        offspring_names.append((full_path, caller, mapper))

    if os.path.isdir(locat) is False:

        try:
            os.system("mkdir -p {}".format(locat))
            os.system("mkdir -p {}/{}".format(locat,sample))
            os.system("mkdir -p {}/{}/vcf".format(locat, sample))
            os.system("mkdir -p {}/{}/genotype".format(locat, sample))
            os.system("mkdir -p {}/{}/match_record".format(locat, sample))

        except:
            print("folder exist.")
    else:
        os.system("mkdir -p {}/{}".format(locat, sample))
        os.system("mkdir -p {}/{}/vcf".format(locat, sample))
        os.system("mkdir -p {}/{}/genotype".format(locat, sample))
        os.system("mkdir -p {}/{}/match_record".format(locat, sample))

    mothers = {}
    mother_combos=[]
    for path, caller, mapper in mother_names:
        print("Now read mother_{}_{}".format(mapper,caller))
        combo = "{}_{}".format(mapper,caller)
        vcf = allel.vcf_to_dataframe(path,fields=LABELS[caller])
        if caller not in ["overlap","noisy"]:
            geno = [x[0] for x in allel.read_vcf(path,fields=["GT"])["calldata/GT"].tolist()]
            vcf["ALT"] = vcf["ALT_1"]
            vcf["GT"] = geno
        else:
            geno = [x.split("/") for x in vcf["GT"].tolist()]
            geno = [[int(x[0]), int(x[1])] for x in geno]
            vcf["GT"] = geno
            vcf["ALT"] = vcf["ALT_1"]
            print("current reading {}".format(combo))
            #print(vcf["GT"])
        vcf["label"]="mother"

        vcf["mapper"] = mapper
        vcf["caller"] = caller
        if caller not in ["overlap","noisy","nanoSV"]:
            vcf.rename(columns={LABELS[caller][8]: "RE"}, inplace=True)
            vcf = vcf[vcf["RE"] >= args.threshold]
        mother_combos.append(combo)
        mothers[combo] = VcfGeno(mapper,caller,vcf)

    offsprings={}
    offspring_combos=[]
    for path, caller, mapper in offspring_names:
        print("Now read offspring_{}_{}".format(mapper, caller))
        combo = "{}_{}".format(mapper, caller)
        vcf = allel.vcf_to_dataframe(path, fields=LABELS[caller])
        if caller not in ["overlap", "noisy"]:
            geno = [x[0] for x in allel.read_vcf(path, fields=["GT"])["calldata/GT"].tolist()]
            vcf["GT"] = geno
            vcf["ALT"] = vcf["ALT_1"]
        else:
            geno = [x.split("/") for x in vcf["GT"].tolist()]
            geno = [[int(x[0]),int(x[1])] for x in geno]
            vcf["GT"] = geno
            vcf["ALT"] = vcf["ALT_1"]
        vcf["label"]="offspring"
        vcf["mapper"] = mapper
        vcf["caller"] = caller
        if caller not in ["overlap","noisy","nanoSV"]:
            vcf.rename(columns={LABELS[caller][8]: "RE"}, inplace=True)
            vcf = vcf[vcf["RE"] >= args.threshold]
        offspring_combos.append(combo)
        offsprings[combo] = VcfGeno(mapper, caller, vcf)

    if len(mother_combos) != len(offspring_combos):
        print("ERROR: unbalanced combos")
        exit()

    else:
        for combo in mother_combos:
            print("Now aligning {}".format(combo))
            mother = mothers[combo].subvcf
            offspring = offsprings[combo].subvcf
            caller = combo.split("_")[1]
            caster = Caster(combo,ALTS[caller].keys(),mother,offspring)
            caster.cast(size=size,eps=eps)
            overlap = caster.output(locat,sample,eps=eps,size=size)
            if overlap is not None:
                print(overlap.index)

if __name__ == "__main__":
    main()