import allel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
from sklearn.cluster import DBSCAN
import os
import argparse
import logging
from datetime import datetime
import csv


LINESTYLE = {"DUP":"solid","DEL":"dashed","INS":"dashdot"}
COLORS={}
CALLERS=["cuteSV","nanoSV","sniffles","svim"]
MAPPERS=["minimap2","ngmlr"]
cnames = ["darkblue","darkgray","darkred","darksalmon","lightblue","lightgray","red","lightsalmon"]
for mapper in MAPPERS:
    for caller in CALLERS:
        COLORS["{}_{}".format(caller,mapper)] = cnames.pop()
scale=[x*10000000 for x in [16,14,12,12,12,12,11,11,10.5,10.5,11,9,8.5,8.5,8.5,8,7,6.5,6.5,7,7,6,5,6.2,4.7,5,4.5,4.5,5,14]]
scale.append(16338)
CHROMOSOME=[str(x) for x in range(1,30)]+["X","MT"]
SCALE_CHR={}
for i in range(len(CHROMOSOME)):
    SCALE_CHR[CHROMOSOME[i]]=scale[i]
LABELS={
    "cuteSV":["CHROM","POS","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "sniffles":["CHROM","POS","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "svim":["CHROM","POS","ID","END","SVTYPE","SVLEN","SUPPORT","samples","GT"],
    "nanoSV":["CHROM","POS","ID","END","SVTYPE","SVLEN","DV","samples","GT"]
}

ALTS={
    "cuteSV": {"INS":"INS","DEL":"DEL","DUP":"DUP","INV":"UNO","BND":"BND"},
    "sniffles": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","INVDUP":"DUP","TRA":"UNO"},
    "svim": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","DUP:TANDEM":"DUP","DUP:INT":"DUP","BND":"BND"},
    "nanoSV": {"DEL":"DEL","DUP":"DUP","BND":"BND","INS":"INS"}
}
CHRS=[""]+["chromosome_{}".format(x) for x in range(1,29)]
DEFAULT=["DEL","INS","DUP"]


class VcfRecord:

    def __init__(self, caller, mapper, vcf: pd.DataFrame):

        self._mapper = mapper
        self._caller = caller
        self._sample_name = "{}_{}".format(mapper,caller)
        self.vcf = vcf
        self.subvcf = None

    def get_caller(self):

        return self._caller

    def get_mapper(self):

        return self._mapper

    def get_combo(self):

        combo = "{}_{}".format(self._caller,self._mapper)

        return combo

    def filter(self,threshold):

        filtered_vcf = None
        if self._caller != "nanoSV":
            try:
                for sp in ["RE", "SUPPORT", "DV"]:
                    if sp in list(self.vcf.columns):
                        filtered_vcf = self.vcf[self.vcf[sp] >= threshold]
            except:
                print("Can't find support count in this vcf: {}".format(self._caller))

        self.vcf = filtered_vcf

    def format_adjust(self, sv, principle="basic"):

        ref = ALTS[self._caller]
        sv_dict = []
        for g in ref.keys():
            if ref[g] == sv:
                sv_dict.append(g)

        if principle == "basic":
            subvcf = self.vcf[self.vcf["SVTYPE"].isin(sv_dict)]
            subvcf.loc[:,"SVTYPE"] = sv
            subvcf.loc[:,"END"] = subvcf["POS"]+abs(subvcf["SVLEN"])
            subvcf.to_csv("./{}_{}_{}.vcf".format(self.get_combo(),"fetal",sv), sep="\t")
            return Subvcf(self.get_caller(), self.get_mapper(),subvcf)
        else:
            subvcf = self.vcf[self.vcf["SVTYPE"].isin(ref.keys())]
            return Subvcf(self.get_caller(), self.get_mapper(),subvcf)

class Subvcf:

    def __init__(self, caller, mapper, vcf: VcfRecord):
        self._mapper = mapper
        self._caller = caller
        self._sample_name = "{}_{}".format(mapper, caller)
        self.vcf = vcf

    def get_combo(self):
        combo = "{}_{}".format(self._caller, self._mapper)

        return combo


def build_joint_alignment(vcfs: dict):

    database=pd.concat([x.vcf for x in vcfs.values()],axis=0)

    return database

class Caster:

    def __init__(self,combos,vcf,sv):
        self.combos = combos
        self.vcf = {}
        self.sv = sv
        self.genome_array = None
        self.subvcf = None
        self.pts = dict()
        for chr in CHROMOSOME:
            self.vcf[chr] = vcf[vcf["CHROM"]==chr]
        self.get_pts()
    def build_genome_array(self):

        for chr in CHROMOSOME:
            self.genome_array[chr] = np.repeat(0,(SCALE_CHR[chr]//5)+1)

    def get_pts(self):

        for chr in CHROMOSOME:
            self.pts[chr] = pd.DataFrame()
            self.pts[chr]['x']=self.vcf[chr]["center"].values
            self.pts[chr]["y"] = self.vcf[chr]["length"].values/2
            self.pts[chr]["combo"] = self.vcf[chr]["combo"].values

    def scan(self,locat,threshold=4):
        self.subvcf = self.vcf
        results = []
        if self.subvcf != None:
            for chr in CHROMOSOME:
                self.subvcf[chr].loc[:,"group"] = self.subvcf[chr].loc[:,"POS"] // 5
                X = self.pts[chr][["x", "y"]].values
                y_pred = DBSCAN(eps=100, min_samples=threshold).fit_predict(X)
                ydf = pd.DataFrame(columns=["group"], data=y_pred)
                self.vcf[chr].loc[:,"group"] = ydf.loc[:,"group"].values
                results.append(self.vcf[chr][self.vcf[chr]["group"] != -1])
            overlap = pd.concat(results, axis=0)
            overlap.drop(["center", "combo", "length"], axis=1, inplace=True)
            print(overlap["group"].value_counts())
            bed = overlap.loc[:,["CHROM","POS","END"]]
            bed.to_csv("./{}/overlap_{}.bed".format(locat,self.sv),index=False,encoding="utf-8",sep='\t',header=False)
            allel.write_vcf("./{}/overlap_{}.vcf".format(locat,self.sv), overlap)
        else:
            print("Here is a empty vcf: {}".format(self.sv))

#            print("DBSCAN view: {}".format(chr),self.pts)
            #plt.scatter(X[:,0],X[:,1],c=y_pred)
            #plt.title(self.sv)
            #fg1 = plt.gcf()
            #fg1.savefig("./{}/{}_{}.pdf".format(locat,chr,self.sv))

#            counts = self.subvcf["group"].value_counts()
#            for pos in range(len(self.genome_array[chr])):
#                self.genome_array[chr][pos] = counts[pos]

def vcfs_eater(file_names):
    record = {}
    vcf_combo = []
    for path, caller, mapper in file_names:
        print(path, caller, mapper)
        df = allel.vcf_to_dataframe(path, fields=LABELS[caller])
        df = df[df["CHROM"].isin(CHROMOSOME)]
        for sv in ALTS[caller].keys():
            df.replace(sv, ALTS[caller][sv], inplace=True)
        #        df["signal"] = ylabel.index(df["combo"].tolist())
        if caller != "nanoSV":
            df.rename(columns={LABELS[caller][6]: "RE"}, inplace=True)
        df.loc[:,"ID"] = "{}_{}".format(caller, mapper) + "_" + df["SVTYPE"] + "_" + df["ID"]
        df = vcf_adjuster(df)
        cleandata = VcfRecord(caller, mapper, df)
        cleandata.filter(1)
        record["{}_{}".format(caller, mapper)] = cleandata
        vcf_combo.append("{}_{}".format(caller, mapper))

    return record,vcf_combo


def vcf_adjuster(df):
    df["combo"] = "{}_{}".format(caller, mapper)
    df.loc[:,"END"] = df["POS"] + abs(df["SVLEN"])
    df["length"] = abs(df["SVLEN"])
    df["center"] = (df["POS"] + df["END"]) // 2
    return df

def main():
    parser = argparse.ArgumentParser(description="test")
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-i', '--input', type=str, help="Input VCF folder.", required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    req_grp.add_argument('-s', '--threshold', type=int, help="filter threshold",default=5)
    req_grp.add_argument('-e', '--eps', type=int, help="clustering eps", default=10)
#    parser.add_argument('-l', '--log', type=str, help="Log file name")
    file_names = []
    args = parser.parse_args()
    locat = args.output.strip('/') + '/'
    vcf_path = args.input.strip('/') + '/'
    vcfs = [x for x in os.listdir(vcf_path) if 'vcf' in x]
    for vcf_file in vcfs:
        print("Now read {}".format(vcf_file))
        full_path = vcf_path + vcf_file
        basename = vcf_file.split(".")[0]
        caller = basename.split("_")[0]
        mapper = basename.split("_")[-1]
        print("mapper: {}, caller: {}".format(mapper, caller))
        file_names.append((full_path, caller, mapper))

    if os.path.isdir(locat) is False:
        try:
            os.system("mkdir {}".format(locat))
        except:
            print("folder exist.")

    """if args.log:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    else:
        out_path = os.path.dirname(args.output)
        logging.basicConfig(filename=(out_path + '/log.log'), level=logging.WARNING)
"""
    logging.info(f"Start time: {datetime.now()}\n")
    logging.info(args)

    record,vcf_combo = vcfs_eater(file_names)

    vcf_database = build_joint_alignment(record)
    casters = {}

    for sv in DEFAULT:
        sv_database = vcf_database[vcf_database["SVTYPE"]==sv]
        casters[sv] = Caster(vcf_combo,sv_database,sv)
        casters[sv].scan(locat)

if __name__ == "__main__":
    main()