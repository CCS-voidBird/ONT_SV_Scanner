import allel
import numpy as np
import matplotlib.pyplot as plt
import BioVenn
import venn
import pandas as pd
pd.options.mode.chained_assignment = None
from sklearn.cluster import DBSCAN
import os
import argparse


INTRO = "SV Scanner -- vcf_cluster was designed for detecting common SVs from different approaches.\nCurrently support Aligners: Minimap2, Ngmlr; SV caller: cuteSV, SVIM, Sniffles, NanoSV."
NOTE={"ALT_1":"SV changes",\
      "SVTYPE":"Type of structural variant",\
      "END":"End position of the variant described in this record",\
      "RE":"Number of read support this record",\
      "SVLEN":"Difference in length between REF and ALT alleles",\
      "GT":"Genotype",\
      "group2": "Cluster name for each clusted SV"}
CALLERS=["cuteSV","nanoSV","sniffles","svim"]
MAPPERS=["minimap2","ngmlr"]
CHROMOSOME=[str(x) for x in range(1,30)]+["X","MT"]
LABELS={
    "cuteSV":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "sniffles":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","RE","samples","GT"],
    "svim":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","SUPPORT","samples","GT"],
    "nanoSV":["CHROM","POS","REF","variants/ALT","ID","END","SVTYPE","SVLEN","QUAL","samples","GT"]
}

ALTS={
    "cuteSV": {"INS":"INS","DEL":"DEL","DUP":"DUP","INV":"UNO","BND":"BND"},
    "sniffles": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","INVDUP":"DUP","TRA":"UNO"},
    "svim": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","DUP:TANDEM":"DUP","DUP:INT":"DUP","BND":"BND"},
    "nanoSV": {"DEL":"DEL","DUP":"DUP","BND":"BND","INS":"INS"}
}
DEFAULT=["DEL","INS","DUP"]


class VcfRecord:
    """
    Object: VCF dataset
    """

    def __init__(self, caller, mapper, vcf: pd.DataFrame):

        self._mapper = mapper
        self._caller = caller
        self.vcf = vcf
        self.subvcf = None

    def get_caller(self):

        return self._caller

    def get_mapper(self):

        return self._mapper

    def get_combo(self):

        combo = "{}_{}".format(self._caller,self._mapper)

        return combo

    """
    ## Unused filtered function
    
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
    """

    """
    ## Unused format change function
    
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
    """
class Subvcf:

    def __init__(self, caller, mapper, vcf: VcfRecord):
        self._mapper = mapper
        self._caller = caller
        self.vcf = vcf

    def get_combo(self):
        combo = "{}_{}".format(self._caller, self._mapper)

        return combo


def build_joint_alignment(vcfs: dict):
    """
    merge vcfs
    :return: a large dataframe
    """

    database=pd.concat([x.vcf for x in vcfs.values()],axis=0)
    return database

class Caster:

    """
    This class deal with clustering works
    """

    def __init__(self,combos,vcf,sv):

        self.combos = combos # combinations
        self.vcf = {} # dict of SVs
        self.sv = sv # SV type
        self.subvcf = None # place holder for subsets
        self.pts = dict()
        for chr in CHROMOSOME:
            self.vcf[chr] = vcf[vcf["CHROM"]==chr]
        self.get_pts()

    def get_pts(self):
        '''
        Transfer SVs to (x,y) pts; X: center postion; Y: half-length
        '''

        for chr in CHROMOSOME:
            self.pts[chr] = pd.DataFrame()
            self.pts[chr]['x']=self.vcf[chr]["center"].values
            self.pts[chr]["y"] = self.vcf[chr]["length"].values/2
            self.pts[chr]["combo"] = self.vcf[chr]["combo"].values

    def scan(self,locat,threshold=2):
        """
        The clustering function
        :param locat: output folder
        :param threshold: Minpts for DBSCAN
        :return: files: overlap vcf; nosie vcf
                 data: clustering records for all SV data
        """
        self.subvcf = self.vcf
        backups=[]
        results = []
        noisy = []
        mapper_stat = {}
        caller_stat = {}
        group_record = []

        sep = LABELS.keys()
        for key in sep:
            caller_stat[key] = []
        for key in ["minimap2","ngmlr"]:
            mapper_stat[key] = []

        if self.subvcf != None:
            for chr in CHROMOSOME:
                X = self.pts[chr][["x", "y"]].values
                if len(X) != 0:
                    o_pred = DBSCAN(eps=10, min_samples=threshold).fit_predict(X)
                    odf = pd.DataFrame(columns=["group2"], data=o_pred)  # group2: clusters

                    backup = [x for x in o_pred if x != -1]
                    backup = [chr+"_"+self.sv+"_"+str(x) for x in backup] # rename clusters across chr.

                    self.vcf[chr].loc[:, "group2"] = odf.loc[:, "group2"].values
                    noisy.append(self.vcf[chr][self.vcf[chr]["group2"] == -1]) # record noise SVs
                    remote = self.vcf[chr][self.vcf[chr]["group2"]==-1]
                    remote["group2"] = remote["combo"]
                    backup_record = self.vcf[chr][self.vcf[chr]["group2"]!=-1]
                    backup_record["group2"] = backup
                    self.vcf[chr] = pd.concat([self.vcf[chr][self.vcf[chr]["group2"]!=-1],remote],axis=0)
                    results.append(self.vcf[chr])
                    backups.append(pd.concat([backup_record,remote],axis=0))

            noisy_set = pd.concat(noisy,axis=0) # a df for noise SVs
            overlap = pd.concat(results, axis=0) # a df for overlapped SVs
            overlap_backup = pd.concat(backups,axis=0) # a backup SV df for restore clustering results

            for combo in self.combos:

                key_mapper = combo.split("_")[1]
                key_caller = combo.split("_")[0]
                cvcf = overlap[overlap["combo"] == combo]
                gs = overlap_backup[overlap_backup["combo"] == combo]["group2"].tolist()
                combo_groups = list(set(cvcf["group2"].tolist()))

                for x in gs:
                    group_record.append([key_caller,key_mapper,x])
                mapper_stat[key_mapper] += combo_groups
                print("{} added {} {} records".format(key_mapper,len(combo_groups),self.sv))
                caller_stat[key_caller] += combo_groups
                print("{} added {} records".format(key_caller,len(combo_groups),self.sv))

            plot_venn(mapper_stat,locat,self.sv)
            print("Making caller venn.")
            plot_venn(caller_stat,locat,self.sv)

            overlap = overlap[overlap["group2"] != -1]
            overlap.drop(["center", "combo", "length"], axis=1, inplace=True)

            gts = [[str(x[0]),str(x[1])] for x in overlap["GT"].tolist()]
            gts_overlap = ["/".join(x) for x in gts]
            overlap["GT"] = gts_overlap
            gts = [[str(x[0]), str(x[1])] for x in noisy_set["GT"].tolist()]
            gts_noisy = ["/".join(x) for x in gts]
            noisy_set["GT"] = gts_noisy
            """
            bed output codes.
            """
            #bed = overlap.loc[:,["CHROM","POS","END"]]
            #bed.to_csv("{}/overlap_{}.bed".format(locat,self.sv),index=False,encoding="utf-8",sep='\t',header=False)
            allel.write_vcf("{}/overlap_{}.vcf".format(locat,self.sv), overlap,description=NOTE)
            allel.write_vcf("{}/noisy_{}.vcf".format(locat,self.sv), noisy_set,description=NOTE)
            r_groups = pd.DataFrame(np.array(group_record),columns=["caller","mapper","group"])
            return r_groups

        else:
            print("Here is a empty vcf: {}".format(self.sv))


def plot_venn(data,locat,sv):

    """
    make venn plots for mapper and caller
    :param data: cluster results
    :param locat: output folder
    :param sv: SV type
    """

    for key in data:
        print(key,"has {} records.".format(len(data[key])))
    if len(data.keys()) == 2:
        print("Doing mapper venn")
        name = "mapper"
        BioVenn.draw_venn(tuple(data["minimap2"]),tuple(data["ngmlr"]),list_z=[],subtitle="",nrtype="pct",
                          xtitle="minimap2",ytitle="ngmlr",title="{} Overlap of Mappers".format(sv),
                          output="png",filename="{}/{}_{}_venn.png".format(locat,name,sv))
    elif len(data.keys()) == 4:
        new_data = {}
        for key in data.keys():
            new_data[key] = set(data[key])
        print("Doing caller venn")
        name = "caller"
        BioVenn.draw_venn(tuple(data["cuteSV"]), tuple(data["sniffles"]), tuple(data["svim"]), subtitle="", nrtype="pct",
                          xtitle="cuteSV", ytitle="sniffles", ztitle="svim", title="{} Overlap of Callers".format(sv),
                          output="png", filename="{}/{}_{}_venn.png".format(locat, name, sv))
        fig = venn.venn(new_data,fmt="{percentage:.2f}%")
        plt.title("Overlap within Callers")
        plt.style.use('seaborn-whitegrid')
        plt.savefig("{}/{}_{}_venn+nanoSV.png".format(locat,name,sv))

def vcfs_eater(file_names,threshold):
    """
    read vcf files by name and filter SVs by threshold
    :param file_names: vcf files
    :param threshold: support counts
    :return: a dict contains SVs from different files, a list contains loaded combinations.
    """

    record = {}
    vcf_combo = []
    for path, caller, mapper in file_names:
        df = allel.vcf_to_dataframe(path, fields=LABELS[caller])
        geno = [x[0] for x in allel.read_vcf(path, fields=["GT"])["calldata/GT"].tolist()]
        df["GT"] = geno
        df = vcf_adjuster(df,caller,mapper,threshold)
        cleandata = VcfRecord(caller, mapper, df)
        record["{}_{}".format(caller, mapper)] = cleandata
        vcf_combo.append("{}_{}".format(caller, mapper))

    return record,vcf_combo


def vcf_adjuster(df,caller,mapper,threshold=5):
    """
    reformat SVs
    :param df: a SV dataframe
    :param caller:  SV caller
    :param mapper:  aligners
    :param threshold: support count
    :return: filtered SV dataframe
    """

    df = df[df["CHROM"].isin(CHROMOSOME)]
    for sv in ALTS[caller].keys():
        df.replace(sv, ALTS[caller][sv], inplace=True)
    if caller != "nanoSV":
        df.rename(columns={LABELS[caller][8]: "RE"}, inplace=True)
        df = df[df["RE"] >= threshold]

    df["caller"] = caller
    df["mapper"] = mapper
    df.loc[:, "ID"] = "{}_{}".format(caller, mapper) + "_" + df["SVTYPE"] + "_" + df["ID"]
    df["combo"] = "{}_{}".format(caller, mapper)
    df.loc[:,"END"] = df["POS"] + abs(df["SVLEN"])
    df["length"] = abs(df["SVLEN"])
    df["center"] = (df["POS"] + df["END"]) // 2
    df["ALT"] = df["ALT_1"]
    return df

def main():
    """
    a parameter reading system
    """
    parser = argparse.ArgumentParser(description=INTRO)
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-i', '--input', type=str, help="Input VCF folder.", required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    req_grp.add_argument('-s', '--threshold', type=int, help="filter threshold",default=5)
    req_grp.add_argument('-e', '--eps', type=int, help="clustering eps", default=10)
    file_names = []
    args = parser.parse_args()
    threshold = int(args.threshold)
    if args.output[0] == "/":
        locat = '/'+ args.output.strip('/') + '/'
    else:
        locat = args.output.strip('/') + '/'
    vcf_path = args.input.strip('/') + '/'

    """
    starting read directory and mkdir for output folder
    """
    vcfs = [x for x in os.listdir(vcf_path) if 'vcf' in x]
    for vcf_file in vcfs:
        print("Now read {}".format(vcf_file))
        full_path = vcf_path + vcf_file
        basename = vcf_file.split(".")[0]
        caller = ""
        mapper = ""
        for i in LABELS.keys():
            if i in basename:
                caller = i
        for m in ["minimap2","ngmlr"]:
            if m in basename:
                mapper = m

        print("mapper: {}, caller: {}".format(mapper, caller))
        file_names.append((full_path, caller, mapper))

    if os.path.isdir(locat) is False:
        try:
            os.system("mkdir {}".format(locat))
        except:
            print("folder exist.")


    record,vcf_combo = vcfs_eater(file_names,threshold)

    vcf_database = build_joint_alignment(record)
    casters = {}
    castbook = []
    for sv in DEFAULT:
        sv_database = vcf_database[vcf_database["SVTYPE"]==sv]
        casters[sv] = Caster(vcf_combo,sv_database,sv)
        castbook.append(casters[sv].scan(locat))
    castbooks = pd.concat(castbook,axis=0)
    castbooks.to_csv("{}/overlap_RECORD.csv".format(locat), index=False, encoding="utf-8", sep='\t')
if __name__ == "__main__":
    main()