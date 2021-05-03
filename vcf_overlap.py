import allel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches
import os
import argparse
import logging
from datetime import datetime

LINESTYLE = {"DUP":"solid","DEL":"dashed","INS":"dashdot"}
COLORS={}
CALLERS=["cuteSV","nanoSV","sniffles","svim"]
MAPPERS=["minimap2","ngmlr"]
cnames = ["darkblue","darkgray","darkred","darksalmon","lightblue","lightgray","red","lightsalmon"]
for mapper in MAPPERS:
    for caller in CALLERS:
        COLORS["{}_{}".format(caller,mapper)] = cnames.pop()
scale=[x*10000000 for x in [16,14,12,12,12,12,11,11,10.5,10.5,11,9,8.5,8.5,8.5,8,7,6.5,6.5,7,7,6,5,6.2,4.7,5,4.5,4.5,5,14]]

CHROMOSOME=[str(x) for x in range(1,30)]+["MT"]
SCALE_CHR={}
for i in range(len(CHROMOSOME)):
    SCALE_CHR[CHROMOSOME[i]]=scale[i]

LABELS={
    "cuteSV":["CHROM","POS","ID","END","SVTYPE","SVLEN","RE","calldata/GT"],
    "sniffles":["CHROM","POS","ID","END","SVTYPE","SVLEN","RE","calldata/GT"],
    "svim":["CHROM","POS","ID","END","SVTYPE","SVLEN","SUPPORT","calldata/GT"],
    "nanoSV":["CHROM","POS","ID","END","SVTYPE","SVLEN","DV","calldata/GT"]
}

ALTS={
    "cuteSV": {"INS":"INS","DEL":"DEL","DUP":"DUP","INV":"UNO","BND":"BND"},
    "sniffles": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","INVDUP":"DUP","TRA":"UNO"},
    "svim": {"DEL":"DEL","INS":"INS","INV":"UNO","DUP":"DUP","DUP:TANDEM":"DUP","DUP:INT":"DUP","BND":"BND"},
    "nanoSV": {"DEL":"DEL","DUP":"DUP","BND":"BND","INS":"INS"}
}
CHRS=[""]+["chromosome_{}".format(x) for x in range(1,29)]
DEFAULT=["DEL","INS","DUP"]
ylabel=[]
for caller in CALLERS:
    for mapper in MAPPERS:
        label = "{}_{}".format(caller,mapper)
        ylabel.append(label)
class VcfRecord:

    def __init__(self, caller, mapper, vcf: pd.DataFrame):

        self._mapper = mapper
        self._caller = caller
        self._sample_name = "{}_{}".format(mapper,caller)
        self.vcf = vcf
        self.subvcf = None
#        self.build_chr_df()

    def get_caller(self):

        return self._caller

    def get_mapper(self):

        return self._mapper

    def get_combo(self):

        combo = "{}_{}".format(self._caller,self._mapper)

        return combo

    def filter(self,threshold):

        filtered_vcf = None
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
            subvcf["SVTYPE"] = sv
            subvcf["END"] = subvcf["POS"]+abs(subvcf["SVLEN"])
            subvcf.to_csv("./{}_{}_{}.vcf".format(self.get_combo(),"fetal",sv), sep="\t")
            return Subvcf(self.get_caller(), self.get_mapper(),subvcf)
        else:
            subvcf = self.vcf[self.vcf["SVTYPE"].isin(ref.keys())]
            return Subvcf(self.get_caller(), self.get_mapper(),subvcf)

"""    def build_chr_df(self):

        self.vcf["scale_pos"] = self.vcf["POS"]/10000000/float(SCALE_CHR[self.vcf["CHROM"]])*100
        self.vcf["scale_svlen"] = self.vcf["SVLEN"]/10000000/float(SCALE_CHR[self.vcf.loc["CHROM"]])*100
        self.vcf["scale_end"] = self.vcf["END"]/10000000/float(SCALE_CHR[self.vcf.loc["CHROM"]])*100"""



class Subvcf:

    def __init__(self, caller,mapper,vcf:VcfRecord):

        self._mapper = mapper
        self._caller = caller
        self._sample_name = "{}_{}".format(mapper, caller)
        self.vcf = vcf

    def get_combo(self):

        combo = "{}_{}".format(self._caller,self._mapper)

        return combo


def build_joint_alignment(vcfs: dict):

    database=pd.concat([x.vcf for x in vcfs.values()],axis=0)

    return database

def holder(x,region):
    return np.where(x[0] in region,x[1],0)

def binary_tree(starts,ends):
    if starts != []:
#        range1 = starts[0]
#        range2 = ends[0]
        ranges = []
#        for i in range(len(starts) - 1):
        for i in range(len(starts)):
            ranges.append([starts[i],ends[i]])
#            if starts[i] > starts[i + 1]:
#                print("This pts list needs to be sorted")
#                break
#            else:
#                if starts[i + 1] <= ends[i]:
#                    range2 = ends[i + 1]
#                else:
#                    ranges.append([range1, range2])
#                    range1 = starts[i + 1]
#                    range2 = ends[i + 1]
#        if range1 != starts[0] and range2 != ends[0]:
#            ranges.append([range1, range2])
        return ranges


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
            df.rename(columns={LABELS[caller][5]: "RE"}, inplace=True)
        df["combo"] = "{}_{}".format(caller, mapper)
        df["ID"] = df["combo"] + "_" + df["SVTYPE"] + "_" + df["ID"]
        df["END"] = df["POS"] + abs(df["SVLEN"])
        print(df.head(2))
        record["{}_{}".format(caller, mapper)] = VcfRecord(caller, mapper, df)
        vcf_combo.append("{}_{}".format(caller, mapper))

    return record,vcf_combo

def build_print_data(df,sv,vcf_combo):
    print("===========")
    data = df[df["SVTYPE"]==sv]
    print_list={}
    for chr in CHROMOSOME:
        chrvcf = data[data["CHROM"] == chr]
        if len(chrvcf) == 0:
            print("Chr {} have no data".format(chr))
            continue
        for combo in vcf_combo:
            print(combo,chr)
            subvcf = chrvcf[chrvcf["combo"]==combo]
            print("scale:",SCALE_CHR[chr])
            #x/SCALE_CHR[chr]*100
            pos_list = [x/SCALE_CHR[chr]*100 for x in subvcf["POS"].tolist()]
            end_list = [x/SCALE_CHR[chr]*100 for x in subvcf["END"].tolist()]
            print(len(pos_list),len(end_list))
            print_list["{}_{}".format(chr,combo)] = binary_tree(pos_list,end_list)
    return print_list


def main():
    parser = argparse.ArgumentParser(description="test")
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-i', '--input', type=str, help="Input VCF folder.", required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    parser.add_argument('-l', '--log', type=str, help="Log file name")
    file_names = []
    args = parser.parse_args()
    locat = args.output.strip('/') + '/'
    vcf_path = args.input.strip('/') + '/'
    # vcf_path = "./vcf_shortcut/"
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

    if args.log:
        logging.basicConfig(filename=args.log, level=logging.INFO)
    else:
        out_path = os.path.dirname(args.output)
        logging.basicConfig(filename=(out_path + '/log.log'), level=logging.WARNING)

    # log = logging.getLogger(__name__)

    logging.info(f"Start time: {datetime.now()}\n")
    logging.info(args)

    record,vcf_combo = vcfs_eater(file_names)
    """
    for path,caller,mapper in file_names:
        print(path,caller,mapper)
        df = allel.vcf_to_dataframe(path, fields=LABELS[caller])
        df = df[df["CHROM"].isin(CHROMOSOME)]
        for sv in ALTS[caller].keys():
            df.replace(sv,ALTS[caller][sv],inplace=True)
#        df["signal"] = ylabel.index(df["combo"].tolist())
        if caller != "nanoSV":
            df.rename(columns={LABELS[caller][5]:"RE"},inplace=True)
        df["combo"] = "{}_{}".format(caller,mapper)
        df["ID"] = df["combo"]+"_"+df["SVTYPE"]+"_"+df["ID"]
        df["END"] = df["POS"] + abs(df["SVLEN"])
        print(df.head(2))
        record["{}_{}".format(caller,mapper)]=VcfRecord(caller,mapper,df)
        vcf_combo.append("{}_{}".format(caller,mapper))
    """

    aln_df = build_joint_alignment(record)
    cnames = [COLORS[x] for x in vcf_combo]
    petchs = [mpatches.Patch(color=color) for color in cnames]
    plt.figure(figsize=(21, 7))
    plt.style.use('dark_background')
    for sv in DEFAULT:
        int_ptx = []
        int_pty = []

#        plt.subplots(figsize=(21,7))
        #fig,ax=plt.subplots()

        print("start aligning: ",sv)
        print_dict = build_print_data(aln_df,sv,vcf_combo)
        #print(print_dict)
        for chr in CHROMOSOME:
            if print_dict is None:
                break
            for combo in vcf_combo:
                key = "{}_{}".format(chr,combo)
                color = COLORS[combo]
                h = CHROMOSOME.index(chr)+(vcf_combo.index(combo)+1)*(1/(len(vcf_combo)+2))
                #print(h)
                if key not in print_dict.keys() or print_dict[key] is None:
                    print("Combo {} in chr {} have no {} data".format(combo,chr,sv))
                    continue
                else:
                    for i in range(len(print_dict[key])):
                        xmin,xmax = print_dict[key][i]
                        int_ptx.append(xmin)
                        int_pty.append(h)
                        #print(xmin,xmax)
                        if sv is not "INS":
                            plt.hlines(h,xmin,xmax,colors=color,linestyles=LINESTYLE[sv],lw=2,zorder=10)
                    if sv is "INS":
                        print("@#$%^&*()#$%^&*()INS34567890-")
                        print(sv,chr,int_ptx,int_pty)
                        plt.scatter(int_ptx,int_pty,marker=2,color="black",s=0.01,zorder=20)
                    else:
                        plt.scatter(int_ptx,int_pty,marker="|",color="black",s=0.01)
        plt.legend(handles=petchs,labels=[x for x in vcf_combo],loc="best")
        plt.title("TEST")
        plt.ylim(0, 32)
        plt.xlim(0,100)
        yl = CHROMOSOME
        plt.yticks(np.arange(30), yl)
        plt.xlabel("Position")
        plt.ylabel("Chromosome")
        fg1=plt.gcf()
        #plt.show()
        fg1.savefig("./{}/{}.pdf".format(locat,sv))
    plt.close()


    print("done",DEFAULT)
if __name__ == "__main__":
    main()
