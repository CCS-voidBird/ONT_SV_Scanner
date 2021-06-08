import allel
import pandas as pd
pd.options.mode.chained_assignment = None
import os
import argparse

INTRO="vcf_stat was designed for collecting numbers from vcf files, include mapper, caller, SV type, absolute length and add sample name to each record."

HEADER = ["SVTYPE","SVLEN""length","caller","mapper","RE","sample"]

CALLERS = ["cuteSV","nanoSV","svim","sniffles","overlap"]

MAPPERS = ["minimap2","ngmlr","overlap"]

def main():
    parser = argparse.ArgumentParser(description=INTRO)
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-i', '--input', type=str, help="Input VCF folder.", required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    req_grp.add_argument('-s', '--sample', type=str, help="sample name",required=True)
    req_grp.add_argument('-r',"--record",type=bool,help="record file (currently not working",default=False)
    file_names = []
    args = parser.parse_args()
    if args.output[0] == "/":
        locat = '/'+ args.output.strip('/') + '/'
    else:
        locat = args.output.strip('/') + '/'
    vcf_path = args.input.strip('/') + '/'
    sample = args.sample
    os.system("mkdir -p {}".format(locat))
    vcfs = [x for x in os.listdir(vcf_path) if 'vcf' in x]
    for vcf_file in vcfs:
        print("Now read {}".format(vcf_file))
        full_path = vcf_path + vcf_file
        basename = vcf_file.split(".")[0]

        file_names.append([full_path,basename])

    stats = [x for x in os.listdir(locat) if 'csv' in x]
    pre_reco = None
    records=[]
    if args.record is True:
        if len(stats) > 1:
            print("Multiple stat records in the output folder!")
            exit()
        if len(stats) ==1:
            reco_path = locat + stats[0]
            pre_reco = pd.read_csv(reco_path,sep="\t")
        elif len(stats) == 0:
            pre_reco = pd.DataFrame()
        records = [pre_reco]

    for vcf in file_names:
        path = vcf[0]
        name = vcf[1]
        record = allel.vcf_to_dataframe(path,fields="*")
        for caller in CALLERS:
            if caller in name:
                record["caller"] = caller
        for mapper in MAPPERS:
            if mapper in name:
                record["mapper"] = mapper
        record["length"] = abs(record["SVLEN"])
        record["sample"] = sample
        if "svim" in record["caller"].values:
            record.rename(columns={"SUPPORT": "RE"}, inplace=True)
        if "nanoSV" in record["caller"].values:
            record["RE"] = record["QUAL"]
        filtered = record[["SVLEN","RE","caller","mapper","SVTYPE","length","sample"]]
        records.append(filtered)
    all = pd.concat(records,axis=0)
    all.to_csv("{}/{}.csv".format(locat,sample),sep="\t",index=False)

if __name__ == "__main__":
    main()



