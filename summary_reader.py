import re
import os
import argparse
import pandas as pd
import numpy as np
CALLER = ['cuteSV', 'sniffles', 'svim', 'nanoSV', 'overlap', 'noisy']
MAPPER = ["minimap2","ngmlr","DEL","INS","DUP","INV"]
HEARERS=["UPSTREAM","DOWNSTREAM","GENE","EXON","INTERGENIC","INTRON","TRANSCRIPT"]
COLS = ["sample","Dtype","caller","mapper"] + HEARERS
def get_information(filepath):
    f = open(filepath,"r",encoding="utf-8")
    html = f.read()
    count_data = []
    pct_data = []
    for head in HEARERS:
        content = re.findall(r'{}(.*?)</tr>'.format(head),html,re.S)
        print(filepath)
        try:
            info = content[0]
            percentage = info.split("%")[0].split(" ")[-1]
            count = info.split("#")[1].split(" ")[1]
            if "," in count:
                count1 = int(count.split(",")[0]) * 1000
                count2 = int(count.split(",")[1]) + count1
            else:
                count2 = int(count)
            count_data.append(count2)
            pct_data.append(percentage)
            print("IN {}, the {} have count {}, {}% of all region.".format(filepath, head, count2, percentage))

        except:
            print("This html do not have {}.".format(head))
            count2 = 0
            percentage = 0
            count_data.append(count2)
            pct_data.append(percentage)
    return count_data, pct_data


def main():
    parser = argparse.ArgumentParser(description="test")
    req_grp = parser.add_argument_group(title='Required')
    req_grp.add_argument('-i', '--input', type=str, help="Input html folder.", required=True)
    req_grp.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    req_grp.add_argument('-s','--sample',type=str,help="input sample name",required=True)
    file_names = []
    args = parser.parse_args()
    if args.output[0] == "/":
        locat = '/'+ args.output.strip('/')
    else:
        locat = args.output.strip('/')
    html_path = args.input.strip('/') + '/'
    os.system("mkdir -p {}".format(locat))
    sample = args.sample
    htmls = [x for x in os.listdir(html_path) if '.html' in x]
    sample_data_count = []
    sample_data_pct = []
    for html in htmls:
        print("Now read {}".format(html))
        full_path = html_path + html
        basename = html.split(".")[0]
        caller = ""
        mapper = ""
        for i in CALLER:
            if i in basename:
                caller = i
        for m in MAPPER:
            if m in basename:
                mapper = m

        print("mapper: {}, caller: {}".format(mapper, caller))
        cdata,pdata = get_information(full_path)
        count_line = [sample,"count",caller,mapper] + cdata
        pct_line = [sample,"pct",caller,mapper] + pdata
        sample_data_count.append(count_line)
        sample_data_pct.append(pct_line)
    out_count = pd.DataFrame(np.array(sample_data_count),columns=COLS)
    out_pct = pd.DataFrame(np.array(sample_data_pct),columns=COLS)
    out_count.to_csv("{}/{}_count_summary.csv".format(locat,sample),sep="\t",index=False)
    out_pct.to_csv("{}/{}_pct_summary.csv".format(locat,sample),sep="\t",index=False)


if __name__ == "__main__":
    main()