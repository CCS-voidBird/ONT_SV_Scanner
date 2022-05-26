# ONT_SV_RANGER

Version 0.5

This Toolkit needs to be maintained.

```
usage: vcf_cluster.py [-h] -i INPUT -o OUTPUT [-s THRESHOLD] [-e EPS]


SV Scanner -- vcf_cluster was designed for detecting common SVs from different
approaches. Currently support Aligners: Minimap2, Ngmlr; SV caller: cuteSV,
SVIM, Sniffles, NanoSV.


optional arguments:
  -h, --help            show this help message and exit


Required:
  -i INPUT, --input INPUT
                        Input VCF folder.
  -o OUTPUT, --output OUTPUT
                        Output directory
  -s THRESHOLD, --threshold THRESHOLD
                        filter threshold
  -e EPS, --eps EPS     clustering eps
```


```
usage: vcf_inheritance_cluster.py [-h] -1 MOTHER -2 OFFSPRING -o OUTPUT
                                  [-s THRESHOLD] [-n SAMPLE] [-e EPS]
                                  [-m MIN_LENGTH]


SV Scanner -- vcf_inheritance_cluster was designed for detecting inherited SVs
from given family pairs. Currently support vcf_cluster outputs and these
Aligners: Minimap2, Ngmlr; SV caller: cuteSV, SVIM, Sniffles, NanoSV.


optional arguments:
  -h, --help            show this help message and exit


Required:
  -1 MOTHER, --mother MOTHER
                        Input the path of mother.
  -2 OFFSPRING, --offspring OFFSPRING
                        Input the path of offspring.
  -o OUTPUT, --output OUTPUT
                        Output directory
  -s THRESHOLD, --threshold THRESHOLD
                        filter threshold
  -n SAMPLE, --sample SAMPLE
                        input the name of group
  -e EPS, --eps EPS     clustering eps
  -m MIN_LENGTH, --min_length MIN_LENGTH
                        min sv length}
```

