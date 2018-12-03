#!/usr/bin/env python3


import os
from sys import argv
from collections import defaultdict

def main():
    infile = argv[1]
    prefix = argv[2]
    loci = set()
    with open(infile, 'r') as mlst:
        for line in mlst:
            loci.add(line.split('\t')[0].split('_')[0])
    num_loci = len(loci)
    list_loci = list(loci)
    list_loci.sort()
    d = defaultdict(dict)
    with open(infile, 'r') as mlst:
        for line in mlst:
            if  line.split('\t')[0].split('_')[0] in d[line.split('\t')[1]]:
                allele_list = d[line.split('\t')[1]][line.split('\t')[0].split('_')[0]]
                allele_list.append(line.strip().split('\t')[0].split('_')[1])
                d[line.split('\t')[1]][line.split('\t')[0].split('_')[0]] = allele_list
            else:    
                d[line.split('\t')[1]][line.split('\t')[0].split('_')[0]] = [line.strip().split('\t')[0].split('_')[1]]
    '''
    for key in d:
        print(key)
        for v in d[key]:
            if len(d[key][v]) > 1:
                print('\t'+ v + ' ' + str(d[key][v]))
    '''
    with open(prefix + '.report.out', 'w') as out:
        for key in d:
            alleles_called = 0
            alleles_multiple = 0
            for v in d[key]:
                alleles_called += 1
                if len(d[key][v]) > 1:
                    alleles_multiple += 1
            out.write(key + '; total: ' + str(alleles_called) + '/' + str(num_loci) + ', multiple: ' + str(alleles_multiple) + '\n')
    
    detailed = open(prefix +'.detailed.tsv', 'w')

    with open(prefix + '.raw.tsv', 'w') as out:
        out.write('\t' +'\t'.join(list_loci) + '\n')
        detailed.write('\t' +'\t'.join(list_loci) + '\n')
        for key in d:
            alleles = [key]
            alleles_detailed = [key]
            for l in list_loci:
                try:
                    a = d[key][l]
                    if len(a) == 1:
                        alleles.append(a[0])
                        alleles_detailed.append(a[0])
                    else:
                        alleles.append('NA')
                        alleles_detailed.append('MULTI')
                except KeyError:
                    alleles.append('NA')
                    alleles_detailed.append('NOT_CALLED')
            out.write('\t'.join(alleles) + '\n')
            detailed.write('\t'.join(alleles_detailed) + '\n')
    detailed.close()

    outfile = open(prefix +'.clean.tsv', 'w')
    dropped = open(prefix +'.dropped.txt', 'w')
    #cleaning up table according to Pasteur recommendations for LM cgMLST, not applicable to 7 gene or other MLST schemes
    with open(prefix + '.raw.tsv', 'r') as infile:
        for line in infile:
            if line.strip().split('\t').count('NA') > 87:
                dropped.write(line.strip().split('\t')[0] + '\n')
            else:
                outfile.write(line)
            
    dropped.close()
    outfile.close()

if __name__ == '__main__':
    main()
