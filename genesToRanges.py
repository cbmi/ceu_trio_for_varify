#!/usr/bin/env python2.7
from cruzdb import Genome

#this takes a list of genes on chromsome20 and gets their transcript coords
g = Genome(db="hg19")

genes = ['AHCY','ARFGEF2','BMP2','DNAJC5','EDN3','GSS','GNAS1','JAG1','PANK2','PRNP','tTG','SALL4','VAPB']

for gene in genes:
    gene_obj = g.refGene.filter_by(name2=gene).first()
    if gene_obj:
        #one based intervals
        #http://gatkforums.broadinstitute.org/discussion/1204/what-input-files-does-the-gatk-accept-require
        print ("{0}:{1}-{2}".format(gene_obj.chrom.replace('chr',''),gene_obj.txStart,gene_obj.txEnd))
