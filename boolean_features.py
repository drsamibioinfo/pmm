"""
This program requires bcftools (accessible from the PATH)
"""
import os
import glob
import pandas as pd
import numpy as np

home_vcf = 'vcfs' # FOLDER CONTAINING VCF FILES
vcf_prefix = 'middle_eastern.sgp' # PREFIX OF YOUR VCF FILE, in my case the files are named: ../../joint.3.chr*.vcf.gz
home_annot = 'sorted' # FOLDER CONTAINING ANNOTATION FILES SORTED BY GENES
annot_prefix = 'middle_eastern.sgp.SortByGene' # PREFIX OF VEP ANNOTATION FILE, in my case the files are named: ../finalAnnot.SortByGene.chr*.txt
gnomAD_ancestry='gnomAD_EAS_AF'  # SELECTED ANCESTRY IN GNOMAD, CHANGE THIS TO YOUR ANCESTRY
exac_ancestry = 'EAS_AF'# SELECTED ANCESTRY IN EXAC, CHANGE THIS TO YOUR ANCESTRY
phenotype = './pmm.pheno.final.csv' # SEE FORMAT IN THE REPOSITORY README FILE

def update_variants_per_gene(data, bool_single, bool_multi, samples, gene):
    gene_with = np.any(data, axis = 0).astype(int)
    bool_single[gene] = {sample:gt for sample, gt in zip(samples, gene_with)}
    gene_with_multi = (np.sum(data, axis = 0) > 1).astype(int)
    bool_multi[gene] = {sample:gt for sample, gt in zip(samples, gene_with_multi)}

def update_gc_per_gene(data, infos,  hetero_homo, bool_aplo, fout, samples, gene):
    n_rows_gene = 0
    n_samples_poli_rare = 0
    if hetero_homo == 'hetero':
        gts_bool = data.astype(bool)
    elif hetero_homo == 'homo':
        gts_bool = (data > 1).astype(bool)
    else:
        raise ValueError('ERROR: wrong value for hetero_homo parameter')
    is_mut = np.any(gts_bool, axis = 0).astype(int)
    poli, poli_inv, counts = np.unique(gts_bool, axis = 1, return_inverse = True, return_counts= True)
    for i_poli in range(poli.shape[1]):
        gc_name = '{}_'.format(gene)
        full_name = gc_name
        flag_variants = False
        for index, variant in infos[poli[:,i_poli]].iterrows():
            gc_name += '{}_{}_{}_'.format(variant['chrom'], variant['start'], variant['alt'])
            full_name += '{}:{}-{}:{}-{}/'.format(variant['chrom'], variant['start'], variant['end'], variant['ref'], variant['alt'])
            flag_variants = True
        if flag_variants:
            gc_name = gc_name[:-1] # to remove the last _
            full_name = full_name[:-1] # to remove the last /
        else:
            gc_name += 'WT' # this is the case where no variants where present
            full_name += 'WT' # this is the case where no variants where present
        row_name = gc_name + ',' + full_name
        if gc_name not in bool_aplo:
            bool_aplo[gc_name] = {sample:0 for sample in samples}
            n_samples_this_poli = 0
        else:
            n_samples_this_poli = sum(bool_aplo[gc_name].values())
        for i_sample, sample in enumerate(samples):
            if poli_inv[i_sample] == i_poli:
                bool_aplo[gc_name][sample] = 1
                n_samples_this_poli += 1
        fout.write('{},{}\n'.format(row_name, n_samples_this_poli))
        fout.flush()
        n_rows_gene += 1

#--- Initialization
bool_single_ultrarare = {}
bool_multi_ultrarare = {}
bool_single_rare = {}
bool_multi_rare = {}
bool_single_common = {}
bool_multi_common = {}
bool_gc_unique_homo = {}
bool_gc_unique_hetero = {}
fout_gc_unique_homo = open('data_gc_unique_homo.txt','wt')
fout_gc_unique_hetero = open('data_gc_unique_hetero.txt','wt')
contigs = ['chr{}'.format(ind) for ind in range(1,23)] + ['chrX', 'chrY']
frequency_ultrarare = 0.001
frequency_rare = 0.01
frequency_common = 0.05
frequency_pathogenic = 0.05
ens2name = pd.read_csv('ensemblid_names.csv')
pheno = pd.read_csv(phenotype)
n_genes = 0

for contig in contigs:
    vcf = '{}/{}.{}.vcf.gz'.format(home_vcf, vcf_prefix, contig)
    if not os.path.exists(vcf):
        raise ValueError('ERROR: missing file {}'.format(vcf))
    if not os.path.exists(vcf+'.tbi'):
        raise ValueError('ERROR: missing file {}'.format(vcf+'.tbi'))
    ann = '{}/{}.{}.txt'.format(home_annot, annot_prefix, contig)
    if not os.path.exists(ann):
        raise ValueError('ERROR: missing file {}'.format(ann))
    cmd = 'bcftools view -h ./{} > header.txt'.format(vcf)
    os.system(cmd)
    samples_contig = None
    with open('./header.txt','rt') as fin:
        for l in fin.readlines():
            l = l.strip().split()
            if l[0] == '#CHROM':
                samples_contig = l[9:]
    if samples_contig is None:
        raise ValueError('ERROR: wrong format in {}'.format(vcf))
    print('Number of samples for contig {}: {}'.format(contig, len(samples_contig)))
    ind_vcf_samples = []
    sex_samples = {}
    samples = []
    for ind_vcf, sample in enumerate(samples_contig):
        if sample in pheno['sample'].values:
            ind_pheno = pheno.index[pheno['sample'] == sample].tolist()
            if len(ind_pheno) != 1:
                raise ValueError('ERROR: {} repetitions for  {} in {}'.format(len(ind_pheno), sample, phenotype))
            ind_pheno = ind_pheno[0]
            ind_vcf_samples.append(ind_vcf)
            if pheno.loc[ind_pheno,'gender'] not in [0,1]:
                raise ValueError('ERROR: wrong sex {} for {}'.format(pheno.loc[ind_pheno,'gender'], sample))
            sex_samples[ind_vcf] = pheno.loc[ind_pheno,'gender']
            samples.append(sample)
    print('Number of samples in {}: {}'.format(phenotype, len(samples)))
    fin_ann = open(ann)
    l = fin_ann.readline()
    gene_old = None
    gene_data = []
    while l:
        l = l.strip()
        if l:
            if l[0] != '#':
                l = l.split()
                idd = l[0]

                #chrom = l[0].split(':')[0]
                #ref = l[0].split(':')[2]
                #alt = l[0].split(':')[3]
                #pos = l[1].split(':')[1]
                #start = int(pos.split('-')[0])
                #if '-' in pos:
                #    end = int(pos.split('-')[1])
                #else:
                #    end = start
                #gene = l[3]
                #feat_type = l[5]
                #cons = l[6]
                #extra = l[13]
                
                # CHANGE: begin
                chrom = l[1].split(':')[0]
                ref = 'R' # I don't see the reference in any position in your file, so I'm using a dummy value here just to avoid changing the rest of the code
                alt = l[2]
                pos = l[1].split(':')[1]
                start = int(pos.split('-')[0])
                if '-' in pos:
                    end = int(pos.split('-')[1])
                else:
                    end = start
                gene = l[3]
                feat_type = l[5]
                cons = l[6]
                extra = l[13] # Could you check is this is correct ? It should be the part starting including the gnomAD annotations
                # CHANGE: end

                if gnomAD_ancestry in extra:
                    ind = extra.index('{}='.format(gnomAD_ancestry))
                    af = float(extra[ind:].split(';')[0].split('=')[1])
                elif exac_ancestry in extra:
                    ind = extra.index('{}='.format(exac_ancestry))
                    af = float(extra[ind:].split(';')[0].split('=')[1])
                else:
                    af = -1.0 # used as a place-hold for missing frequency
                switch = False
                if af >= 0.5: # here switch ref./alt.
                    af = 1 - af
                    dummy = alt
                    alt = ref
                    ref = dummy
                    switch = True
                if 'CLIN_SIG' in extra:
                    ind = extra.index('CLIN_SIG=')
                    clin_sig = extra[ind:].split(';')[0].split('=')[1]
                else:
                    clin_sig = ''
                if (gene != gene_old):
                    if (gene_old is not None) and ('ENSG' in gene_old):
                        gene_name = np.unique(ens2name['Gene name'].loc[(ens2name['Gene stable ID'] == gene_old)].values)[0]
                        info = pd.DataFrame(gene_data, columns = ['id', 'chrom', 'start', 'end', 'ref', 'alt', 'gene', 'feat_type', 'cons', 'af', 'switch', 'clin_sig'])
                        info.drop_duplicates(subset = 'id', inplace = True)
                        info.sort_values(by = 'id', axis = 0, inplace = True)
                        with open('ids.txt','wt') as fout:
                            for i in info['id']:
                                fout.write(i+'\n')
                        #cmd = 'bcftools filter -r {}:{}-{} -Ou ./{} | bcftools view -i "%ID=@./ids.txt" -Ou | bcftools query -f "%ID [ %GT] \n" > ./variants_{}.txt'.format(chrom, min(info['start']), max(info['start']), vcf, gene_name)
                        cmd = 'bcftools filter -r {}:{}-{} -Ou ./{} | bcftools view -i "%ID=@./ids.txt" -Ou | bcftools query -f "%ID [ %GT] \n" > ./variants.txt'.format(chrom, min(info['start']), max(info['start']), vcf)
                        os.system(cmd)
                        gts = []
                        #with open('variants_{}.txt'.format(gene_name),'rt') as fin:
                        with open('variants.txt'.format(gene_name),'rt') as fin:
                            for l in fin.readlines():
                                l = l.strip()
                                if l:
                                    l = l.split()
                                    id_variant = l[0] # id of the variant
                                    row = [id_variant,]
                                    switch = info.loc[info['id'] == id_variant, 'switch'].values
                                    if len(switch) != 1:
                                        raise ValueError('ERROR: not-unique variant id')
                                    switch = switch[0]
                                    for ind, gt in enumerate(l[1:]):
                                        if ind in ind_vcf_samples:
                                            if not switch:
                                                if (gt == '1/1') or (gt == '1|1'):
                                                    row.append(2)
                                                elif (contig == 'chrX') and (sex_samples[ind] == 0) and ((gt == '0/1') or (gt == '0|1') or (gt == '1/0') or (gt == '1|0')):
                                                    row.append(2)
                                                elif (gt == '0/1') or (gt == '0|1') or (gt == '1/0') or (gt == '1|0'):
                                                    row.append(1)
                                                elif (gt == '0/0') or (gt == '0|0') or (gt == './.') or (gt == '.|.') or (gt == '.'):
                                                    row.append(0)
                                                else:
                                                    raise ValueError('ERROR: unexpected genotype {}'.format(gt))
                                            else: # in this case ref. and alt. are switched
                                                if (gt == '0/0') or (gt == '0|0'):
                                                    row.append(2)
                                                elif (contig == 'chrX') and (sex_samples[ind] == 0) and ((gt == '0/1') or (gt == '0|1') or (gt == '1/0') or (gt == '1|0')):
                                                    row.append(2)
                                                elif (gt == '0/1') or (gt == '0|1') or (gt == '1/0') or (gt == '1|0'):
                                                    row.append(1)
                                                elif (gt == '1/1') or (gt == '1|1') or (gt == './.') or (gt == '.|.') or (gt == '.'):
                                                    row.append(0)
                                                else:
                                                    raise ValueError('ERROR: unexpected genotype {}'.format(gt))
                                    gts.append(row)
                        gts = pd.DataFrame(gts, columns = ['id',]+samples)
                        gts.sort_values(by = 'id', axis = 0, inplace = True)
                        #gts.to_csv('variants_genotyped_{}.csv'.format(gene_name))
                        if gts.shape[0] != info.shape[0]:
                            #raise ValueError('ERROR: wrong dimensions in info {} or gts {}'.format(info.shape, gts.shape))
                            print('WARNING: wrong dimensions in info {} or gts {} for gene {}'.format(info.shape, gts.shape, gene_old))
                        else:
                            data = info.merge(gts, left_on = 'id', right_on = 'id')
                            #data.to_csv('variants_annotated_{}.csv'.format(gene_name))
                            keep = np.logical_or.reduce( [data['cons'].str.contains(variant_type) for variant_type in ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification', 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant']] )
                            #same_len = data['ref'].str.len() == data['alt'].str.len()
                            same_len = data['start'] == data['end'] # CHANGE: as you don't have the ref. I'm using the length to separate SNPs from INDELs
                            unique = np.sum(data[samples].astype(int), axis = 1).astype(int) == 1
                            pat = data['clin_sig'].str.contains('athogenic')
                            var1 = np.logical_and.reduce( ( keep, same_len, np.logical_or(data['af'] < frequency_ultrarare, data['af'] == -1.0) ) ) # SNP to keep with frequency < frequency_ultrarare
                            var2 = np.logical_and.reduce( ( keep, np.logical_not(same_len), np.logical_or(data['af'] < frequency_ultrarare, np.logical_and(data['af'] == -1.0, unique)) ) ) # indels to keep with frequency < frequency_ultrarare (indels that are not-unique with missing frequency are discarded)
                            var3 = np.logical_and( pat, data['af'] < frequency_pathogenic ) # pathogenic with frequency < frequency_pathogenic
                            var4 = np.logical_and.reduce( ( keep, same_len, np.logical_or(data['af'] < frequency_rare, data['af'] == -1.0) ) )
                            var5 = np.logical_and.reduce( ( keep, np.logical_not(same_len), np.logical_or(data['af'] < frequency_rare, np.logical_and(data['af'] == -1.0, unique)) ) )
                            var6 = np.logical_and.reduce( ( keep, same_len, np.logical_or(data['af'] < frequency_common, data['af'] == -1.0) ) )
                            var7 = np.logical_and.reduce( ( keep, np.logical_not(same_len), np.logical_or(data['af'] < frequency_common, np.logical_and(data['af'] == -1.0, unique)) ) )
                            var8 = np.logical_and.reduce( ( keep, same_len, np.logical_or(data['af'] >= frequency_common, data['af'] == -1.0) ) )
                            var9 = np.logical_and.reduce( ( keep, np.logical_not(same_len), np.logical_or(data['af'] >= frequency_common, np.logical_and(data['af'] == -1.0, unique)) ) )
                            ultrarare = np.logical_or.reduce( (var1, var2, var3) )
                            rare = np.logical_and( np.logical_or.reduce( (var4, var5, var3) ), np.logical_not(ultrarare) )
                            common = np.logical_and( np.logical_or.reduce( (var6, var7, var3) ), np.logical_not(rare) )
                            combination = np.logical_or.reduce( ( var8, var9 ) )
                            if np.any(ultrarare):
                                update_variants_per_gene(data[samples].loc[ultrarare], bool_single_ultrarare, bool_multi_ultrarare, samples, gene_name)
                            if np.any(rare):
                                update_variants_per_gene(data[samples].loc[rare], bool_single_rare, bool_multi_rare, samples, gene_name)
                            if np.any(common):
                                update_variants_per_gene(data[samples].loc[common], bool_single_common, bool_multi_common, samples, gene_name)
                            if np.any(combination):
                                update_gc_per_gene(data[samples].loc[combination], data[['id','chrom','ref','alt','start','end']].loc[combination], 'homo', bool_gc_unique_homo, fout_gc_unique_homo, samples, gene_name)
                                update_gc_per_gene(data[samples].loc[combination], data[['id','chrom','ref','alt','start','end']].loc[combination], 'hetero', bool_gc_unique_hetero, fout_gc_unique_hetero, samples, gene_name)
                            n_genes += 1
                            print('Gene {}; number of genes: {}'.format(gene_old, n_genes), flush = True)
                    gene_data = []
                    gene_old = gene
                gene_data.append([idd, chrom, start, end, ref, alt, gene, feat_type, cons, af, switch, clin_sig])
        l = fin_ann.readline()

data = pd.DataFrame(bool_single_ultrarare)
data.to_csv('data_al1_ultrarare.csv')
data = pd.DataFrame(bool_multi_ultrarare)
data.to_csv('data_al2_ultrarare.csv')
data = pd.DataFrame(bool_single_rare)
data.to_csv('data_al1_rare.csv')
data = pd.DataFrame(bool_multi_rare)
data.to_csv('data_al2_rare.csv')
data = pd.DataFrame(bool_single_common)
data.to_csv('data_al1_common.csv')
data = pd.DataFrame(bool_multi_common)
data.to_csv('data_al2_common.csv')
data = pd.DataFrame(bool_gc_unique_homo)
data.to_csv('data_gc_unique_homo.csv')
data = pd.DataFrame(bool_gc_unique_hetero)
data.to_csv('data_gc_unique_hetero.csv')

fout_gc_unique_homo.close()
fout_gc_unique_hetero.close()
