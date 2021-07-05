#!/usr/bin/env python3
import os
from executor import execute


def main():
    prefix = "middle_eastern.sgp"
    vep_script = "/home/snouto/mendel/bin/vep/vep"
    out_dir = "/home/snouto/mendel/out"
    sorted_dir = os.path.join(out_dir, "sorted")
    if not os.path.exists(sorted_dir):
        os.makedirs(sorted_dir)
    vcf_file = "/home/snouto/mendel/data/middle_eastern478_finalgatk4_newcovid19_recal.vcf.gz"
    with open("/home/snouto/mendel/data/chromosomes.txt") as cfile:
        chrs = cfile.readlines()
        chrs = [x.replace("\n", "") for x in chrs]

    for chr in chrs:
        out_file = os.path.join(out_dir, f"{prefix}.{chr}.vcf")
        print(f"Splitting VCF for Chromosome: {chr}")
        split_cmd = f"bcftools view --regions {chr} {vcf_file} > {out_file}"
        print(f"Running Command: {split_cmd}")
        execute(split_cmd)
        print(f"Annotating: {out_file}")
        vep_cmd = f"vep -i {out_file} \
            --offline \
            --fasta /root/.vep/homo_sapiens/102_GRCh38/ \
            --plugin dbNSFP,/home/snouto/mendel/data/nsfp/dbNSFP4.2a_variant.{chr}.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
            -everything \
            --buffer_size 100000 \
            --force_overwrite \
            --dir_cache /root/.vep/ \
            --dir_plugins /root/.vep/Plugins \
            --cache -o {out_file}.txt"
        execute(vep_cmd)
        vep_out = f"{out_file}.txt"
        sorted_file_name = f"{prefix}.SortByGene.{chr}.txt"
        execute(f"grep -v \# {vep_out} | sort -k 4 > {sorted_dir}/{sorted_file_name}")

    print("All Finished")


if __name__ == '__main__':
    main()
