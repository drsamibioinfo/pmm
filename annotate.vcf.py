#!/usr/bin/env python
import os, argparse, sys


print_msg = lambda msg, error=False: print(f"[*] Message: {msg}") if not error else print(f"[*] Error: {msg}")

execute = lambda cmd: os.system(cmd)

def setup():
    p = argparse.ArgumentParser(usage="This script is used to annotate multiple vcf files in a folder through VEP")
    p.add_argument("--input", "-i", help="Input Directory containing VCF files to use", required=True)
    p.add_argument("--output", "-o", help="Output Directory to save annotated VCF Files", required=True)
    p.add_argument("--tabular", "-t", const=True, nargs="?", default=False,
                   help="The output Files should be in Tabular Format")
    p.add_argument("--vcf", "-c", const=True, nargs="?", default=False, help="The output files should be in VCF Format")
    p.add_argument("--fasta", "-f", required=True, help="Homo Sapiens (hg38) Assembly reference Fasta File")
    p.add_argument("--vep", "-p", required=True, help="VEP root directory")
    p.add_argument("--db", "-d", help="dbNSFP Database gz file which contains all chromosomes", required=True)
    return p


def annotate(args):
    input = args.input
    output = args.output
    tabular = args.tabular
    fasta = args.fasta
    vep_dir = args.vep
    dbnsfp_file = args.db
    print_msg("Reading input directory")
    if not os.path.isdir(input):
        print_msg("Input Parameter should be a directory.", error=True)
        raise Exception("Input Parameter should be a directory")
    files = os.listdir(input)
    for file in files:
        # ignore index file
        if ".tbi" in file:
            continue
        if not ".gz" in file:
            continue
        file_path = os.path.join(input, file)
        file_name, _ = os.path.splitext(file)
        output_file = os.path.join(output, f"{file_name}.norm.vcf.gz")
        normalize_cmd = f"bcftools norm -f {fasta} {file_path} | bgzip -c > {output_file}"
        print_msg(msg=f"Normalizing VCF File: {file_path}")
        execute(normalize_cmd)
        print_msg(msg=f"Creating an index for the output File: {output_file}")
        execute(f"tabix {output_file}")
        vep_fasta = os.path.join(vep_dir, "homo_sapiens/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz")
        vep_cmd = f"vep -i {output_file} \
        --offline \
        --fasta {vep_fasta} \
        --plugin dbNSFP,{dbnsfp_file},Ensembl_transcriptid,Uniprot_acc \
        --everything \
        --buffer_size 100000 \
        --force_overwrite \
        --dir_cache {vep_dir} \
        --dir_plugins {os.path.join(vep_dir, 'Plugins')} \
        --symbol \
        --pick \
        --hgvs \
        --terms SO \
        --cache -o {output_file}.txt"
        print_msg(msg=f"Annotating VCF File: {output_file}")
        execute(vep_cmd)
        break
    print(f"All Annotation Done, Files saved into {output}")


def main():
    parser = setup()
    if len(sys.argv) <= 1:
        parser.print_help()
        return
    args = parser.parse_args()
    annotate(args)



if __name__ == '__main__':
    main()
