#!/usr/bin/env python3
"""
Конвертирует merged VCF с числовыми GT в таблицу с REF|ALT
и пересчитывает AF по всей популяции.
"""

import argparse
from cyvcf2 import VCF
import pandas as pd

def vcf_to_refalt_af_table(input_vcf, output_file):
    """
    Преобразует VCF в TSV с GT в формате REF|ALT и пересчитанной AF.
    """
    vcf = VCF(input_vcf)
    samples = vcf.samples
    rows = []

    for var in vcf:
        chrom = var.CHROM
        pos = var.POS
        ref = var.REF
        alt = var.ALT[0]  # Берём первый альтернативный аллель

        # Пересчёт AF по всем образцам
        alt_count = 0
        total_alleles = 0
        genotypes = []
        for g in var.genotypes:  # g = [allele1, allele2, phased]
            a1, a2 = g[0], g[1]
            total_alleles += 2
            if a1 == 1:
                alt_count += 1
            if a2 == 1:
                alt_count += 1
            g1 = ref if a1 == 0 else alt
            g2 = ref if a2 == 0 else alt
            genotypes.append(f"{g1}|{g2}")

        af = alt_count / total_alleles if total_alleles > 0 else 0

        rows.append([chrom, pos, ref, alt, round(af,6)] + genotypes)

    # Создаём DataFrame
    df = pd.DataFrame(rows, columns=["CHROM","POS","REF","ALT","AF"] + samples)
    df.to_csv(output_file, sep="\t", index=False)
    print(f"[INFO] Успешно создан файл с пересчитанной AF: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert merged VCF to REF|ALT Beagle-style table with recalculated AF")
    parser.add_argument("--input", "-i", required=True, help="Входной VCF файл (merged, с числовыми GT)")
    parser.add_argument("--output", "-o", required=True, help="Выходной TSV файл с REF|ALT и пересчитанной AF")
    args = parser.parse_args()

    vcf_to_refalt_af_table(args.input, args.output)

if __name__ == "__main__":
    main()
