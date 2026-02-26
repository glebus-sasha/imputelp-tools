#!/usr/bin/env python3
"""
vcf_to_refalt_table.py

Конвертирует merged VCF с числовыми GT в Beagle-подобную таблицу с REF|ALT.
"""

import argparse
from cyvcf2 import VCF
import pandas as pd

def vcf_to_refalt_table(input_vcf, output_file, info_field="AF"):
    """
    Преобразует VCF в TSV с GT в формате REF|ALT
    """
    vcf = VCF(input_vcf)
    samples = vcf.samples
    rows = []

    for var in vcf:
        chrom = var.CHROM
        pos = var.POS
        ref = var.REF
        alt = var.ALT[0]  # Берем первый альтернативный аллель
        info = var.INFO.get(info_field, ".")  # Используем заданное поле INFO
        
        genotypes = []
        for g in var.genotypes:  # [allele1, allele2, phased]
            g1 = ref if g[0] == 0 else alt
            g2 = ref if g[1] == 0 else alt
            genotypes.append(f"{g1}|{g2}")
        
        rows.append([chrom, pos, ref, alt, info] + genotypes)
    
    df = pd.DataFrame(rows, columns=["CHROM","POS","REF","ALT","INFO"] + samples)
    df.to_csv(output_file, sep="\t", index=False)
    print(f"[INFO] Успешно создан Beagle-подобный файл: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Convert merged VCF to REF|ALT Beagle-style table")
    parser.add_argument("--input", "-i", required=True, help="Входной VCF файл (можно gzipped)")
    parser.add_argument("--output", "-o", required=True, help="Выходной TSV файл")
    parser.add_argument("--info-field", "-f", default="AF", help="Поле INFO для использования в столбце INFO (по умолчанию AF)")

    args = parser.parse_args()
    vcf_to_refalt_table(args.input, args.output, args.info_field)

if __name__ == "__main__":
    main()
