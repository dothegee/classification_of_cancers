# convert_external_to_csv
import os
import gzip
import shutil
from tqdm import tqdm
import glob

import pandas as pd


workspace_path = 'data/'

# 폴더명만 리스트로 가져오기
folder_list = [folder for folder in os.listdir(workspace_path) if os.path.isdir(os.path.join(workspace_path, folder))]

train = pd.read_csv("train_data/train.csv")

columns = [
        "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome",
        "Start_Position", "End_Position", "Strand", "Variant_Classification",
        "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2", "dbSNP_RS", "dbSNP_Val_Status",
        "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
        "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2",
        "Tumor_Validation_Allele1", "Tumor_Validation_Allele2",
        "Match_Norm_Validation_Allele1", "Match_Norm_Validation_Allele2",
        "Verification_Status", "Validation_Status", "Mutation_Status",
        "Sequencing_Phase", "Sequence_Source", "Validation_Method",
        "Score", "BAM_File", "Sequencer", "Tumor_Sample_UUID",
        "Matched_Norm_Sample_UUID", "HGVSc", "HGVSp", "HGVSp_Short",
        "Transcript_ID", "Exon_Number", "t_depth", "t_ref_count",
        "t_alt_count", "n_depth", "n_ref_count", "n_alt_count",
        "all_effects", "Allele", "Gene", "Feature", "Feature_type",
        "One_Consequence", "Consequence", "cDNA_position", "CDS_position",
        "Protein_position", "Amino_acids", "Codons", "Existing_variation",
        "DISTANCE", "TRANSCRIPT_STRAND", "SYMBOL", "SYMBOL_SOURCE",
        "HGNC_ID", "BIOTYPE", "CANONICAL", "CCDS", "ENSP",
        "SWISSPROT", "TREMBL", "UNIPARC", "UNIPROT_ISOFORM",
        "RefSeq", "MANE", "APPRIS", "FLAGS", "SIFT", "PolyPhen",
        "EXON", "INTRON", "DOMAINS", "1000G_AF", "1000G_AFR_AF",
        "1000G_AMR_AF", "1000G_EAS_AF", "1000G_EUR_AF", "1000G_SAS_AF",
        "ESP_AA_AF", "ESP_EA_AF", "gnomAD_AF", "gnomAD_AFR_AF",
        "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF", "gnomAD_FIN_AF",
        "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF", "MAX_AF",
        "MAX_AF_POPS", "gnomAD_non_cancer_AF", "gnomAD_non_cancer_AFR_AF",
        "gnomAD_non_cancer_AMI_AF", "gnomAD_non_cancer_AMR_AF",
        "gnomAD_non_cancer_ASJ_AF", "gnomAD_non_cancer_EAS_AF",
        "gnomAD_non_cancer_FIN_AF", "gnomAD_non_cancer_MID_AF",
        "gnomAD_non_cancer_NFE_AF", "gnomAD_non_cancer_OTH_AF",
        "gnomAD_non_cancer_SAS_AF", "gnomAD_non_cancer_MAX_AF_adj",
        "gnomAD_non_cancer_MAX_AF_POPS_adj", "CLIN_SIG", "SOMATIC",
        "PUBMED", "TRANSCRIPTION_FACTORS", "MOTIF_NAME", "MOTIF_POS",
        "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "miRNA", "IMPACT",
        "PICK", "VARIANT_CLASS", "TSL", "HGVS_OFFSET", "PHENO",
        "GENE_PHENO", "CONTEXT", "tumor_bam_uuid", "normal_bam_uuid",
        "case_id", "GDC_FILTER", "COSMIC", "hotspot", "RNA_Support",
        "RNA_depth", "RNA_ref_count", "RNA_alt_count", "callers"
    ]



# maf 파일 읽고 DataFrame으로 받기
def read_maf_to_dataframe(file_path, columns):
    # 파일에서 주석 줄을 무시하고 데이터 읽기
    #ann = []
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # 주석 줄 무시
                #ann.append(line.strip().split('\t'))
            data.append(line.strip().split('\t'))  # 탭으로 분리하여 리스트에 추가
    # DataFrame 생성
    df = pd.DataFrame(data, columns=columns)
    return df


# 모든 데이터 합치기
combined_df = pd.DataFrame(columns=train.columns)

for x in tqdm(folder_list):
# 폴더 경로 설정
    folder_path = f'C:\\workspace\\cancer_project\\data\\{x}'

    # .maf 확장자를 가진 파일들만 리스트로 만들기
    maf_files = glob.glob(f"{folder_path}/*.maf")


    #print(maf_files)
    df = pd.DataFrame(columns=columns)

    for maf in tqdm(maf_files):
            maf_df = read_maf_to_dataframe(maf, columns)
            # 칼럼명이 데이터 프레임 0번째에 들어감
            maf_df = maf_df.iloc[1:]
            df = pd.concat([df, maf_df], axis=0)

    # 데이터 전처리

    # index 초기화
    df.reset_index(inplace=True)
    # 단백질을 나타내는 p. 지우기
    df['HGVSp_Short'] = df['HGVSp_Short'].str.replace('p.', '', regex=False)
    # silent 돌연변이 기호 바꾸기 (train data와 같이)
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: x.replace('=', x[0]) if '=' in x else x)
    # fs 뒤에 지우기
    df['HGVSp_Short'] = df['HGVSp_Short'].apply(lambda x: x.split('fs')[0] + 'fs' if 'fs' in x else x)
    select_col = ['case_id','Tumor_Sample_Barcode','Hugo_Symbol','HGVSp_Short','Variant_Classification','Variant_Type','SIFT', 'PolyPhen','CLIN_SIG','COSMIC']
    # 필요한 col로 filter
    filtered_df = df[select_col]
    del df
    # pivot df 만들기 (index='Tumor_Sample_Barcode', columns='Hugo_Symbol', values='HGVSp_Short')
    pivot_df = filtered_df.pivot_table(index='Tumor_Sample_Barcode', columns='Hugo_Symbol', values='HGVSp_Short', aggfunc='first')
    pivot_df.columns.name = None  # 'Hugo_Symbol' 이름을 제거
    # 인덱스 제거
    pivot_df = pivot_df.reset_index()
    # 칼럼명을 train과 같은 이름으로 변경
    pivot_df.rename(columns={'Tumor_Sample_Barcode':'ID'}, inplace=True)
    # SUBCLASS 추가
    pivot_df['SUBCLASS'] = x
    final_df = pivot_df.reindex(columns=train.columns)
    # 외부 데이터를 하나의 csv파일로 생성
    del pivot_df
    combined_df = pd.concat([combined_df, final_df], ignore_index=True)
    # 암종별 csv 생성 및 저장
    final_df.to_csv(f'train_data/{x}.csv', index=False)
    del final_df
combined_df.to_csv(f'train_data/external_data.csv', index=False)

