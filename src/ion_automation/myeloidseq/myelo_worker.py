"""myelo_worker.py: Automated workflow for oncomine myeloid assay."""
__author__      = "Kelsey Zhu"
__copyright__   = "Copyright 2022, Langone Pathlab"
__version__ = "1.0.2"

import requests
import os
import glob
import pandas as pd
import numpy as np
import re
import ast
from time import time
import logging
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import configparser
from ion_automation.myeloidseq.dropout_worker import dropout

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("ion_reporter")
RESULTS = list()

class myeloseq(object):

    def __init__(self, conf_file):
        config = configparser.ConfigParser()
        config.read(conf_file)
        self.HOST = config['DEFAULT']['HOST']
        self.TOKEN = config['DEFAULT']['TOKEN']
        self.UID = config['DEFAULT']['UID']
        self.MYELOSEQ_HOME = config['MYELOSEQ']['MYELOSEQ_HOME']
        self.VAR_HOME = config['MYELOSEQ']['VAR_HOME']
        self.WORK_DIR = config['MYELOSEQ']['WORK_DIR']
        self.DEST_PATH = config['MYELOSEQ']['DEST_PATH']
        self.AA_CODES = config['MYELOSEQ']['AA_CODES']
        self.MYELOSEQ_GENES = config['MYELOSEQ']['MYELOSEQ_GENES']
        self.INCL_FUNCS = config['MYELOSEQ']['INCL_FUNCS'].split(",")
        self.EXCL_CALLS = config['MYELOSEQ']['EXCL_CALLS'].split(",")
        self.LOCATIONS = config['MYELOSEQ']['LOCATIONS'].split(",")
        self.STYLE_CSS = config['MYELOSEQ']['STYLE_CSS']
        self.VAR_TYPES = config['SOLID']['VAR_TYPES'].split(";")

        self.codon_df = pd.read_csv(self.AA_CODES, sep=',')
        self.myeloseq_68genes = pd.read_excel(self.MYELOSEQ_GENES, engine='openpyxl', sheet_name='68-gene', skiprows=0)
        self.nyu_myeloseq_50genes = pd.read_excel(self.MYELOSEQ_GENES, engine='openpyxl', sheet_name='50-gene', skiprows=0)

        self._workbook = None
        self._dropout = dropout(conf_file)

    @property
    def workbook(self):
        return self._workbook

    @workbook.setter
    def workbook(self, value):
        self._workbook = value

    def __del__(self):
        self._workbook = None

    def get_coverage(self, row):
        try:
            m = re.search(r'.*;DP=([0-9]+);.*', row['INFO'])
            return m.group(1)
        except:
            return None

    def get_AF(self, row):
        try:
            m = re.search(r'AF=(.+);AO=.*;TYPE=(.*);VARB=(.*);HS;.*', row['INFO'])
            af_list = m.group(1).split(",")
            if len(af_list) > 1:
                idx = af_list.index(max(af_list))
                return "%s:%s:%s:%s" %(row['ALT'].split(",")[idx], max(af_list),
                                       m.group(2).split(",")[idx], m.group(3).split(",")[idx])
            else:
                return "%s:%s:%s:%s" %(row['ALT'],m.group(1),m.group(2), m.group(3))
        except:
            return None

    def get_func_row(self, row):
        try:
            m = re.search(r'AF=.*;HS;FUNC=(.*)', row['INFO'])
            func = ast.literal_eval(m.group(1))[0]
            return "%s:%s:%s" %(func['transcript'], func['gene'],func['exon'])
        except:
            return None

    def get_tumor_AF(self, row):
        try:
            alt_allele = row['genotype'].split("/")[1]
            if alt_allele in row['allele_frequency_%']:
                alt_allele = row['genotype'].split("/")[1]
                m = re.search(r'(\d+.\d+)(.*)',row['allele_frequency_%'].split("%s="%alt_allele)[1])
                return m.group(1)
            else:
                return row['allele_frequency_%']
        except:
            return row['allele_frequency_%']

    def is_artifact(self, row):
        if row['type'] == 'SNV':
            ref, alt = tuple(row['genotype'].split("/"))
            if ref != alt and row['ref'] not in row['genotype']:
                return True
        return False

    def get_ExAC_info(self, row):
        try:
            m = re.search(r'AMAF=(.+):GMAF=(.+):EMAF=(.+)', row['5000Exomes'])
            return "%s:%s:%s"%(m.group(1), m.group(2), m.group(3))
        except:
            return "NA:NA:NA"

    def get_read_depth(self, row):
        try:
            return row['allele_coverage'].split(",")[1].split("=")[1]
        except:
            return None

    def get_alt_maf(self, row):
        try:
            return min(row['maf'].split(":"))
        except:
            return row['maf']

    def is_hotspot(self, row):
        try:
            return "Hotspot" in row['Oncomine Variant Annotator v3.2']
        except:
            return False

    def is_mnv(self, row):
        try:
            return "[" in row['protein'] or "," in row['protein'] or ";" in row['protein']
        except:
            return False

    def get_download_link(self, sample):
        try:
            headers = {'Content-Type': "application/x-www-form-urlencoded", 'Authorization': self.TOKEN}
            i = 1
            res = None
            params = {'name': "%s_v%s" % (sample,i)}
            r = requests.get('https://%s/api/v1/getvcf' % self.HOST, headers=headers, params=params, verify=False)
            while str(r) == "<Response [200]>":
                res = r
                i += 1
                params = {'name': "%s_v%s" % (sample, i)}
                r = requests.get('https://%s/api/v1/getvcf' % self.HOST, headers=headers, params=params, verify=False)
            return str(res.json()[0]['data_links']).split("filePath=")[1]
        except:
            return None

    def get_tsv_file(self, sample):
        try:
            download_link = self.get_download_link(sample)
            if download_link:
                cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, download_link,
                                                            os.path.join(self.DEST_PATH,"downloads"))
                os.system(cp_cmd)
                zip_file = os.path.basename(download_link)
                unzip_cmd = 'echo %s | sudo unzip -q -o %s'%(self.UID, os.path.join(self.DEST_PATH,
                                                                "downloads",zip_file))
                os.system(unzip_cmd)
                sample_pair = os.path.dirname(download_link).split("/")[6]
                file_path = os.path.join(self.VAR_HOME,sample_pair)
                return glob.glob(os.path.join(file_path, "%s*-full.tsv" % sample_pair))[0],\
                       glob.glob(os.path.join(file_path, "%s*_Non-Filtered_*.vcf" %sample_pair))[0], \
                       glob.glob(os.path.join(file_path, "%s*_Non-Filtered_*-oncomine.tsv" % sample_pair))[0]
        except:
            return None, None, None

    def clean_up(self):
        rm_downloads_cmd = 'echo %s | sudo -S rm -f %s/*.zip' % (self.UID, os.path.join(self.DEST_PATH, "downloads"))
        os.system(rm_downloads_cmd)

    def copy_files(self, run_id):
        logger.info("Generating QC plots...")
        QC_plot_cmd = "Rscript MyeloSeq-QC-plot.R"
        os.system(QC_plot_cmd)
        logger.info("Copying the QC plots over to the Z drive...")
        mkdir_cmd = 'echo %s | sudo -S mkdir -p %s' % (self.UID, os.path.join(self.DEST_PATH,
                                                                    "reports/%s" % run_id))
        os.system(mkdir_cmd)
        plot_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "*.pdf",
                                            os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(plot_cp_cmd)
        csv_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "*.csv",
                                            os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(csv_cp_cmd)
        logger.info("Copying the report over to the Z drive")
        cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "%s.xlsx" % run_id,
                                                os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(cp_cmd)
        if os.path.exists("%s-dropouts.html" % run_id):
            html_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "%s-dropouts.html" % run_id,
                                                os.path.join(self.DEST_PATH, "reports/%s" % run_id))
            os.system(html_cmd)
        sc_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "sc_filtered_variants.tsv",
                                        os.path.join(self.DEST_PATH, "reports/%s/%s_SC_Variants.tsv" % (
                                                            run_id, run_id)))
        os.system(sc_cp_cmd)
        if os.path.exists("sc2_filtered_variants.tsv"):
            sc2_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "sc2_filtered_variants.tsv",
                                            os.path.join(self.DEST_PATH, "reports/%s/%s_SC2_Variants.tsv" % (
                                                                run_id, run_id)))
            os.system(sc2_cp_cmd)
        # copy varaint tsv/vcf files over to the Z drive
        logger.info("Copying variant files over to the Z drive")
        if os.path.exists(os.path.join(self.DEST_PATH, "downloads/%s" % run_id)):
            rm_cmd = 'echo %s | sudo -S rm -rf %s' % (self.UID,
                                            os.path.join(self.DEST_PATH, "downloads/%s" % run_id))
            os.system(rm_cmd)
        cp_var_cmd = 'echo %s | sudo -S cp -r %s %s' % (self.UID, self.VAR_HOME,
                    os.path.join(self.DEST_PATH, "downloads/%s" % run_id))
        os.system(cp_var_cmd)
        rm_var_cmd = 'echo %s | sudo -S rm -rf %s/*'%(self.UID, self.VAR_HOME)
        os.system(rm_var_cmd)
        rm_outfiles_cmd = 'echo %s | sudo -S rm -f %s/*.csv %s/*.pdf %s/*.xlsx %s/*.tsv %s/*.html'%(self.UID,
                self.MYELOSEQ_HOME,self.MYELOSEQ_HOME, self.MYELOSEQ_HOME, self.MYELOSEQ_HOME, self.MYELOSEQ_HOME)
        os.system(rm_outfiles_cmd)
        self.clean_up()

    def write_to_excel(self, df):
        try:
            df = df.fillna("")
            run_id = df.iloc[0,0]
            logger.info("run ID: %s"%run_id)
            # Create a Pandas XlsxWriter engine.
            writer = pd.ExcelWriter("%s.xlsx"%run_id, engine='xlsxwriter')
            df.to_excel(writer, sheet_name='DNA', index=False)
            workbook = writer.book

            # Light red fill with dark red text.
            format11 = workbook.add_format({'bg_color': '#FFC7CE',
                                           'font_color': '#9C0006'})

            # Light yellow fill with dark yellow text.
            format22 = workbook.add_format({'bg_color': '#FFEB9C',
                                           'font_color': '#9C6500'})

            # Green fill with dark green text.
            format33 = workbook.add_format({'bg_color': '#C6EFCE',
                                           'font_color': '#006100'})

            worksheet = writer.sheets['DNA']

            # Add some cell formats.
            format1 = workbook.add_format({'num_format': '#,##0.00'})
            format2 = workbook.add_format({'num_format': '0%'})

            # Define the formats
            format6 = workbook.add_format({'bg_color': '#B9D3EE', 'border_color': '#BFBFBF', 'border': 1})  # dark blue
            format4 = workbook.add_format({'bg_color': '#FFFFFF', 'border_color': '#BFBFBF', 'border': 1})  # white
            format5 = workbook.add_format({'bg_color': '#E5F8F3', 'border_color': '#BFBFBF', 'border': 1})  # light green
            format3 = workbook.add_format({'bg_color': '#DBE8F0', 'border_color': '#BFBFBF', 'border': 1})  # light blue

            for row in range(1, df.shape[0] + 1):
                worksheet.write(row, 0, df.iloc[row-1, 0], format5)

            print (list(self.nyu_myeloseq_50genes['NYU.Myeloid']))
            print(list(self.myeloseq_68genes['Oncomine.Myeloid']))

            # Set the color for the starting row
            current_color = 'Dark'
            # Format the 1st row

            for column in range(1, df.shape[1]):  # format the first 2 columns
                worksheet.write(1, column, df.iloc[0, column], format3)

            # Start formatting from the 2nd row until the end of the df
            for row in range(2,df.shape[0]+1):
                # if the id of the row is the same as the id of the previous row
                if df.iloc[row - 1, 1] == df.iloc[row - 2, 1]:
                    if current_color == 'Dark':
                        sample_format = format3
                    elif current_color == 'Light':
                        sample_format = format4
                    for column in range(1, df.shape[1]):  # format the first 2 columns
                        worksheet.write(row, column, df.iloc[row - 1, column], sample_format)
                # if it's different than that of the previous row switch the colors
                else:
                    if current_color == 'Dark':
                        current_color = 'Light'
                    elif current_color == 'Light':
                        current_color = 'Dark'
                    for column in range(1, df.shape[1]):  # format the first 2 columns
                        worksheet.write(row, column, df.iloc[row - 1, column],
                                        format3 if current_color == 'Dark' else format4)
            # Set the column width and format.
            worksheet.set_column('B:B', 22, format1)
            worksheet.set_column('D:D', 18, format1)
            worksheet.set_column('N:N', 18, format1)
            worksheet.set_column('P:P', 18, format1)

            # Set the format but not the column width.
            worksheet.set_column('L:L', None, format2)

            for row in range(1, df.shape[0] + 1):
                if df.iloc[row-1, 5] != 'FUSION' and df.iloc[row-1, 4] != 'NA':
                    print(df.iloc[row - 1, 4], df.iloc[row - 1, 5])
                    if df.iloc[row-1, 4] in list(self.nyu_myeloseq_50genes['NYU.Myeloid']):
                        worksheet.write(row, 4, df.iloc[row - 1, 4], format33)
                    elif df.iloc[row-1, 4] in list(self.myeloseq_68genes['Oncomine.Myeloid']):
                        worksheet.write(row, 4, df.iloc[row-1, 4], format22)
                    else:
                        worksheet.write(row, 4, df.iloc[row-1, 4], format11)

            # Close the Pandas Excel writer and output the Excel file.
            writer.save()
        except:
            raise

    def get_vcf_fusion_key(self, row):
        try:
            m = re.search(r'(\w*)\[*\]*(chr\d+):(\d+)\[*\]*(\w)*',row['ALT'])
            return row['ID'].split("_")[0]
        except Exception as e:
            return None

    def get_tsv_fusion_key(self, row):
        try:
            if row['type'] == 'RNAExonVariant':
                if row['gene'] == 'EGFR|EGFR':
                    return 'EGFR-EGFR.E1E8.DelPositive.1'
                elif row['gene'] == 'MET|MET':
                    return 'MET-MET.M13M15'
            else:
                return row['# locus'].split("_")[1]
        except:
            return row['# locus']

    def get_fusion_read_counts(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*);READ_COUNT=(.+);GENE_NAME=(.+);RPM=(.+);NORM_COUNT=(.+)',
                            row['INFO']);
            return m.group(2)
        except:
            return 0

    def get_fusion_RPM(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*);READ_COUNT=(.+);GENE_NAME=(.+);RPM=(.+);NORM_COUNT=(.+)',
                          row['INFO']);
            return m.group(4)
        except:
            return 0

    def get_codon_letter(self, code):
        try:
            idx = self.codon_df.loc[self.codon_df['Codon'] == code].index
            return self.codon_df.iloc[idx]['Letter'].values[0]
        except:
            return None

    def get_codon_code(self, AA_change):
        try:
            return AA_change.replace("Ala", "A") \
                             .replace("Arg", "R") \
                             .replace("Asn", "N") \
                             .replace("Asp", "D") \
                             .replace("Cys", "C") \
                             .replace("Gln", "Q") \
                             .replace("Glu", "E") \
                             .replace("Gly", "G") \
                             .replace("His", "H") \
                             .replace("Ile", "I") \
                             .replace("Leu", "L") \
                             .replace("Lys", "K") \
                             .replace("Met", "M") \
                             .replace("Phe", "F") \
                             .replace("Pro", "P") \
                             .replace("Ser", "S") \
                             .replace("Thr", "T") \
                             .replace("Trp", "W") \
                             .replace("Tyr", "Y") \
                             .replace("Val", "V") \
                             .replace("Ter", "X")
        except:
            return AA_change

    def get_AA_Change(self, row):
        if isinstance(row['Amino Acid Change'], str):
            aa_change = row['Amino Acid Change'].split("|")[0]
            return self.get_codon_code(aa_change)
        else:
            return 'NA'

    def empty_row(self,row):
        try:
            if row['Type'] in ['SNV', 'INDEL'] and row['% Frequency'].astype(float) > 0:
                return False
            elif row['Type'] == 'RNAExonVariant':
                if row['Read Counts'] == 'NA' and row['Read/M'] == 'NA':
                    return True
                elif row['Genes'] in (['BRAF','EGFR']):
                    return True
                elif row['Read Counts'].isnull().values.any() or \
                    row['Read/M'].isnull().values.any():
                    return True
            elif row['Type'] in ['SNV', 'INDEL'] and str(row['Variant Effect']) == 'nan' \
                    and str(row['Transcript']) == 'nan' and str(row['Coding']) == 'nan':
                return True
            else:
                return False
        except:
            return False

    def select_gene(self, row):
        try:
            if row['type'] == 'FUSION':
                genes = row['gene'].split("|")
                exons = row['exon'].split("|")
                return "%s(%s) - %s(%s)" %(genes[0], exons[0], genes[1], exons[1])
            if row['gene'].endswith("|"):
                return row['gene'].replace("|","")
            elif not "|" in row['gene']:
                return row['gene']
            else:
                row_genes = row['gene'].strip().split("|")
                sel_gene = list(set(row_genes) & set(list(self.myeloseq_68genes['Oncomine.Myeloid'])))
                return sel_gene[0]
        except:
            return row['gene']

    def remove_function_bar(self, row):
        try:
            if row['function'].startswith("|") or row['function'].endswith("|"):
                return row['function'].replace("|", "")
            else:
                return row['function']
        except:
            return row['function']

    def get_location(self, row):
        try:
            return row['location'].split(":")[1]
        except:
            return row['location']

    def process_sample(self, args):
        sample, run_id, bar_code, logger = tuple(args)
        logger.info("start processing %s from %s"%(sample,run_id))
        filtered_tsv, filtered_vcf, oncomine_tsv = self.get_tsv_file(sample)
        if filtered_tsv:
            if filtered_vcf:
                try:
                    ion_fusions = pd.read_csv(filtered_vcf, sep="\t", skiprows=180)
                    ion_fusions = ion_fusions.loc[ion_fusions['FILTER'].isin(["PASS","."])]
                    ion_fusions['fusion_key'] = ion_fusions.apply(lambda x: self.get_vcf_fusion_key(x), axis=1)
                    ion_fusions['Read Counts'] = ion_fusions.apply(lambda x: self.get_fusion_read_counts(x), axis=1)
                    ion_fusions['Read/M'] = ion_fusions.apply(lambda x: self.get_fusion_RPM(x), axis=1)
                    ion_fusions['Read/M'] = ion_fusions['Read/M'].astype(float)
                    ion_fusions['Read/M'] = ion_fusions['Read/M'].apply(np.int64)
                    ion_fusions = ion_fusions[['fusion_key',"Read Counts","Read/M"]]
                    ion_fusions = ion_fusions.loc[ion_fusions["Read Counts"].astype(int) >= 15]
                    print (ion_fusions.to_string())
                except Exception as e:
                    logger.error(str(e))
                    ion_fusions = None
            ion_variants = pd.read_csv(filtered_tsv, sep="\t", skiprows=2)

            # filters that are applied
            ion_variants = ion_variants.loc[(ion_variants['filter'].isin(['PASS','GAIN','.','LOSS'])) &
                                                            (~ion_variants['type'].isin(self.EXCL_CALLS))]
            ion_variants['tumor_AF'] = ion_variants.apply(lambda x: self.get_tumor_AF(x), axis=1)
            ion_variants['ExAC_info'] = ion_variants.apply(lambda x: self.get_ExAC_info(x), axis=1)
            ion_variants['DP'] = ion_variants.apply(lambda x: self.get_read_depth(x), axis=1)
            ion_variants['MAF'] = ion_variants.apply(lambda x: self.get_alt_maf(x), axis=1)
            ion_variants['HS'] = ion_variants.apply(lambda x: "yes" if self.is_hotspot(x) else "", axis=1)
            ion_variants['artifact'] = ion_variants.apply(lambda x: self.is_artifact(x), axis=1)
            ion_variants['function'] = ion_variants.apply(lambda x: self.remove_function_bar(x), axis=1)
            ion_variants['splice_site'] = ion_variants.apply(lambda x: self.get_location(x), axis=1)
            ion_variants['MNV'] = ion_variants.apply(lambda x: "MNV" if self.is_mnv(x) else "", axis=1)

            ion_variants_snv = ion_variants.loc[ (ion_variants['HS'] == 'yes') | (ion_variants['type'] == 'FLT3ITD') |
                                    (ion_variants['MNV'] == 'MNV') |
                                    ((ion_variants['type'].isin(self.VAR_TYPES))
                                    & (ion_variants['function'].isin(self.INCL_FUNCS) | ion_variants['splice_site'].isin(self.LOCATIONS))
                                    & ((ion_variants['MAF'].isnull()) |(ion_variants['MAF'].astype(float) <= 0.01))
                                    & (~ion_variants['artifact']))]

            ion_variants_snv = ion_variants_snv.loc[ ((is_numeric_dtype(ion_variants_snv['tumor_AF']))
                                                          & (ion_variants_snv['tumor_AF'].astype(float) >= 0.03)) |
                                                     (is_string_dtype(ion_variants_snv['tumor_AF'])) |
                                                     (ion_variants_snv['type'] == 'FLT3ITD')]
            ion_variants_fusion = ion_variants.loc[(ion_variants['type'].isin(['FUSION','RNAExonVariant']))]
            ion_variants = pd.concat([ion_variants_snv, ion_variants_fusion], ignore_index=True)
            ion_variants['fusion_key'] = ion_variants.apply(lambda x: self.get_tsv_fusion_key(x), axis=1)
            #select oncomine genes
            ion_variants['gene'] = ion_variants.apply(lambda x: self.select_gene(x), axis=1)
            ion_variants['# locus'] = ion_variants.apply(lambda x: self.rename_locus(x), axis=1)
            ion_variants['anno_key'] = ion_variants.apply(lambda x: "%s:%s"%(x['gene'],x['# locus']),axis=1)
            #update annotations
            if oncomine_tsv:
                onco_anno = pd.read_csv(oncomine_tsv, sep="\t", skiprows=2)
                onco_anno = onco_anno.loc[(onco_anno['call'] == 'POS')]
                try:
                    onco_anno['anno_key'] = onco_anno.apply(
                            lambda x: "%s:%s:%s" % (x['FUNC1.gene'], x['CHROM'], x['POS']) if x['FUNC1.gene'] != 'CSDE1' else
                                        "%s:%s:%s" % (x['FUNC2.gene'], x['CHROM'], x['POS']), axis=1)
                    NRAS_row = onco_anno.loc[onco_anno['FUNC2.gene'] == 'NRAS']
                    onco_anno.loc[onco_anno['FUNC2.gene'] == 'NRAS', 'FUNC1.transcript'] = NRAS_row.iloc[0]['FUNC2.transcript']
                    onco_anno.loc[onco_anno['FUNC2.gene'] == 'NRAS', 'FUNC1.protein'] = NRAS_row.iloc[0][
                        'FUNC2.protein']
                    onco_anno.loc[onco_anno['FUNC2.gene'] == 'NRAS', 'FUNC1.coding'] = NRAS_row.iloc[0][
                        'FUNC2.coding']
                    onco_anno.loc[onco_anno['FUNC2.gene'] == 'NRAS', 'FUNC1.exon'] = NRAS_row.iloc[0][
                        'FUNC2.exon']
                except:
                    onco_anno['anno_key'] = onco_anno.apply(
                            lambda x: "%s:%s:%s" % (x['FUNC1.gene'], x['CHROM'], x['POS']), axis=1)
                try:
                    onco_anno = onco_anno[['anno_key', 'FUNC1.coding', 'FUNC1.exon', 'FUNC1.function', 'FUNC1.transcript', 'FUNC1.protein']]
                except:
                    onco_anno['FUNC1.coding'] = 'NA'
                    onco_anno['FUNC1.function'] = 'NA'
                    onco_anno['FUNC1.transcript'] = 'NA'
                    onco_anno['FUNC1.protein'] = 'NA'
                    onco_anno = onco_anno[['anno_key', 'FUNC1.coding', 'FUNC1.exon', 'FUNC1.function', 'FUNC1.transcript', 'FUNC1.protein']]
                logger.info(onco_anno.head(3).to_string())

                ion_variants = ion_variants.merge(onco_anno, left_on="anno_key", right_on="anno_key", how="left")
                print(ion_variants.shape)
            # filter out artifacts
            if not sample.startswith("NC-DNA") or not sample.startswith("NTC-DNA") or not sample.startswith("SC-DNA"):
                ion_variants = ion_variants.loc[~((ion_variants['gene'] == 'CBL') &
                                                  (ion_variants['# locus'] == 'chr11:119149355')
                                                  & (ion_variants['type'] == "INDEL") &
                                                  (ion_variants['genotype'] == "TATGATGATGATGATGATGA/TATGATGATGATGATGA") &
                                                 (ion_variants['tumor_AF'].astype(float) >= 0.03) )]
                ion_variants = ion_variants.loc[~((ion_variants['gene'] == 'SH2B3') &
                                                  (ion_variants['# locus'].isin(['chr12:111885351','chr12:111885350'])))]
                ion_variants = ion_variants.loc[~( (ion_variants['type'].isin(["SNV","FLT3ITD"])) &(ion_variants['genotype'] == "C/C") & ( ion_variants['gene'].isin(['ASXL1','FLT3'])))]
            if ion_variants.empty:
                return pd.DataFrame({"Run": run_id, "Sample": sample, "Barcode": bar_code,
                                                   "Locus": 'negative', "Genes": 'NA', "Type": 'NA',
                                                   "Exon": 'NA', "Transcript": 'NA', "Coding": 'NA',
                                                   "Variant Effect": 'NA', 'Genotype': 'NA',
                                                   "% Frequency": "NA", "ExAC_AF": 'NA', "Amino Acid Change": 'NA',
                                                   "Read Counts": 'NA', "Read/M": 'NA', "AMAF": 'NA', "GMAF": 'NA',
                                                    "EMAF": 'NA',"HS":"NA","Length":'NA','Coverage':'NA'},
                                                  index=[0])
            new = ion_variants['ExAC_info'].str.split(":", n=2, expand=True)
            ion_variants.insert(0, "Run", run_id)
            ion_variants.insert(1, "Sample", sample)
            ion_variants.insert(2, "Barcode", bar_code)

            print("Before filters are applied")
            print(ion_variants.head(3).to_string())
            if not ion_fusions.empty:
                ion_variants = ion_variants.merge(ion_fusions, left_on="fusion_key", right_on="fusion_key",how="left")
            else:
                ion_variants['Read Counts'] = 'NA'
                ion_variants['Read/M'] = 'NA'
            ion_variants.drop(columns=['ExAC_info', 'go', '5000Exomes', 'hrun', 'drugbank',
                                               'fusion_presence', 'ratio_to_wild_type',
                                               'norm_count_within_gene', 'filter',
                                               'allele_coverage', 'allele_ratio', 'pvalue','dgv',
                                               'allele_frequency_%', 'MyVariantDefaultDb_hg19',
                                               'phylop', 'pfam', 'location', 'maf',
                                               'sift', 'polyphen', 'grantham', 'normalizedAlt',
                                               'NamedVariants', 'function', 'coding', 'transcript', 'protein', 'exon'], inplace=True)
            ion_variants.rename(columns={'# locus': 'Locus', 'gene': 'Genes', 'FUNC1.exon': 'Exon',
                                         'FUNC1.transcript': 'Transcript', 'genotype': 'Genotype',
                                         'type': 'Type', 'FUNC1.coding': 'Coding', 'FUNC1.function': 'Variant Effect',
                                         'ref': 'Ref', 'tumor_AF': '% Frequency','length':'Length',
                                         'FUNC1.protein': 'Amino Acid Change', 'coverage': 'Coverage',
                                         'exac': 'ExAC', 'MAF':'ExAC_AF'},
                                        inplace=True)

            ion_variants['AMAF'] = new[0]
            ion_variants["GMAF"] = new[1]
            ion_variants["EMAF"] = new[2]
            ion_variants = ion_variants[
                ["Run", "Sample", "Barcode", "Locus", "Genes", "Type", "Exon", "Transcript", "Coding", "Variant Effect", "Genotype",
                 "% Frequency", "ExAC_AF", "Amino Acid Change", "Read Counts", "Read/M", "AMAF", "GMAF", "EMAF","HS","Length","Coverage"]]
            logger.info("%s processed"%sample)
            logger.info("After filters are applied")
            ion_variants = ion_variants.drop_duplicates(subset=['Sample','Barcode','Locus','Genes','Type'])
            ion_variants['bad'] = ion_variants.apply(lambda x: self.empty_row(x), axis=1)
            ion_variants = ion_variants.loc[ion_variants['bad'] != 1]
            logger.info(ion_variants.head(3).to_string())
            ion_variants.drop(columns=['bad'], inplace=True)
            if ion_variants.empty:
                return pd.DataFrame({"Run": run_id, "Sample": sample, "Barcode": bar_code,
                                 "Locus": 'negative', "Genes": 'NA', "Type": 'NA',
                                 "Exon": 'NA', "Transcript": 'NA', "Coding": 'NA',
                                 "Variant Effect": 'NA', 'Genotype': 'NA',
                                 "% Frequency": "NA", "ExAC_AF": 'NA', "Amino Acid Change": 'NA',
                                 "Read Counts": 'NA', "Read/M": 'NA', "AMAF": 'NA', "GMAF": 'NA', "EMAF": 'NA',
                                 "HS": "NA","Length":"NA","Coverage":"NA"},
                                index=[0])
            else:
                return ion_variants
        else:
            logger.warning("%s tsv file is not found or errors occur while processing the tsv"%sample)

    def rename_genes(self, row):
        try:
            if row['Type'] == 'FUSION':
                genes = row['Genes'].split("|")
                exons = row['Exon'].split(",")
                return "%s(%s) - %s(%s)" %(genes[0],exons[0],genes[1],exons[1])
            else:
                if row['Genes'].endswith("|"):
                    return row['Genes'].replace("|", "")
                else:
                    return row['Genes'].replace("|", ",")
        except:
            return row['Genes']

    def rename_locus(self, row):
        try:
            if row['type'] == 'FUSION':
                if "_" in row['# locus']:
                    return row['# locus'].split("_")[0]
            return row['# locus']
        except Exception as e:
            logger.error(str(e))
            return row['# locus']

    def reform(self, colname,row):
        try:
            if row[colname].endswith("|"):
                return row[colname].replace("|", "")
            else:
                return row[colname].replace("|", ",")
        except:
            return row[colname]

    def start(self):
        self.clean_up()
        RESULTS = list()
        try:
            ts = time()
            # read oncomine sample sheet
            logger.info(self.workbook)
            header = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', nrows=1, skiprows=2)
            run_id = list(header.columns.values)[2].replace("-DNA","")
            sample_sheet = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', skiprows=5)
            sample_sheet['sample_id'] = sample_sheet['Accession #'] + "-" + sample_sheet['DNA #']
            logger.info(sample_sheet.to_string())  # shows headers with top 5 rows

            sc_sample_name = list(sample_sheet['sample_id'])[0]
            sc2_sample_name = sc_sample_name.replace("m-Seraseq","m2-Seraseq")
            logger.info("SC sample: %s" %sc_sample_name)
            for sample, barcode in zip(list(sample_sheet['sample_id']), list(sample_sheet['Bar code'])):
                try:
                    if sample == "" or sample == None or str(sample) == 'nan': continue
                    RESULTS.append(self.process_sample([sample,run_id,barcode,logging.getLogger(sample)]))
                except:
                    pass
            if RESULTS and len(RESULTS) > 0:
                df = pd.concat(RESULTS)
                sc_df = df.loc[df['Sample'] == sc_sample_name]
                sc_df['Read Counts'] = sc_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read Counts'], axis=1)
                sc_df['Read/M'] = sc_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read/M'], axis=1)
                sc_df.to_csv("sc_filtered_variants.tsv", index = False, sep = "\t")
                sc2_df = df.loc[df['Sample'] == sc2_sample_name]
                if not sc2_df.empty:
                    sc2_df['Read Counts'] = sc2_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read Counts'], axis=1)
                    sc2_df['Read/M'] = sc2_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read/M'], axis=1)
                    sc2_df.to_csv("sc2_filtered_variants.tsv", index = False, sep = "\t")
                sample_df = df.loc[~(df['Sample'].isin([sc_sample_name,sc2_sample_name]))]
                logger.info(sample_df.to_string())
                self.write_to_excel(sample_df)
                add_df = df.loc[~(df['Sample'].str.contains("SC-DNA|NC-DNA"))]
                add_df = add_df.drop(columns=['Run','AMAF','GMAF','EMAF','Read Counts','Read/M','ExAC_AF'])
                add_df['AA'] = add_df.apply(lambda x: self.get_AA_Change(x), axis=1)
                add_df = add_df.rename(columns={'Barcode':'BC','Variant Effect':'Variant.Effect','Amino Acid Change':'Amino.Acid.Change',
                                       '% Frequency':'Frequency','HS':'Info'})
                add_df = add_df[['Sample','BC','Locus','Genes','Type','Exon','Transcript','Coding','Variant.Effect',
                                 'Genotype','Info','Length','Frequency','Amino.Acid.Change','AA','Coverage']]
                print(add_df.head(3).to_string())
                add_df.to_csv("%s.csv" % run_id, index=False, sep=",")
                self._dropout.workbook = self.workbook
                self._dropout.start()
                # Copy files over to the share drive
                self.copy_files(run_id)
                logger.info('Took %s seconds to process samples', time() - ts)

        except Exception as e:
             logger.error(str(e))

