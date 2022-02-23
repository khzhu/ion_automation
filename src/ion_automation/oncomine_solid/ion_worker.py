import requests
import os
import glob
import pandas as pd
import re
import ast
from time import time
import logging
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import configparser

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("ion_reporter")

class oncomine_solid(object):

    def __init__(self, conf_file):
        config = configparser.ConfigParser()
        config.read(conf_file)
        self.HOST = config['DEFAULT']['HOST']
        self.TOKEN = config['DEFAULT']['TOKEN']
        self.UID = config['DEFAULT']['UID']
        self.HOME_DIR = config['SOLID']['HOME_DIR']
        self.VAR_DIR = config['SOLID']['VAR_DIR']
        self.WORK_DIR = config['SOLID']['WORK_DIR']
        self.DEST_PATH = config['SOLID']['DEST_PATH']
        self.INCL_FUNCS = config['SOLID']['INCL_FUNCS'].split(",")
        self.LOCATIONS = config['SOLID']['LOCATIONS'].split(",")
        self.EXCL_CALLS = config['SOLID']['EXCL_CALLS'].split(",")
        self.VAR_TYPES = config['SOLID']['VAR_TYPES'].split(";")
        self.CNV_GENES = config['SOLID']['CNV_GENES'].split(",")
        self.MUT_GENES = config['SOLID']['MUT_GENES'].split(",")
        self.FUSION_GENES = config['SOLID']['FUSION_GENES'].split(",")
        self._workbook = None

    @property
    def workbook(self):
        return self._workbook

    @workbook.setter
    def workbook(self, value):
        self._workbook = value

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
                #CTTTTTTTTTTTTGA=6.61,CTTTTTTTTTTTTGAT=4.57,CTTTTTTTTTTTTTGA=75.23,CTTTTTTTTTTTTTGAT=13.59
                alt_allele = row['genotype'].split("/")[1]
                m = re.search(r'(\d+.\d+)(.*)',row['allele_frequency_%'].split("%s="%alt_allele)[1])
                return m.group(1)
            else:
                return row['allele_frequency_%']
        except:
            print(row['allele_frequency_%'])
            return row['allele_frequency_%']

    def is_artifact(self, row):
        if row['type'] == 'SNV' and row['ref'] not in row['genotype']:
            return True
        else:
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

    def get_copy_number(self, row):
        try:
            if row['type'] == 'CNV' and row['gene'] in self.CNV_GENES:
                return row['iscn'].split("x")[1]
        except:
            return None

    def get_CI_copy_number(self, row):
        try:
            if row['type'] == 'CNV' and row['gene'] in self.CNV_GENES:
                m = re.search(r'5%:(.+),95%:(.+)', row['confidence'])
                return m.group(1)
        except:
            return None

    def oncomine_in(self, row):
        try:
            return row['Oncomine Variant Annotator v3.2'].str.len() > 0
        except:
            return False

    def is_hotspot(self, row):
        #Oncomine=[transcript=NM_020975.5;normalizedRef=C;normalizedPos=43617397;normalizedAlt=T;type=Hotspot;...]
        try:
            return "Hotspot" in row['Oncomine Variant Annotator v3.2']
        except:
            return False

    def get_download_link(self, sample):
        try:
            params = {'name': "%s_v1"%sample}
            headers = { 'Content-Type':"application/x-www-form-urlencoded", 'Authorization' : self.TOKEN}
            r = requests.get('https://%s/api/v1/getvcf'%self.HOST, headers=headers, params=params, verify=False)
            return str(r.json()[0]['data_links']).split("filePath=")[1]
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
                unzip_cmd = 'echo %s | sudo unzip -q -o %s'%(self.UID, os.path.join(self.DEST_PATH,"downloads",zip_file))
                os.system(unzip_cmd)
                sample_pair = os.path.dirname(download_link).split("/")[6]
                file_path = os.path.join(self.VAR_DIR,sample_pair)
                return glob.glob(os.path.join(file_path, "%s*-full.tsv" % sample_pair))[0],\
                       glob.glob(os.path.join(file_path, "%s*_Non-Filtered_*.vcf" %sample_pair))[0]
        except:
            return None, None

    def copy_files(self, run_id):
        logger.info("Generating QC plots...")
        QC_plot_cmd = "Rscript Oncosolid_QC_plot.R"
        os.system(QC_plot_cmd)
        logger.info("Copying the QC plots over to the Z drive...")
        mkdir_cmd = 'echo %s | sudo -S mkdir -p %s' % (self.UID, os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(mkdir_cmd)
        plot_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "*.pdf", os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(plot_cp_cmd)
        csv_cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "*.csv", os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(csv_cp_cmd)
        logger.info("Copying the report over to the Z drive")
        cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "%s.xlsx" % run_id,
                                                            os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(cp_cmd)
        cp_cmd = 'echo %s | sudo -S cp -f %s %s' % (self.UID, "sc_filtered_variants.tsv",
                                                            os.path.join(self.DEST_PATH, "reports/%s/%s_SC_Variants.tsv" % (
                                                            run_id, run_id)))
        os.system(cp_cmd)
        # copy varaint tsv/vcf files over to the Z drive
        logger.info("Copying variant files over to the Z drive")
        cp_var_cmd = 'echo %s | sudo -S cp -r %s %s' % (self.UID, self.VAR_DIR,
            os.path.join(self.DEST_PATH, "downloads/%s" % run_id))
        os.system(cp_var_cmd)
        rm_var_cmd = 'echo %s | sudo -S rm -rf %s/*'%(self.UID, self.VAR_DIR)
        os.system(rm_var_cmd)
        rm_outfiles_cmd = 'echo %s | sudo -S rm -f %s/*.csv %s/*.pdf %s/*.xlsx %s/*.tsv'%(self.UID, self.HOME_DIR,
                                                                    self.HOME_DIR,self.HOME_DIR, self.HOME_DIR)
        os.system(rm_outfiles_cmd)
        rm_downloads_cmd = 'echo %s | sudo -S rm -f %s/*.zip' % (self.UID, os.path.join(self.DEST_PATH, "downloads"))
        os.system(rm_downloads_cmd)

    def write_to_excel(self, df):
        try:
            df = df.fillna("")
            run_id = df.iloc[0,0]
            logger.info("run ID: %s"%run_id)
            # Create a Pandas XlsxWriter engine.
            writer = pd.ExcelWriter("%s.xlsx"%run_id, engine='xlsxwriter')
            df.to_excel(writer, sheet_name='DNA', index=False)
            workbook = writer.book
            #format workbook rows/cells
            red_format = workbook.add_format({'bg_color': 'red'})
            green_format = workbook.add_format({'bg_color': 'green'})
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
            # Close the Pandas Excel writer and output the Excel file.
            writer.save()

            # Copy files over to the share drive
            self.copy_files(run_id)
        except:
            raise

    def get_vcf_fusin_key(self, row):
        try:
            if "_" in row['ID']:
                return row['ID'].split("_")[0]
            else:
                return row['ID']
        except Exception as e:
            return None

    def get_location(self, row):
        try:
            return row['location'].split(":")[1]
        except:
            return row['location']

    def get_tsv_fusin_key(self, row):
        try:
            #chr1:156100564-chr1:156844698_LMNA-NTRK1.L2N11
            if row['type'] == 'RNAExonVariant':
                if row['gene'] == 'EGFR|EGFR':
                    return 'EGFR-EGFR.E1E8.DelPositive.1'
                elif row['gene'] == 'MET|MET':
                    return 'MET-MET.M13M15'
                elif row['gene'] == 'MET':
                    return 'MET.M13M14M15.WT'
            else:
                return row['# locus'].split("_")[1]
        except:
            return None

    def get_fusion_read_counts(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*)(.*);READ_COUNT=(\d+);(.*)',row['INFO'])
            return m.group(3)
        except:
            return 0

    def get_fusion_RPM(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*)(.*);RPM=(\d+.\d+);(.*)',row['INFO'])
            return m.group(3)
        except:
            return 0

    def process_sample(self, args):
        sample, run_id, bar_code, logger = tuple(args)
        logger.info("start processing %s from %s"%(sample,run_id))
        filtered_tsv, filtered_vcf = self.get_tsv_file(sample)
        if filtered_tsv:
            if filtered_vcf:
                try:
                    ion_fusions = pd.read_csv(filtered_vcf, sep="\t", skiprows=210)
                    ion_fusions = ion_fusions.loc[(ion_fusions['FILTER'].isin(['PASS','.'])) &
                                                  (ion_fusions['ALT'] != "<CNV>")]
                    ion_fusions['fusion_key'] = ion_fusions.apply(lambda x: self.get_vcf_fusin_key(x), axis=1)
                    ion_fusions['Read Counts'] = ion_fusions.apply(lambda x: self.get_fusion_read_counts(x), axis=1)
                    ion_fusions['Read/M'] = ion_fusions.apply(lambda x: self.get_fusion_RPM(x), axis=1)
                    ion_fusions = ion_fusions[['fusion_key',"Read Counts","Read/M"]]
                    ion_fusions = ion_fusions.loc[(ion_fusions["Read Counts"].astype(int) >= 15)]
                    print (ion_fusions.to_string())
                except Exception as e:
                    logger.error(str(e))
                    ion_fusions = None

            ion_variants = pd.read_csv(filtered_tsv, sep="\t", skiprows=2)

            # filters that are applied
            ion_variants = ion_variants.loc[(ion_variants['filter'].isin(['PASS','GAIN','LOSS','.'])) &
                                                            (~ion_variants['type'].isin(self.EXCL_CALLS))]
            ion_variants['tumor_AF'] = ion_variants.apply(lambda x: self.get_tumor_AF(x), axis=1)
            ion_variants['ExAC_info'] = ion_variants.apply(lambda x: self.get_ExAC_info(x), axis=1)
            ion_variants['DP'] = ion_variants.apply(lambda x: self.get_read_depth(x), axis=1)
            ion_variants['Copy Number'] = ion_variants.apply(lambda x: self.get_copy_number(x), axis=1)
            ion_variants['5% CI'] = ion_variants.apply(lambda x: self.get_CI_copy_number(x), axis=1)
            ion_variants['MAF'] = ion_variants.apply(lambda x: self.get_alt_maf(x), axis=1)
            ion_variants['HS'] = ion_variants.apply(lambda x: "yes" if self.is_hotspot(x) else "", axis=1)
            ion_variants['ONCOIN'] = ion_variants.apply(lambda x: "in" if self.oncomine_in(x) else "", axis=1)
            ion_variants['artifact'] = ion_variants.apply(lambda x: self.is_artifact(x), axis=1)
            ion_variants['splice_site'] = ion_variants.apply(lambda x: self.get_location(x), axis=1)

            ion_variants_snv = ion_variants.loc[ (ion_variants['HS'] == 'yes') | (ion_variants['ONCOIN'] == 'in') |
                                    ((ion_variants['type'].isin(self.VAR_TYPES))
                                    & (ion_variants['function'].isin(self.INCL_FUNCS) | ion_variants['splice_site'].isin(self.LOCATIONS))
                                    & ((ion_variants['MAF'].isnull()) |(ion_variants['MAF'].astype(float) <= 0.01))
                                    & (~ion_variants['artifact']))]

            ion_variants_snv = ion_variants_snv.loc[ (ion_variants_snv['tumor_AF'].astype(float) >= 3) &
                                                     ion_variants_snv['gene'].isin(self.MUT_GENES)]
            ion_variants_cnv = ion_variants.loc[(ion_variants['type'] == 'CNV') &
                                                (ion_variants['gene'].isin(self.CNV_GENES)) &
                                                                (ion_variants['Copy Number'].astype(float) >= 5) &
                                                                (ion_variants['5% CI'].astype(float) >= 4)]
            ion_variants_fusion = ion_variants.loc[(ion_variants['type'].isin(['FUSION','RNAExonVariant']))]
            ion_variants = pd.concat([ion_variants_snv, ion_variants_cnv, ion_variants_fusion],
                                                ignore_index=True)
            new = ion_variants['ExAC_info'].str.split(":", n=2, expand=True)
            ion_variants.insert(0, "Run", run_id)
            ion_variants.insert(1, "Sample", sample)
            ion_variants.insert(2, "Barcode", bar_code)
            ion_variants['fusion_key'] = ion_variants.apply(lambda x: self.get_tsv_fusin_key(x), axis=1)
            print("Before filtering")
            print(ion_variants.head(3).to_string())
            if not ion_fusions.empty:
                ion_variants = ion_variants.merge(ion_fusions, left_on="fusion_key", right_on="fusion_key",how="left")
                ion_variants.drop(ion_variants[(ion_variants['type'] == 'RNAExonVariant') &
                                                (ion_variants['Read Counts'].isnull())].index, inplace = True)
            else:
                ion_variants['Read Counts'] = 'NA'
                ion_variants['Read/M'] = 'NA'
            if ion_variants.empty:
                return pd.DataFrame({"Run": run_id, "Sample": sample, "Barcode": bar_code,
                                                   "Locus": 'negative', "Genes": 'NA', "Type": 'NA',
                                                   "Exon": 'NA', "Transcript": 'NA', "Coding": 'NA',
                                                   "Variant Effect": 'NA', 'Genotype': 'NA',
                                                   "% Frequency": "NA", "ExAC_AF": 'NA', "Amino Acid Change": 'NA',
                                                   "Read Counts": 'NA', "Read/M": 'NA', "Copy Number": 'NA',
                                                   "CNV Confidence": 'NA', "AMAF": 'NA', "GMAF": 'NA', "EMAF": 'NA',"HS":""},
                                                  index=[0])
            ion_variants.drop(columns=['ExAC_info', 'go', '5000Exomes', 'hrun', 'drugbank',
                                               'fusion_presence', 'ratio_to_wild_type', 'length',
                                               'norm_count_within_gene', 'iscn', 'filter',
                                               'allele_coverage', 'allele_ratio', 'pvalue',
                                               'precision', 'Subtype', 'dgv', 'cnv_pvalue',
                                               'allele_frequency_%', 'MyVariantDefaultDb_hg19',
                                               'phylop', 'pfam', 'location', 'maf',
                                               'sift', 'polyphen', 'grantham', 'normalizedAlt',
                                               'NamedVariants'], inplace=True)
            ion_variants.rename(columns={'confidence': 'CNV Confidence',
                                                 '# locus': 'Locus', 'gene': 'Genes', 'exon': 'Exon',
                                                 'transcript': 'Transcript', 'genotype': 'Genotype',
                                                 'type': 'Type', 'coding': 'Coding', 'function': 'Variant Effect',
                                                 'ref': 'Ref', 'tumor_AF': '% Frequency',
                                                 'protein': 'Amino Acid Change', 'coverage': 'Coverage',
                                                 'exac': 'ExAC', 'MAF':'ExAC_AF'},
                                        inplace=True)

            ion_variants['AMAF'] = new[0]
            ion_variants["GMAF"] = new[1]
            ion_variants["EMAF"] = new[2]
            ion_variants = ion_variants[
                ["Run", "Sample", "Barcode", "Locus", "Genes", "Type", "Exon", "Transcript", "Coding", "Variant Effect", "Genotype",
                 "% Frequency", "ExAC_AF", "Amino Acid Change", "Read Counts", "Read/M", "Copy Number", "CNV Confidence", "AMAF", "GMAF", "EMAF","HS"]]
            logger.info("%s processed"%sample)
            print("After filters applied")
            print(ion_variants.head(3).to_string())

            return ion_variants.drop_duplicates(subset=['Sample','Barcode','Locus','Genes','Type'])
            #return ion_variants
        else:
            logger.warning("%s tsv file not found"%sample)

    def run(self):
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
            logger.info("SC sample: %s" %sc_sample_name)
            for sample, barcode in zip(list(sample_sheet['sample_id']), list(sample_sheet['Bar code'])):
                try:
                    if barcode == "" or barcode == None: continue
                    RESULTS.append(self.process_sample([sample,run_id,barcode,logging.getLogger(sample)]))
                except:
                    continue
            if RESULTS and len(RESULTS) > 0:
                df = pd.concat(RESULTS)
                sc_df = df.loc[df['Sample'] == sc_sample_name]
                sc_df['Read Counts'] = sc_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read Counts'], axis=1)
                sc_df['Read/M'] = sc_df.apply(lambda x: 'NA' if (x['Type'] == "SNV" or x['Type'] == "INDEL") else x['Read/M'], axis=1)
                sc_df.to_csv("sc_filtered_variants.tsv", index = False, sep = "\t")
                sample_df = df.loc[df['Sample'] != sc_sample_name]
                logger.info(df.to_string())
                self.write_to_excel(sample_df)
                logger.info('Took %s seconds to process samples', time() - ts)
        except Exception as e:
             logger.error(str(e))
