import pandas as pd
import numpy as np
import argparse
import os

def classify_sex(row, sex_coverage_threshold):
    row = row.astype(float)
    if row['AMELY'] < sex_coverage_threshold and row['AMELX'] < sex_coverage_threshold:
        return 'INSUFFICIENT_COVERAGE'
    elif row['AMELY'] > 0:
        return 'Male' if row['AMELX'] / row['AMELY'] < 2 else 'Female'
    else:
        return 'Female'
    
def load_and_merge_vars(high_conf_snvs, high_conf_indels):
    snvs = pd.read_csv(high_conf_snvs, sep = '\t')
    indels = pd.read_csv(high_conf_indels, sep = '\t')
    var_call_df = pd.concat([snvs, indels])

    return var_call_df
    
def prepare_data_for_report(var_call_df, configuration_df):

    # initialise a dictionary to hold processed data
    data_for_report = {
        'summary_stats': {},
        'samples_with_mutations': [],
        'samples_without_mutations': configuration_df['id'].tolist(),  # all samples are initially assumed to have no mutations
        'qc_metrics': {}
    }

    # summary statistics
    data_for_report['summary_stats']['total_samples'] = len(data_for_report['samples_without_mutations'])
    data_for_report['summary_stats']['total_samples'] = configuration_df['id'].nunique()
    data_for_report['summary_stats']['total_variants'] = len(var_call_df)
    data_for_report['summary_stats']['avg_variants_per_sample'] = var_call_df.groupby('sample_ID').size().mean()

    # detailed data for each sample
    for sample_id, sample_data in var_call_df.groupby('sample_ID'):
        # change '.' in sample_data keys to '_'
        sample_data.columns = sample_data.columns.str.replace('.', '_')
        sample_details = {
            'sample_id': sample_id,
            'total_variants': len(sample_data),
            'variants': sample_data[['chr', 'pos', 'ref', 'alt', 'gene', 'protein', 'variant_type', 'cadd_scaled', 'allele_frequency', 'SSCS_non_ref_counts', 'SSCS_total_depth', 'SSCS_allele_frequency']].to_dict(orient='records')
        }
        sample_details['variants'] = [{**variant, 'genomic_change': f"{variant['chr'].replace('chr', '')}:g.{variant['pos']}{variant['ref']}>{variant['alt']}"} for variant in sample_details['variants']]
        # extract the coding change from the protein field by regex for ':' and 'c.'
        sample_details['variants'] = [{**variant, 'coding_change': variant['protein'].split(':')[-2] if pd.notnull(variant['protein']) else 'NA'} for variant in sample_details['variants']]
        # extract protein change from the protein field by regex for ':' and 'p.' 
        sample_details['variants'] = [{**variant, 'protein_change': variant['protein'].split(':')[-1] if pd.notnull(variant['protein']) else 'NA'} for variant in sample_details['variants']]
        data_for_report['samples_with_mutations'].append(sample_details)
        # remove samples with mutations from samples_without_mutations
        if sample_id in data_for_report['samples_without_mutations']:
            data_for_report['samples_without_mutations'].remove(sample_id)

    return data_for_report

def filter_variants(var_call_df):
    var_call_df = var_call_df[var_call_df['SSCS.allele.frequency'] >= 0.01]
    var_call_df = var_call_df[var_call_df['SSCS.total.depth'] >= 300]
    print(var_call_df)
    return var_call_df

def format_g_dot(variant_dict):
    return f"{variant_dict['gene']}:{variant_dict['chr']}:g.{variant_dict['pos']}{variant_dict['ref']}>{variant_dict['alt']}"

def main(cleaned_bams_dir, high_conf_snvs, high_conf_indels, configuration_csv, sex_coverage_threshold, output_csv):

    data = pd.DataFrame()
    for f in os.listdir(cleaned_bams_dir):
        if f.endswith("_UMI_usage_per_smMIP.txt"):
            sample = '_'.join(f.split("_")[:-4]) # everything but the '_UMI_usage_per_smMIP.txt' part
            df = pd.read_csv(os.path.join(cleaned_bams_dir, f), sep="\t")
            extracted_col  = df["Total_unique_UMI_pairs_(molecules)"]
            data.insert(len(data.columns), sample, extracted_col)         
    data = data.T
    data.columns = df["smMIP"]
    data = data.sort_index()
    data = data.loc[:, ~data.columns.str.contains('SampleID')] # drop columns contnaining 'SampleID' in the name

    tmp_df = pd.DataFrame() # a temporary dataframe to hold the sum of coverage for each gene/smMIP region
    for col in data.columns:
        region = col.split("_")[0] # gene/smMIP region
        if region not in tmp_df.columns:
            tmp_df[region] = data[col]
        else:
            tmp_df[region] += data[col]
    
    # classify the sex of samples based on AMELX and AMELY coverage            
    sex = tmp_df.apply(lambda row: classify_sex(row, sex_coverage_threshold), axis=1)

    # apply coverage thresholds to the smMIPs/gene regions
    output_df = pd.DataFrame()
    # for col in tmp_df.columns: except AMELX and AMELY
    for col in tmp_df.columns.difference(['AMELX', 'AMELY']):
        output_df[col] = np.where(tmp_df[col] > 300, 'PASS', np.where(tmp_df[col] < 250, 'FAIL', 'WARN'))
    output_df.index = tmp_df.index
    # add AMELX and AMELY columns with values from tmp_df
    output_df['AMELX'] = tmp_df['AMELX']
    output_df['AMELY'] = tmp_df['AMELY']

    output_df['sex'] = sex

    configuration_df = pd.read_csv(configuration_csv, header=0)

    var_call_df = load_and_merge_vars(high_conf_snvs, high_conf_indels)
    var_call_df = filter_variants(var_call_df)
    data_for_report = prepare_data_for_report(var_call_df, configuration_df)

    unique_genes = list(var_call_df['gene'].unique())
    unique_genes = [gene for gene in unique_genes if pd.notnull(gene)] # remove nan from unique_genes


    for sample in output_df.index.tolist():
        # if sample is in data_for_report['samples_without_mutations'], then the GENE_VAR columns shuould be 'NMD' for all genes
        if sample in data_for_report['samples_without_mutations']:
            for gene in unique_genes:
                output_df.at[sample, gene + '_VAR'] = 'NMD'
        else:
            sample_data = [sample_data for sample_data in data_for_report['samples_with_mutations'] if sample_data['sample_id'] == sample][0]
            for gene in unique_genes:
                gene_variants = [variant for variant in sample_data['variants'] if variant['gene'] == gene]

                if len(gene_variants) > 0:
                    for mut in gene_variants:
                        print(mut)
                        # output semi-colon separated list of variants as the 'genomic_change' found in the variant
                        output_df.at[sample, gene + '_VAR'] = '; '.join([f"{format_g_dot(mut)} @ {round(mut['SSCS_allele_frequency'] * 100, 1)}%" if pd.isnull(mut['protein']) else f"{mut['protein'].split(',')[-1]} @ {round(mut['SSCS_allele_frequency'] * 100, 1)}%"])
                else:
                    output_df.at[sample, gene + '_VAR'] = 'NMD'
    
    print(output_df)

    # save var_call_df to csv
    var_call_df.to_csv('var_call_df.csv', index=False)
    output_df.to_csv(output_csv, index_label='CMD_ID')

parser = argparse.ArgumentParser(description="Create a heatmap of coverage per smMIP across samples.")
parser.add_argument("--cleaned-bams-dir", type=str, help="`cleaned_bams` directory that contains the `*_UMI_usage_per_smMIP.txt` files", required=True)
parser.add_argument("--high-conf-snvs", type=str, help="High confidence SNVs file", required=True)
parser.add_argument("--high-conf-indels", type=str, help="High confidence indels file", required=True)
parser.add_argument('--configuration-csv', type=str, help='Path to the CSV file containing configuration data.', default='configuration.csv')
parser.add_argument("--sex-coverage-threshold", type=int, help="Coverage threshold to call sample sex from AMELX and AMELY", default=50)
parser.add_argument("--output-csv", type=str, help="Output file.", default="OUTPUT_SUMMARY.csv")

if __name__ == "__main__":

    args = parser.parse_args()
    
    main(args.cleaned_bams_dir, args.high_conf_snvs, args.high_conf_indels, args.configuration_csv, args.sex_coverage_threshold, args.output_csv)