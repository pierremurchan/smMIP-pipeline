import pandas as pd
import argparse
from jinja2 import Environment, FileSystemLoader
#from bs4 import BeautifulSoup

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process variant calling results and generate an HTML report.")
    parser.add_argument('--called-mutations', type=str, help='Path to the input TXT file containing variant calling results.')
    parser.add_argument('--template-path', type=str, help='Path to the directory containing the HTML template.', default='./')
    parser.add_argument('--coverage-heatmap', type=str, help='Path to the HTML file containing coverage heatmap.', default='coverage_heatmap.html')
    #parser.add_argument('--multiqc-report', type=str, help='Path to the HTML file containing MultiQC report.', default='multiqc_report.html')
    #parser.add_argument('--nextflow-execution-report', type=str, help='Path to the HTML file containing Nextflow execution report.', default='pipeline_info/execution_report_2024-02-21_09-47-31.html')
    #parser.add_argument('--nextflow-timeline-report', type=str, help='Path to the HTML file containing Nextflow timeline report.', default='pipeline_info/execution_timeline_2024-02-21_09-47-31.html')
    parser.add_argument('--configuration-csv', type=str, help='Path to the CSV file containing configuration data.', default='configuration.csv')
    parser.add_argument('--output-file', type=str, help='Path to the output HTML report.', default='dev_variant_report.html')
    return parser.parse_args()

def read_variant_data(input_file):
    return pd.read_csv(input_file, sep='\t')

def read_configuration_data(configuration_csv):
    config_df = pd.read_csv(configuration_csv, sep=',')
    return config_df[config_df['type'] != 'control']

def read_and_sanitize_html(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()

    # sanitize HTML content
    #soup = BeautifulSoup(html_content, 'html.parser')
    #sanitized_content = soup.prettify()  # Adjust this line as needed for your sanitization logic
    
    #return sanitized_content
    return html_content

def check_execution_report(nextflow_execution_report):
    try:
        with open(nextflow_execution_report) as f:
            pass
    except FileNotFoundError:
        print(f"Execution report not found at {nextflow_execution_report}. Please check the path and try again.")
        exit(1)

def prepare_data_for_report(var_call_df, configuration_df):

    # initialise a dictionary to hold processed data
    data_for_report = {
        'summary_stats': {},
        'samples_with_mutations': [],
        'samples_without_mutations': configuration_df['id'].tolist(),
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
        # extract protein change from the protein field by regex for ':' and 'p.' 
        sample_details['variants'] = [{**variant, 'protein_change': variant['protein'].split(':')[-1] if pd.notnull(variant['protein']) else 'NA'} for variant in sample_details['variants']]
        data_for_report['samples_with_mutations'].append(sample_details)
        if sample_id in data_for_report['samples_without_mutations']:
            data_for_report['samples_without_mutations'].remove(sample_id)

    # QC metrics ?
    return data_for_report

def report_mutations_in_genes(data_for_report):
    # initialize a dictionary to hold the information about the mutations in each gene
    data_for_report['mutations_in_genes'] = {}
    CALR_alle_frequnecy_threshold = 0.01
    JAK2_MPL_alle_frequnecy_threshold = 0.01

    # check if there are samples with mutations
    if len(data_for_report['samples_with_mutations']) > 0:
        # iterate over the samples with mutations
        for sample in data_for_report['samples_with_mutations']:
            # iterate over the mutations in each sample
            for mutation in sample['variants']:
                # check if the mutation is in CALR
                if mutation['gene'] == 'CALR' and mutation['SSCS_allele_frequency'] > CALR_alle_frequnecy_threshold:
                    # if the gene is not in the dictionary, add it with an empty list
                    if 'CALR' not in data_for_report['mutations_in_genes']:
                        data_for_report['mutations_in_genes']['CALR'] = []
                    # add the sample to the list of samples with mutations in CALR
                    data_for_report['mutations_in_genes']['CALR'].append(sample['sample_id'])
                # check if the mutation is in JAK2
                elif mutation['gene'] == 'JAK2' and mutation['SSCS_allele_frequency'] > JAK2_MPL_alle_frequnecy_threshold:
                    # if the gene is not in the dictionary, add it with an empty list
                    if 'JAK2' not in data_for_report['mutations_in_genes']:
                        data_for_report['mutations_in_genes']['JAK2'] = []
                    # add the sample to the list of samples with mutations in JAK2
                    data_for_report['mutations_in_genes']['JAK2'].append(sample['sample_id'])
                # check if the mutation is in MPL
                elif mutation['gene'] == 'MPL' and mutation['SSCS_allele_frequency'] > JAK2_MPL_alle_frequnecy_threshold:
                    # if the gene is not in the dictionary, add it with an empty list
                    if 'MPL' not in data_for_report['mutations_in_genes']:
                        data_for_report['mutations_in_genes']['MPL'] = []
                    # add the sample to the list of samples with mutations in MPL
                    data_for_report['mutations_in_genes']['MPL'].append(sample['sample_id'])
        # check if there are samples with mutations in CALR
        if 'CALR' not in data_for_report['mutations_in_genes']:
            # if there are no samples with mutations in CALR, add the gene
            data_for_report['mutations_in_genes']['CALR'] = 'No mutation detected in gene CALR'
        # check if there are samples with mutations in JAK2
        if 'JAK2' not in data_for_report['mutations_in_genes']:
            # if there are no samples with mutations in JAK2, add the gene
            data_for_report['mutations_in_genes']['JAK2'] = 'No mutation detected in gene JAK2'
        # check if there are samples with mutations in MPL
        if 'MPL' not in data_for_report['mutations_in_genes']:
            # if there are no samples with mutations in MPL, add the gene
            data_for_report['mutations_in_genes']['MPL'] = 'No mutation detected in gene MPL'
    else:
        # if there are no samples with mutations, add the genes with the value 'No mutation detected in gene X'
        data_for_report['mutations_in_genes']['CALR'] = 'No mutation detected in gene CALR'
        data_for_report['mutations_in_genes']['JAK2'] = 'No mutation detected in gene JAK2'
        data_for_report['mutations_in_genes']['MPL'] = 'No mutation detected in gene MPL'
    return data_for_report

def generate_html_report(data, template_path, output_file):
    env = Environment(loader=FileSystemLoader(template_path))
    template = env.get_template('variant_report_template.html')
    html_output = template.render(data=data)
    with open(output_file, 'w') as f:
        f.write(html_output)

def main():
    args = parse_arguments()
    #check_execution_report(args.nextflow_execution_report)
    var_call_df = read_variant_data(args.called_mutations)
    configuration_df = read_configuration_data(args.configuration_csv)
    data_for_report = prepare_data_for_report(var_call_df, configuration_df)
    data_for_report = report_mutations_in_genes(data_for_report)
    data_for_report['coverage_heatmap'] = args.coverage_heatmap
    #data_for_report['multiqc_report'] = read_and_sanitize_html(args.multiqc_report)
    #data_for_report['nextflow_execution_report'] = read_and_sanitize_html(args.nextflow_execution_report)
    #data_for_report['nextflow_timeline_report'] = read_and_sanitize_html(args.nextflow_timeline_report)
    generate_html_report(data_for_report, args.template_path, args.output_file)

if __name__ == "__main__":
    main()

