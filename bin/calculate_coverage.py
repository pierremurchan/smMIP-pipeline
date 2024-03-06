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

def main(input_dir, output, sex_coverage_threshold):

    data = pd.DataFrame()
    for f in os.listdir(input_dir):
        if f.endswith("_UMI_usage_per_smMIP.txt"):
            sample = '_'.join(f.split("_")[:-4]) # everything but the '_UMI_usage_per_smMIP.txt' part
            df = pd.read_csv(os.path.join(input_dir, f), sep="\t")
            extracted_col  = df["Total_unique_UMI_pairs_(molecules)"]
            data.insert(len(data.columns), sample, extracted_col)         
    data = data.T
    data.columns = df["smMIP"]
    data = data.sort_index()
    data = data.loc[:, ~data.columns.str.contains('SampleID')] # drop columns contnaining 'SampleID' in the name

    tmp_df = pd.DataFrame()
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

    output_df.to_csv(output, index_label='CMD_ID')

parser = argparse.ArgumentParser(description="Create a heatmap of coverage per smMIP across samples.")
parser.add_argument("--input-dir", type=str, help="Input directory with `*_UMI_usage_per_smMIP.txt` files", required=True)
parser.add_argument("--sex-coverage-threshold", type=int, help="Coverage threshold to call sample sex from AMELX and AMELY", default=50)
parser.add_argument("--output", type=str, help="Output file.", default="coverage.csv")

if __name__ == "__main__":

    args = parser.parse_args()
    
    main(args.input_dir, args.output, args.sex_coverage_threshold)
