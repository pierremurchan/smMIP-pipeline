import pandas as pd
import plotly.express as px
import numpy as np
import argparse
import os

def discrete_colorscale(bvals, colors):
    """
    Genreates a discrete colorscale for plotly
    -----

        bvals - list of values bounding intervals/ranges of interest
        colors - list of rgb or hex colorcodes for values in [bvals[k], bvals[k+1]],0<=k < len(bvals)-1
        returns the plotly  discrete colorscale
    """
    if len(bvals) != len(colors)+1:
        raise ValueError('len(boundary values) should be equal to  len(colors)+1')
    bvals = sorted(bvals)
    nvals = [(v-bvals[0])/(bvals[-1]-bvals[0]) for v in bvals]  # normalised values
    
    dcolorscale = []
    for k in range(len(colors)):
        dcolorscale.extend([[nvals[k], colors[k]], [nvals[k+1], colors[k]]])
    return dcolorscale

def main(input_dir, output):

    data = pd.DataFrame()
    for f in os.listdir(input_dir):
        if f.endswith("_UMI_usage_per_smMIP.txt"):
            sample = '_'.join(f.split("_")[:-4])
            df = pd.read_csv(os.path.join(input_dir, f), sep="\t")
            extracted_col  = df["Total_unique_UMI_pairs_(molecules)"]
            data.insert(len(data.columns), sample, extracted_col)            
    data = data.T


    bvals = [0, 250, 300, 100000] # defining thresholds for flagging
    colors = ['#d9534f', '#fee391', '#5cb85c']

    dcolorsc = discrete_colorscale(bvals, colors)

    fig = px.imshow(
        data, 
        labels=dict(x="Sample", y="smMIP", color="Coverage"),
        x=data.columns,
        y=data.index,
        aspect="auto",
        color_continuous_scale=dcolorsc
    )

    fig.update_layout(
        title="Heatmap of Coverage per smMIP across Samples",
        xaxis_title="smMIP",
        yaxis_title="Sample",
        xaxis={'side': 'bottom'},
        yaxis_nticks=len(data.index),
    )

    fig.update_traces(hovertemplate='Sample: %{x}<br>smMIP: %{y}<br>Coverage: %{z}<extra></extra>')
    fig.write_html(output)

parser = argparse.ArgumentParser(description="Create a heatmap of coverage per smMIP across samples.")
parser.add_argument("--input-dir", type=str, help="Input directory with `*_UMI_usage_per_smMIP.txt` files", required=True)
parser.add_argument("--output", type=str, help="Output file.", default="coverage_heatmap.html")

if __name__ == "__main__":

    args = parser.parse_args()
    
    main(args.input_dir, args.output)