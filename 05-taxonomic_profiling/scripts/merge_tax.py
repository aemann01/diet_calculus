import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--taxonomy')
parser.add_argument('-i', '--input')
parser.add_argument('-o', '--output')
args = parser.parse_args()

taxonomy = pd.read_csv(args.taxonomy, sep="\t")
k1 = pd.read_csv(args.input, sep="\t")
merged = pd.DataFrame.merge(taxonomy, k1, how='right', left_on="taxid", right_on="taxid")
merged.to_csv(args.output, sep="\t", index=False)