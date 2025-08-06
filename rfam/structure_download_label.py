import sys
import os

project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if project_root not in sys.path:
    sys.path.append(project_root)

from utils import pdb_utils

import pandas as pd

# TODO: Rename the columns
# Label
rfam_pdb_file = pd.read_csv('rfam/Rfam.pdb',sep='\t')
rfam_label = pd.read_csv('rfam/rfam_seed_header.csv',sep='\t')

merged_df = pd.merge(rfam_pdb_file,rfam_label[['AC','TP','WK','CL']],how='left', left_on='rfam_acc', right_on='AC')
merged_df.dropna()

merged_df.to_csv('rfam/rfam_pdb_labeled.csv',sep='\t',index=False)

# TODO: process the structures