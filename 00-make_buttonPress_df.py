import pandas as pd
import numpy as np
import os
from os.path import join, isdir, dirname, basename, exists
import itertools
from glob import glob
import argparse
# load the data from the .dat (json) file into a dictionary  
from json import load


eCON= '/data/bswift-1/Pessoa_Lab/eCON'
DATPAT = join(eCON,'onsetdir/{subj}/subj{subj}_run{run}.dat')
yoked = pd.read_excel(join(eCON,'onsetdir/CON_yoked_table.xlsx'))
yoked = yoked.query('use == 1')

def _get_button_presses():
    # inirialize and empty dataframe
    df = pd.DataFrame()
    # Loop thru every unctrol participant
    for i, row in yoked.iterrows():
        subj = row['uncontrol']
        # Filter out runs
        runs = np.arange(6)[row.loc['run0':'run5'].astype(bool)]
        # loop thru runs
        for j, run in enumerate(runs):
            # load run .dat file
            path = DATPAT.format(subj=subj, run=run)
            with open(path, 'r') as f:
                data = load(f)
            
            # number of button presses by control
            numEsc = len(list(itertools.chain(*data['Escapes'])))
            # number of button presses by the uncontrol yoke
            numNonEsc = len(data['nonEscapes'])
            tmp_df_uncon = pd.DataFrame(['P{:02d}'.format(i),'uncontrol',subj,j,numNonEsc],
                                        index=['Pair','Group','Subject','run','buttPress']).T
            tmp_df_con = pd.DataFrame(['P{:02d}'.format(i),'control',row['control'],j,numEsc],
                                      index=['Pair','Group','Subject','run','buttPress']).T
            df = pd.concat([df,tmp_df_con,tmp_df_uncon],axis=0,ignore_index=True)

    # Sum up button presses across runs
    summed_df = df.groupby(['Pair','Group','Subject'])['buttPress'].sum().reset_index()
    # make 'Pair' the index
    summed_df.set_index('Pair',inplace=True)
            
    return summed_df
            
if __name__ == '__main__':
    '''
    Extract total number of button-presses
    '''
    
    parser = argparse.ArgumentParser(description='Extract total number of button-presses')
    parser.add_argument('-o','--output',type=str,
                        help='output path to save the button press df',
                        default='data/behavioral')
    parser.add_argument('--overwrite',type=int,default=0,
                        help='overwrite the existing file? (0 or 1)')
    
    args = parser.parse_args()
    
    if exists(args.output+'/button_presses.txt') and not args.overwrite:
        raise OSError(args.output+"/button_presses.txt: File exists. Use '-overwrite 1' to overwrite")
    else:
        os.makedirs(args.output,exist_ok=True)
        df = _get_button_presses()
        df.to_csv(args.output+'/button_presses.txt',sep='\t')
    
    