from os.path import join, exists
import os
import numpy as np
import pandas as pd
import argparse

eCON = '/data/bswift-1/Pessoa_Lab/eCON'

def _get_STAI():
    # load full STAI scores and re-organize them
    scores = pd.read_excel(join(eCON,'STAIscores/scores.xlsx'))
    scores.rename(columns={'SCORING':'TRAIT'},inplace=True)
    scores.drop(columns=['PAIR SCORE'],inplace=True)
    scores['Subject'] = scores['SubID'].apply(lambda name: ''.join(name.strip().split("_")))
    scores.drop('SubID',axis=1,inplace=True)
    return scores


if __name__ == '__main__':
    '''
    This script gets every participants state and trait scores
    '''
    
    parser = argparse.ArgumentParser(description='gets every participants state and trait scores')
    parser.add_argument('-o','--output',type=str,
                        help='output path to save the button press df',
                        default='data/behavioral')
    parser.add_argument('--overwrite',type=int,default=0,
                        help='overwrite the existing file? (0 or 1)')
    
    args = parser.parse_args()
    
    if exists(args.output+'/STAI_scores.txt') and not args.overwrite:
        raise OSError(args.output+"/STAI_scores.txt: File exists. Use '-overwrite 1' to overwrite")
    else:
        os.makedirs(args.output,exist_ok=True)
        df = _get_STAI()
        df.to_csv(args.output+'/STAI_scores.txt',sep='\t',index=False)