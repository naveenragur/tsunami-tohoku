"""
To filter  earhtquale events from the dtopo list to remove unrealistic high magnitude events

"""
from curses.ascii import NUL
import os
import numpy as np
import pandas as pd

EQdtopo=pd.read_csv("gis/dtopo_TMOD/tohoku_dtopo_list.csv")

TEST = None
for row in len(EQdtopo):
    if EQdtopo.loc[row,"EveMw"] == 9 and EQdtopo.loc[row,"EveDepth"] > 11 or \
    EQdtopo.loc[row,"EveMw"] == 8.5 and EQdtopo.loc[row,"EveDepth"] > 20 or \
        EQdtopo.loc[row,"EveMw"] == 8 and EQdtopo.loc[row,"EveDepth"] >  30:
        TEST=True
else:
    TEST=False
    
        EQdtopo.drop(row,inplace=True)

