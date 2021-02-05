import numpy as np
import pandas as pd

faults = pd.read_csv('faults_173a.csv')
print(faults.iloc[:10])

faults2 = faults.drop_duplicates(subset = ["Strike","Dip","Rake"])
print(faults2.iloc[:10])
