import pandas as pd

strike_rang = [*range(-180, 180, 30)]
dip_rang = [*range(0,90,30)]
rake_rang = [*range(-100,100,30)]

strike_ls=[];dip_ls=[];rake_ls=[]
for st in strike_rang:
    for dp in dip_rang:
        for rk in rake_rang:
            strike_ls.append(st); dip_ls.append(dp); rake_ls.append(rk)

data = {
        'Strike': strike_ls,
        'Dip': dip_ls,
        'Rake': rake_ls,
        }

df = pd.DataFrame(data, columns = ['Strike', 'Dip', 'Rake'])
df.to_csv('trialdataframe_large.csv')
