import pandas as pd, numpy as np
p='/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv'
df=pd.read_csv(p); df['cd_min']=df[['cd_ref','cd_min_alt']].min(axis=1)
AD='autodock_mgltools_exh128_gnina'; DD='diffdock_smina'; EB='equibind_unguided_gnina'; arms=[AD,DD,EB]
prots=sorted(set(df[df.method==DD].protein)&set(df[df.method==AD].protein)&set(df[df.method==EB].protein))
for col,cd,lab in [('rmsd','cd_ref','REF'),('rmsd_nearest_copy','cd_min','NEAREST')]:
    s=df[(df.method.isin(arms))&(df.eff_rank<=15)&(df.protein.isin(prots))&(df.pb_valid)&(df[cd]<=8)]
    print(lab,'near-site pool',len(s),'clipped (inplace>5 or bestfit>5):',int(((s[col]>5)|(s.bestfit_rmsd>5)).sum()),'clipped (inplace>5 only):',int((s[col]>5).sum()))
# cost arithmetic from printed Table 8 cells
print('CPU total from AD row',1752.8*311,'from Vina row',1805.0*302,'GPU total',118.9*311)
print('AD CPU/442',545112/442,'Vina CPU/431',545112/431,'GPU/442',36972/442,'GPU/success 214',36972/214)
print('DD CPU',2.93*1622/2068,'GPU',49.9*1622/2068,'EB CPU',165.3*466/483,'GPU',1.4*466/483)
# :197 footnote basis: 33.7% under DiffDock confidence ranking; candidates
s=df[(df.method==DD)&(df.eff_rank==1)]
for n in [303,300,297]:
    print('33.7% of',n,'=',0.337*n)
print('DD rank1 NN ref',int((s.rmsd<=2).sum()),'NN&valid',int(((s.rmsd<=2)&s.pb_valid).sum()),'| nearest NN',int((s.rmsd_nearest_copy<=2).sum()),'NN&valid',int(((s.rmsd_nearest_copy<=2)&s.pb_valid).sum()))
# Kurzfassung character counts
for a,b in [('36,6\\% für AutoDock, 34,0\\% für DiffDock und 18,2\\%','49,2\\% für AutoDock, 43,6\\% für DiffDock und 18,8\\%'),('\\ensuremath{-}4,2 bis +9,4','\\ensuremath{-}1,7 bis +12,8'),('etwa sechs von zehn','etwa sechs bis sieben von zehn')]:
    print(len(a),len(b),len(b)-len(a))
