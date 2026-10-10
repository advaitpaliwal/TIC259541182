# Build TRICERATOPS target for TIC 259541182 (glue from hunt5b/trc/run_fpp.py target block) + TIC<->Gaia exclusion map (from deep/375939912/fpp/scripts/xmatch.py).
import sys, pickle, os, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
O='work/hunt7/fpp/259541182/'; H='work/hunt7/'
sys.path.insert(0,O+'scripts'); import gaia_patch
import triceratops.triceratops as tr
tic=259541182; secs=[4,5,31,32,97,98]; pk=O+'trc/out/target_259541182.pkl'
if not os.path.exists(pk):
    tg=tr.target(ID=tic,sectors=np.array(secs)); aps=[]
    for i,s in enumerate(secs):
        fn=H+f'lc/{tic}/spoc_s{s}.fits'; ap=None
        if os.path.exists(fn):
            with fits.open(fn) as h:
                a=h['APERTURE'].data.astype(int); c0=h['APERTURE'].header['CRVAL1P']; r0=h['APERTURE'].header['CRVAL2P']
            rows,cols=np.nonzero(np.bitwise_and(a,2))
            if len(rows): ap=np.column_stack([cols+c0,rows+r0])
        if ap is None:
            tp=np.round(tg.pix_coords[i][0]); ap=np.array([np.repeat(np.arange(tp[0]-2,tp[0]+3),5),np.tile(np.arange(tp[1]-2,tp[1]+3),5)]).T
        aps.append(ap); print('sector',s,'ap pix',len(ap),'spoc' if os.path.exists(fn) else 'qlp5x5')
    tg._aps=aps; pickle.dump(tg,open(pk,'wb'))
tg=pickle.load(open(pk,'rb')); st=tg.stars
print(st.iloc[0][['ID','Tmag','Teff','rad','mass']].to_dict())
L=H+'localize/259541182/'
g=pd.read_csv(L+'gaia_field.csv',dtype={'Source':str,'source':str}); g=g.rename(columns={'source':'Source'}) if 'Source' not in g else g
d=pd.read_csv(L+'neighbour_sums.csv',dtype={'source':str}); g=g.merge(d[['source','delta_chi2_sum']],left_on='Source',right_on='source',how='left')
rac=[c for c in ('RA_ICRS','ra') if c in g][0]; dec=[c for c in ('DE_ICRS','dec') if c in g][0]
pmr=g['pmRA'] if 'pmRA' in g else g.get('pmra',0); pmd=g['pmDE'] if 'pmDE' in g else g.get('pmdec',0)
d0=np.deg2rad(g[dec]); ra=g[rac]-16*pd.Series(pmr).fillna(0)/3.6e6/np.cos(d0); de=g[dec]-16*pd.Series(pmd).fillna(0)/3.6e6
rows=[]
for _,s in st.iterrows():
    sep=np.hypot((ra-s.ra)*np.cos(np.deg2rad(s.dec)),de-s.dec)*3600; j=int(np.argmin(sep.values)); m=sep.iloc[j]<1.5; r=g.iloc[j]
    rows.append(dict(tic=str(s.ID),Tmag=s.Tmag,trc_sep=s['sep (arcsec)'],gaia=r.Source if m else '',match_sep=round(float(sep.iloc[j]),3),delta_chi2_sum=r.delta_chi2_sum if m else np.nan))
o=pd.DataFrame(rows); o['excluded']=(o.delta_chi2_sum>=25)&(o.tic!=str(tic)); o.to_csv(O+'exclusions.csv',index=False)
print(o[o.gaia!=''].to_string()); print('excluded',int(o.excluded.sum()),'of',len(o)-1,'; unmatched localized:',sorted(set(d.source)-set(o.gaia)))
