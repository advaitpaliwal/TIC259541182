# Build fitdata.npz for TIC 259541182 from native-cadence FITS in hunt7/lc (same sectors/products TLS used: SPOC 4,5,31,32 PDCSAP; QLP 97,98 DET_FLUX).
# quality==0, wotan biweight (transit masked), keep |phase| < 0.12 d. Exposure from each file's TIMEDEL.
import numpy as np, json
from astropy.io import fits
from wotan import flatten, transit_mask
HERE='work/deep/259541182/fit'; LC='work/hunt7/lc/259541182'
R=json.load(open('work/hunt7/res/259541182.json')); P0,T00=R['P'],R['T0']
D0=json.load(open('work/hunt7/raw_gate_259541182.json'))['duration_days']
FILES=[('spoc_s4','PDCSAP_FLUX','PDCSAP_FLUX_ERR'),('spoc_s5','PDCSAP_FLUX','PDCSAP_FLUX_ERR'),('spoc_s31','PDCSAP_FLUX','PDCSAP_FLUX_ERR'),
       ('spoc_s32','PDCSAP_FLUX','PDCSAP_FLUX_ERR'),('qlp_s97','DET_FLUX','DET_FLUX_ERR'),('qlp_s98','DET_FLUX','DET_FLUX_ERR')]
T,F,E,S,X=[],[],[],[],[]; log={}
for nm,fc,ec in FILES:
    h=fits.open(f'{LC}/{nm}.fits'); d=h[1].data; hd=h[1].header; s=h[0].header['SECTOR']
    exp=int(round(hd['TIMEDEL']*86400))
    noerr=not np.isfinite(d[ec]).any()   # QLP s97 DET_FLUX_ERR is all NaN -> constant error = OOT MAD (jitter term free per exposure group)
    m=(d['QUALITY']==0)&np.isfinite(d[fc])&np.isfinite(d['TIME'])&(noerr|np.isfinite(d[ec]))
    t,f=d['TIME'][m].astype(float),d[fc][m].astype(float); e=np.ones_like(f) if noerr else d[ec][m].astype(float); e=e/np.nanmedian(f); f=f/np.nanmedian(f)
    msk=transit_mask(time=t,period=P0,duration=3*D0,T0=T00)
    tr=flatten(t[~msk],f[~msk],method='biweight',window_length=0.5,return_trend=True)[1]
    g=np.isfinite(tr); fl=np.interp(t,t[~msk][g],tr[g]); f,e=f/fl,e/fl
    r=f-1; sig=1.4826*np.median(np.abs(r[~msk]-np.median(r[~msk])))
    if noerr: e=np.full_like(f,sig)
    ph=((t-T00+0.5*P0)%P0)-0.5*P0; w=(r<5*sig)&(r>-10*sig)&(np.abs(ph)<0.12)
    T+=[t[w]];F+=[f[w]];E+=[e[w]];S+=[np.full(w.sum(),s)];X+=[np.full(w.sum(),exp)]
    log[nm]=dict(sector=s,exp_s=exp,timesys=hd.get('TIMESYS'),timeref=hd.get('TIMEREF'),kept=int(w.sum()),oot_mad_ppm=round(sig*1e6),med_err_ppm=round(float(np.median(e))*1e6),err_column=('none: set to OOT MAD' if noerr else ec),flux_column=fc)
    print(nm,log[nm],flush=True)
np.savez(HERE+'/fitdata.npz',t=np.concatenate(T),f=np.concatenate(F),e=np.concatenate(E),s=np.concatenate(S),x=np.concatenate(X))
json.dump(dict(P0=P0,T0_0=T00,dur0_d=D0,files=log),open(HERE+'/prep_log.json','w'),indent=1)
