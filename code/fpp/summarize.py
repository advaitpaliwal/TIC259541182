# mean +/- sample SD over loc_molusc runs (from deep/375939912/fpp/scripts/summarize.py).
import json,glob,pandas as pd
F='work/hunt7/fpp/259541182/trc/out/'
d=pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(F+'res_loc_molusc_*.json'))])
s=pd.concat([pd.read_csv(F+f'probs_loc_molusc_{k}.csv').groupby('scenario')['prob'].sum() for k in d.run],axis=1).fillna(0).mean(axis=1)
print(f"n={len(d)} FPP {d.FPP.mean():.5f} +/- {d.FPP.std(ddof=1):.5f} [{d.FPP.min():.5f},{d.FPP.max():.5f}] NFPP {d.NFPP.mean():.3e} +/- {d.NFPP.std(ddof=1):.3e} n_stars {d.n_stars[0]} depth {d.depth_ppm[0]:.1f}ppm err {d.err_ppm[0]:.1f} nbins {d.nbins[0]} dropped {d.dropped[0]}")
print(' '.join(f'{k}={v:.5f}' for k,v in s.items() if v>0))
