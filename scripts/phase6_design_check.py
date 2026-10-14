"""Phase 6 design check of O1 (threshold donors only): null SD vs a TCGA-direction shift.

For three block weightings (binomial, + co-methylation, WZ), each threshold donor is scored
leave-one-out. Prints the f at which a shift p0 + f * delta moves the unnormalised score by
2.5 null SDs. Reads no eval donor, cancer or spike record.

Usage: uv run python scripts/phase6_design_check.py
"""
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings('ignore')
spec=importlib.util.spec_from_file_location("s6","scripts/phase6_score.py"); s6=importlib.util.module_from_spec(spec); spec.loader.exec_module(s6)
p5s,p6=s6.p5s,s6.p6
data=p5s.data_dir()/"phase5"
smp=pd.read_csv('pipeline/samples/phase5_samples.tsv',sep='\t'); thr=smp.run[smp.role=='threshold'].tolist()
bl=p5s.load_blocks(data/"blocks_cpg.tsv"); mk=pd.read_csv('pipeline/markers/phase5_markers.tsv',sep='\t')
rho={'atlas':0.173575,'tumour':0.0821275}  # Phase 5 results/phase5/run_info.tsv
for s in ['atlas','tumour']:
    b=bl[bl.set==s].reset_index(drop=True); nb=len(b)
    T=mk[mk.set==s].groupby('block').tumour_mean.mean().reindex(b.block).to_numpy()
    S=[p6.summarise(p5s.load_cache(str(data/"records"/f"{d}.dedup.npz"),[s])[0][s],nb,rho[s]) for d in thr]
    out={}
    for name in ['binom','comethyl','wz']:
        sc=[];sig=[]
        for i in range(14):
            m=p6.fit(S[:i]+S[i+1:],T)
            M,N,G=[p6._get(S[i],f) for f in ['M', 'N', 'G']]
            ok=(N>0)&np.isfinite(m.p0); q=np.clip(m.p0,.01,.99)
            if name=='binom': v=q*(1-q)/np.maximum(N,1)
            elif name=='comethyl': v=q*(1-q)*G/np.maximum(N,1)**2
            else: v=q*(1-q)*(G/np.maximum(N,1)**2+m.phi)
            D=m.delta
            sc.append(np.where(ok,D*(M/np.maximum(N,1)-m.p0)/v,0).sum()); sig.append(np.where(ok,D**2/v,0).sum())
        sc=np.array(sc);sig=np.array(sig)
        # f needed for mean shift = 2.5 null SD (TCGA-direction signal)
        print(s,name,'null mean/sd (units of signal per f)',round(sc.mean()/sig.mean(),4),round(sc.std(ddof=1)/sig.mean(),4),' f for 2.5 SD',round(2.5*sc.std(ddof=1)/sig.mean(),4))
