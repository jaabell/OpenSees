#!/usr/bin/env python3
# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2026, Jose Antonio Abell. All rights reserved.
"""Generate and run a matched explicit legacy/combined-element comparison.

Uses the SAME executable, mesh, four-rank partition, dt, damping and output.
All decks/results are written to the required external --output directory.
Example:
  python compare_split_elastic.py --binary /path/to/OpenSeesMP --ranks 4 \
      --output /tmp/split-comparison --steps 4000 --repeats 2
Requires NumPy. This is a validation/throughput example, not a PML benchmark.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np


def generate(output, custom, ranks, steps, dx, warmup):
    nx,ny=round(4/dx),round(2/dx)
    assert abs(nx*dx-4)<1.e-10 and abs(ny*dx-2)<1.e-10
    coords={1+j*(nx+1)+i:(-2+i*dx,-2+j*dx) for j in range(ny+1) for i in range(nx+1)}
    tol=dx/100
    shared={n:abs(x)<=1+tol and y>=-1-tol for n,(x,y) in coords.items()}
    duplicate=max(coords)+1
    cells=[]; stars={n:set() for n in coords}
    volume={n:0. for n in coords}
    def material(y): return 1 if y>=-.5 else 2 if y>=-1.2 else 3
    def depth(x,y): return max(abs(x)-1,-1-y,0.)
    factor={n:max(0.,1-depth(x,y)) for n,(x,y) in coords.items()}
    for j in range(ny):
        for i in range(nx):
            nodes=[1+j*(nx+1)+i,2+j*(nx+1)+i,2+(j+1)*(nx+1)+i,1+(j+1)*(nx+1)+i]
            x,y=-2+(i+.5)*dx,-2+(j+.5)*dx
            exterior=abs(x)>1+tol or y<-1-tol
            part=min(ranks-1,i*ranks//nx)
            tag=1+j*nx+i
            # Smooth element-layer coordinate starts at zero on first layer.
            r=max(0.,min(1.,(depth(x,y)-dx/2)/(1-dx)))
            cells.append((tag,nodes,exterior,part,material(y),3*r*r-2*r*r*r))
            for n in nodes:
                stars[n].add(part)
                if exterior:volume[n]+=dx*dx/4
    controls=[]
    for wanted in [(0.,-.5),(-.5,-.5),(.5,-.5),(0.,-.9)]:
        n=min(coords,key=lambda n:sum((a-b)**2 for a,b in zip(coords[n],wanted)))
        if n not in controls:controls.append(n)
    out=[]
    out+=['wipe','set rank [getPID]',f'if {{[getNP] != {ranks}}} {{error "rank count mismatch"}}']
    out+=['model BasicBuilder -ndm 2 -ndf 2']
    E=[600.,450.,300.]; rho=2.; nu=.2
    for m,e in enumerate(E,1):
        out.append(f'nDMaterial ElasticIsotropic {m} {e:.17g} {nu} {rho}')
        out.append(f'nDMaterial ElasticIsotropic {10+m} {e/2:.17g} {nu} {rho/2}')
    counts=[]
    for p in range(ranks):
        local=[n for n in coords if p in stars[n]]
        quads=0;links=0
        out.append(f'if {{$rank == {p}}} {{')
        out.append('model BasicBuilder -ndm 2 -ndf 2')
        for n in local:
            if not custom or shared[n]:
                x,y=coords[n];out.append(f'node {n} {x:.17g} {y:.17g}')
                if not custom and not shared[n]:out.append(f'node {n+duplicate} {x:.17g} {y:.17g}')
        if custom:
            out.append('model BasicBuilder -ndm 2 -ndf 4')
            for n in local:
                if not shared[n]:
                    x,y=coords[n];out.append(f'node {n} {x:.17g} {y:.17g}')
        # Ordinary quad factories still inspect the current builder NDF.
        out.append('model BasicBuilder -ndm 2 -ndf 2')
        for n in local:
            x,y=coords[n]
            outer=abs(x)>=2-tol or y<=-2+tol
            if outer:
                out.append(f'fix {n} '+('1 1 0 0' if custom else '1 1'))
        for tag,nodes,ext,part,m,g in cells:
            if part!=p:continue
            conn=' '.join(map(str,nodes))
            if ext and custom:
                factors=' '.join(f'{factor[n]:.17g}' for n in nodes)
                out.append(f'element SplitElastic {tag} {conn} {E[m-1]} {nu} {rho} -coupling 10 -couplingFactors {factors}')
                out.append(f'region {tag} -ele {tag} -rayleigh {g*.4:.17g} {g*.0002:.17g} 0 0')
                quads+=1
            else:
                out.append(f'element quad {tag} {conn} 1 PlaneStrain {10+m if ext else m}')
                quads+=1
                if ext:
                    other=tag+nx*ny; second=' '.join(str(n if shared[n] else n+duplicate) for n in nodes)
                    out.append(f'element quad {other} {second} 1 PlaneStrain {10+m}')
                    out.append(f'region {tag} -ele {tag} {other} -rayleigh {g*.4:.17g} {g*.0002:.17g} 0 0')
                    quads+=1
        if not custom:
            for n in local:
                if shared[n] or min(stars[n])!=p or factor[n]==0:continue
                c=10*factor[n]*volume[n]
                out.append(f'uniaxialMaterial Viscous {n} {c:.17g} 1')
                out.append(f'element zeroLength {n+2*nx*ny} {n} {n+duplicate} -mat {n} {n} -dir 1 2')
                links+=1
        counts.append({'rank':p,'solids':quads,'connectors':links})
        out.append('}')
    dt=.0005
    values=[]
    for k in range(warmup+steps+2):
        a=math.pi*8*(k*dt-.12);values.append(.01*(1-2*a*a)*math.exp(-a*a))
    out.append('timeSeries Path 1 -dt '+str(dt)+' -values '+' '.join(f'{v:.17g}' for v in values))
    out.append('pattern Plain 1 1 {')
    for n,(x,y) in coords.items():
        if shared[n] and abs(y)<tol and abs(x)<=.5+tol:
            amp=math.cos(math.pi*x)
            parts=' || '.join(f'$rank == {p}' for p in sorted(stars[n]))
            out.append(f' if {{{parts}}} {{sp {n} 2 {amp:.17g}}}')
    out.append('}')
    for p in range(ranks):
        nodes=[n for n in controls if min(stars[n])==p]
        if nodes:
            out.append(f'if {{$rank == {p}}} {{recorder Node -file controls-{p}.txt -precision 16 -time -dT .005 -node '+ ' '.join(map(str,nodes))+' -dof 1 2 disp}')
    out+=['constraints Transformation','numberer ParallelPlain' if ranks>1 else 'numberer Plain',
          'system MPIDiagonal' if ranks>1 else 'system Diagonal','algorithm Linear','integrator ExplicitDifference','analysis Transient',
          f'if {{[analyze {warmup} {dt}] != 0}} {{error "warmup failed"}}','barrier' if ranks>1 else '# serial',
          'set start [clock microseconds]',f'if {{[analyze {steps} {dt}] != 0}} {{error "analysis failed"}}',
          'set elapsed [expr {([clock microseconds]-$start)/1.e6}]',
          'set status [open /proc/self/status r]; set rss 0',
          'foreach line [split [read $status] "\\n"] {if {[regexp {^VmHWM:\\s+(\\d+)} $line -> value]} {set rss $value}}',
          'close $status','puts "METRIC $rank $elapsed $rss"',
          'set f [open "state-$rank.txt" w]']
    for p in range(ranks):
        out.append(f'if {{$rank == {p}}} {{')
        for n in coords:
            if shared[n] and min(stars[n])==p:
                out.append(f'puts $f "{n} [nodeDisp {n} 1] [nodeDisp {n} 2] [nodeVel {n} 1] [nodeVel {n} 2]"')
        out.append('}')
    out+=['close $f','remove recorders','wipe','exit']
    output.mkdir(parents=True,exist_ok=True)
    (output/'model.tcl').write_text('\n'.join(out)+'\n')
    return {'counts':counts,'dt':dt,'timed_steps':steps,'warmup_steps':warmup,'dx':dx}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--ranks',type=int,default=1)
    parser.add_argument('--steps',type=int,default=4000)
    parser.add_argument('--warmup',type=int,default=100)
    parser.add_argument('--dx',type=float,default=.05)
    parser.add_argument('--repeats',type=int,default=2)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output=args.output.resolve();args.binary=args.binary.resolve()
    root=Path(__file__).resolve().parents[2]
    if args.output==root or root in args.output.parents:
        parser.error('--output must be outside the source tree')
    if args.ranks<1 or args.repeats<1 or args.steps<1 or args.warmup<1 or args.dx<=0 or args.dx>=1:
        parser.error('positive ranks/repeats/steps/warmup and 0 < dx < 1 required')
    runs=[];history={};state={}
    for repeat in range(args.repeats):
        # Alternating order: legacy/custom/custom/legacy for two repeats.
        for custom in ([False,True] if repeat%2==0 else [True,False]):
            name='combined' if custom else 'legacy'
            path=args.output/f'{repeat}-{name}'
            metadata=generate(path,custom,args.ranks,args.steps,args.dx,args.warmup)
            cmd=[str(args.binary),'model.tcl']
            if args.ranks>1:cmd=['mpirun','--bind-to','core','-np',str(args.ranks),*cmd]
            start=time.monotonic()
            result=subprocess.run(cmd,cwd=path,capture_output=True,text=True,check=False)
            (path/'run.log').write_text(result.stdout+'\n'+result.stderr)
            if result.returncode:raise RuntimeError(f'run failed; see {path}/run.log')
            metrics=[]
            for line in (result.stdout+'\n'+result.stderr).splitlines():
                if line.startswith('METRIC '):
                    _,rank,seconds,rss=line.split();metrics.append((int(rank),float(seconds),int(rss)))
            if len(metrics)!=args.ranks:raise RuntimeError(f'missing rank metrics; see {path}/run.log')
            run={'method':name,'repeat':repeat,'analysis_seconds':max(m[1] for m in metrics),
                 'aggregate_process_peak_rss_MiB':sum(m[2] for m in metrics)/1024,
                 'launch_build_analysis_seconds':time.monotonic()-start,**metadata}
            runs.append(run);print(json.dumps(run),flush=True)
            history[(repeat,name)]={p.name:np.loadtxt(p) for p in path.glob('controls-*.txt')}
            state[(repeat,name)]={p.name:np.loadtxt(p) for p in path.glob('state-*.txt') if p.stat().st_size}
    errors=[]
    for repeat in range(args.repeats):
        worst=0.
        for collection in [history,state]:
            a,b=collection[(repeat,'legacy')],collection[(repeat,'combined')]
            assert a.keys()==b.keys()
            for key in a:
                np.testing.assert_allclose(b[key],a[key],rtol=2.e-8,atol=2.e-11)
                worst=max(worst,float(np.max(np.abs(a[key]-b[key]))))
        errors.append(worst)
    medians={name:{key:float(np.median([r[key] for r in runs if r['method']==name]))
                    for key in ['analysis_seconds','aggregate_process_peak_rss_MiB']}
             for name in ['legacy','combined']}
    summary={'binary':str(args.binary),'host':os.uname().nodename,'ranks':args.ranks,
             'runs':runs,'medians':medians,'max_absolute_history_or_state_difference':errors,
             'note':'Same executable/partition/integration/output; local process peaks include baseline libraries. No PML inference.'}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
