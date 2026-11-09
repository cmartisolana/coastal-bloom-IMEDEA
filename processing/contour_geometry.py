"""Contour-distance preselection of META4.0 eddies for the occupancy calculation."""
from pathlib import Path
import os, sys, json
import numpy as np
from scipy.spatial import cKDTree
os.chdir(Path(__file__).resolve().parents[1])
Path('results/rebuilt').mkdir(parents=True,exist_ok=True)
tag=sys.argv[1]
x0,x1,y0,y1=json.load(open('config/regions.json'))[tag]['box']
R=6371.0088
def xyz(lon,lat):
    lo,la=np.radians(lon),np.radians(lat)
    return np.stack([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)],-1)
cz=np.load('data/processed/ne10m_coast_filtered_2km.npz')
clon=((cz['coast_lon']+180)%360)-180;clat=cz['coast_lat']
nr=(clon>x0-8)&(clon<x1+8)&(clat>y0-8)&(clat<y1+8)
tree=cKDTree(xyz(clon[nr],clat[nr]))
def distance(lo,la):return R*2*np.arcsin(np.clip(tree.query(xyz(lo,la))[0]/2,0,1))
for pol in ['anticyclonic','cyclonic']:
    z=np.load(f'data/intermediate/meta4_{tag}_{pol}.npz')
    dc=distance(z['lon'],z['lat'])
    L,A=z['contour_lon'].astype(float),z['contour_lat'].astype(float)
    Ld=np.concatenate([L,.5*(L+np.roll(L,-1,1))],1)
    Ad=np.concatenate([A,.5*(A+np.roll(A,-1,1))],1)
    near=dc<500;dv=np.full(Ld.shape,np.inf)
    dv[near]=distance(Ld[near].ravel(),Ad[near].ravel()).reshape(-1,Ld.shape[1])
    inside=(Ld>=x0)&(Ld<=x1)&(Ad>=y0)&(Ad<=y1)
    np.savez_compressed(f'results/rebuilt/s2geom_{tag}_{pol}.npz',dmin_contour_inbox_km=np.where(inside,dv,np.inf).min(1))
