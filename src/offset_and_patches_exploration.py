"""Tonight's checks: (1) EMIT<->S2 co-registration search, (2) candidate dense-vegetation patches for visual mangrove check.
Usage: python3 -I night_check.py S2_DIR EMIT_DIR OUT_DIR"""
import tifffile,glob,sys,json,os,numpy as np
from pyproj import Transformer
from scipy import ndimage as ndi
S2,EM,OUT=sys.argv[1:4];os.makedirs(OUT,exist_ok=True)
bands=['B02','B03','B04','B05','B06','B07','B08','B8A','B11']
g=lambda b:tifffile.imread(glob.glob(f'{S2}/*_{b}_(Raw)*')[0]).astype(float)
B={b:g(b) for b in bands}
tg=tifffile.TiffFile(glob.glob(f'{S2}/*_B04_(Raw)*')[0]).pages[0].tags
sc=tg['ModelPixelScaleTag'].value;tp=tg['ModelTiepointTag'].value
x0,y0,dx,dy=tp[3],tp[4],sc[0],sc[1]
valid=np.all([B[b][...,1]>0 for b in bands],0)
red,nir=B['B04'][...,0],B['B08'][...,0]
ndvi=(nir-red)/(nir+red+1e-9)
H,W=valid.shape
X,Y=np.meshgrid(x0+(np.arange(W)+.5)*dx,y0-(np.arange(H)+.5)*dy)
e=tifffile.imread(f'{EM}/emit.tif');t=np.loadtxt(f'{EM}/emit.tfw');px,py,ox,oy=t[0],t[3],t[4],t[5]
ev=e[...,0]>-9000
en=(e[...,65]-e[...,38])/(e[...,65]+e[...,38]+1e-9)
nr,nc=ev.shape
to=Transformer.from_crs(32639,4326,always_xy=True)
def agg(sE,sN):
    lon,lat=to.transform(X+sE,Y+sN)
    c=np.floor((lon-(ox-px/2))/px).astype(int);r=np.floor(((oy+py/2)-lat)/(-py)).astype(int)
    ok=valid&(r>=0)&(r<nr)&(c>=0)&(c<nc)
    idx=(r*nc+c)[ok];n=np.bincount(idx,minlength=nr*nc);s=np.bincount(idx,weights=ndvi[ok],minlength=nr*nc)
    with np.errstate(invalid='ignore',divide='ignore'):m=(s/n).reshape(nr,nc)
    return m,n.reshape(nr,nc)
res=[]
for sE in range(-120,121,15):
    for sN in range(-120,121,15):
        m,n=agg(sE,sN);k=ev&(n>=20)&np.isfinite(m)
        if k.sum()<400:continue
        a,b=en[k],m[k];res.append((sE,sN,int(k.sum()),np.corrcoef(a,b)[0,1],np.sqrt(((a-b)**2).mean())))
res.sort(key=lambda z:-z[3])
import csv
with open(f'{OUT}/offset_search.csv','w',newline='') as f:
    w=csv.writer(f);w.writerow(['shift_east_m(S2 moved)','shift_north_m','n_px','r','rmse']);w.writerows([[a,b,c,round(d,4),round(r_,4)] for a,b,c,d,r_ in res])
base=[z for z in res if z[0]==0 and z[1]==0][0]
print('baseline',base);print('top5',[tuple(round(v,3) for v in z) for z in res[:5]])
# surface of r for plot
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
E_=sorted({z[0] for z in res});N_=sorted({z[1] for z in res});Z=np.full((len(N_),len(E_)),np.nan)
for a,b,c,d,r_ in res:Z[N_.index(b),E_.index(a)]=d
fig,ax=plt.subplots(1,2,figsize=(13,5))
im=ax[0].imshow(Z,origin='lower',extent=[E_[0]-7.5,E_[-1]+7.5,N_[0]-7.5,N_[-1]+7.5],cmap='viridis');plt.colorbar(im,ax=ax[0],label='r');ax[0].plot(0,0,'r+');ax[0].set_title('Correlation vs shift of S2 (m)');ax[0].set_xlabel('east');ax[0].set_ylabel('north')
# (2) candidate patches: dense canopy proxy
dense=valid&(ndvi>0.6)&(B['B08'][...,0]>0.3)
dense=ndi.binary_opening(dense,iterations=1)
lab,nl=ndi.label(dense);sizes=ndi.sum(dense,lab,range(1,nl+1))
keep=[i+1 for i,s in enumerate(sizes) if s*dx*dy>=1000]  # >=0.1 ha
ax[1].imshow(np.where(valid,ndvi,np.nan),cmap='RdYlGn',vmin=-.2,vmax=.8)
feats=[];kml=[]
inv=Transformer.from_crs(32639,4326,always_xy=True)
for n_,i in enumerate(sorted(keep,key=lambda i:-sizes[i-1]),1):
    mk=lab==i;ys,xs=np.nonzero(mk);cy,cx=ys.mean(),xs.mean();ax[1].text(cx,cy,str(n_),color='k',fontsize=8,ha='center',weight='bold')
    ax[1].contour(mk,[0.5],colors='k',linewidths=.6)
    # polygon from pixel union via shapely
    from shapely.geometry import box,mapping;from shapely.ops import unary_union
    polys=unary_union([box(x0+x*dx,y0-(y+1)*dy,x0+(x+1)*dx,y0-y*dy) for y,x in zip(ys,xs)]).simplify(1)
    from shapely.ops import transform as tr_
    p84=tr_(lambda a,b,z=None:inv.transform(a,b),polys)
    mean_ndvi=float(ndvi[mk].mean());ha=float(sizes[i-1]*dx*dy/1e4)
    feats.append({'type':'Feature','properties':{'id':n_,'area_ha':round(ha,2),'mean_ndvi':round(mean_ndvi,2),'mean_nir':round(float(B['B08'][...,0][mk].mean()),2),'mean_swir':round(float(B['B11'][...,0][mk].mean()),2),'decision':''},'geometry':mapping(p84)})
json.dump({'type':'FeatureCollection','features':feats},open(f'{OUT}/candidate_dense_vegetation.geojson','w'))
ax[1].set_title(f'Candidate dense-vegetation patches ({len(feats)})');ax[1].axis('off')
plt.savefig(f'{OUT}/night_check.png',dpi=110,bbox_inches='tight')
print('patches',[(f['properties']['id'],f['properties']['area_ha'],f['properties']['mean_ndvi']) for f in feats][:15])
