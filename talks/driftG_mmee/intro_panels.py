"""
Opening-slide visual (two panels): individual reproductive variance -> fluctuating abundance.

  TOP    reproduction diagram: identical individuals leave different numbers of offspring.
  BOTTOM one population's abundance over time -- it wanders (demographic stochasticity)
         but persists (extinction is saved for the next slide). Same abundance model as
         the extinction/fixation slide, so it reads as one of those same populations.

Deck palette; on-brand teal; transparent background. Panels rendered separately (each
tuned on its own) then composited with a soft down-arrow, so the causal jump is explicit.
"""
import numpy as np, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image, ImageDraw

INK='#1b1f24'; TEAL='#1b8c98'; LINK='#b2bcc6'; LABEL='#5f6b76'
AXIS='#b7c0c9'; TICK='#5f6b76'; MEANL='#cbd4dc'
mpl.rcParams.update({'font.family':'DejaVu Sans','text.color':INK,'axes.labelcolor':INK,
    'xtick.color':TICK,'ytick.color':TICK,'axes.edgecolor':AXIS,'axes.linewidth':1.1,'svg.fonttype':'none'})

def trim(path, m=24):
    im=Image.open(path); a=np.array(im)[...,3]; ys,xs=np.where(a>0)
    box=(max(xs.min()-m,0),max(ys.min()-m,0),min(xs.max()+m,im.width),min(ys.max()+m,im.height))
    im.crop(box).save(path); return Image.open(path)

# ---------- TOP: reproduction diagram ----------
def render_top(fname):
    counts=[2,0,4,1,0,3,1]; GAP=1.15; fams=[]; x=0.0
    for c in counts:
        if c==0: fams.append((x+0.5,[])); x+=1+GAP
        else:
            blk=[x+i+0.5 for i in range(c)]; fams.append((x+c/2.0,blk)); x+=c+GAP
    xR=x-GAP; YP,YO=1.7,0.0
    fig,ax=plt.subplots(figsize=(10.5,3.6))
    for px,blk in fams:
        for ox in blk: ax.plot([px,ox],[YP,YO],color=LINK,lw=1.4,alpha=0.9,solid_capstyle='round',zorder=1)
    off=[o for f in fams for o in f[1]]
    ax.scatter(off,[YO]*len(off),s=470,color=TEAL,edgecolors='white',linewidths=1.6,zorder=3)
    ax.scatter([f[0] for f in fams],[YP]*len(fams),s=830,color=TEAL,edgecolors='white',linewidths=1.8,zorder=3)
    xL=-1.4
    ax.text(xL,YP,'parents',ha='right',va='center',color=LABEL,fontsize=15)
    ax.text(xL,YO,'offspring',ha='right',va='center',color=LABEL,fontsize=15)
    ax.set_xlim(xL-3.0,xR+0.6); ax.set_ylim(-0.55,2.25); ax.axis('off')
    fig.tight_layout(pad=0.3); fig.savefig(fname,dpi=200,transparent=True,bbox_inches='tight'); plt.close(fig)
    return trim(fname)

# ---------- BOTTOM: single persisting, fluctuating abundance path ----------
def render_bottom(fname, seed=12, v=2.5):
    K,r,T,dt=30.0,0.6,110.0,0.03; c=r/K; rng=np.random.default_rng(seed); steps=int(T/dt)
    t=np.linspace(0,T,steps+1); n=np.empty(steps+1); n[0]=K
    for k in range(steps):
        n[k+1]=max(n[k]+(r-c*n[k])*n[k]*dt+np.sqrt(max(n[k],0)*v)*np.sqrt(dt)*rng.standard_normal(),0)
    fig,ax=plt.subplots(figsize=(10.06,3.0))
    ax.axhline(K,color=MEANL,lw=1.2,ls=(0,(5,4)),zorder=1)      # deterministic level
    ax.axhline(0,color=AXIS,lw=1.1,zorder=2)                     # zero (no extinction here)
    ax.plot(t,n,color=TEAL,lw=2.0,solid_capstyle='round',zorder=3)
    ax.set_xlim(0,T); ax.set_ylim(-2.5,max(n)*1.10)
    ax.set_xlabel('time  \u2192',labelpad=2); ax.set_ylabel('population size,  $n$')
    ax.set_xticks([]); ax.set_yticks([0,K]); ax.set_yticklabels(['$0$','$n_0$'])
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.tick_params(length=3)
    fig.tight_layout(pad=0.3); fig.savefig(fname,dpi=200,transparent=True,bbox_inches='tight'); plt.close(fig)
    print(f'  abundance path: min={n.min():.0f} max={n.max():.0f} std={n.std():.1f} (persists, never 0)')
    return trim(fname)

top=render_top('repro_variance.png')
bot=render_bottom('abundance_single.png')

# ---------- composite with a soft down-arrow ----------
GAP=118; W=max(top.width,bot.width); H=top.height+GAP+bot.height
cv=Image.new('RGBA',(W,H),(0,0,0,0))
cv.paste(top,((W-top.width)//2,0),top)
cv.paste(bot,((W-bot.width)//2,top.height+GAP),bot)
d=ImageDraw.Draw(cv); cx=W//2; y0=top.height+24; y1=top.height+GAP-20; col=(178,190,202,255)
d.line([(cx,y0),(cx,y1-14)],fill=col,width=5)
d.polygon([(cx-14,y1-16),(cx+14,y1-16),(cx,y1)],fill=col)
cv.save('intro_repro_abundance.png')
print('composite',cv.size,' (top',top.size,'/ bottom',bot.size,')')
