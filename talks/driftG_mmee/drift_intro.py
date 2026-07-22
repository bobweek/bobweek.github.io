"""
Collapsed opening slide: variation of types + reproductive variance -> random genetic drift.

  TOP    reproduction diagram with TWO heritable types (cyan / pink alleles). Each parent
         leaves a different number of offspring, and offspring inherit the parent's type.
  BOTTOM the fixation/loss figure (Wright-Fisher drift): the type frequency wanders until
         one type fixes (cyan paths -> 1) and the other is lost (pink paths -> 0).

Same two colours run through both panels, so the jump reads directly: cyan & pink
individuals -> their frequencies drift -> fixation of one, loss of the other.
"""
import numpy as np, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
from drift_figs import make_fixation           # reuse the fixation/loss panel (cyan/pink)

CYAN='#007c91'; PINK='#c2185b'; LINK='#b2bcc6'; LABEL='#5f6b76'
mpl.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})

def trim(path, m=24):
    im=Image.open(path); a=np.array(im)[...,3]; ys,xs=np.where(a>0)
    box=(max(xs.min()-m,0),max(ys.min()-m,0),min(xs.max()+m,im.width),min(ys.max()+m,im.height))
    im.crop(box).save(path); return Image.open(path)

def render_top(fname):
    # two heritable types; both have a "jackpot" (3) and a "failure" (0) -> reads as chance
    types  = ['C','P','C','C','P','P','C']
    counts = [ 3 , 0 , 1 , 0 , 2 , 3 , 1 ]     # cyan offspring 5, pink offspring 5
    cmap={'C':CYAN,'P':PINK}
    GAP=1.15; fams=[]; x=0.0
    for t,c in zip(types,counts):
        if c==0:
            fams.append((x+0.5, [], t)); x+=1+GAP
        else:
            blk=[x+i+0.5 for i in range(c)]; fams.append((x+c/2.0, blk, t)); x+=c+GAP
    xR=x-GAP; YP,YO=1.7,0.0
    fig,ax=plt.subplots(figsize=(11.0,3.6))
    for px,blk,t in fams:                       # links behind
        for ox in blk:
            ax.plot([px,ox],[YP,YO],color=LINK,lw=1.4,alpha=0.9,solid_capstyle='round',zorder=1)
    for t in ['C','P']:                         # circles, coloured by type
        ox=[o for (px,blk,tt) in fams if tt==t for o in blk]
        pxs=[px for (px,blk,tt) in fams if tt==t]
        ax.scatter(ox,[YO]*len(ox),  s=470, color=cmap[t], edgecolors='white', linewidths=1.6, zorder=3)
        ax.scatter(pxs,[YP]*len(pxs),s=830, color=cmap[t], edgecolors='white', linewidths=1.8, zorder=3)
    xL=-1.4
    ax.text(xL,YP,'parents',   ha='right',va='center',color=LABEL,fontsize=15)
    ax.text(xL,YO,'offspring', ha='right',va='center',color=LABEL,fontsize=15)
    ax.set_xlim(xL-3.0,xR+0.6); ax.set_ylim(-0.55,2.25); ax.axis('off')
    fig.tight_layout(pad=0.3); fig.savefig(fname,dpi=200,transparent=True,bbox_inches='tight'); plt.close(fig)
    return trim(fname)

top = render_top('repro_variance_types.png')
make_fixation('fixation.png', seed=0)           # cyan fixes / pink fixes / gray segregating
bot = trim('fixation.png')
if bot.width != top.width:                      # match widths for clean stacking
    nw=top.width; nh=int(round(bot.height*nw/bot.width)); bot=bot.resize((nw,nh), Image.LANCZOS)

GAP=118; W=max(top.width,bot.width); H=top.height+GAP+bot.height
cv=Image.new('RGBA',(W,H),(0,0,0,0))
cv.paste(top,((W-top.width)//2,0),top)
cv.paste(bot,((W-bot.width)//2,top.height+GAP),bot)
d=ImageDraw.Draw(cv); cx=W//2; y0=top.height+24; y1=top.height+GAP-20; col=(178,190,202,255)
d.line([(cx,y0),(cx,y1-14)],fill=col,width=5)
d.polygon([(cx-14,y1-16),(cx+14,y1-16),(cx,y1)],fill=col)
cv.save('drift_intro.png')
print('top',top.size,' bottom',bot.size,' composite',cv.size)
