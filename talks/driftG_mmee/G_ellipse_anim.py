#
# G-matrix drift with a simple (illustrative) recombination term.
#
#   dG = -(v/n) G dt + sqrt(v/n) * sqrt(2) * sqrtG X sqrtG * sqrt(dt)      [drift, paper eq.21]
#        - r_rec * offdiag(G) dt                                          [recombination, illustrative]
#
# Drift drives the genetic correlation to +/-1 (ellipse collapses to a line, random
# orientation across replicates -- the Phillips picture).  The recombination term
# relaxes the *covariance* (off-diagonal) back toward zero, i.e. pushes rho -> 0, so
# with recombination the ellipses stay round.  This is only a cartoon of the
# drift-recombination balance, not the full multilocus model.
#
import numpy as np, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.animation import FuncAnimation, PillowWriter

# ---- deck palette (daytime.css: white slides) ----
INK='#1b1f24'; GUIDE='#dfe4ea'; SUB='#8a97a6'
REPCOL=['#1b8c98','#c2185b','#5a4cc7','#3c8500','#007c91']   # acid, pink, violet, lime, cyan

def sqrtm_sym(G):
    w,V=np.linalg.eigh(G); w=np.clip(w,1e-12,None); return (V*np.sqrt(w))@V.T

def simulate(G0,v,n0,T,N,R,r_rec=0.0,seed=0):
    rng=np.random.default_rng(seed); dt=T/N; delta=v/n0
    hist=np.zeros((R,N+1,2,2))
    for r in range(R):
        G=G0.copy(); hist[r,0]=G
        for k in range(1,N+1):
            sG=sqrtm_sym(G)
            a,b=rng.standard_normal(),rng.standard_normal(); c=rng.standard_normal()/np.sqrt(2)
            X=np.array([[a,c],[c,b]])
            stoch=np.sqrt(2)*(sG@X@sG)
            off=G-np.diag(np.diag(G))                       # off-diagonal (covariance) part
            G=G-delta*G*dt-r_rec*off*dt+np.sqrt(delta)*np.sqrt(dt)*stoch
            G=0.5*(G+G.T)
            w,V=np.linalg.eigh(G); G=(V*np.clip(w,1e-9,None))@V.T
            hist[r,k]=G
    return hist

def corr(G):
    d=np.sqrt(np.diag(G)); return G/np.outer(d,d)

def rho(G): return G[0,1]/np.sqrt(G[0,0]*G[1,1])

# ---- validate the illustration: drift -> |rho|~1 ; drift+recomb -> rho~0 ----
Hd=simulate(np.eye(2),1.0,80.0,150.0,1500,60,r_rec=0.0,seed=1)
for rr in [0.05,0.10,0.20]:
    Hr=simulate(np.eye(2),1.0,80.0,150.0,1500,60,r_rec=rr,seed=1)
    rd=np.array([abs(rho(Hd[i,-1])) for i in range(60)])
    rrho=np.array([abs(rho(Hr[i,-1])) for i in range(60)])
    print(f"r_rec={rr:.2f}:  drift-only mean|rho|={rd.mean():.2f} (frac>0.9={np.mean(rd>0.9):.2f}) | +recomb mean|rho|={rrho.mean():.2f} (frac>0.9={np.mean(rrho>0.9):.2f})")

# ============================ RENDERING ============================
def ell_params(M, fixed_major=None):
    w,V=np.linalg.eigh(M); w=np.clip(w,0,None)
    i1,i0=(1,0) if w[1]>=w[0] else (0,1)
    major=2*np.sqrt(w[i1]); minor=2*np.sqrt(w[i0])
    ang=np.degrees(np.arctan2(V[1,i1],V[0,i1]))
    if fixed_major and major>1e-9:
        s=fixed_major/major; major*=s; minor*=s
    return major,minor,ang

def setup_cell(ax, guide='circle'):
    ax.set_xlim(-1.7,1.7); ax.set_ylim(-1.7,1.7); ax.set_aspect('equal')
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    if guide=='circle':
        ax.add_patch(Ellipse((0,0),2,2,fill=False,ec=GUIDE,lw=1.1,ls=(0,(4,4))))
    ax.plot([-1.7,1.7],[0,0],color=GUIDE,lw=0.8,zorder=0)
    ax.plot([0,0],[-1.7,1.7],color=GUIDE,lw=0.8,zorder=0)

# ---------- GIF A: drift-only G-matrix orientation (opening, next to Phillips) ----------
def make_orient_gif(fname, R=5, seed=6):
    H=simulate(np.eye(2),1.0,60.0,260.0,1700,R,r_rec=0.0,seed=seed)
    Nst=H.shape[1]-1; frames=np.linspace(0,Nst,90).astype(int)
    fig,axes=plt.subplots(1,R,figsize=(11,2.5)); fig.patch.set_alpha(0)
    ells=[]
    for c in range(R):
        setup_cell(axes[c], guide='circle')
        e=Ellipse((0,0),1,1,fc=REPCOL[c],ec=REPCOL[c],alpha=0.32,lw=2.4); axes[c].add_patch(e); ells.append(e)
    fig.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.02,wspace=0.08)
    def upd(fi):
        k=frames[fi]
        for c in range(R):
            mj,mn,an=ell_params(H[c,k],fixed_major=2.4)
            ells[c].width,ells[c].height,ells[c].angle=mj,mn,an
        return ells
    FuncAnimation(fig,upd,frames=len(frames),blit=False).save(fname,writer=PillowWriter(fps=14),savefig_kwargs={'transparent':True})
    plt.close(fig)
    print('saved',fname,'final |rho| per rep:',[f'{abs(rho(H[c,-1])):.2f}' for c in range(R)])

# ---------- GIF B: drift vs drift+recombination (correlation ellipses) ----------
def make_compare_gif(fname, R=5, seed=6, r_rec=0.10):
    Hd=simulate(np.eye(2),1.0,60.0,260.0,1700,R,r_rec=0.0,seed=seed)
    Hr=simulate(np.eye(2),1.0,60.0,260.0,1700,R,r_rec=r_rec,seed=seed)
    Nst=Hd.shape[1]-1; frames=np.linspace(0,Nst,90).astype(int)
    fig,axes=plt.subplots(2,R,figsize=(11,4.7)); fig.patch.set_alpha(0)
    ells=[]
    for row,(H,col) in enumerate([(Hd,'#c2185b'),(Hr,'#1b8c98')]):
        for c in range(R):
            setup_cell(axes[row,c], guide='circle')
            e=Ellipse((0,0),1,1,fc=col,ec=col,alpha=0.32,lw=2.4); axes[row,c].add_patch(e); ells.append((e,H))
    axes[0,0].set_ylabel('drift only',fontsize=15,color=INK,labelpad=8); axes[0,0].yaxis.label.set_visible(True)
    axes[1,0].set_ylabel('drift +\nrecombination',fontsize=15,color=INK,labelpad=8); axes[1,0].yaxis.label.set_visible(True)
    fig.subplots_adjust(left=0.11,right=0.99,top=0.99,bottom=0.02,wspace=0.08,hspace=0.10)
    def upd(fi):
        k=frames[fi]; idx=0; arts=[]
        for row in range(2):
            for c in range(R):
                e,H=ells[idx]
                mj,mn,an=ell_params(corr(H[c,k]))
                e.width,e.height,e.angle=mj,mn,an; arts.append(e); idx+=1
        return arts
    FuncAnimation(fig,upd,frames=len(frames),blit=False).save(fname,writer=PillowWriter(fps=14),savefig_kwargs={'transparent':True})
    plt.close(fig)
    print('saved',fname)
    print('  drift-only final |rho|:',[f'{abs(rho(Hd[c,-1])):.2f}' for c in range(R)])
    print('  +recomb   final |rho|:',[f'{abs(rho(Hr[c,-1])):.2f}' for c in range(R)])

make_orient_gif('G_drift_orient.gif', R=5, seed=6)
make_compare_gif('G_recomb_compare.gif', R=5, seed=6, r_rec=0.10)
