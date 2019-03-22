# regular spiking params, lower threshold of vth=20
import math
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import time

from mpl_toolkits.axes_grid1.inset_locator import inset_axes,zoomed_inset_axes, mark_inset

from matplotlib import rc
from scipy.interpolate import interp1d

from xppcall import xpprun, read_numerics, read_pars, read_inits

fontsize=15


matplotlib.rcParams['text.latex.preamble'] = [r' \usepackage{bm} \usepackage{xcolor} \setlength{\parindent}{0pt}']
matplotlib.rcParams.update({'figure.autolayout': True})

rc('text', usetex=True)
rc('font', family='serif', serif=['Computer Modern Roman'])

sizeOfFont = 20
fontProperties = {'weight' : 'bold', 'size' : sizeOfFont}



vdat = np.loadtxt('izk_lc_v_30.dat')
udat = np.loadtxt('izk_lc_u_30.dat')

# first point is just after jump, last point is just before jump
v_plus = vdat[0]
u_plus = udat[0]

v_minus = vdat[-1]
u_minus = udat[-1]

period =  vdat[-1,0]

a=0.02;b=0.2;I=10.


def v_interp(t):
    t = np.mod(t,period)
    return interp1d(vdat[:,0],vdat[:,1])(t)

def u_interp(t):
    t = np.mod(t,period)
    return interp1d(udat[:,0],udat[:,1])(t)

def phase_reset(phi, dx=0., dy=0., steps_per_cycle = 500000,
                num_cycles = 5, return_intermediates=False):

    # total period
    T = period

    steps_before = int(phi/(2*math.pi) * steps_per_cycle) + 1
    t1 = np.linspace(0, phi/(2*math.pi) * T, steps_before)
    #vals1 = odeint(D2Attractor,[LCinit(0,fixedpts),0],t1, args=((a1,b1),(a2,b2),(a3,b3),(a4,b4)))
    vals1 = sim1(t1,[vdat[0,1],udat[0,1]])

    t2 = np.linspace(phi/(2*math.pi) * T, T * num_cycles,steps_per_cycle * num_cycles - steps_before)
    #vals2 = odeint(D2Attractor,list(vals1[-1,:] + np.array([dx, dy])),t2, args=((a1,b1),(a2,b2),(a3,b3),(a4,b4)))
    print vals1[-1,:],t1[-1]
    vals2 = sim1(t2,list(vals1[-1,:] + np.array([dx, dy])))
    print vals2[0,:],t2[0]
    
    """
    crossings = ((vals2[:-1,0] > 0) * (vals2[1:,0] <= 0) 
    * (vals2[1:,1] > 0))
    
    if len(crossings) == 0:
    raise RuntimeError("No complete cycles after the perturbation")
    crossing_fs = ( (vals2[1:,0][crossings] - 0)
    / (vals2[1:,0][crossings]-vals2[:-1,0][crossings]) )
    """
        
    if False:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(t1,vals1[:,0])
        ax.plot(t2,vals2[:,0])
        plt.show()
        
    crossings = ((vals2[:-1,0] <= 0) * (vals2[1:,0] > 0)* (vals2[1:,1] < 0))
    
    print 'crossings', sum(crossings)
    if len(crossings) == 0:
        raise RuntimeError("No complete cycles after the perturbation")
    crossing_fs = ( (vals2[1:,0][crossings] - 0)
                    / (vals2[1:,0][crossings]-vals2[:-1,0][crossings]) )
    crossing_times = (crossing_fs * t2[:-1][crossings]
                      + (1-crossing_fs) * t2[1:][crossings])
    crossing_phases = np.fmod(crossing_times, T)/T * 2 * math.pi
    crossing_phases[crossing_phases > math.pi] -= 2*math.pi
    
    if return_intermediates:
        return dict(t1=t1, vals1=vals1, t2=t2, vals2=vals2,
                    crossings=crossings,
                    crossing_times=crossing_times,
                    crossing_phases=crossing_phases)
    else:
        return -crossing_phases[-1]
    

def izk1(y,t):
    v = y[0];u = y[1]
    return np.array([0.04*v**2. + 5.*v + 140. - u + I,
            a*(b*v-u)])

def adjoint(y,t):
    return np.dot(-np.transpose(A(t)),y)

def LC(t):
    # interpolation of limit cycle
    t = np.mod(t,period)
    return np.array([v_interp(t),u_interp(t)])

# saltation matrix
F1m,F2m = izk1([vdat[-1,1],udat[-1,1]],0)
F1p,F2p = izk1([vdat[0,1],udat[0,1]],0)

print 'F1m,F2m',F1m,F2m
print 'F1p,F2p',F1p,F2p

def S():    
    return np.array([ [F1p, 0.],[F2p-F2m, F1m] ])/F1m
    #return np.array([ [F1p, 0.],[F2p-F2m, 0] ])/F1m


def A(t):
    return np.array([ [2*(0.04)*v_interp(t)+5, -1.],
                      [a*b, -a] ])


def sim1(t,y0):
    #T = 200;dt = 0.01
    #t = np.linspace(0,T,int(T/dt))
    
    

    c=-65;d=2
    
    y = np.zeros((len(t),2))
    y[0,0] = y0[0]
    y[0,1] = y0[1]
    
    for i in range(1,len(t)):
        if y[i-1,0] >= 30:
            y[i,0] = c
            y[i,1] += d
        else:
            dt = t[i]-t[i-1]
            y[i,:] = y[i-1,:] + dt*izk1(y[i-1,:],t[i-1])
    

    return y


def sim2():
    T = -200;dt = -0.0005
    t = np.linspace(0,T,int(T/dt))
    s = 0

    c=-65;d=2
    
    y = np.zeros((len(t),2))
    y[0,0] = .1
    y[0,1] = -.1

    a = np.dot(y[0,:],izk1(LC(t[0]),t[0]))

    y[0,:] /= a

    
    for i in range(1,len(t)):        
        
        if np.abs(s)>=period and (s<dt):
            
            y[i,:] = np.dot(np.transpose(S()),y[i-1,:])
            #y[i,:] = np.dot(S(),y[i-1,:])

            print t[i-1],y[i-1,:],y[i,:],a, (t[i-1]<dt)
            
            s = 0

        else:
            y[i,:] = y[i-1,:] + dt*adjoint(y[i-1,:],t[i])
            

            #if i%1000 == 0:
        s += dt


    if False:
        fig2 = plt.figure()
        ax5 = fig2.add_subplot(111)
        ax5.plot(t,y[:,0])
        ax5.plot(t,y[:,1])
        plt.show()

    # extract 1 period near end and normalize
    jump_pos = (np.abs(y[1:,0]-y[:-1,0])>.1)
    
    idx = np.arange(len(y[1:]))
    start_idx = idx[jump_pos][-2]
    end_idx = idx[jump_pos][-1]
    adjoint_v = y[start_idx:end_idx,0][::-1]
    adjoint_u = y[start_idx:end_idx,1][::-1]
    
    # normalize adjoint ( use oone point because some points are near zero)
    a = np.dot([adjoint_v[0],adjoint_u[0]],izk1(LC(t[0]),t[0]))
    adjoint_v /= a
    adjoint_u /= a

    #for i in range(len(adjoint_v)):
    #    #print adjoint_v[i],izk1(LC(t[i]),t[i])
    #    
    #    #print a
        

    
    
    if False:
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        tt1 = t[start_idx:end_idx]
        tt = np.linspace(0,period,len(t[start_idx:end_idx]))

        ax1.plot(tt,adjoint_v)
        ax2.plot(tt,adjoint_u)

        ax1.set_title('non-normalized voltage iprc (adjoint)')
        ax2.set_title('non-normalized u iprc (adjoint)')
        #ax.plot(t,y[:,0])
        #ax.plot(t,y[:,1])
        #ax.plot(t,LC(t)[0])
        #ax.plot(y[:,0],y[:,1])

        plt.show()

    return adjoint_v,adjoint_u

prc_x = np.loadtxt('matlabprc_dx.csv')
prc_y = np.loadtxt('matlabprc_dy.csv')

def H(t):
    x = np.linspace(0,period,len(prc_x))
    t = np.mod(t,period)
    return interp1d(x,prc_x[::-1])(t)/period

def H2odd(phi,t):
    return H(-phi)-H(phi)

def phase2var(phase):
    # assuming phase in [0,T)
    return v_interp(phase),u_interp(phase)    
    
def phase_vs_full():

    fname = 'izk2.ode'

    pars = read_pars(fname)
    eps = 0.1

    T = 2000;dt = 0.02
    t = np.linspace(0,T,int(T/dt))
    
    y = np.zeros(len(t))

    # three different initial conditions
    
    # solution values of each phase


    phase_inits = [11, # quarter phase
                   -16, # negative quarter phase
                   20 # near antiphase
    ]

    colors = ['tab:blue','tab:orange','tab:green']

    # negative quarter phase (a bit closer to antiphase)

    # almost antiphase
    fig = plt.figure(figsize=(4,3))
    ax1 = fig.add_subplot(111)


    for j in range(len(phase_inits)):

        if j == 0:
            label1 = 'Theory'
            label2 = 'Numerics'
        else:
            label1 = None
            label2 = None

        y[0] = -phase_inits[j]

        for i in range(1,len(t)):
            y[i] = y[i-1] + dt*H2odd(y[i-1],t[i])        

        ax1.plot(t/eps,np.mod(y/period,1),label=label1,color=colors[j])
        

        v1,u1 = phase2var(phase_inits[j])
        v2,u2 = phase2var(0)

        #print v1,u1,v2,u2,type(v1)
        
        npa, vn = xpprun(fname,
                         xppname='xppaut',
                         inits={'v1':float(v1),'u1':float(u1),
                                'v2':float(v2),'u2':float(u2)},
                         parameters={'eps':eps,'total':t[-1]/eps},
                         clean_after=True)

        t_full = npa[:,0]
        sv = npa[:,1:]

        skip = 10000
        ax1.scatter(t_full[::skip][1:],np.mod(sv[:,vn.index('psi')][::skip][1:],1),facecolors='none',edgecolors=colors[j],label=label2)

    ax1.plot([0,T/eps],[0,0],color='k')
    ax1.plot([0,T/eps],[1,1],color='k')
    #ax1.plot([0,T/eps],[-.5,-.5],color='gray',ls='--')
    ax1.plot([0,T/eps],[.5,.5],color='gray',ls='--')

    ax1.set_yticks([0,.5,1])#ax1.set_yticks([-.5,0,.5])
    ax1.set_yticklabels([r'$0$',r'$\pi$',r'$2\pi$'])#ax1.set_yticklabels([r'$-T/2$',r'$0$',r'$T/2$'])
    #ax1.yticks([-.5,0,.5],[r'$-T/2$',r'$0$',r'$T/2$'])


    ax1.set_ylim(-.1,1.1)
    ax1.set_xlim(0,T/eps)
    ax1.set_xlabel(r'$t$',fontsize=fontsize)
    ax1.set_ylabel(r'$\phi$',fontsize=fontsize)

    ax1.legend()
    
    #plt.show()

def plot_H():
    """
    plot H functions and phase RHS
    """

        
    adjoint_v,adjoint_u = sim2()
    x_adjoint = np.linspace(0,period,len(adjoint_v))

    prc_x = np.loadtxt('matlabprc_dx.csv')
    prc_y = np.loadtxt('matlabprc_dy.csv')

    x = np.linspace(0,period,len(prc_x))
    

    fig = plt.figure(figsize=(6,5))
    ax11 = fig.add_subplot(221)
    ax12 = fig.add_subplot(222)
    ax21 = fig.add_subplot(223)
    ax22 = fig.add_subplot(224)
    
    ax11.plot(x_adjoint,adjoint_v,lw=2,color='k')
    ax11.scatter(x[::3],prc_x[::3],facecolors='none',edgecolors='tab:blue')

    H = prc_x[::-1]/period
    ax12.plot(x,H,color='black')

    
    ax22.plot([0,period],[0,0],color='gray',lw=.5)
    ax22.plot(x,H[::-1]-H,lw=2,color='k')

    
    axins = inset_axes(ax22, width=.7, height=0.6,loc='lower right')
    axins.plot([0,period],[0,0],color='gray',lw=.5)
    axins.plot(x,H[::-1]-H,lw=2,color='k')
    axins.set_xlim(18,27)
    axins.set_ylim(-.0002,.0002)
    axins.set_xticks([])
    axins.set_yticks([])

    axins.annotate('',
                   xy=(23, 0), 
                   xytext=(25, .00), 
                   arrowprops=dict(arrowstyle="-|>, head_width=.25, head_length=.5",facecolor='k'))

    axins.annotate('',
                   xy=(22, 0), 
                   xytext=(20, .00), 
                   arrowprops=dict(arrowstyle="-|>, head_width=.25, head_length=.5",facecolor='k'))

    
    mark_inset(ax22, axins, loc1=1, loc2=3, fc="none", ec="0.5")

    ax22.annotate('',
                  xy=(5, 0), xycoords='data',
                  xytext=(2, 0), textcoords='offset points',
                  arrowprops=dict(arrowstyle="-|>, head_width=.4, head_length=1",facecolor='k'))

    ax22.annotate('',
                  xy=(38, 0), 
                  xytext=(35, .00), 
                  arrowprops=dict(arrowstyle="-|>, head_width=.4, head_length=1",facecolor='k'))

    
    
    ax21.plot(x_adjoint,adjoint_u,lw=2,color='k')
    ax21.scatter(x[::3],prc_y[::3],facecolors='none',edgecolors='tab:blue')

    ax21.set_xlabel(r'$\phi$',fontsize=fontsize)
    ax22.set_xlabel(r'$\phi$',fontsize=fontsize)

    ax11.set_title(r'\textbf{A}',fontsize=fontsize,x=0)
    ax12.set_title(r'\textbf{B}',fontsize=fontsize,x=0)
    ax21.set_title(r'\textbf{C}',fontsize=fontsize,x=0)
    ax22.set_title(r'\textbf{D}',fontsize=fontsize,x=0)
    

    ax11.set_ylabel(r'$Z_v(\phi)$',fontsize=fontsize)
    ax21.set_ylabel(r'$Z_u(\phi)$',fontsize=fontsize)
    ax12.set_ylabel(r'$H(\phi)$',fontsize=fontsize)
    ax22.set_ylabel(r'$H(-\phi)-H(\phi)$',fontsize=fontsize)

    
    ax11.set_xticks([0,period/2,period])
    ax12.set_xticks([0,period/2,period])
    ax21.set_xticks([0,period/2,period])
    ax22.set_xticks([0,period/2,period])
    
    
    ax11.set_xticklabels([r'$0$',r'$T/2$',r'$T$'])
    ax12.set_xticklabels([r'$0$',r'$T/2$',r'$T$'])
    ax21.set_xticklabels([r'$0$',r'$T/2$',r'$T$'])
    ax22.set_xticklabels([r'$0$',r'$T/2$',r'$T$'])

    ax11.set_xlim(0,period)
    ax12.set_xlim(0,period)
    ax21.set_xlim(0,period)
    ax22.set_xlim(0,period)

    ax11.set_xticks([0,period/2,period])
    ax12.set_xticks([0,period/2,period])
    ax21.set_xticks([0,period/2,period])
    ax22.set_xticks([0,period/2,period])
    
    ax11.set_xticklabels([r'$0$',r'$\pi$',r'$2\pi$'])
    ax12.set_xticklabels([r'$0$',r'$\pi$',r'$2\pi$'])
    ax21.set_xticklabels([r'$0$',r'$\pi$',r'$2\pi$'])
    ax22.set_xticklabels([r'$0$',r'$\pi$',r'$2\pi$'])
    
    #ax22.set_title('rhs of weak coupling',fontsize=fontsize)
    
    #plt.show()

    
def main():

    phase_vs_full()

    plot_H()

    plt.show()

    # numerical iPRC not working in python due to discontinuity and need for high tolerance/accuracy in numerics. see matlab files for numerical iPRC estimation.
    if False:
        # total perturbations
        total_perts = 30

        # phase values where perturbations will be applied
        phis = np.linspace(0,2*math.pi,total_perts)

        pert = 1e-1 # keep perturbation to single variable for now
        print "Calculating iPRC via direct method for perturbations in the x direction..."
        x_prc = np.array([
            phase_reset(phi, dx=pert, dy=0)
            for phi in phis
        ])
        print " "

        print "Calculating iPRC via direct method for perturbations in the y direction..."
        y_prc = np.array([
            phase_reset(phi, dx=0, dy=pert)
            for phi in phis
        ])
        print " "

        #fig2 = plt.figure(figsize=(6,6))
        fig2 = plt.figure()
        ax1 = fig2.add_subplot(211)
        ax2 = fig2.add_subplot(212)

        #axes2 = fig2.add_axes([.1,.1,.8,.8])
        #mp.figure()
        ax1.title('iPRC of voltage (blue) and u (gray) (direct perturbations)')
        #mp.title('iPRC of 2D Glass Network (Inhibitory Feedback): Analytical vs Direct')
        #p2, = axes2.plot(np.linspace(0,1,len(adjoint_solx)),adjoint_soly,color='.5',lw=2) # analytic
        #p1, = axes2.plot(np.linspace(0,1,len(adjoint_soly)),adjoint_solx,color='b',lw=2)
        ax1.plot(np.linspace(0,1,total_perts),y_prc/pert/(2*math.pi),marker='o',linestyle='None',color='.5',markeredgecolor='.5') # direct
        ax2.plot(np.linspace(0,1,total_perts),x_prc/pert/(2*math.pi),'bo')

        #mp.legend([p1,p2,p3,p4],['Adjoint x', 'Adjoint y', 'Direct x','Direct y'])
        #mp.xlim((0,1))
        #mp.ylim((0,1))

        plt.show()

if __name__ == "__main__":
    main()
