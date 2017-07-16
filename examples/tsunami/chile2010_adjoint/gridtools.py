
from __future__ import print_function
import numpy as np

def grid_eval(X, Y, Q, xout, yout, return_ma=True):
    from scipy.interpolate import RegularGridInterpolator
    from numpy import ma  # for masked arrays
    
    Qdim = Q.ndim
    if Qdim == 2:
        # change to 3d array of shape (1, Q.shape[0], Q.shape[1]):
        Q = np.array([Q])
        
    nvars = Q.shape[0]  # number of arrays to interpolate

    ndim_out = len(xout.shape)
    xout1 = np.ravel(xout)
    yout1 = np.ravel(yout)
    x1 = X[:,0]
    y1 = Y[0,:]
    dx = x1[1] - x1[0]
    dy = y1[1] - y1[0]

    # augment Q with border of values on all 4 sides:
    x1 = np.hstack((x1[0]-0.501*dx, x1, x1[-1]+0.501*dx))
    y1 = np.hstack((y1[0]-0.501*dy, y1, y1[-1]+0.501*dy))
    Q1 = np.empty((nvars,len(x1),len(y1)))
    Q1[:,1:-1, 1:-1] = Q   # center portion
    Q1[:,1:-1,0] = Q[:,:,0]
    Q1[:,1:-1,-1] = Q[:,:,-1]
    Q1[:,0,1:-1] = Q[:,0,:]
    Q1[:,-1,1:-1] = Q[:,-1,:]
    # corners:
    Q1[:,0,0] = Q[:,0,0]
    Q1[:,0,-1] = Q[:,0,-1]
    Q1[:,-1,-1] = Q[:,-1,-1]
    Q1[:,-1,0] = Q[:,-1,0]

    #print(x1)
    #print(y1)
    qout = np.empty([nvars]+list(xout.shape))
    for k in range(nvars):
        evalfunc = RegularGridInterpolator((x1,y1), Q1[k,:,:], method='nearest',
                bounds_error=False, fill_value=np.nan)
        xyout = np.vstack((xout1,yout1)).T
        qout_k = evalfunc(xyout)
        #print(qout1)
        if ndim_out == 2:
            qout_k = np.reshape(qout_k, xout.shape)   
            qout[k,:,:] = qout_k
        else:
            qout[k,:] = qout_k

    if Qdim==2 and ndim_out==2:
        qout = qout[0,:,:]  # revert back to 2d array
    if Qdim==2 and ndim_out==1:
        qout = qout[0,:]  # revert back to 1d array

    if return_ma:
        # convert from an array with nan's to a masked array:
        qout = ma.masked_where(qout != qout, qout)

    #print('type is %s' % type(qout))

    return qout
    
def grid_output(framesoln, out_var, xout, yout, levels='all', return_ma=True):
    from numpy import ma  # for masked arrays
    if levels == 'all':
        levels = range(1,100)  # more levels than will ever use

    qout = np.empty(xout.shape)
    qout[:] = np.nan
    xmin = xout.min()
    xmax = xout.max()
    ymin = yout.min()
    ymax = yout.max()
    for stateno,state in enumerate(framesoln.states):
        state = framesoln.states[stateno]
        patch = state.patch

        if patch.level not in levels:
            # skip this patch
            continue

        if (xmin > state.grid.x.upper) or (xmax < state.grid.x.lower) \
                or (ymin > state.grid.y.upper) or (ymax < state.grid.y.lower):
            # no overlap
            continue
            
        #print('overlap at level %i' % patch.level)
        Xc,Yc = state.grid.c_centers
        xc = Xc[:,0]
        yc = Yc[0,:]

        if type(out_var) == int:
            Q = state.q[out_var, :, :]
        else:
            Q = out_var(state.q)
        #import pdb; pdb.set_trace()

        qout1 = grid_eval(Xc, Yc, Q, xout, yout, return_ma=False)
        #nononan = np.count_nonzero(qout1==qout1)  # number of non-nan elements
        #print('%i new non-nans' % nononan)
        #print(qout1)
        qout = np.where(np.isnan(qout1), qout, qout1)
        if return_ma:
            # convert from an array with nan's to a masked array:
            qout = ma.masked_where(qout != qout, qout)
        #print(qout)
    return qout

        
def sea_surface(q):
    h = q[0,:,:]
    eta = q[3,:,:]
    surf = np.where(h>0, eta, np.nan)
    return surf

def load_frame(frameno, outdir='_output'):
    from clawpack.visclaw.data import ClawPlotData
    plotdata = ClawPlotData()
    plotdata.outdir = outdir
    framesoln = plotdata.getframe(frameno, plotdata.outdir)
    return framesoln

def adjoint_ip_chile2010(frameno_f, frameno_a, outdir_f='_output',
               outdir_a='adjoint/_output'):

    from numpy import ma  # for masked arrays

    #x = np.linspace(-110,-65,181)
    #y = np.linspace(-60,0,241)
    x = np.linspace(-110,-65,91)
    y = np.linspace(-60,0,121)
    xout,yout = np.meshgrid(x,y,indexing='ij')

    frame_f = load_frame(frameno_f, outdir_f)
    q_f = grid_output(frame_f, lambda q:q, xout, yout)
    eta_f = ma.masked_where(q_f[0,:,:]<1e-3, q_f[3,:,:])

    frame_a = load_frame(frameno_a, outdir_a)
    q_a = grid_output(frame_a, lambda q:q, xout, yout)
    eta_a = ma.masked_where(q_a[0,:,:]<1e-3, q_a[3,:,:])

    inner_product = abs(eta_f*eta_a + q_f[1,:,:]*q_a[1,:,:] + \
                    q_f[2,:,:]*q_a[2,:,:])

    t_f = frame_f.t / 3600.
    t_a = frame_a.t / 3600.
    #return eta_f, t_f, eta_a, t_a, inner_product

    #def test_adjoint_ip_chile2010(frameno_f=10, frameno_a=6):
    import matplotlib.pyplot as plt
    from clawpack.visclaw import geoplot, colormaps


    #eta_f, t_f, eta_a, t_a, inner_product = \
    #    adjoint_ip(frameno_f, frameno_a, xout, yout)

    plt.figure(1,figsize=(13,5))
    plt.clf()
    plt.subplot(131)
    plt.pcolor(xout,yout,eta_f,cmap=geoplot.tsunami_colormap)
    plt.clim(-.2,.2)
    plt.plot([-86.392], [-17.975],'ko')
    land = ma.masked_where(q_f[0,:,:]>0, q_f[3,:,:])
    plt.contourf(xout, yout, land, [0,10000], colors=['g'])
    plt.title('Forward solution at t = %g' % t_f)
    plt.gca().set_aspect(1./np.cos(-30*np.pi/180.))

    plt.subplot(132)
    plt.pcolor(xout,yout,inner_product,cmap=colormaps.white_red)
    plt.clim(0,0.0001)
    plt.plot([-86.392], [-17.975],'ko')
    plt.contourf(xout, yout, land, [0,10000], colors=['g'])
    plt.title('Inner product')
    plt.gca().set_aspect(1./np.cos(-30*np.pi/180.))

    plt.subplot(133)
    plt.pcolor(xout,yout,eta_a,cmap=geoplot.tsunami_colormap)
    plt.clim(-.002,.002)
    plt.plot([-86.392], [-17.975],'ko')
    plt.contourf(xout, yout, land, [0,10000], colors=['g'])
    plt.title('Adjoint solution at t = -%g' % t_a)
    plt.gca().set_aspect(1./np.cos(-30*np.pi/180.))


