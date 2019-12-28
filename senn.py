# import modules
import numpy as np


def write_ode(node_num, first_nl_node=-1, last_nl_node=-1):
    """
    'write_ode' writes the system of ODE including Frankenhaeuser and Huxley equations

    Parameters
    ----------
    node_num (int):
        number of nodes in the fiber
    first_nl_node (int):
        first nonlinear node (zero-based numeration --> 0 = first)
    last_nl_node (int):
        last nonlinear node (zero-based numeration --> node_num - 1 = last)

    Returns
    -------
    the file 'eqdiff.py' is written for later import and use

    """

    # index for nonlinear nodes: (Frankenhaeuser and Huxley equations)
    ifh = range(first_nl_node, last_nl_node + 1)
    # additional 4 equations for each nonlinear node
    iextra = range(node_num, node_num + (last_nl_node - first_nl_node + 1) * 4)

    # open files where ode system is defined
    fid = open('eqdiff.py', 'w')

    # import modules
    fid.write("import numpy as np\n\n\n")

    # define constant
    fid.write('def eqdiff(t, y, Ga, Gm, Cm, ve, d, l, Vr):\n')
    fid.write('\n')
    fid.write('    # Parameters\n')
    fid.write('    PNa = 8*1e-3*1e-2 # (m/s)\n')
    fid.write('    PK = 1.2*1e-3*1e-2 # (m/s)\n')
    fid.write('    Pp = 0.54*1e-3*1e-2 # (m/s)\n')
    fid.write('    gl = 30.3*1e-3/1e-4 # (1/ohm/m^2)\n')
    fid.write('    Vl = 0.026e-3 # (V)\n')
    fid.write('    Nao = 114.5 # (mM) = (m mol /lit) = (m mol /dm^3) = (mol/m^3)\n')
    fid.write('    Nai = 13.7 # (mM) = (m mol /lit) = (m mol /dm^3) = (mol/m^3)\n')
    fid.write('    Ko = 2.5 # (mM) = (m mol /lit) = (m mol /dm^3) = (mol/m^3)\n')
    fid.write('    Ki = 120 # (mM) = (m mol /lit) = (m mol /dm^3) = (mol/m^3)\n')
    fid.write('    F = 96514 # (C/mol)\n')
    fid.write('    R = 8.3144 # (J/K/mol)\n')
    fid.write('    T = 295.18 # (K)\n')
    fid.write('    # additional parameters\n')
    fid.write('    H = F/R/T\n\n')

    # write system
    fid.write('    dydt = [\n')
    # set counter
    cnt = 1
    for k in range(node_num):
        # check if the node is nonlinear
        if (k >= first_nl_node) & (k <= last_nl_node):
            if k == 0:   # first node
                fid.write('     1/Cm*(Ga*(       -  y[{}] +  y[{}]             -  ve(t,{}) + ve(t,{}))\n'.format(k,k+1,k,k+1))
                fid.write('                -np.pi*d*l*((F*H*(y[{}]+Vr)/(1-np.exp((y[{}]+Vr)*H)))*((PNa*y[{}]*y[{}]**2 + Pp*y[{}]**2)*(Nao - Nai*np.exp((y[{}]+Vr)*H)) + PK*y[{}]**2*(Ko-Ki*np.exp((y[{}]+Vr)*H)))\n'.format(k, k, iextra[cnt * 4 - 3], iextra[cnt * 4 - 4], iextra[cnt * 4 - 1], k, iextra[cnt * 4 - 2], k))
                fid.write('                + gl*(y[{}]+Vl))),\n'.format(k))

            elif (k == node_num - 1):    # last node node
                fid.write('     1/Cm*(Ga*( y[{}] -  y[{}]           + ve(t,{}) -  ve(t,{}))\n'.format(k-1,k,k-1,k))
                fid.write('                -np.pi*d*l*((F*H*(y[{}]+Vr)/(1-np.exp((y[{}]+Vr)*H)))*((PNa*y[{}]*y[{}]**2 + Pp*y[{}]**2)*(Nao - Nai*np.exp((y[{}]+Vr)*H)) + PK*y[{}]**2*(Ko-Ki*np.exp((y[{}]+Vr)*H)))\n'.format(k, k, iextra[cnt * 4 - 3], iextra[cnt * 4 - 4], iextra[cnt * 4 - 1], k, iextra[cnt * 4 - 2], k))
                fid.write('                + gl*(y[{}]+Vl))),\n'.format(k))

            else:   # inner node
                fid.write('     1/Cm*(Ga*( y[{}] -  2*y[{}] +  y[{}] + ve(t,{}) -  2*ve(t,{}) + ve(t,{}))\n'.format(k-1,k,k+1,k-1,k,k+1))
                fid.write('                -np.pi*d*l*((F*H*(y[{}]+Vr)/(1-np.exp((y[{}]+Vr)*H)))*((PNa*y[{}]*y[{}]**2 + Pp*y[{}]**2)*(Nao - Nai*np.exp((y[{}]+Vr)*H)) + PK*y[{}]**2*(Ko-Ki*np.exp((y[{}]+Vr)*H)))\n'.format(k, k, iextra[cnt * 4 - 3], iextra[cnt * 4 - 4], iextra[cnt * 4 - 1], k, iextra[cnt * 4 - 2], k))
                fid.write('                + gl*(y[{}]+Vl))),\n'.format(k))

            # update counter
            cnt += 1

        else:   # linear node
            if k == 0:   # first node
                fid.write('     1/Cm*(Ga*(       -  y[{}] +  y[{}] +           -  ve(t,{}) + ve(t,{}))   -Gm*y[{}]),\n'.format(k, k + 1, k, k + 1, k))

            elif (k == node_num - 1) & (first_nl_node == -1) & (last_nl_node == -1):   # last node and zero nonlinear node
                fid.write('     1/Cm*(Ga*(y[{}] -  y[{}]           + ve(t,{}) -  ve(t,{}))                -Gm*y[{}])]\n'.format(k-1, k , k-1, k, k))

            elif (k == node_num - 1):   # last node
                fid.write('     1/Cm*(Ga*(y[{}] -  y[{}]           + ve(t,{}) -  ve(t,{}))                -Gm*y[{}]),\n'.format(k-1, k , k-1, k, k))

            else:   # inner node
                fid.write('     1/Cm*(Ga*(y[{}] -  2*y[{}] +  y[{}] + ve(t,{}) -  2*ve(t,{}) + ve(t,{}))  -Gm*y[{}]),\n'.format(k - 1, k, k + 1, k-1 , k, k + 1, k))

    # write additional 4 equations for each nonlinear node
    if (first_nl_node != -1) & (last_nl_node != -1):
        fid.write('     # additional FH equations\n')
        cnt = 1
        for k in range(len(ifh)):
            fid.write('     0.36e6*(y[{}]-22e-3)/(1-np.exp((22e-3-y[{}])/3e-3))*(1-y[{}]) - 0.4e6*(13e-3-y[{}])/(1-np.exp((y[{}]-13e-3)/20e-3))*y[{}],\n'.format(ifh[k],ifh[k],iextra[cnt*4-4],ifh[k],ifh[k],iextra[cnt*4-4]))
            fid.write('     0.1e6*(-10e-3-y[{}])/(1-np.exp((y[{}]+10e-3)/6e-3))*(1-y[{}]) - 4.5e3/(1+np.exp((45e-3-y[{}])/10e-3))*y[{}],\n'.format(ifh[k],ifh[k],iextra[cnt*4-3],ifh[k],iextra[cnt*4-3]))
            fid.write('     0.02e6*(y[{}]-35e-3)/(1-np.exp((35e-3-y[{}])/10e-3))*(1-y[{}]) - 0.05e6*(10e-3-y[{}])/(1-np.exp((y[{}]-10e-3)/10e-3))*y[{}],\n'.format(ifh[k],ifh[k],iextra[cnt*4-2],ifh[k],ifh[k],iextra[cnt*4-2]))

            if cnt == len(ifh):
                fid.write('     0.006e6*(y[{}]-40e-3)/(1-np.exp((40e-3-y[{}])/10e-3))*(1-y[{}]) - 0.09e6*(-25e-3-y[{}])/(1-np.exp((y[{}]+25e-3)/20e-3))*y[{}]]\n'.format(ifh[k],ifh[k],iextra[cnt*4-1],ifh[k],ifh[k],iextra[cnt*4-1]))
            else:
                fid.write('     0.006e6*(y[{}]-40e-3)/(1-np.exp((40e-3-y[{}])/10e-3))*(1-y[{}]) - 0.09e6*(-25e-3-y[{}])/(1-np.exp((y[{}]+25e-3)/20e-3))*y[{}],\n'.format(ifh[k],ifh[k],iextra[cnt*4-1],ifh[k],ifh[k],iextra[cnt*4-1]))
            cnt += 1

    fid.write('    return dydt')
    # close file
    fid.close()


def define_stimulus(kind, **kwargs):
    """
    'define_stimulus' generates the lambda function to estimate the external nodal voltages

    Parameters
    ----------
    kind (str):
        define the stimulus: 'efield', 'monoph', 'biph', 'sine'
        'efield' requires: length, node_num, tp, magnitude
        'monoph' requires: length, node_num, tp, magnitude, rhoe
        'biph'   requires: length, node_num, tp, magnitude, rhoe, tp1
        'sine'   requires: length, node_num, tp, magnitude, rhoe, freq

    The following parameters have to be passed through **kwargs:
        magnitude (float):
            magnitude of the stimulus
        tp (float):
            duration of the stimulus
        node_num (int):
            number of node in the fiber
        length (float):
            length of the fiber
        rhoe (float):
            medium resistivity
        tp1 (float):
            required only for 'biph' stimulus
            tp1 --> high part of the stimulus
            (consequently: tp - tp1 --> low part of the stimulus)
        freq (float):
            frequency of the sine wave

    Returns
    -------
    ve (function):
        lambda function to estimate the external nodal voltages
    inod (ndarray):
        index of node to be plotted
        (the ones where the action potential is originated)
    leg (str):
        legend to be used in the plot of the action potential

    """

    if kind.lower() == 'efield':
        vext_norm = kwargs['length'] * np.arange(kwargs['node_num'] - 1, -1, -1)
        ve = lambda t, k: (t <= kwargs['tp']) * kwargs['magnitude'] * vext_norm[k]
        # plot variables
        inod = np.arange(0, 6, dtype=int)
        leg = ['node #' + str(k) for k in inod]

    else:
        # Stimulation via electrode in central position
        ye = 2e-3
        x = (np.linspace(0, kwargs['node_num'] - 1, kwargs['node_num']) - (kwargs['node_num'] - 1) / 2) * kwargs['length']
        r = np.sqrt(x ** 2 + ye ** 2)

    if kind.lower() == 'monoph':
        # Electrical stimulation
        ve = lambda t, k: (t <= kwargs['tp']) * kwargs['rhoe'] * kwargs['magnitude'] / (4 * np.pi * r[k])
        # plot variables
        inod = np.arange((kwargs['node_num'] - 1) / 2, (kwargs['node_num'] - 1) / 2 + 6, dtype=int)
        leg = ['node #' + str(k) for k in inod]

    elif kind.lower() == 'biph':
        # Electrical stimulation
        ve = lambda t, k: (kwargs['rhoe'] * kwargs['magnitude']) / (4 * np.pi) * ((t <= kwargs['tp1']) / r[k] - ((t > kwargs['tp1']) & (t <= kwargs['tp'])) / r[k])
        # plot variables
        inod = np.arange((kwargs['node_num'] - 1) / 2, (kwargs['node_num'] - 1) / 2 + 6, dtype=int)
        leg = ['node #' + str(k) for k in inod]

    elif kind.lower() == 'sine':
        # Electrical stimulation
        ve = lambda t, k: (t <= kwargs['tp']) * kwargs['rhoe'] * kwargs['magnitude'] * np.sin(2 * np.pi * kwargs['freq'] * t) / (4 * np.pi * r[k])
        # plot variables
        inod = np.arange((kwargs['node_num'] - 1) / 2, (kwargs['node_num'] - 1) / 2 + 6, dtype=int)
        leg = ['node #' + str(k) for k in inod]

    return ve, inod, leg


if __name__ == "__main__":
    from scipy.integrate import ode
    import matplotlib.pylab as plt
    from timeit import default_timer as timer
    plt.ion()

    # FIXED PARAMETERS
    D = 20e-6           # fiber diamater
    d = 0.7*D           # axon diamater
    L = 100 * D           # internodal length
    rhoi = 110*1e-2     # axoplasm resistivity
    rhoe = 300 * 1e-2   # external medium resistivity
    gm = 30.4*10        # membrane conductance per unit length
    l = 2.5e-6          # nodal gap width
    cm = 2e-2           # membrane conductance per unit area
    Vr = -70e-3         # rest potential

    # COMPUTED PARAMETERS
    Ga = np.pi * d**2 / (4 * rhoi * L)
    Gm = gm * np.pi * d * l
    Cm = cm * np.pi * d * l

    # write & import equation
    N = 21       # total number of node
    nlin1 = 0    # first non-linear node
    nlin2 = 20   # last non-linear node
    write_ode(N, nlin1, nlin2)
    from eqdiff import eqdiff

    # callback to save solution at each iteration of the integration
    def solout(t, y):
        time.append(t)
        sol.append(y.copy())

    # define stimulus
    kind = 'sine'
    if kind == 'monoph':
        tp = 100e-6
        NTp = 10
        param = {'length':L, 'node_num':N, 'tp':tp, 'magnitude':0, 'rhoe':rhoe}
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'biph':
        tp = 100e-6
        tp1 = 50e-6
        NTp = 10
        param = {'length': L, 'node_num': N, 'tp': tp, 'magnitude': 0, 'rhoe': rhoe, 'tp1':tp1}
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'efield':
        tp = 100e-6
        NTp = 10
        param = {'length': L, 'node_num': N, 'tp': tp, 'magnitude': 0}
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'sine':
        tp = 20e-6
        NTp = 25
        param = {'length': L, 'node_num': N, 'tp': tp, 'magnitude': 0, 'rhoe':rhoe, 'freq':(1/tp)}
        ktime = 1e6
        umes_time = "($\mu$s)"



    # define simulation time range and initial condition
    icond = [0] * N + [0.0005, 0.8249, 0.0268, 0.0049] * (nlin2 - nlin1 + 1)

    # Identification of bound values
    Isub = 0
    Isup = 0.1e-3
    flag = 1

    print('\n------------------------------------------------------')
    print('Identifying sub and supra threshold values...')
    print('------------------------------------------------------')
    ts = timer()
    while flag:
        # define stumulus
        param['magnitude'] = -Isup
        ve, inod, leg = define_stimulus(kind, **param)

        # initialize solution variable
        time = []
        sol = []
        # define integrator
        r = ode(eqdiff).set_integrator('dopri5')
        # set initial conditions
        r.set_initial_value(icond, 0).set_f_params(Ga, Gm, Cm, ve, d, l, Vr)
        # store solution at each iteration step
        r.set_solout(solout)
        # integrate
        r.integrate(NTp * tp)
        # get complete solution
        t = np.array(time)
        x = np.array(sol)

        # get number of nodes with voltage > 80 mV
        N80 = (np.max(x[:, 0:N], axis=0) > 80e-3).sum()
        if N80 > 3:
            Isim = (Isub + Isup) / 2
            flag = 0
        else:
            Isub = Isup
            Isup = 2 * Isup
    te = timer()
    print('...done. (sub, sup) = ({},{})'.format(Isub, Isup))
    print('\n elapsed time: {:3f} ms'.format(te-ts))


    # Identification of the thresholds
    toll = 0.5
    err = 100
    Iold = -Isim
    Inew = -Isim
    cnt = 1
    print('\n------------------------------------------------------')
    print('Looking for threshold...')
    print('------------------------------------------------------')
    ts = timer()
    while (err > toll):
        # define stumulus
        param['magnitude'] = Inew
        ve, inod, leg = define_stimulus(kind, **param)

        # initialize solution variable
        time = []
        sol = []
        # define integrator
        r = ode(eqdiff).set_integrator('dopri5')
        # set initial conditions
        r.set_initial_value(icond, 0).set_f_params(Ga, Gm, Cm, ve, d, l, Vr)
        # store solution at each iteration step
        r.set_solout(solout)
        # integrate
        r.integrate(NTp * tp)
        # get complete solution
        t = np.array(time)
        x = np.array(sol)

        # get number of nodes with voltage > 80 mV
        N80 = (np.max(x[:, 0:N], axis=0) > 80e-3).sum()
        if N80 > 3:
            Isup = -Inew
        else:
            Isub = -Inew

        Iold = -Inew
        Inew = -(Isub + Isup) / 2
        err = np.abs(Inew + Iold) / Iold * 100
        print('{}) I = {:.6f} A \t Err = {} % \t N80 = {}'.format(cnt, Iold, err, N80))
        cnt += 1
    te = timer()
    print('\n elapsed time: {:3f} ms'.format(te - ts))

    # one shot simulation
    # -------------------
    # define stumulus
    param['magnitude'] = Inew
    ve, inod, leg = define_stimulus(kind, **param)

    # initialize solution variable
    time = []
    sol = []
    # define integrator
    r = ode(eqdiff).set_integrator('dopri5')
    # set initial conditions
    r.set_initial_value(icond, 0).set_f_params(Ga, Gm, Cm, ve, d, l, Vr)
    # store solution at each iteration step
    r.set_solout(solout)
    # integrate
    r.integrate(NTp * tp)
    # get complete solution
    t = np.array(time)
    sol = np.array(sol)

    # print max. values
    vmax = sol[:, 0:N].max()
    i, j = np.unravel_index(np.argmax(sol[:, 0:N]), sol[:, 0:N].shape)
    tmax = t[i]
    print('\n------------------------------------------------------')
    print('Stimulus: ' + kind)
    print('Strength: {}'.format(Inew))
    print('Duration: {}'.format(tp * ktime) + umes_time.replace('$', '').replace('\\', '').replace('mu', 'u'))
    print('max. membrane voltage: {:.3f} mV'.format(vmax*1e3))
    print('max. reached at: {:.3f} '.format(tmax * ktime) + umes_time.replace('$','').replace('\\','').replace('mu','u'))
    print('max. stimulation at node number: {} (0 = first node)'.format(j))

    # action potential check
    N80 = (np.max(sol[:, 0:N], axis=0) > 80e-3).sum()
    if N80 > 3:
        print('action potential started: {} nodes with V > 80 mv'.format(N80))
    else:
        print('action potential not started: {} nodes with V > 80 mV'.format(N80))

    # plot stimulus
    plt.figure(facecolor='w')
    plt.plot(t * ktime, -ve(t, inod[0]) * 1e3, linewidth=2)
    plt.xlabel('time ' + umes_time,fontsize=16)
    plt.ylabel('external potential (mV)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()

    # plot nodes
    plt.figure(facecolor='w')
    plt.plot(t * ktime, sol[:, inod] * 1e3, linewidth=2)
    plt.xlabel('time ' + umes_time, fontsize=16)
    plt.ylabel('external potential (mV)', fontsize=16)
    plt.legend(leg,loc='best')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()