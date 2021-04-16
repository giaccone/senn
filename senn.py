# import modules
import numpy as np
from timeit import default_timer as timer
from scipy.integrate import ode


# classes
class AxonModel:
    """
    Class that defines the Axon Model according to McNeal 1976 and Reilly 1985.

    Attributes
    ----------
    D (float): fiber diameter
    rhoi (float): axoplasm resistivity
    rhoe (float): external medium resistivity
    gm (float): membrane conductance per unit length
    l (float): nodal gap width
    cm (float): membrane capacitance per unit area
    Vr (float): rest potential
    node_num (int): total number of node
    inl1 (int): first non-linear node
    inl2 (int): second non-linear node
    icond (list): initial conditions for ODE
    d (float): axon diameter
    L (float): internodal length
    Ga (float): internonal conductance
    Gm (float): transmembrane conductance
    Cm  (float): transmembrane capacitance
    """

    def __repr__(self):
        return "AxonModel: {} nodes, {} nonlinear nodes".format(self.node_num, self.inl2 - self.inl1 + 1)

    def __str__(self):

        umes = {'D':' (m)',
                'rhoi':'(ohm m)',
                'rhoe':'(ohm m)',
                'gm': '(S / m**2)',
                'l':'(m)',
                'cm':'(F / m**2)',
                'Vr':'(V)',
                'node_num':'',
                'inl1':'',
                'inl2':'',
                'icond':'',
                'd':'(m)',
                'L':'(m)',
                'Ga':'(S)',
                'Gm':'(S)',
                'Cm':'(F)'}


        msg = "AxonModel:\n"
        for attr, value in self.__dict__.items():
            if attr == 'node_num':
                msg += "\t{}\t{} {}\n".format(attr, value, umes[attr])
            elif attr == 'icond':
                msg += "\t{}\t\t[{} {} ... {} {}] {}\n".format(attr, value[0], value[1], value[-2], value[-1], umes[attr])
            else:
                msg += "\t{}\t\t{} {}\n".format(attr, value, umes[attr])
        return msg

    def __init__(self, D=20e-6, rhoi=110e-2, rhoe=300e-2, gm=304, l=2.5e-6, cm=2e-2, Vr=-70e-3, node_num=21, inl1=0, inl2=20):
        """

        Parameters
        ----------
        D (float): fiber diameter
        rhoi (float): axoplasm resistivity
        rhoe (float): external medium resistivity
        gm (float): membrane conductance per unit area
        l (float): nodal gap width
        cm (float): membrane conductance per unit area
        Vr (float): rest potential
        node_num (int): total number of node
        inl1 (int): first non-linear node
        inl2 (int): second non-linear node
        """

        # given parameters
        self.D = D          # fiber diameter
        self.rhoi = rhoi    # axoplasm resistivity
        self.rhoe = rhoe    # external medium resistivity
        self.gm = gm        # membrane conductance per unit area
        self.l = l          # nodal gap width
        self.cm = cm        # mebrane conductance per unit area
        self.Vr = Vr        # rest potential

        self.node_num = node_num    # total number of node
        self.inl1 = inl1            # first non-linear node
        self.inl2 = inl2            # second non-linear node

        # computed parameters
        self.icond = [0] * node_num + [0.0005, 0.8249, 0.0268, 0.0049] * (inl2 - inl1 + 1)  # initial condition

        self.d = 0.7 * self.D   # axon diameter
        self.L = 100 * self.D   # internodal length

        self.Ga = np.pi * self.d ** 2 / (4 * self.rhoi * self.L)
        self.Gm = self.gm * np.pi * self.d * self.l
        self.Cm = self.cm * np.pi * self.d * self.l


class StimulusModel:
    """
        Class that defines the Stimulus.

        Attributes
        ----------
        kind (str): kind of stimulus ('efield', 'monoph', 'biph', 'sine')
        magnitude (float): magnitude of the stimulus
        tp (float): duration of the strimulus
        tp1 (float): additional duration (for 'biph')
            tp1 --> high part of the stimulus
            (consequently: tp - tp1 --> low part of the stimulus)
        freq (float): frequency (for 'sine')
        waveform (function): lambda function that generates the stimulus over time
        voltage_ext (function): lambda function that generates the exitation voltage over time at a given point
        ye (ndarray): depth of the electrode
        x (ndarray): location of Ranvier's nodes
        r (ndarray): distance of the Ranvier's nodes fron the electrode
        inod (list): index of node to be plotted (the ones where the action potential is originated)
        leg (list): legend to be used in the plot of the action potential
        """

    def __repr__(self):
        return "StimulusModel: {}".format(self.kind)

    def __str__(self):
        msg = "StimulusModel:\n"
        msg += "\tkind:\t\t{}\n".format(self.kind)
        msg += "\tmagnitude:\t{}\n".format(self.magnitude)
        msg += "\tpulse duration:\t{} (s)\n".format(self.tp)
        msg += "\ttotal duration:\t{} (s)\n".format(self.tend)
        if self.kind.lower() == 'biph':
            msg += "\tother duration\t{} (s)\n".format(self.tp1)
        if self.kind.lower() == 'sine':
            msg += "\tfrequency:\t{} (Hz)\n".format(self.freq)

        return msg

    def __init__(self, axon, kind, magnitude, tp, tend=None, tp1=None, freq=None):
        """

        Parameters
        ----------
        axon (AxonModel): axon model
        kind (str): kind of stimulus
        magnitude (float): magnitude of the stimulus
        tp (float): pulse duration
        tend (float): end of the stimulus (it corresponds to the total simulation time)
        tp1 (float): additional duration (for 'biph')
            tp1 --> high part of the stimulus
            (consequently: tp - tp1 --> low part of the stimulus)
        freq (float): frequency (for 'sine')
        """
        # given attributes
        self.kind = kind
        self.magnitude = magnitude
        self.tp = tp
        if tend is None:
            self.tend = 10 * self.tp
        else:
            self.tend = tend
        self.tp1 = tp1
        self.freq = freq

        # initialization of attributes defined later with 'build_voltage_ext' and 'build_waveform'
        self.voltage_ext = None
        self.waveform = None
        self.ye = None
        self.x = None
        self.r = None
        self.inod = None
        self.leg = None

        # update voltage_ext, waveform, ye, x, r, inod, leg
        self.update_stimulus(axon)

    def update_stimulus(self, axon):
        """
        'update_stimulus' generates the lambda function to define excitation voltage over time at a given point
        and the lambda function to define waveform over time

        Parameters
        ----------
        axon (AxonModel): axon model

        """
        if self.kind.lower() == 'efield':
            normalized_waveform = axon.L * np.arange(axon.node_num - 1, -1, -1)
            self.voltage_ext = lambda t, k: (t <= self.tp) * self.magnitude * normalized_waveform[k]
            self.waveform = lambda t: (t <= self.tp) * self.magnitude
            self.inod = np.arange(0, 6, dtype=int)
        else:
            self.ye = 2e-3
            self.x = (np.linspace(0, axon.node_num - 1, axon.node_num) - (axon.node_num - 1) / 2) * axon.L
            self.r = np.sqrt(self.x ** 2 + self.ye ** 2)

        if self.kind.lower() == 'monoph':
            # Electrical stimulation
            self.voltage_ext = lambda t, k: (t <= self.tp) * axon.rhoe * self.magnitude / (4 * np.pi * self.r[k])
            self.waveform = lambda t: (t <= self.tp) * self.magnitude
            # plot variables
            self.inod = np.arange((axon.node_num - 1) / 2, (axon.node_num - 1) / 2 + 6, dtype=int)

        elif self.kind.lower() == 'biph':
            # Electrical stimulation
            self.voltage_ext = lambda t, k: (axon.rhoe * self.magnitude) / (4 * np.pi) * ((t <= self.tp1) / self.r[k] - ((t > self.tp1) & (t <= self.tp)) / self.r[k])
            self.waveform = lambda t: self.magnitude * ((t <= self.tp1) - ((t > self.tp1) & (t <= self.tp)))
            # plot variables
            self.inod = np.arange((axon.node_num - 1) / 2, (axon.node_num - 1) / 2 + 6, dtype=int)

        elif self.kind.lower() == 'sine':
            # Electrical stimulation
            self.voltage_ext = lambda t, k: (t <= self.tp) * axon.rhoe * self.magnitude * np.sin(2 * np.pi * self.freq * t) / (4 * np.pi * self.r[k])
            self.waveform = lambda t: (t <= self.tp) * self.magnitude * np.sin(2 * np.pi * self.freq * t)
            # plot variables
            self.inod = np.arange((axon.node_num - 1) / 2, (axon.node_num - 1) / 2 + 6, dtype=int)

        self.leg = ['node #' + str(k) for k in self.inod]


def write_ode(node_num, first_nl_node=-1, last_nl_node=-2):
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

            elif (k == node_num - 1) & (first_nl_node == -1) & (last_nl_node == -2):   # last node and zero nonlinear node
                fid.write('     1/Cm*(Ga*(y[{}] -  y[{}]           + ve(t,{}) -  ve(t,{}))                -Gm*y[{}])]\n'.format(k-1, k , k-1, k, k))

            elif (k == node_num - 1):   # last node
                fid.write('     1/Cm*(Ga*(y[{}] -  y[{}]           + ve(t,{}) -  ve(t,{}))                -Gm*y[{}]),\n'.format(k-1, k , k-1, k, k))

            else:   # inner node
                fid.write('     1/Cm*(Ga*(y[{}] -  2*y[{}] +  y[{}] + ve(t,{}) -  2*ve(t,{}) + ve(t,{}))  -Gm*y[{}]),\n'.format(k - 1, k, k + 1, k-1 , k, k + 1, k))

    # write additional 4 equations for each nonlinear node
    if (first_nl_node != -1) & (last_nl_node != -2):
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


def find_sub_supra(axon, stimulus, eqdiff, sub_value=0, sup_value=0.1e-3):
    """
    'find_sub_supra' computes boundary values for the bisection method (used to identify the threeshold)

    Parameters
    ----------
    axon (AxonModel): axon model
    stimulus (StimulusModel): stimulus model
    eqdiff (function): function that defines the ODE system
    sub_value (float): initial guess of sub-threshold value (default is 0)
    sup_value (float): initial guess of supra-threshold value (default is 0.1e-3)

    Returns
    -------
    sub_value (float): sub-threshold value
    sup_value (float): supra-threshold value

    """

    # Identification of bound values
    flag = 1

    print('\n------------------------------------------------------')
    print('Identifying sub and supra threshold values...')
    print('------------------------------------------------------')
    ts = timer()
    while flag:
        # update stimulus
        stimulus.magnitude = -sup_value
        stimulus.update_stimulus(axon)

        # callback to save solution at each iteration of the integration
        def solout(t, y):
            time.append(t)
            sol.append(y.copy())

        # initialize solution variable
        time = []
        sol = []
        # define integrator
        r = ode(eqdiff).set_integrator('dopri5')
        # set initial conditions
        r.set_initial_value(axon.icond, 0).set_f_params(axon.Ga, axon.Gm, axon.Cm, stimulus.voltage_ext, axon.d, axon.l, axon.Vr)
        # store solution at each iteration step
        r.set_solout(solout)
        # integrate
        r.integrate(stimulus.tend)
        # get complete solution
        x = np.array(sol)

        # get number of nodes with voltage > 80 mV
        N80 = (np.max(x[:, 0:axon.node_num], axis=0) > 80e-3).sum()
        if N80 > 3:
            flag = 0
        else:
            sub_value = 1*sup_value
            sup_value = 2 * sup_value
    te = timer()
    print('...done. (sub, sup) = ({},{})'.format(sub_value, sup_value))
    print('\n    elapsed time: {:3f} ms'.format(te - ts))
    return sub_value, sup_value


def find_threshold(axon, stimulus, eqdiff, sub_value, sup_value, toll=0.5):
    """
    'find_threshold' computes the threshold up to a given tolerance.

    Parameters
    ----------
    axon (AxonModel): axon model
    stimulus (StimulusModel): stimulus model
    eqdiff (function): function that defines the ODE system
    sub_value (float): sub-threshold value
    sup_value (float): supra-threshold value
    toll (float): tolerance to identify the threshold (default is 0.5 %)

    Returns
    -------
    threshold (float): threshold value up to a given tolerance 'toll'
    stimulus (StimulusModel): stimulus model with magnitude = threshold
    """
    # Identification of the thresholds
    err = 100
    N80 = 0
    sim_value = 0.5 * (sub_value + sup_value)
    old_value = -sim_value
    new_value = -sim_value
    cnt = 1

    # callback to save solution at each iteration of the integration
    def solout(t, y):
        time.append(t)
        sol.append(y.copy())

    print('\n------------------------------------------------------')
    print('Looking for threshold...')
    print('------------------------------------------------------')
    ts = timer()
    while ((err > toll) or (N80 < 3)):
        # update stimulus
        stimulus.magnitude = new_value
        stimulus.update_stimulus(axon)

        # initialize solution variable
        time = []
        sol = []
        # define integrator
        r = ode(eqdiff).set_integrator('dopri5')
        # set initial conditions
        r.set_initial_value(axon.icond, 0).set_f_params(axon.Ga, axon.Gm, axon.Cm, stimulus.voltage_ext, axon.d, axon.l, axon.Vr)
        # store solution at each iteration step
        r.set_solout(solout)
        # integrate
        r.integrate(stimulus.tend)
        # get complete solution
        x = np.array(sol)

        # get number of nodes with voltage > 80 mV
        N80 = (np.max(x[:, 0:axon.node_num], axis=0) > 80e-3).sum()
        if N80 > 3:
            sup_value = -new_value
        else:
            sub_value = -new_value

        old_value = -new_value
        new_value = -(sub_value + sup_value) / 2
        err = np.abs(new_value + old_value) / old_value * 100
        print('{}) I = {:.6f} A \t Err = {} % \t N80 = {}'.format(cnt, old_value, err, N80))
        cnt += 1

    te = timer()
    print('\n    elapsed time: {:3f} ms'.format(te - ts))
    print('    threshold: {}  (tolerance {} %%)'.format(old_value, toll))



    # update stimulus
    stimulus.magnitude = -old_value
    stimulus.update_stimulus(axon)
    return old_value, stimulus


def evaluate_senn(axon, stimulus, eqdiff):
    """
    'evaluate_senn' computes the solution of the SENN model for a given stimulus

    Parameters
    ----------
    axon (AxonModel): axon model
    stimulus (StimulusModel): stimulus model
    eqdiff (function): function that defines the ODE system

    Returns
    -------
    t (ndarray): time values at which the solution is defined
    sol (ndarray):  solution
    """

    # callback to save solution at each iteration of the integration
    def solout(t, y):
        time.append(t)
        sol.append(y.copy())

    # initialize solution variable
    time = []
    sol = []
    # define integrator
    r = ode(eqdiff).set_integrator('dopri5')
    # set initial conditions
    r.set_initial_value(axon.icond, 0).set_f_params(axon.Ga, axon.Gm, axon.Cm, stimulus.voltage_ext, axon.d, axon.l, axon.Vr)
    # store solution at each iteration step
    r.set_solout(solout)
    # integrate
    r.integrate(stimulus.tend)
    # get complete solution
    t = np.array(time)
    sol = np.array(sol)

    return t, sol


def analyze_solution(time, solution, axon, stimulus):
    """
    'analyze_solution' print a synthetic summary of the output

    Parameters
    ----------
    time (ndarray): simulation time
    solution (ndarray): solution
    axon (AxonModel): axon model
    stimulus (StimulusModel): stimulus model

    """
    # print max. values
    vmax = solution[:, 0:axon.node_num].max()
    i, j = np.unravel_index(np.argmax(solution[:, 0:axon.node_num]), solution[:, 0:axon.node_num].shape)
    tmax = time[i]
    print('\n------------------------------------------------------')
    print('Stimulus: ' + stimulus.kind)
    print('Magnitude: {:g}'.format(stimulus.magnitude))
    print('Duration: {:g} (s)'.format(stimulus.tp))
    print('max. membrane voltage: {:g} (V)'.format(vmax))
    print('max. reached at: {:g} (s)'.format(tmax))
    print('max. stimulation at node number: {:d} (0 = first node)'.format(j))

    # action potential check
    N80 = (np.max(solution[:, 0:axon.node_num], axis=0) > 80e-3).sum()
    if N80 >= 3:
        print('action potential started: {:d} nodes with V > 80 mv'.format(N80))
    else:
        print('action potential not started: {:d} nodes with V > 80 mV'.format(N80))


if __name__ == "__main__":
    from scipy.integrate import ode
    import matplotlib.pylab as plt
    from timeit import default_timer as timer
    plt.ion()

    # Define axon model
    D = 20e-6           # fiber diamater
    d = 0.7*D           # axon diamater
    L = 100 * D           # internodal length
    rhoi = 110*1e-2     # axoplasm resistivity
    rhoe = 300 * 1e-2   # external medium resistivity
    gm = 30.4*10        # membrane conductance per unit length
    l = 2.5e-6          # nodal gap width
    cm = 2e-2           # membrane conductance per unit area
    Vr = -70e-3         # rest potential
    node_num = 21  # total number of node
    inl1 = 0  # first non-linear node
    inl2 = 20  # last non-linear node

    axon = AxonModel(D=D, rhoi=rhoi, rhoe=rhoe, gm=gm, l=l, cm=cm, Vr=Vr, node_num=node_num, inl1=inl1, inl2=inl2)

    # write & import equation
    write_ode(axon.node_num, axon.inl1, axon.inl2)
    from eqdiff import eqdiff

    # define stimulus
    kind = 'sine'
    if kind == 'monoph':
        tp = 100e-6
        tend = 10 * tp
        magnitude = 0
        stimulus = StimulusModel(axon, kind, magnitude, tp, tend)
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'biph':
        tp = 100e-6
        tp1 = 50e-6
        tend = 10 * tp
        magnitude = 0
        stimulus = StimulusModel(axon, kind, magnitude, tp, tend, tp1=tp1)
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'efield':
        tp = 100e-6
        tend = 10 * tp
        magnitude = 0
        stimulus = StimulusModel(axon, kind, magnitude, tp)
        ktime = 1e6
        umes_time = "($\mu$s)"
    elif kind == 'sine':
        magnitude = 0
        tp = 20e-6
        freq = 1 / tp
        tend = 25 * tp
        stimulus = StimulusModel(axon, kind, magnitude, tp, tend, freq=freq)
        ktime = 1e6
        umes_time = "($\mu$s)"

    # Identification of boundaries for bisection method
    Isub, Isup = find_sub_supra(axon, stimulus, eqdiff)

    # Identification of the thresholds
    It, stimulus = find_threshold(axon, stimulus, eqdiff, Isub, Isup, toll=0.5)

    # one shot simulation
    t, sol = evaluate_senn(axon, stimulus, eqdiff)

    # analyze output
    analyze_solution(t, sol, axon, stimulus)

    # plot stimulus
    plt.figure(facecolor='w')
    plt.plot(t * ktime, -stimulus.waveform(t) * 1e3, linewidth=2)
    plt.xlabel('time ' + umes_time, fontsize=16)
    plt.ylabel('stimulus (mA)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()

    # plot stimulus
    plt.figure(facecolor='w')
    plt.plot(t * ktime, -stimulus.voltage_ext(t, stimulus.inod[0]) * 1e3, linewidth=2)
    plt.xlabel('time ' + umes_time,fontsize=16)
    plt.ylabel('external potential (mV)', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()

    # plot nodes
    plt.figure(facecolor='w')
    plt.plot(t * ktime, sol[:, stimulus.inod] * 1e3, linewidth=2)
    plt.xlabel('time ' + umes_time, fontsize=16)
    plt.ylabel('external potential (mV)', fontsize=16)
    plt.legend(stimulus.leg,loc='best')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
