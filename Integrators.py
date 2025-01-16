# Contains numerical ODE integrators

import numpy as np

# ==============================================================================

class RK45Integrator:
  """General embedded Runge-Kutta 4th/5h-order solver with adaptive stepsize control"""

  def __init__(self, method="Fehlberg"):
    self.solver_initialized = False
    self.solver = None
    self.nstages = None
    self.init_solver(method)

  # ----------------------------------------------------------------------------

  def init_solver(self, method):
    # Note: the C coefficients are for the 5th-order solution, while
    # the D coeficients for the embedded 4th-order one
    if method.lower() in ["fehlberg", "fehl", "fb", "rkf"]:
      self.solver = "Fehlberg"
      self.nstages = 6
      A = [1/4, 3/8, 12/13, 1, 1/2]
      B2 = [1/4]
      B3 = [3/32, 9/32]
      B4 = [1932/2197, -7200/2197, 7296/2197]
      B5 = [439/216, -8, 3680/513,-845/4104]
      B6 = [-8/27, 2, -3544/2565, 1859/4104, -11/40]
      C = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55]
      D = [25/216, 0, 1408/2565, 2197/4101, -1/5, 0]
    elif method.lower() in ["cashkarp", "ck", "rkck"]:
      self.solver = "Cash-Karp"
      self.nstages = 6
      A = [1/5, 3/10, 3/5, 1, 7/8]
      B2 = [1/5]
      B3 = [3/40, 9/40]
      B4 = [3/10, -9/10, 6/5]
      B5 = [-11/54, 5/2, -70/27, 35/27]
      B6 = [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096]
      C = [37/378, 0, 250/621, 125/594, 0, 512/1771]
      D = [825/27648, 0, 18575/48384, 13525/55296, 277/1433, 1/4]
    elif method.lower() in ["dormandprince", "dormand-price", "dp", "rkdp"]:
      self.solver = "Dormand-Prince"
      self.nstages = 7
      A = [1/5, 3/10, 4/5, 8/9, 1, 1]
      B2 = [1/5]
      B3 = [3/40, 9/40]
      B4 = [44/45, -56/15, 32/9]
      B5 = [19372/6561, -25360/2187, 64448/6561, -212/729]
      B6 = [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656]
      B7 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84]  #not a typo XD
      C = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
      D = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]
    else:
      print(f"Solver {method} not recognized!")
      sys.exit()
    self.a2, self.a3, self.a4, self.a5, self.a6 = A
    self.b21 = B2
    self.b31, self.b32 = B3
    self.b41, self.b42, self.b43 = B4
    self.b51, self.b52, self.b53, self.b54 = B5
    self.b61, self.b62, self.b63, self.b64, self.b65 = B6
    self.c1, self.c2, self.c3, self.c4, self.c5, self.c6 = C
    self.d1, self.d2, self.d3, self.d4, self.d5, self.d6 = D
    self.solver_initialized = True

  # ----------------------------------------------------------------------------

  def solve_IVP(self, f, y0, t0, tf, h0=1.0, dtout="all", max_step=None, method="Fehlberg", atol=1e-6, rtol=1e-3, ret_steps=False, stop=None, verbose=False):
    """Solves an IVP using a RK45 integrator

    This function receives an IVP specified by:
      f: derivatives of the dependent variables, a function f(t, y), where t
         is the independent variable and y is the vector of dependent variables,
         that returns a vector of derivatives
      y0: vector of start values of the dependent variables, y = (y1, y2, ...)
      t0, tf: start and end values of the independent variable
      h0: initial time step to use (will be adjusted)

    The supported options are:
      dtout: the interval of the dependent variable at which values should be
        saved and returned at the end of integration. Specify 'all' to return
        all computed values, with variable stepsize (this is the default)
      method: specifies the coefficients to use for the RK45 solver, must be
        one of "Fehlberg" (default), "CashKarp"
      atol: absolute tolerance of the error (default is 1e-6)
      rtol: relative tolerance of the error (default is 1e-3)
      ret_steps: also return the stepsizes
      stop: stopping criterion, a function f(t, y) called after every step
        that takes the current simulation state and returns whether to stop

    Returns:
      (ts, ys) if ret_steps is False (the default)
        or
      (ts, ys, hs) otherwise
    where:
      ts: list of values of the indep var at which the solution is evaluated
      ys: list of corresponding values of the dependent variables
      hs: list of corresponding actual setpsizes
    """
    self.init_solver(method)
    t = t0;
    y = y0
    ts = [t0]; ys = [y0]; hs = [np.nan]
    if dtout != "all": nextout = t0 + dtout
    h = h0
    next_rep = (tf-t0)/20
    while t < tf:

      hmax = np.inf
      if dtout != "all":
        hmax = min(tf - t, dtout)
      if max_step is not None:
        hmax = min(hmax, max_step)

      tnew, ynew, hnext = self.step(t, y, f, h, atol, rtol, hmax)

      hlast = tnew - t
      t = tnew; y = ynew; h = hnext

      if dtout == "all" or t >= nextout:
        ts.append(t); ys.append(y)
        if ret_steps: hs.append(hlast)
        if dtout != "all": nextout += dtout
      # print(f"{t:.4f}", f"{hlast:.3e}", y); input()

      if verbose and t >= next_rep:
        print(f"{t:.4f}", f"{hlast:.3e}", y)
        next_rep += (tf-t0)/20

      if stop is not None:
        if stop(t, y): break

    ts = np.array(ts)
    ys = np.array(ys)
    hs = np.array(hs)

    if ret_steps: return ts, ys, hs
    else: return ts, ys

  # ----------------------------------------------------------------------------

  def step(self, t, y, f, h, atol=1e-6, rtol=1e-3, hmax=np.inf):
    """Does a single explicit Runge-Kutta step with adaptive stepsize control

    Options are the same as solve_IVP above, except:
      hmax: maximum step size (defaults to infinity, so no maximum)

    Returns:
      (tf, yf, hf)
    where:
      tf: stepped independent variable
      yf: vector of stepped dependent variables
      hf: next step size
    """
    if not self.solver_initialized:
      self.init_solver()
    h = min(h, hmax)
    while True:

      # print("\ntrying h=", h, ", hmax=", hmax)

      s = self
      k1 = h * f(t, y)
      k2 = h * f(t + s.a2*h, y + s.b21*k1)
      k3 = h * f(t + s.a3*h, y + s.b31*k1 + s.b32*k2)
      k4 = h * f(t + s.a4*h, y + s.b41*k1 + s.b42*k2 + s.b43*k3)
      k5 = h * f(t + s.a5*h, y + s.b51*k1 + s.b52*k2 + s.b53*k3 + s.b54*k4)
      k6 = h * f(t + s.a6*h, y + s.b61*k1 + s.b62*k2 + s.b63*k3 + s.b64*k4 + s.b65*k5)
      y5 = y + s.c1*k1 + s.c2*k2 + s.c3*k3 + s.c4*k4 + s.c5*k5 + s.c6*k6
      y4 = y + s.d1*k1 + s.d2*k2 + s.d3*k3 + s.d4*k4 + s.d5*k5 + s.d6*k6

      # Error estimation, scaled per equation
      yerr = y5 - y4
      yscal = atol + np.maximum(np.abs(y), np.abs(y5)) * rtol
      errs = np.abs(yerr/yscal)
      errmax = np.max(errs)
      # print("yerr=", yerr)
      # print("yscal=", yscal)
      # print("errs=", errs)
      # print("errmax=", errmax)

      if errmax <= 1.0:
        # print("accepted"); input()
        success = True
        tnew = t + h
      else:
        # print("rejected")
        success = False
        # Decrease stepsize
        hnew = 0.9*h*errmax**(-0.25)
        # Decrease no more than a factor of 10
        h = max(hnew, 0.1*h) if h >= 0 else min(hnew, 0.1*h)
        # print("new h=", h); input()
        # h = min(h, hmax)
        # If proposed h is ~zero, raise exception
        if (t + h == t):
          raise(Exception("stepsize underflow in RKF45"))

      if success: break

    # Ajust stepsize for next iteration; increase no more than a factor of 5
    hnext = 0.9*h*errmax**(-0.2)
    hnext = min(hnext, 5*h) if h >= 0 else max(hnext, 5*h)

    return tnew, y5, hnext

# ==============================================================================

# Test unit: Harmonic oscillator
def harmonic_oscillator():

  from math import sqrt, pi, sin, cos
  import matplotlib.pyplot as plt
  import matplotlib.ticker as mticker

  # Harmonic oscillator
  k = 1.0
  m = 1.0
  w = sqrt(k/m)
  x0 = 1.0
  v0 = 0
  T = 2*pi/w
  A = sqrt(x0**2 + (v0/w)**2)
  tfin = 10*T

  def derivs(t, yvars):
    x, v = yvars
    return np.array([v, -k*x])

  y0 = np.array([x0, v0])
  RK45 = RK45Integrator("Fehlberg")
  ts1, ys1, hs1 = RK45.solve_IVP(y0, derivs, 0, tfin, dtout=0.1, tol=1e-6, steps=True)

  xs1 = ys1[:, 0]; vs1 = ys1[:, 1]
  Es1 = 0.5*m*vs1**2 + 0.5*k*xs1**2

  tse = np.linspace(0, tfin, 500)
  xse = A*np.cos(w*tse)
  vse = -A*w*np.sin(w*tse)
  Ese = 0.5*m*vse**2 + 0.5*k*xse**2
  E0 = 0.5*m*v0**2 + 0.5*k*x0**2

  E_err1 = (Es1 - E0) / E0

  plt.figure(figsize=(12, 8))

  plt.suptitle("Runge-Kutta-Fehlberg with adaptive stepsize: simple harmonic oscillator (m=1, k=1, x0=1, v0=0)")

  plt.subplot(4, 1, 1)
  plt.plot(tse, xse, "-", color="C0", label="Exact")
  plt.plot(ts1, xs1, "o", mfc="none", color="C0", label="RKF45")
  plt.title("")
  plt.grid(ls=":")
  plt.legend(loc=3, title="Position")

  plt.subplot(4, 1, 2)
  plt.plot(tse, vse, "-", color="C1", label="Exact")
  plt.plot(ts1, vs1, "o", mfc="none", color="C1", label="RKF45")
  plt.grid(ls=":")
  plt.legend(loc=2, title="Speed")

  plt.subplot(4, 1, 3)
  plt.plot(ts1, E_err1, color="C3", label="Energy error (relative)")
  plt.grid(ls=":")
  plt.legend(loc="upper center")
  plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%0.0e'))

  plt.subplot(4, 1, 4)
  plt.plot(ts1, hs1, color="C4", label="Stepsize")
  plt.grid(ls=":")
  plt.legend(loc="upper center")
  plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%f'))


  plt.tight_layout()
  plt.subplots_adjust(hspace=0)

  plt.show()

# ------------------------------------------------------------------------------

# Test unit: Halley's comet
def Halleys_comet():

  from math import sqrt, pi, sin, cos
  import matplotlib.pyplot as plt
  import matplotlib.patches as mpatches

  AU = 1.495978e11
  G = 6.67430e-11
  MU = 1.3271244e20
  YR = 86400*365.2425

  rap = 35.082 * AU
  rpe = 0.586 * AU
  a = (rpe+rap)/2
  x0 = 0
  y0 = rap
  vx0 = sqrt(MU*(2/rap-1/a))
  vy0 = 0
  vars0 = np.array([x0, y0, vx0, vy0])
  tfin = 76 * YR
  atol = np.array([1000, 1000, 10, 10])
  rtol = 1e-4

  E0 = 0.5*(vx0**2+vy0**2) - MU/np.sqrt(x0**2 + y0**2)

  def derivs(t, vars):
    x, y, vx, vy = vars
    r3 = (x**2 + y**2)**(3/2)
    ax = -MU/r3 * x
    ay = -MU/r3 * y
    return np.array([vx, vy, ax, ay])

  RK45 = RK45Integrator()
  ts, vars, hs = RK45.solve_IVP(derivs, vars0, 0, tfin, dtout=86400, atol=atol, rtol=rtol, ret_steps=True)

  xs = vars[:, 0]
  ys = vars[:, 1]
  vxs = vars[:, 2]
  vys = vars[:, 3]

  Es = 0.5*(vxs**2+vys**2) - MU/np.sqrt(xs**2 + ys**2)
  E_err = (Es - E0)/E0

  from scipy.integrate import solve_ivp

  result = solve_ivp(derivs, (0, tfin), vars0, max_step=0.01*YR)
  ts1 = result.t
  xs1 = result.y[0, :]
  ys1 = result.y[1, :]
  vxs1 = result.y[2, :]
  vys1 = result.y[3, :]
  Es1 = 0.5*(vxs1**2+vys1**2) - MU/np.sqrt(xs1**2 + ys1**2)
  E_err1 = (Es1 - E0)/E0

  plt.figure(figsize=(4,12))
  plt.plot(xs/AU, ys/AU)
  plt.plot(xs1/AU, ys1/AU)
  plt.scatter([0], [0], marker="x", color="k")
  plt.gca().add_artist(mpatches.Circle((0,0), 1.0, fill=None, ls=":"))
  plt.gca().set_aspect("equal")
  plt.xlabel("x [au]")
  plt.ylabel("y [au]")
  plt.tight_layout()

  plt.figure()
  plt.plot(ts/YR, E_err)
  plt.plot(ts1/YR, E_err1)
  plt.title("Energy error (relative)")
  plt.xlabel("t [yr]")
  plt.ylabel("(E-E0)/E")
  plt.tight_layout()

  plt.figure()
  print(hs)
  plt.plot(ts/YR, hs/3600)
  plt.title("Stepsizes")
  plt.xlabel("t [yr]")
  plt.ylabel("stepsize [hours]")
  plt.tight_layout()

  plt.show()

# ------------------------------------------------------------------------------

if __name__ == "__main__":

  # harmonic_oscillator()
  Halleys_comet()
