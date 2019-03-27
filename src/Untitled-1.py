propagator.propagate_adaptive(rtol=1e-7, atol=1e-12)

def propagate_adaptive(self, h=1e-6, atol=1e-12, safety=0.9,
                           rtol=1e-5, cwise=False):
    uftemp, error = self.model.step(self.h, self.uf)
    self.uf = uftemp
    return ifft(self.uf)

def step(self, h, Efin):
        """ one step with error determined by step-doubling """
        self.setup_disp(h)
        Ec = self.one_step(h, Efin)
        self.setup_disp(h/2.0)
        Ef1 = self.one_step(h/2.0, Efin)
        Ef2 = self.one_step(h/2.0, Ef1)
        # TODO: check this error code and Ef
        error = (Ef2 - Ec)/(2.0**4 - 1.0)
        Ef = Ef2 + error
        return Ef, abs(error)

def one_step(self, h, Efin):
    """ RK4IP step after J. Hult, JLT 25,12 p.3770 (2007) """
    self.D = sp.exp(self.system.L*h/2.0)
    self.uf = self.model.one_step(self.h, self.uf)
    Efin = self.uf
    Ei = self.D*Efin
    k1 = h*self.D*self.nl_polar(Efin)
    k2 = h*self.nl_polar(Ei + k1/2.0)
    k3 = h*self.nl_polar(Ei + k2/2.0)
    k4 = h*self.nl_polar(self.D*(Ei + k3))
    res = self.D*(Ei + (k1/6.0 + k2/3.0 + k3/3.0)) + k4/6.0
    res[self.system.invalid_ws_i] = 0.0
    return res

def nl_polar(self, Efin):
    Et = ifft(Efin)
    Et2 = Et**2
    rk_polar = self.kerr_polar(Et, Et2, Efin)
    nl_polar = fft(rk_polar)
    self.nl_polar_coef = -1j*self.grid.W/(2*c*self.neff)
    return self.system.nl_polar_coef*nl_polar

def kerr_polar(self, Et, Et2, Efin):
    self.kpolar_coef = (1.0 - self.fr)*self.chi3
    ret = self.system.kpolar_coef*Et*Et2
    return ret


pl = fftw3.Plan(i, o, self.direction, flags=('measure',))
ut = ifft(self.uf)