import ising

class rqe:
    def __init__(self, problem_hamiltonian:ising.ising, primary_shadow_qubit_coupling:str = "1:1"):
        self.Np, self.J, self.kappa, self.lattice, self.disorder = problem_hamiltonian.get_params()