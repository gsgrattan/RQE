
#TODO: Implement to store density matrix and expectation values to be returned out of algo.simulate() from both RQE and VQE implementations
import cirq
import numpy as np
import helpers
import ising

class result:
    def __init__(self, sim_result:cirq.SimulationTrialResult, N:int):

        self.pq = self.__extract_primary_qubits(sim_result.qubit_map)
        self.N = len(self.pq)

        self.cirq_result = sim_result

        #initialize the variables this class calculates
        self.negativity = 0
        self.expectation = 0

        
    
    def __extract_primary_qubits(self, qubit_map):
        primary_qubits = []
        for qubit in list(qubit_map):
            if (qubit.name[0] == "p"):
                primary_qubits.append(qubit)

        return primary_qubits

    def __generate_pauli_strings(self,lattice):
        ZZ = []
        X = []

        #Do the Ising ZZ interaction for a general case connectivity
        for pair in lattice:
            q1, q2 = pair
            q1 = self.pq[q1]
            q2 = self.pq[q2]
            ZZ.append(cirq.ZZ.on(q1, q2))
        for qubit in self.pq:
            X.append(cirq.X.on(qubit))

        return  cirq.PauliString(ZZ), cirq.PauliString(X)
    
    def expect(self, problem_hamiltonian:ising.ising):
        self.Np, self.J, self.kappa, lattice, self.disorder = problem_hamiltonian.get_params()

        pauli_ZZ, pauli_X = self.__generate_pauli_strings(lattice)
        final_state = self.cirq_result.state_vector()
        qubit_map = self.cirq_result.state_vector()

        ZZ_ev = pauli_ZZ.expectation_value_from_state_vector(final_state, qubit_map = qubit_map)*self.J

        X_ev = pauli_X.expectation_value_from_state_vector(final_state, qubit_map = qubit_map)*self.kappa
        return ZZ_ev + X_ev



    

    def negativity(self):
         #Get the Density matrix cooresponding to the primary qubits
        self.rho_p =self.cirq_result.density_matrix_of(self.pq)

        #calculate the negativity
        self.negativity = helpers.negativity(self.rho_p)

    def get_final_density_matrix(self):
        return self.rho_p

    def get_final_state_vector(self):
        return self.cirq_result.final_state_vector()

    

        



    

    
