
#TODO: Implement to store density matrix and expectation values to be returned out of algo.simulate() from both RQE and VQE implementations
import cirq
try:
    import qutip
except ImportError:
    import pip 
    pip.main(["install", "--user", "qutip"])
    import qutip
    
import models

import numpy as np

class result:
    """
    Object to store results and calculate evaluation metrics

    Attributes:
    -----------

    Methods:
    --------
        expect(problem_hamiltonian):
            problem_hamiltonian:ising - The problem hamiltonian of interest
    """
    def __init__(self, sim_result:cirq.SimulationTrialResult):

        self.pq = self.__extract_primary_qubits(sim_result.qubit_map)
        self.N = len(self.pq)

        self.cirq_result = sim_result

    def __extract_primary_qubits(self, qubit_map):
        primary_qubits = []
        for qubit in list(qubit_map):
            if (qubit.name[0] == "p"):
                primary_qubits.append(qubit)

        return primary_qubits

    def __generate_ising_pauli_strings(self,lattice):
        ZZ = []
        X = []

        #Do the Ising ZZ interaction for a general case connectivity
        for pair in lattice:
            q1, q2 = pair
            q1 = self.pq[q1]
            q2 = self.pq[q2]
            ZZ.append(cirq.PauliString(cirq.Z(q1), cirq.Z(q2)))
        for q1 in self.pq:
            X.append(cirq.PauliString(cirq.X(q1)))
        ZZ = cirq.PauliSum.from_pauli_strings(ZZ)
        X = cirq.PauliSum.from_pauli_strings(X)
        return ZZ, X

    def __generate_heisenberg_pauli_strings(self,lattice):
        heisenberg_interaction_gates = []
        
        for pair in lattice:
            q1, q2 = pair
            q1 = self.pq[q1]
            q2 = self.pq[q2]
            #Generate the XX, YY, ZZ interaction for each edge in the lattice
            heisenberg_interaction_gates.append(cirq.PauliString(cirq.X(q1), cirq.X(q2)))
            heisenberg_interaction_gates.append(cirq.PauliString(cirq.Y(q1), cirq.Y(q2)))
            heisenberg_interaction_gates.append(cirq.PauliString(cirq.Z(q1), cirq.Z(q2)))

        return cirq.PauliSum.from_pauli_strings(heisenberg_interaction_gates)

    #TODO: Verify my correctness
    def expect(self, problem_hamiltonian:models.Model):
        """
        Calculates the expectation value of the transverse field ising model using the connectivity in problem_hamiltonian

        Parameters
        ----------
            problem_hamiltonian : ising 
                ising object defined in another 
        """

        expectation_value = 0
        self.Np = problem_hamiltonian.get_Np()
        self.lattice = problem_hamiltonian.get_lattice()

        interactions = problem_hamiltonian.get_interactions()
        #If we are dealing with an ising model

        final_state = self.cirq_result.state_vector()
        qubit_map = self.cirq_result.qubit_map


        if (type(problem_hamiltonian) == models.Model.ising):
            J = interactions[2][0]
            kappa = interactions[1][0]
            pauli_ZZ, pauli_X = self.__generate_ising_pauli_strings(self.lattice)
            ZZ_ev = -1*J*pauli_ZZ.expectation_from_state_vector(final_state, qubit_map = qubit_map)


            X_ev = -1*kappa*pauli_X.expectation_from_state_vector(final_state, qubit_map = qubit_map)

            expectation_value += np.real(ZZ_ev + X_ev)
        
        #If we are dealing with a heisenberg model
        elif (type(problem_hamiltonian) == models.Model.heisenberg):
            J = interactions[2][0]

            pauli_XXYYZZ = self.__generate_heisenberg_pauli_strings(self.lattice)

            XXYYZZ_ev = J*pauli_XXYYZZ.expectation_from_state_vector(final_state, qubit_map = qubit_map)

            expectation_value += np.real(XXYYZZ_ev)



       
        return expectation_value


    def negativity(self, density_data):
        """
        Calculates the negativity over a splitting partition for the given density matrix data
        NOTE: assumes ring topology of Ising Model
        
        Parameters
        ----------
            density_data : array_like
                2^(Np) x 2^(Np) numpy array that stores the data for the density matrix for the primary system

        """
        N = int(np.log2(len(density_data)))

        rho = qutip.Qobj(inpt = density_data, dims=[[2]*N, [2]*N])
        splitting_subspace = ([1]* (N//2)) +  ([0] * (N//2))

        if N %2 == 1:
            splitting_subspace.append(1)

        eigenvalues = qutip.partial_transpose(rho, mask = splitting_subspace).eigenenergies()
        
        return self.__calc_negativity(eigenvalues)

    def __calc_negativity(self,eig):
        """
        Sums the negative eigenvalues
        """
        neg = 0
        for val in np.real(eig):
            if val < 0:
                neg += abs(val)

        return neg


    def get_final_density_matrix(self):
        return self.cirq_result.density_matrix_of(self.pq)

    def get_final_state_vector(self):
        return self.cirq_result.final_state_vector()

    

        



    

    
