import ising
import result

import numpy as np
from numpy import random

import cirq
import qsimcirq


#TODO: Verify my functionality
class rqe:
    def __init__(self, problem_hamiltonian:ising.ising, thermometry:bool= False, thermometry_error:bool = False):
        """
        Initializes the RQE object with the Problem Hamiltonian of choice
        """
        self.Np, self.J, self.kappa, self.lattice, self.disorder = problem_hamiltonian.get_params()
        self.Hp = problem_hamiltonian
        self.thermometry = thermometry
        self.thermometry_error = thermometry_error



        self.pq = cirq.NamedQubit("p").range(self.Np, prefix="p")

    def set_circuit_parameters(self, dt:float=0.0667, N_lps:int=20, N_r:int=10, p2:float=0.0) -> None:
        """
        Set the parameters for the RQE circuit
        """
        self.dt = dt
        self.p2 = p2
        self.p1 = p2/10
        self.N_lps = N_lps
        self.N_r = N_r


        self.phase_sum = np.pi/4

        self.ramping_function_scaling_factor = self.__normalize_ramping_function(N_lps)
    
    def set_shadow_qubit_parameters(self, primary_shadow_coupling:str="1:1", B0:float=5, Bf:float=0.75, t1:float=0.75) -> None:
        """
        Set the parameters for the shadow qubits
        """
        if (primary_shadow_coupling == "1:1"):
            self.Ns = self.Np
        else:
            self.Ns = self.Np//2
        
        self.sq = cirq.NamedQubit("s").range(self.Ns, prefix="s")
        self.B0 = B0
        self.Bf = Bf
        self.t1 = t1
        self.sq_slope = (Bf - B0)/t1

    
    def simulate(self, i:int =1):
        """
        Simulates the RQE Circuit
        """
        self.__construct_circuit()

        simulator = qsimcirq.QSimSimulator()

        results = result.result(simulator.simulate(self.circuit))

        return results

    def __construct_circuit(self):
        """
        Constructs the RQE circuit based on the provided parameters
        """
        self.circuit = cirq.Circuit()

        for r in range(self.N_r):
            #If it is the last reset and we are analyzing thermometry, do a thermometry cycle
            if ((r == self.N_r -1) and self.thermometry):
                self.circuit.append(self.__thermometry_cycle())
            else:
                self.circuit.append(self.__rqe_cycle())

    def __algorithmic_layer(self, t):
        yield self.__primary_layer(t)
        yield self.__shadow_layer(t)
        yield self.__primary_shadow_layer(t)

    def __rqe_cycle(self):
        for k in range(self.N_lps):
            t = (k+1)/(self.N_lps +1)
            yield self.__primary_layer(t)
            yield self.__shadow_layer(t)
            yield self.__primary_shadow_layer(t)
        yield self.__reset_layer()

           

    def __thermometry_cycle(self):
        for k in range(self.N_lps):
            t = (k+1)/(self.N_lps+1)
            yield self.__primary_layer(t)
            yield self.__thermometry_layer(t)
            yield self.__primary_shadow_layer(t)
        yield cirq.measure(self.sq, key="s")



    def __primary_layer(self, t):
        """
        Implements the Trotter layer for the Primary system
        """
        #initialize the power for the cirq.ZZ gate (Ising Interaction)
        power = -2*self.dt
        for pair in self.lattice:
            qubit0 = self.pq[pair[0]]
            qubit1 = self.pq[pair[1]]

            #Add the ising interaction and a potential error
            yield ((cirq.ZZ**(self.J*power/np.pi)).on(qubit0, qubit1))
            yield self.__double_qubit_error(qubit0, qubit1)


        #Add the disorder
        if (self.disorder):
            gates = []
            errors = []
            h_vals = self.Hp.get_disorder_vals()
            for i, qubit in enumerate(self.pq):
                h = h_vals[i]
                gates.append(cirq.rz(h*power).on(qubit))
                errors.append(self.__single_qubit_error(qubit))

            yield cirq.Moment(gates)
            yield errors

        #Transverse Field
        gates = []
        errors = []

        for i, qubit in enumerate(self.pq):
            gates.append(cirq.rx(self.kappa*power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))

        yield cirq.Moment(gates)
        yield errors

    def __primary_shadow_layer(self, t):
        """
        Implements the Trotter layer for the primary shadow interaction
        """
        #Initialize the power of the YY gate
        power = -2*self.dt * self.__ramping_function(t) * self.ramping_function_scaling_factor
        primary_shadow_gate = cirq.YY**(power/np.pi)

        #If we have a 1:1 coupling
        if (self.Np == self.Ns):

            for i in range(self.Np):
                primary_shadow_gate = cirq.YY**(power/np.pi)

                qubit0 = self.pq[i]
                qubit1 = self.sq[i]

                yield primary_shadow_gate.on(qubit0, qubit1)
                yield self.__double_qubit_error(qubit0, qubit1)

        #If we have a 2:1 coupling
        else:
            for i in range(self.Ns):
                primary_shadow_gate = cirq.YY**(power/np.pi)

                qubit0 = self.pq[2*i]
                qubit1 = self.sq[i]

                yield primary_shadow_gate.on(qubit0, qubit1)
                yield self.__double_qubit_error(qubit0, qubit1)

    def __shadow_layer(self, t):
        """
        Implements the Trotter layer for the shadow qubit system
        """
        gates = []
        errors = []
        #Iterate through the shadow qubits
        for index, qubit in enumerate(self.sq):
            energy = self.__sq_energy_sweep(t)

            power = -2*self.dt*energy/2
            gates.append(cirq.rz(power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))

        yield cirq.Moment(gates)
        yield errors

    def __thermometry_layer(self,t):
        gates = []
        errors = []
        energies = self.__sq_energy_sample(t)

        for index, qubit in enumerate(self.sq):
            energy = energies[index]

            power = -2*self.dt*energy/2
            gates.append(cirq.rz(power).on(qubit))
            if (self.thermometry_error):
                errors.append(self.__single_qubit_error(qubit))

        yield cirq.Moment(gates)
        yield errors       


    def __reset_layer(self):
        """
        Adds a reset to all the shadow qubits
        """
        gates = []
        for qubit in self.sq:
            gates.append(cirq.reset(qubit))
        yield cirq.Moment(gates)

    def __double_qubit_error(self, qubit1, qubit2):
        error_qubit = random.choice([qubit1, qubit2])

        yield cirq.X.on(error_qubit).with_probability(self.p2/3)
        yield cirq.Y.on(error_qubit).with_probability(self.p2/3)
        yield cirq.Z.on(error_qubit).with_probability(self.p2/3)

    def __single_qubit_error(self, qubit):
        yield cirq.X.on(qubit).with_probability(self.p1/3)
        yield cirq.Y.on(qubit).with_probability(self.p1/3)
        yield cirq.Z.on(qubit).with_probability(self.p1/3)


    def __sq_energy_sweep(self, t):
        """
        Implements the shadow qubit sweeping function from the Berg paper
        """
        if (t > self.t1):
            return self.Bf
        else: 
            return self.B0 + self.sq_slope*t
    def __sq_energy_sample(self, t):
        if (abs(t - 1/(self.N_lps + 1)) < 0.001):
            self.sq_energy = np.random.uniform(low=1,high=6, size=self.Ns)
        return self.sq_energy

    def __ramping_function(self, t):
        return 4*t*(1-t)

    def __normalize_ramping_function(self, N_lps):
        "Cumulative Sum of the Phases"
        cum_sum = 0
        for k in range(N_lps):
            t = (k + 1) / (N_lps + 1)
            cum_sum += self.dt * self.__ramping_function(t)

        return self.phase_sum / cum_sum