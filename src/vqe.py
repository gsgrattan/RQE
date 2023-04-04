import cirq
import qsimcirq
import random as rand

import numpy as np

import ising
import result

#TODO: Verify my functionality
class vqe:
    def __init__(self, problem_hamiltonian:ising.ising):

        self.Np, self.J, self.kappa, self.lattice, self.disorder = problem_hamiltonian.get_params()

        self.q = cirq.NamedQubit("p").range(self.Np, prefix="p")
    
    
    #Set the Parameters of the 
    def set_parameters(self, dt, p2, N_layers):
        self.dt = dt
        self.p2 = p2
        self.p1 = p2/10

        self.N_layers = N_layers


    #Callthe simulator
    def simulate(self,i:int =1):
        self.__construct_circuit()

        simulator = qsimcirq.QSimSimulator()
        results = result.result(simulator.simulate(self.circuit))

        return results

        
    
    def __construct_circuit(self):
        self.circuit = cirq.Circuit()
        self.circuit.append(self.__initial_state())
        for k in range(self.N_layers):
            t = (k+1)/(self.N_layers + 1)
            self.circuit.append(self.__algorithmic_layer(t))



    def __algorithmic_layer(self, t):

        yield self.__H0_layer(t)
        yield self.__Hp_layer(t)

    def __Hp_layer(self,t):
        power = -2*self.dt*self.__f(t)

        for pair in self.lattice:
            qubit0 = self.q[pair[0]]
            qubit1 = self.q[pair[1]]
            
            yield (cirq.ZZ**(self.J*power/np.pi)).on(qubit0, qubit1)
            yield self.__double_qubit_error(qubit0, qubit1)

        transverse_field = []
        errors = []
        for qubit in self.q:
            transverse_field.append(cirq.rz(self.kappa*power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))

        yield cirq.Moment(transverse_field)
        yield errors

    def __H0_layer(self, t):
        power = -2*self.dt*(1-self.__f(t))
        gates = []
        errors = []
        for qubit in self.q:
            gates.append(cirq.rx(power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))
            
        yield cirq.Moment(gates)
        yield errors

    def __initial_state(self):
        gates = []
        for q in self.q:
            gates.append(cirq.H.on(q))
        yield cirq.Moment(gates)
    def __f(self, t):
        return (1-(1-t)**2)

    def __double_qubit_error(self, qubit1, qubit2):
        "Choose a qubit to have an error"
        error_qubit = rand.choice([qubit1, qubit2])

        "Choose a random pauli operator to be the error (X,Y,or Z)"
        yield cirq.X.on(error_qubit).with_probability(self.p2/3)
        yield cirq.Y.on(error_qubit).with_probability(self.p2/3)
        yield cirq.Z.on(error_qubit).with_probability(self.p2/3)

    def __single_qubit_error(self, qubit):
        """If the random number in [0,1) is less than p1 An error has occured
        Choose a random pauli operator to be the error (X,Y,or Z)"""
        yield cirq.X.on(qubit).with_probability(self.p1/3)
        yield cirq.Y.on(qubit).with_probability(self.p1/3)
        yield cirq.Z.on(qubit).with_probability(self.p1/3)


    