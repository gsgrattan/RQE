import cirq
import qsimcirq
import random as rand

import numpy as np

from problem import problem


class VQE:
    def __init__(self, problem_instance:problem):
        self.Np = problem_instance.Np
        self.J = problem_instance.J
        self.kappa = problem_instance.kappa
        self.ZZ_couplings = problem_instance.ZZ

        self.q = cirq.NamedQubit("q").range(self.Np, prefix="q")
    
    
    #Set the Parameters of the object

    def set_parameters(self, dt, p2, N_layers):
        self.dt = dt
        self.p2 = p2
        self.N_layers = N_layers


    #Callthe simulator
    def simulate(self,i):
        self.__construct_circuit()

        simulator = qsimcirq.QSimSimulator()
        results = simulator.simulate(self.circuit)

        
    
    def __construct_circuit(self):
        self.circuit = cirq.Circuit()
        for k in range(self.N_layers):
            t = (k+1)/(self.N_layers + 1)
            self.circuit.append(self.__algorithmic_layer(t))



    def __algorithmic_layer(self, t):

        yield self.__Hp_layer(t)
        yield self.__H0_layer(t)

    def __Hp_layer(self,t):
        power = -2*self.dt*self.__f(t)

        for pair in self.ZZ_couplings:
            qubit0 = self.q[pair[0]]
            qubit1 = self.q[pair[1]]
            yield (cirq.ZZ**(power/np.pi)).on(qubit0, qubit1)
            yield self.__double_qubit_error(qubit0, qubit1)

        transverse_field = []
        errors = []
        for qubit in self.q:
            transverse_field.append(cirq.rz(power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))

        yield [cirq.Moment(transverse_field), cirq.Moment(errors)]

    def __H0_layer(self, t):
        power = -2*self.dt*(1-self.__f(t))
        gates = []
        errors = []
        for qubit in self.q:
            gates.append(cirq.rx(power).on(qubit))
            errors.append(self.__single_qubit_error(qubit))
            
        yield [cirq.Moment(gates), cirq.Moment(errors)]


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


    