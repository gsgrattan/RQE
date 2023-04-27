
import random
import networkx as nx
import cirq
class ising:
    """
    A class to represent a transverse field Ising model.

    Attributes: 
    -----------
    Np : int 
        number of spins/qubits in the model
    J : float
        strength of the ising interaction
    kappa : float
        strength of the transverse field
    disorder : bool
        True if we want a longitufdional disorder on each of the spins
    periodic_bc : bool
        True if the lattice has periodic boundry conditions

    Methods:
    --------
    get_params():
        Returns the parameters needed to simulate the hamiltonian.

    Methods defined here:

        __init__(self, Np, J, kappa, disorder, periodic_bc)
            Constructs the ising model with all the necessary parameters.
        
            Parameters
            ----------
                Np : int
                    The number of spins/qubits in the model
                J : float = 1 
                    strength of the ising interaction
                kappa : float = 0.5
                    strength of the transverse field
                disorder : bool = False
                    True if we want longitudional disorder on each of the spins
                periodic_bc : bool = True
                    True if the lattice has periodic boundry conditions

        get_params(self)
                Returns the parameters needed to simulate the hamiltonian.

            Returns
            -------
                Np : int
                    The number of spins/qubits in the model
                J : float = 1 
                    strength of the ising interaction
                kappa : float = 0.5
                    strength of the transverse field
                lattice : array( tuple(int,int) )
                    coupling of sites in the ising model
                disorder : bool
                    True if we want longitudional on each of the spins
    """
   
    def __init__(self, Np:int, J: float = 1.0, kappa: float = 0.5, disorder: bool = False, periodic_bc:bool = True):
        # Store the necessary m 
        self.Np = Np

      
        self.J = 1*J

        self.kappa = kappa
        self.periodic = periodic_bc
        self.disorder = disorder 

        self.interactions = {1:[cirq.X], 2: [cirq.ZZ]}


        
        self.lattice = self.__create_lattice(Np, periodic_bc)
        self.disorder_vals = self.__create_h_vals(Np, disorder)    

    def __create_lattice(self, N, periodic):
        
        #Creates a 1D lattice of qubits, with or without periodic boundary conditions

        edges = []
        
        for i in range(N - 1):
            edges.append((i, i+1))
       
        if periodic:
            edges.append((0, N - 1))
       
        return edges


    def __create_h_vals(self, N, is_disorder):

        #Strengths of a varying longitudional field that creates disorder
         
        hvals = [0.0]*N
        if (is_disorder):
            for i in range(N):
                hvals[i] = (2*self.J - -2*self.J) * random.random() - 2*self.J
        return hvals
    
    #Getters and Setters
    def get_params(self):
        """
        Returns the parameters needed to simulate the hamiltonian.
        """
        return self.Np, self.J, self.kappa, self.lattice, self.disorder
    
    def get_disorder_vals(self):
        """
        Returns the disorder values
        """
        return self.disorder_vals

    def get_interactions(self):
        return self.interactions

    def set_lattice(self, lattice_edges):
        self.lattice = lattice_edges
       


class heisenberg: #Insert breaking bad reference here

    def __init__(self, Np:int, J:float, disorder:bool = False, periodic_bc:bool =True):
        self.Np = Np
        self.J = J
        self.disorder = disorder
        self.periodic = periodic_bc

        self.interactions = {2: [cirq.XX, cirq.YY, cirq.ZZ]}

        self.lattice = self.__create_lattice(N=Np, periodic=periodic_bc)
        self.disorder_vals = self.__create_h_vals(N=Np, is_disorder=disorder)

   
    def __create_lattice(self, N, periodic):
        
        #Creates a 1D lattice of qubits, with or without periodic boundary conditions

        edges = []
        
        for i in range(N - 1):
            edges.append((i, i+1))
       
        if periodic:
            edges.append((0, N - 1))
       
        return edges

    def __create_h_vals(self, N, is_disorder):

        #Strengths of a varying longitudional field that creates disorder
         
        hvals = [0.0]*N
        if (is_disorder):
            for i in range(N):
                hvals[i] = (2*self.J - -2*self.J) * random.random() - 2*self.J
        return hvals


    #Getters and Setters
    def set_lattice(self, lattice_edges):
        self.lattice = lattice_edges
        
    def get_params(self):
        return self.Np, self.J, self.lattice, self.disorder
    
    def get_disorder_vals(self):
        return self.disorder_vals

    def get_interactions(self):
        return self.interactions

