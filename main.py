#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np #importing mathematical library.
import math     #importing real numbers library
import cmath #importing complex numbers library


NumTypes=(int, float, complex, np.int64) 
IntTypes=(int, np.int64)

class reg():
    """Base class for registers.

     Based on numpy.ndarray 1*2^n sized array.
     """
     #TODO: add empty register feature
    def __init__(self,data):#constructor of register with data
  
        
        """
        Constructor.
        Takes only one external parameter, register initial state as numpy.ndarray.
        The parameter is checked by reg.datacheck() method, normalized by reg.normalize and
        written into reg.data variable. Length of the state vector is stored in reg.size variable.

        :param data: initial state
        :return:
        """
        #else:
        self.datacheck(data) #here we send the data to be checked.
        self.data=self.normalize(data);
        self.size=len(data);

    def datacheck(self,data):
        """
        Data checking function.
        Takes one external parameter, checks if it is a numpy.ndarray,
        if it has length equal to an integer order of 2, and checks all the elements to be numbers. In case of
        any non-accordance, raises error. No returned data.

        :param data:
        :return:
        """
        if not isinstance(data, np.ndarray):
            raise NameError('Wrong input data type') #if our input data is not numpy array, generate error
        if math.log(len(data),2)%1!=0:
            raise NameError('Wrong size of register. Not order of 2.')  #Checks size of register # TODO: Implement padding to make it of order 2?
        for elem in data:
            if not isinstance(elem, NumTypes):
                raise NameError('Values of register must be numbers') #Checks if values in array are numbers
    def normalize(self,data):
        """
        Normalizes any entered register state.
        Takes numpy.ndarray as external parameter, finds norm of the
        state vector and divides the array by the norm. Returns normalized array.

        :param data:
        :return:
        """
        norm=np.dot(abs(data),abs(data))
        if norm==1: #checks for normalization
            return data #if normalized, returns data
        else: 
            return (complex(1)/math.sqrt(norm))*data #normalizes data
    def read(self):
        """
         returns the register state -- reg.data.
        :return:
        """
        return self.data;
    def getsize(self):
        """
        returns the register size -- reg.size.
        :return:
        """
        return self.size;   
    def rewrite(self,data):
        """
         Checks external parameter using reg.datacheck(), checks length of
        the array to be equal the register’s state length, rewrite reg.data with the new data. No
        return.

        :param data:
        :return:
        """
        self.datacheck(data) #send the new value to be checked
        if self.size!=len(data):
            raise NameError('Size of the new data is not equal to size of the register') # Check size of new data
        self.data=data; #rewrite the value after
    #def edit(self, state_number, new_value):
    #    pass
    def measure(self): #simulates wavefunction collapse and measurement
        """
        Reduces the register’s quantum state, returns number of basis state.

        :return: number of state
        """
        n=self.size # get size to define range of integers to choose from
        probs=abs(self.data)**2  # use the state amplitudes to create distribution to sample from              
        k=np.random.choice(n,1,p=probs) # choose random integer that corresponds to measured state

        return k
     
    def outer(self): #takes the matrix outer product of 2 registers to make a gate
        """
         Returns gate object made from the register’s state outer product matrix.
        :return:
        """
        return gate(np.outer(self.read(),self.read()))
    def mult(self,num): #allows for multiplication of registers by scalars
        """
        Does multiplication of the register by scalar
        :param num: scalar
        :return:
        """
        if not isinstance(num, NumTypes):
                raise NameError('Multiplication factor must be a number') # Checks if it is scalar
        return reg(num*self.data)
           
class gate(): #base class for gates
 
    
    #TODO: add empty gate feature
    #Size is a single parameter, because we consider the matrix to be square
  
    def __init__(self, data):#constructor of gate with data
        '''
         Constructor. Takes external parameter numpy.ndarray data,
        checks it with gate.data check, writes down data and size of the matrix as the class pa-
        rameters.

        :param data:
        :return:
        '''
        self.data_check(data);
        #nparray of matrix elements
        self.data=data; #after the checks, we can finally save the data
        #size of the matrix
        self.size=len(data); #change the size
    def data_check(self, data):
        '''
         Data checking function. Checks external parameter to be
        square numpy.ndarray, check every element to a number .

        :param data:
        :return:
        '''
        if not isinstance(data, np.ndarray):
            raise NameError('Wrong input data type') #if our input data is not numpy array, generate error
        for row in data:
            if len(row)!=len(data):
                raise NameError('The data is not square') #we check each row to have the same length as the number of rows. If not, generate error
            for elem in row:
                if not isinstance(elem, NumTypes):
                    raise NameError('Elements of the matrix must be numbers') # Checks element type, need numbers
    def act(self, register):    #apply the gate on a register
        '''
         Applies the gate on a register. Implemented by numpy.ndarray
        multiplication. Returns changed state register.

        :param register:
        :return:
        '''
        self.regcheck(register);
        k=np.dot(self.data,register.read()) #matrix multiplication by numpy library
        return reg(k);  #rewrite the new value into the register
    def regcheck(self,register):
        '''
        Checks possibility of applying the gate on a register. Raises an error in case of problems.

        :param register:
        :return:
        '''
        if not isinstance(register,reg):
            raise NameError('The register to apply must be of the register class')  #Needs to be register defined using register calss
            raise NameError('Size of the gate does not match size of the register') # Sizes need to match
    def edit(self,x,y,value):   #change value of a matrix element (x,y)
        '''
        Changes value of a matrix element (x,y).

        :param x: coordinate
        :param y: coordinate
        :param value:
        :return: no return
        '''
        self.data[x,y]=value
        return self
    def read(self):
        '''
         Returns numpy.ndarray of matrix elements.

        :return: numpy.ndarray
        '''
        return self.data;
    def getsize(self):
        '''
         returns size of according quantum state.

        :return: size (int)
        '''
        return self.size
    def mult(self,num):
        '''
         multiplies all the gate matrix elements by a scalar.
         

        :param num: scalar factor
        :return:
        '''
        if not isinstance(num, NumTypes):
                raise NameError('Multiplication factor must be a number') # Checks if it is scalar
        return gate(num*self.data)


def tensor(gates): #takes in an array of gates and returns the tensor product of all of them

    for elem in gates: #check that every element is a gate
            if not isinstance(elem, gate):
                raise NameError('Inputs must be gates')                                
    n=len(gates)
    newgate=gates[0].read() #create new array to hold tensor product
    for i in range(1,n): #iteratively create tensor product
        newgate=np.kron(newgate,gates[i].read()); #numpy's implementation of tensor product
    return gate(newgate);

def tensorReg(regs): #exact copy of tensor function but for registers instead of gates
    '''
    exact copy of tensor function but for registers instead of gates
    :param regs:
    :return:
    '''
    for elem in regs: #check that every element is a register
            if not isinstance(elem, reg):
                raise NameError('Inputs must be registers')                                
    n=len(regs)
    newreg=regs[0].read() #create new array to hold tensor product
    for i in range(1,n): #iteratively create tensor product
        newreg=np.kron(newreg,regs[i].read()); #numpy's implementation of tensor product
    return reg(newreg);



#---------------------------------------------------------------------------------------------------
#Now define basic gates and registers necessary for quantum computation including Grover's algorithm

Had=np.ones((2,2),dtype=np.complex)
Had[1,1]=-1
Had=(1/math.sqrt(2))*Had
H=gate(Had) #1-qubit Hadamard gate
Id=np.ones((2,2),dtype=np.complex)
Id[0,1]=0
Id[1,0]=0
I=gate(Id) #1-qubit identity gate
pi=math.pi

up=np.zeros(2,dtype=np.complex)
up[0]=1
U=reg(up) #1-qubit register in the |0> state

down=np.zeros(2,dtype=np.complex)
down[1]=1
D=reg(down) #1-qubit register in the |1> state

def phase(phi): #1-qubit phase gate with inputted choice of phi
    """
    1-qubit phase gate with inputted choice of phi
    :param phi: phase
    :return: gate
    """
    mat=np.zeros((2,2),dtype=np.complex)
    mat[0,0]=1
    mat[1,1]=cmath.exp(phi*1j).real+(cmath.exp(phi*1j).imag)*1j
    return gate(mat)