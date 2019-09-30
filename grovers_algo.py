#!/usr/bin/python3
# -*- coding: utf-8 -*-

from main import * #importing from main.py
import sys
import numpy as np #importing mathematical library
import math     #importing real numbers library


#Create the gates needed to build the pieces of Grover's algorithm, and runs it
def build(n):

    N=2**n; #Total number of states

    H_array=[H] #Creates n-qubit Hadamard gate
    for i in range(1,n):
        H_array.append(H)
    bigH=tensor(H_array)

    buildI=[I] #Creates n-qubit identity gate
    for i in range(1,n):
        buildI.append(I)
    bigI=tensor(buildI)

    init_array=[U]; #Create the state |00....0>, i.e. state with 1 as the first index and 0 everywhere else
    for i in range(1,n):
        init_array.append(U)
    register=tensorReg(init_array); 
    register=bigH.act(register) #Act the Hadamard gate on |00...0>, creating the initial Grover state
    print("The initial state of the register is:")
    print(register.read()) #Print out initial state as a check. Should be equiprobable state
    print("\n")

#Now build Grover's diffusion operator
    grover=bigI.mult(-1).data+register.outer().mult(2).data #Makes gate representing (I-2|psi><psi|) for n-qubits
    grover=gate(grover)

#Now build oracle that randomly selects state to flip
    oracle=tensor(buildI) #Start with n-qubit identity gate

    k=np.random.randint(0,high=N, dtype=int) # pick random number for the state flipped by oracle
    oracle.edit(k,k,-1) #edit oracle entry to flip
 
    print("The oracle flipped state  %d" %k)

    returned = algorithm(N,n,register, oracle,grover,k) #Run grover's algorithm
    
    print("Grover's algorithm retrieved state %d " %returned)
    return k, returned

#Grover's algorithm
def algorithm(N,n,register, oracle,grover,m):
    itnum=int(pi*math.sqrt(N)/4.) #Define the number of iterations to be pi*sqrt(N)/4
    starter=register
    for i in range(0,itnum): #Apply O and then G to the initial input register itnum times
        starter=grover.act(oracle.act(starter))
    return starter.measure() #Measure the state according to the modified probability densities, should retrieve state m


if __name__ == '__main__': # This was relevant when being imported to larger project
    num_q = int(input("Enter number of qubits: ")); #User inputs number of qubits to simulate
    if num_q==0 or num_q<0:
        print("The number of qubits is physically impossible. Terminating the execution.")
        sys.exit()


    build(num_q)
