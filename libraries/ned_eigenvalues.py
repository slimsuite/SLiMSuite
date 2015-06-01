#!/usr/bin/python

# Modified N. Davey RLC module
# Copyright (C) 2009 Norman E. Davey & Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       ned_eigenvalues
Description:  Modified N. Davey Relative Local Conservation module
Version:      1.0
Last Edit:    03/09/09
Copyright (C) 2009 Norman E. Davey & Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is not for standalone running and has no commandline options (including 'help'). All options are handled
    by the parent module.

Uses general modules: operator, math, random
"""
#########################################################################################################################
import operator, math, random

def isList(l):
	for method in ['__getitem__', '__setitem__']:
		if method not in dir(l):
			return 0
	return 1

def dot(v1,v2): 
	summer = 0
	for i in range(len(v1)):
		summer += v1[i]*v2[i]
	
	return summer
	
def list_approx_equaltozero(l):
	for v in l:
		if abs(v) > 0.000001:
			return False
			
	return True

def approx_equalto(z):
	return abs(z) < 0.000001

def getconj(z):
	try:
		return z.conjugate()
	except AttributeError:
		return z
		
def norm(table):  
	return math.sqrt(abs(dot(table,map(getconj,table)) ))

def normalise(table): 
	return divide(table,norm(table))

def divide(mat,val):
	res = []
	for x in mat:
		if val != 0:
			try:	
				res.append(x/val)
			except:
				print x,val
		else:
			res.append(0.0)
		
	return res
	
def equal(x):
	equal_matrix = []
	for rows in range(x):
		row_val = []
		for cols in range(x):
			row_val.append(rows == cols)
		equal_matrix.append(row_val)
	
	return equal_matrix
	
def add(a,b):
	min_matrix = []
	for rows in range(len(a)):
		row_val = []
		for cols in range(len(a[rows])):
			row_val.append(a[rows][cols] + b[rows][cols])
		min_matrix.append(row_val)
	return min_matrix
	
	
def minus(a,b):
	min_matrix = []
	for rows in range(len(a)):
		row_val = []
		for cols in range(len(a[rows])):
			row_val.append(a[rows][cols] - b[rows][cols])
		min_matrix.append(row_val)
	return min_matrix
	
def multiply(mat,val):
	return [x*val for x in mat]

	
def convert_to_Matrix(matrix):
	new_matrix = []
	for row in 	matrix:
		new_matrix.append(row[0])
	
	return new_matrix
	
def mmul_vect(a,b):
	nrows = len(a)
	w = [None] * nrows
	for row in range(nrows):
		w[row] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, a[row], b))
	return w
	

def mmul_int(a,b):
	mul_matrix = []
	
	for row in a:
		row_val = []
		mul_matrix.append(multiply(row,b))
		
	if len(mul_matrix) == 1:
		mul_matrix = mul_matrix[0]
		
	return  mul_matrix
	
def mmul(a,b):
	mul_matrix = []
	
	for row in a:
		row_val = []
		for col in zip(*b):
			summer = 0
			for i, j in zip(row, col):
				summer += i*j
			
			row_val.append(summer)

		mul_matrix.append(row_val)
		
	return  mul_matrix

def outer(v1,v2):
	new_matrix =[]
	for a in v1:
		
		row = []
		for b in v2:
			row.append(a*b)
		
		new_matrix.append(row)
		
	
	return new_matrix



def read_scoring_matrix(sm_file):
    '''Read in a scoring matrix from a file, e.g., blosum80.bla, and return it as an array.'''
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
    # dictionary to map from amino acid to its row/column in a similarity matrix
    aa_to_index = {}
    for i, aa in enumerate(amino_acids): aa_to_index[aa] = i


    aa_index = 0
    first_line = 1
    row = []
    list_sm = [] # hold the matrix in list form

    try:
		matrix_file = open(sm_file, 'r')

		for line in matrix_file:

			if line[0] != '#' and first_line:
				first_line = 0
			if len(amino_acids) == 0:
				for c in line.split():
					aa_to_index[string.lower(c)] = aa_index
					amino_acids.append(string.lower(c))
					aa_index += 1
			elif line[0] != '#' and first_line == 0:
				if len(line) > 1:
					row = line.split()
					list_sm.append(row)

    except IOError, e:
		print "Could not load similarity matrix: %s. Using identity matrix..." % sm_file
		return identity(20)
	
    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
		for i in range(0,19):
		    for j in range(i+1, 20): list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
		for j in range(len(list_sm[i])):
			try: list_sm[i][j] = float(list_sm[i][j])
			except: print'Fuck',i,j,list_sm[i][j]

    return list_sm
    #sim_matrix = array(list_sm, type=Float32)
    #return sim_matrix

		
class Eigenvalues:

	def forwardSub(self,matrix, b ):
		b = b[0] 
		x = []
		for i in range(len(matrix)):
			try:
				x.append(( b[i] - dot(x,matrix[i][:i])) / matrix[i][i] )
			except:	
				x.append(b[i] )
		return x

	
	def qr_algorithm(self,matrix):
		m = len(matrix)
		n = len(matrix[0])
		R = matrix
		for i in range(min(m,n)):
			v, beta = self.householder_Vector(zip(*R)[i],i)
			R = minus(R,outer(v,mmul_int([mmul_vect(zip(*R),v)],beta)))
		
		for i in range(1,min(n,m)): 
			R[i][:i] = [0] * i
	
		R = R[:n]	
		Q_temp = zip(self.iterate_Solution(zip(*R),zip(*matrix)))
		
		Q = []
		for q in Q_temp:
			Q.append(q[0])
			
		return Q,R
		
		
	def iterate_Solution(self,matrix,b):
		
		if len(b) > 1:
			Q = []
			for vec in zip(*b):
				sol = zip(self.iterate_Solution(matrix,[vec]))
				Q.append(sol[0][0])
				
			return Q
		else:
			x = self.forwardSub(matrix,b)
			diff = minus(b,[mmul_vect(matrix,x)])
			maxdiff = dot(diff[0],diff[0])
			
			i = 0
			while i < 10:
				i += 1
				if isList(x[0]):
					add1 = x
				else:
					add1 = [x]
					
				add2 = self.forwardSub(matrix, diff)
				
				
				if isList(add2[0]):
					pass
				else:
					add2 = [add2]
				
				xnew = add(add1,add2)
			
				diffnew = minus(b,[mmul_vect(matrix,xnew[0])])
				
				maxdiffnew = dot(diffnew[0],diffnew[0])
				if approx_equalto(maxdiffnew - maxdiff): 
					i = 10
					
				x, diff, maxdiff = xnew, diffnew, maxdiffnew
	
	
		return x
		
	def householder_Vector(self,vector, index ):
		
		v = normalise([0.0]*index + list(vector[index:]))
	
		t = v[index]
		sigma = 1.0 - t**2
		
		if sigma != 0.0:
			t = v[index] = t<=0 and t-1.0 or -sigma / (t + 1.0)
			v = divide(v,t)
			
		return v, 2.0 * t**2 / (sigma + t**2)
		
	def hessenberg_matrix(self,matrix ):
	
		for i in range(len(matrix)-2):
			v,beta =  self.householder_Vector(zip(*matrix)[i],i+1)
			
			matrix = minus(matrix,outer(v,mmul_int([mmul_vect(zip(*matrix),v)],beta)))
			matrix = minus(matrix,outer(mmul_vect(matrix,v),mmul_int([v],beta)))
			
		return matrix
	
	def get_eigenvalues(self,matrix):
		matrix = self.hessenberg_matrix((matrix))
		
		eigvalues = []
		for i in range(len(matrix)-1,0,-1):
			iter = 0
			while not list_approx_equaltozero(matrix[i][:i]):
				iter += 1
				if iter > 100:
					break
					
				shift = mmul_int(equal(i+1) , matrix[i][i])
				q, r = self.qr_algorithm(minus(matrix,shift))
				
				matrix = add(mmul(r,q),shift)
				
				
			eigvalues.append(matrix[i][i] )
			
			matrix = [matrix[r][:i] for r in range(i)] 
		eigvalues.append( matrix[0][0] )
		
		
		return eigvalues
	

