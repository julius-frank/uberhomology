
import numpy as np

def main():
	# For test purposes: reads in 1-sphere.
	# X=[{frozenset([1]),frozenset([2]),frozenset([3])},
	#    {frozenset([1,2]),frozenset([1,3]),frozenset([3,2])},
	#    {}]
	## reads a finite simplicial complex through input:
	X=input_simplex()  
	dim=len(X)-1 # The dimension of the complex
	vert=len(X[0]) # The number of vertices
	
	
	## For test purposes: colors the vertices 1,2 black
	#eps=set([1,2]) 
	## Read a coloring through input:
	#eps=input_coloring(vert) 
	
	uberhomology(X)


def uberhomology(Y): 
	"""
	Takes a simplicial complex Y and computes its uberhomology as a tri-graded F_2-module.
	Y should be given as a list, Y[n] is the set of n-simplices.
	Every element of Y[n] is a size-n+1-subset of the vertices 1,...,len(Y[0])+1
	It is assumed that Y is closed under taking faces.

	The überhomology of Y is graded by (i,n,w), where 
	 * i is the homological degree,
	 * n is the simplex dimension,
	 * w is the simplex weight.
	"""
	

	###
	#first, every set of simplices Y[n] is made into a (randomly) ordered list X[n]: 
	#This allows us to identify simplices by their index, instead of carrying around the set.
	X=[]
	for n in range(len(Y)):
		X.append(list(Y[n]))
	dim=len(X)-1 #The dimension of the complex
	vert=len(X[0]) #The number of vertices
	###

	###
	#Now, we write up all the possible colorings into a list, according to their number of black vertices: 
	col=[[frozenset()]] #we start with one empty coloring, which has no black vertices (i.e., lives in col[0])
	for i in range(vert):
		col=new_vertex_coloring(col) #new_vertex_coloring takes a list of colorings for i vertices, and gives all the colorings for i+1 vertices
	#note that col has length dim+1.
	###

	###
	#Now, we calculate uberhomology. 
	#This is tri-graded: We'll use the variables (n,w) for the grading by dimension and weight of the simplices, 
	#and the variable i for the homological degree, which corresponds to the length of the coloring. (i.e., no. of black vertices)
	#(n,w) are done at once for each horizontal homology (maybe optimize later?)
	#i needs to be done only iteratively: To get H^i, we only need the colorings of lengths i-1,i,i+1.

	#We'll do the degrees first, and implement explicit bases at a later time.
	horiz_prev={} #this is the horizontal homology of degree i-1: For i=0, this will be 0.
	horiz_deg=[{} for _ in range(vert+1)] #this will be the list of graded degrees of horizontal homology for each coloring
	horiz_degs=[np.zeros([dim+1,dim+2], int) for _ in range(vert+1)] #this will be the list of total graded degrees of horizontal homologies for all colorings of a given length
	
	#Start with calculation for the 1-coloring:
	eps_1=col[vert][0] #The all-1-coloring
	horiz_prev[eps_1]=horizontal_homology(X,eps_1) #get horizontal homology (can simplify, this is just the chain complex!)
	horiz_deg[vert][eps_1]=horiz_prev[eps_1][2] #write the degree matrix into the degree list
	horiz_degs[vert]+=horiz_deg[vert][eps_1] #add to the matrix of cumulative degrees for the given coloring length
	print('The total degrees of horizontal homology for colorings of weight '+str(vert)+' (i.e., '+str(vert)+' black vertices) is:')
	print(horiz_degs[vert])
	Dprev=[[np.zeros([0,horiz_degs[vert][n,w]], int) for w in range(dim+2)] for n in range(dim+1)] #sets up the differential D^vert as (0,deg(H_n,w))-matrices.

	horiz={frozenset(eps_1):horizontal_homology(X,frozenset())} #horiz is a dictionary: It maps a coloring eps to the horizontal homology with that coloring.
	#horiz is initialized with the all-1-coloring.

	# The results will be recorded in terms of degree and basis:
	# For each tri-grading (i,n,w), uber_rank records the rank of that uberhomology.
	uber_rank=np.zeros([vert+1,dim+1,dim+2], int)
	# uber_basis: Expresses the homology in grade (i,n,w) in terms of coloring and simplices
	uber_basis=[[[None for _ in range(dim+2)]for _ in range(dim+1)] for _ in range(vert+1)] #B gives a basis for the homologies in terms of simplices
	#uber_cycles=[[[None for _ in range(dim+2)]for _ in range(dim+1)] for _ in range(vert+1)] #E goes the other way: It expresses cycles in terms of the homology basis
	r=np.zeros([vert+1,dim+1,dim+2], int) #rank of differential starting at (i,n,w)
	c=np.zeros([vert+1,dim+1, dim+2], int) #rank of ubercomplex in degree (i,n,w)
	###c is hust horiz_degs... take one
	for n in range(dim+1):
		for w in range(dim+2):
			c[vert,n,w]=horiz_degs[vert][n,w] #write in complex degrees for i=vert
	Uprev=[[np.identity(horiz_degs[vert][n,w],int) for w in range(dim+2)] for n in range(dim+1)] #Calculation aid 
	Uiprev=[[np.identity(horiz_degs[vert][n,w],int) for w in range(dim+2)] for n in range(dim+1)] #Calculation aid 


	for i in range(vert,0,-1): #calculates H^i for i ranging from vert to 1. (H^0 seperately after.)
		###
		#First, calculate horizontal homology for all colorings of length i-1: (This is the slow part!)
		
		for eps in col[i-1]:
			horiz[eps]=horizontal_homology(X,eps) 
			horiz_deg[i-1][eps]=horiz[eps][2] #write the degree matrix into the degree list
			horiz_degs[i-1]+=horiz_deg[i-1][eps] #add to the matrix of cumulative degrees for the given coloring length
		# print('\r\n We will now enter the loop i='+str(i)+'.')
		# print('We will compute horizontal homologies for colorings with '+str(i-1)+' black vertices.')
		# print('The differential into degree '+str(i)+' comes from those.')
		print('\r\nThe total degrees of horizontal homology for colorings of weight '+str(i-1)+' is:')
		print(horiz_degs[i-1])
		###



		###
		#Next: For each grading (n,w), construct the differential and compute the überhomology:
		D=[[np.zeros([horiz_degs[i][n,w],horiz_degs[i-1][n,w]], int) for w in range(dim+2)] for n in range(dim+1)] #sets up a 0-matrix of the correct size for all (n,w)

		for n in range(dim+1):
			for w in range(dim+2): 
				if not(0 in np.shape(D[n][w])): #Consider only differentials of positive shape
					#print('The shape of D^'+str(i-1)+'_('+str(n)+','+str(w)+') is '+str(np.shape(D[n][w])))
					#To fill D^i+1_n,w, we need to find the basis elements of C^i+1, express them as simplices, apply the differential, and express this in the basis of C^i again.
					b=0 #column coordinate for D^i+1_n,w
					for eps in col[i-1]:
						if horiz_deg[i-1][eps][n][w]>0: #If the horizontal homology is nonzero for this coloring in degree (n,w)
							B=horiz[eps][0][n][w] #this is the basis matrix for H^i+1
							S=[s for s in X[n] if w==n+1-len(set(s).intersection(set(eps)))] #Lists all n-simplices of weight w
							# print('(i-1,n,w)=('+str(i-1)+','+str(n)+','+str(w)+'). eps='+str(eps)+'.')
							# print(str(S))
							# print(B)
							a=0 # Row coordinate for D^i-1_n,w
							for zet in col[i]: #now, consider the higher-degree colorings:
								if horiz_deg[i][zet][n][w]>0:#If the horizontal homology is nonzero for this coloring in degree (n,w)								    
									T=[s for s in S if w==n+1-len(set(s).intersection(set(zet)))]  #Lists all n-simplices that have weight w in eps AND zet (OPTIMIZABLE)
									if len(T)>0 and eps.issubset(zet):#If these exist AND zet just adds one vertex to eps (OPTIMIZABLE)
										#Then there might be nonzero entries
										U=[s for s in X[n] if w==n+1-len(set(s).intersection(set(zet)))]  #Lists all n-simplices with weight w in zet. 
										# print(zet)
										E=horiz_prev[zet][1][n][w] #E expresses cycles in terms of the basis of Hzet
										########
										#This list should only be made once? Or is it quicker like this, because we don't need to store it for every coloring?
										for j in range(horiz_deg[i-1][eps][n][w]): #go through basis of H eps
											#Go though B[:,j], look up simplices, check whether they live in zeta, reconstruct sum, express with E.
											Dj=np.zeros([len(U),1], int) #Dj is the differential of j in terms of the simplices
											for k,s in enumerate(S):
												if B[k,j]!=0 and (s in T): #i.e., s is a summand of the j-th basis element, and retains its weight in zet
													Dj[U.index(s)]=1 #Then this deserves a nonzero entry in Dj
													#But remember: Dj still only goes to the chains, not to homology!
											dj=E@Dj #this expresses Dj in terms of the chosen basis of Hzet
											for k,v in enumerate(dj):
												D[n][w][a+k,b+j]=v #The corresponding entries in the b+j-th column of D_n,w get made
											# print('\r\n Heres a differential: \r\n')
											# print(D[n][w][a:a+len(dj),b+j])
									a+=horiz_deg[i][zet][n][w] #increase row degree 
							b+=horiz_deg[i-1][eps][n][w] #increase column degree 
					


					# print('the indices are (i,n,w)='+str((i,n,w)))
					# print('The differential D^'+str(i-1)+'_('+str(n)+','+str(w)+') in degree '+str(i-1)+ ' is:')
					# print(D[n][w])
					# print('The basis change matrix is ')
					# print(Uprev[n][w])
					# print('r'+str((i,n,w))+'='+str(r[i,n,w]))
					# print('the corresponding differential is: ')
					# print(Dprev[n][w])

					#breakpoint()

				# Now, to homology: Since we have the differentials, this is now the same algorithm as before.
				# First, apply the basis change to D_i-1, and project to the kernel
				UD=(Uprev[n][w]@D[n][w])[r[i,n,w]:,:] 
				SNF=smith_normal_form(UD)
				 # The rank of the differential D^i-1:
				r[i-1,n,w]=SNF[4]
				# Can directly write U into Uprev, we used the previous one already
				Uprev[n][w]=SNF[0] 
				V=SNF[2]
				# Now record the uberhomology rank and basis:
				uber_rank[i,n,w]=horiz_degs[i][n,w]-r[i,n,w]-r[i-1,n,w]
				uber_basis[i][n][w]=Uiprev[n][w]@np.vstack((np.zeros([r[i,n,w],horiz_degs[i][n,w]-r[i,n,w]], int),V))[:,r[i-1,n,w]:] 
				#B[n][w]=Uiprev@np.vstack((np.zeros([r[n],k[n]-r[n]], int),V))[:,r[n+1]:] #This is a composition: First take V, then include into C_n, then act by Ui. 
				#Of this composite, the basis for H_n is then given by the columns after r_n+1

				#Now remember the inverse of U for next round
				Uiprev[n][w]=SNF[1]


				#else: #If the differential D^i-1 is 0, 
					#we could still have a nonzero kernel for D^i:

					#we still need the correct dimensions for Uprev and Uiprev:
				#	Uprev[n][w]=np.identity(horiz_degs[i-1][n,w], int) 
				#	Uiprev[n][w]=np.identity(horiz_degs[i-1][n,w], int)

									

				


		horiz_prev=horiz #shift H_i to H_i-1 
	
	
	#Finally, we need to consider i=0, where D^-1 is trivial. Hence, we just need ker D^0:
	for n in range(dim+1):
		for w in range(dim+2):
			uber_rank[0,n,w]=horiz_degs[0][n,w]-r[0,n,w] # The degree is just the total degree minus the rank
			uber_basis[0][n][w]=Uiprev[n][w][r[0,n,w]:,:] # The Basis is just the last columns of Uinverse

	###
	print('\r\n ***********************************************************\r\n')

	print('\r\n Calculation complete. \r\n')
	print('The uberhomology of the simplicial complex is the following:')
	for i in range(vert+1):
		print('\r\nIn homological degree '+str(i)+' (i.e., for colorings with this many black vertices),')
		print('the ranks of graded uberhomology are given as follows:')
		if i==0: print('(The matrix entry (n,w) corresponds to dimension n, weight w)')
		print(uber_rank[i])





def horizontal_homology(X,eps): #takes a simplicial complex X and a coloring eps, computes horizontal homology 
	#(Specifically, it returns lists of matrices B, E to compute the bases for horizontal homology, and a matrix h of degrees)
	dim=len(X)-1 #The dimension of the complex
	vert=len(X[0]) #The number of vertices
	#print('dim(X)='+str(dim)+', number of vertices='+str(vert))
	C=[[[] for _ in range(dim+2)] for _ in range(dim+1)] #creates a list (we'll need the ordering!) that's bigraded with respect to dimension and weight of the simplices
	#print('C has length (i.e., how far n goes) '+str(len(C))+' and width (i.e., how far w goes)' + str(len(C[0])))

	for n in range(len(X)): #Fill C with simplices 
		for s in X[n]: 
			w=n+1-len(set(s).intersection(eps)) #the weight of s is the number of 0-colored simplices, i.e. n+1 minus the number of 1-colored simplices
			C[n][w].append(frozenset(s)) #then add s to the simplices of dimension n, and weight w.
			#print('The simplex '+str(s)+' has dimension ' +str(n) +' and weight ' +str(w))


	diff=[[0 for _ in range(dim+2)] for _ in range(dim+1)] #creates an empty list for the differentials: There will be one starting at each C[n][w], 
	for n in range(1,dim+1): #Differentials from dimension 0 are 0
		for w in range(dim+2): #Now, write the differentials: These will be a matrix M:C[n][w]->C[n-1][w].
			b=len(C[n][w])
			a=len(C[n-1][w]) #Then M is an a*b-matrix
			M=np.zeros([a,b], int) #Create M with 0-entries
			for i in range(len(C[n][w])): #Go through simplices of C[n][w] (except for 0-simplices)
				s=set(C[n][w][i]) 
				for j in set(s).intersection(eps):
					t=set(s)
					#print('j: '+str(j))
					#print('t: '+str(t))
					t.remove(j)
					#print('t-j'+str(t))
					#print('C('+str(n-1)+')('+str(w)+') = '+str(C[n-1][w]))
					M[C[n-1][w].index(frozenset(t)),i]+=1#Find the j-face of s in the list of n-1-simplices, and increase the corresponding matrix entry.
			#print('The differential M('+str(n)+','+str(w)+') is '+ str(M)+'.')
			diff[n][w]=M

	#Now, compute homology: This is bigraded, and we can look at different weights separately.
	#These three variables will collect the results: 
	h=np.zeros([dim+1,dim+2], int)#h[n,w] is the degree of the n-th homology in weight w.
	B=[[None for _ in range(dim+2)]for _ in range(dim+1)] #B gives a basis for the homologies in terms of simplices
	E=[[None for _ in range(dim+2)]for _ in range(dim+1)] #E goes the other way: It expresses cycles in terms of the homology basis

	for w in range(dim+2): #Homology in weight-grading w
		#The crucial part is that we need the homology with a basis represented in chains. We'll achieve this by calculating an iterated Smith normal form:
		#That is, the differential D_n can be written as a product VDU, with V,U invertible and D diagonal with 1s in the top left, and 0s else.
		#Then the inverse of U gives us a basis for the kernel, by taking the last columns (where D is 0)
		#Also, D_n+1 goes to precisely those summands: We can replace D_n+1 by the projection of VD_n+1, and get a map to the kernel. 
		#Taking SNF for this modified D_n+1, we get a basis for H_n as the last columns of inv(U_n)*[[0],[V_n+1]].
		k=[len(C[n][w]) for n in range(dim+1)]# k[n] is the number of n-simplices of weight w, i.e., the dimension of C[n,w]
		r=[0] # r[n] will be the rank of D_n. The rank of D_0 is 0.
		Uprev=np.identity(k[0], int) #the previous U: For the first step, this is the identity, since all C_0 is the kernel of D_0.
		Uiprev=np.identity(k[0], int) #likewise for the previous U-inverse

		for n in range(0,dim): #Iteration that calculates H_n. For that, we need the SNF of D_n+1. (We kept the one for D_n) 
			#Also note that the highest dim is done seperately.
			# print('H_'+str(n)+','+str(w)+'(X)')
			# print('Uprev:')
			# print(Uprev)
			# print('D_'+str(n+1)+':')
			# print(diff[n+1][w])
			UD=(Uprev@diff[n+1][w])%2
			# print('This is UD, minus the zero rows:')
			D=UD[r[n]:,:] #Take every row from zeroD, starting from r_n (Might be inefficient)
			# print('This is supposed to be the remaining rows:')
			# print(D)

			SNF=smith_normal_form(D)#Calculates SNF for D: Get D=V diag_r U like above.
			#SNF is a tuple: 0 is U, 1 is inv(U), 2 is V, and 3 is the rank of D
			U=SNF[0]
			Ui=SNF[1]
			V=SNF[2]
			Vi=SNF[3]
			r.append(SNF[4]) #D_n+1 has rank r
			# print('U, V for D_'+str(n+1)+' are:')
			# print(U)
			# print(V)
			# print('The rank of D_'+str(n+1)+' is '+str(r[n+1]))

			#For the explicit basis of H_n, we need to chase it up to C_n:
			B[n][w]=Uiprev@np.vstack((np.zeros([r[n],k[n]-r[n]], int),V))[:,r[n+1]:] #This is a composition: First take V, then include into C_n, then act by Ui. 
			#Of this composite, the basis for H_n is then given by the columns after r_n+1
			#For the other direction, we take U and project C_n to the kernel of D_n, then take Vi and project to the cokernel of D_n.
			
			# print('V, V_i for D_'+str(n+1)+' are:')
			# print(V)
			# print(Vi)
			# print(V@Vi)

			# print('Uprev has dimensions '+str(Uprev.shape))#
			# print('Uprev[r[n]:,:] then has dimensions '+str(Uprev[r[n]:,:].shape))

			E[n][w]=(Vi@(Uprev[r[n]:,:]))[r[n+1]:,:] #The projections in the appropiate bases are then just selections of rows.
			
			#print('In degree ('+str(n)+','+str(w)+'), we have B and E as follows:')
			#print(B[n][w])

			#print(E[n][w])
			Uprev=np.copy(U) #Copies U and Ui into the previous variables
			Uiprev=np.copy(Ui)

		#For H_dim, just need the kernel of D_n:
		r.append(0) #D_dim+1 has degree 0
		# print('weight='+str(w))
		# print(Uiprev)
		B[dim][w]=Uiprev[:,r[dim]:] #The basis for H_n is given by the last k_n-r_n columns of the inverse of U
		E[dim][w]= Uprev[r[dim]:,:]#likewise, the other direction is obtained by taking the last k_n-r_n rows of U
		# print('In degree ('+str(dim)+','+str(w)+'), we have B and E as follows:')
		# print(B[dim][w])
		# print(E[dim][w])
		h[:,w]=np.array([k[n]-r[n]-r[n+1] for n in range(dim+1)]) #puts the degree of the homologies into h

	# print('The horizontal homology for the coloring '+str(eps)+' has the following ranks: ')
	# print(str(h))
	# print('(The row index is the dimension grading, the column index the weight grading.)')
	# print('The non-zero bases are given by:')
	# for i in range(len(B)):
	# 	for j in range(len(B[i])):
	# 		if B[i][j].shape[0]!=0 and B[i][j].shape[1]!=0:
	# 			print('The basis B_'+str(i)+','+str(j)+' is:')
	# 			print(B[i][j])
	return (B,E,h)

	#Should return h and a list of basis elements for nonzero homology: E.g.
	#B=[((0,0), ), ####->The homology in degree 0,0 has generators 
	#   ((1,3), )]

def input_coloring(v): #reads a coloring of v vertices through input
	print('Please specify a coloring of the '+str(v)+' vertices.')
	print('This should be in the form of a string of 0s and 1s. E.g., "011" colors three vertices with 0, 1, and 1.')
	i_eps=input()
	eps=set()
	for i in range(0,v): #eps is the set of vertices that are colored with 1.
		if i_eps[i]=='1':
			eps.add(i+1)
	print('The 1-colored vertices are '+ str(eps)+'.')
	return eps


def input_simplex():
	print('Warning: At the moment, there is no error catches. Please behave.')
	print('How many vertices does your simplicial complex have?')

	i_vertices=input()
	vertices=int(i_vertices) #number of vertices

	print('Your simplicial complex has the vertices 1,...,' + i_vertices + '.')

	print('Now we will record the simplices: To add a simplex, type in the list of vertices spanning it. ')
	print('For example, the input "1,3,5" adds a 2-simplex between these vertices, together with all its faces.')
	print('If you want to add no more simplices, just hit enter.')

	i_s=input() #input simplex

	i_simplices=set() #the set of all input simplices (not necessarily including faces)

	while i_s != '':
		t=set([int(x) for x in i_s.split(",")]) #make the new simplex into a set
		i_simplices.add(frozenset(t))
		print('The new simplex has dimension ' +  str(len(t)-1) + '. Insert next simplex, or hit enter to continue.')
		i_s=input()

	print('\r\n ***********************************************************\r\n')

	print('You have entered a total of '+ str(len(i_simplices)) + ' simplices.')

	simplices=[set() for x in range(vertices)] # The list of simplices, ordered by dimension: simplices[n] will be the set of n-simplices. Max dimension is the number of vertices - 1.

	dim=-1 #dimension of the simplicial set
	for s in i_simplices: #read in the input simplices
		n=len(s)-1 #n is the dimension of s, i.e., the number of vertices - 1.
		dim=max(dim, n) #update dimension
		simplices[n].add(frozenset(s)) #adds s as an n-simplex

	while len(simplices)>dim+1:
		simplices.pop() #remove the list of simplices larger than the given dimension

	simplices=add_faces(simplices, vertices) #add in all the faces

	print('All faces have been added. In total, you have:')
	for i in range(0,dim+1):
		print(str(len(simplices[i]))+ ' simplices of dimension '+ str(i)+ ', namely:')
		for j in simplices[i]:
			print(str(list(j)))

	print('\r\n ***********************************************************\r\n')


	return simplices

def smith_normal_form(M): #takes a matrix M, and calculates its smith normal form over F_2.
	#specifically, finds matrices U, V and D, such that M=VDU, U,V are invertible, and D is diagonal with the first r entries 1, and the rest 0.
	#The function returns the matrices U, V, their inverses, and the rank r as an integer. (We need the matrices for the homology algorithm, and to calculate bases)
	
	m=M.shape[0]
	n=M.shape[1]
	V=np.identity(m, int) #U, V and their inverses Ui, Vi start as identity matrices.
	Vi=np.identity(m,int)
	U=np.identity(n, int)
	Ui=np.identity(n,int)
	D=np.copy(M) #D starts out as M
	r=0 #r starts out as 0, gets increased in the loop

	
	#First, we do Gauß-Jordan elimination, to bring D into row echolon form and calculate U and Ui
	h=0 #pivot row
	k=0 #pivot column

	while h < m and k <n:
		#first, find nonzero elements in column k:
		nz=np.nonzero(D[h:,k])[0]+h #nz is the array of indices of nonzero elements starting from pivot index h:
		# print('nonzeros in column '+str(k)+' after the pivot:')
		# print(nz)
		# print(D)
		if len(nz)==0:
			k+=1 #no new pivot, look in the next column
		else:
			if nz[0]!=h:
				# print('swap rows '+str(nz[0])+' and '+str(h))
				D[[nz[0],h],:]=D[[h,nz[0]],:]#swap rows nz(0) and h in D	
				Vi[[nz[0],h],:]=Vi[[h,nz[0]],:]#swap rows nz(0) and h in Vi			
				V[:,[nz[0],h]]=V[:,[h,nz[0]]]#swap columns nz(0) and h in V

			for i in nz[1:]:#for all nonzero rows below h, i.e., the other entries of nz, clear out the row in D, and the column in V
				# print('add row'+str(h)+' to '+str(i))
				# print(nz)
				# #Note: For fields other than F_2, you need to add in a factor!
				# print('D,U,Ui before:')
				# print(D)
				# print(Ui)
				# print(U)

				D[i,:]=(D[i,:]+D[h,:])%2
				Vi[i,:]=(Vi[i,:]+Vi[h,:])%2
				V[:,h]=(V[:,h]+V[:,i])%2 #Note that we modify row h here! 
				# print('D,U,Ui after:')
				# print(D)
				# print(Ui)
				# print(U)

			#Now that we've cleared out the rows below, we'll clear out the columns to the right.
			nzr=np.nonzero(D[h,:])[0] #the indices of nonzero elements in row h
			for i in nzr[1:]: #the first one is our pivot k, all other need to be eliminated
				# print('k is '+str(k)+' and i is '+str(i))
				# print('D,V before:')
				# print(D)
				# print(V)
				D[h,i]=0 #Add column k to column i in D.
				#(As the k-column of D now only has one entry, we only have to consider the row h.)
				Ui[:,i]=(Ui[:,i]+Ui[:,k])%2 # Add column k to i in Ui.
				U[k,:]=(U[i,:]+U[k,:])%2 #subtract row i from row k in U.
				# print('D,V after:')
				# print(D)
				# print(V)

			#finally, swap the pivot column k to h:
			if h!=k:
				D[:,[h,k]]=D[:,[k,h]]
				Ui[:,[h,k]]=Ui[:,[k,h]]
				U[[h,k],:]=U[[k,h],:]
			h+=1
			k+=1

	#Note that h is now precisely the rank of our matrix M!

	# print('V:')
	# print(V)
	# print('D:')
	# print(D)
	# print('U:')
	# print(U)
	# print('VDU:')
	# print((V@D@U)%2) #should always be M
	# print('M:')
	# print(M)
	# print('UUi:')
	# print((U@Ui)%2)
	return (U, Ui, V, Vi, h)

def add_faces(initial_simplices, vertices): #takes a set of simplices, the maximal dimension and number of vertices, and adds in all the faces.
	simplices=initial_simplices
	dim=len(simplices)-1 #dimension of complex
	for n in range(dim, 1, -1): # start at maximal-dimension simplices, then go down the dimension
		#print('n='+str(n))
		for s in simplices[n]:#go through n-simplices
			for i in s: #adds all faces of s as n-1-simplices
				d=set(s)
				d.remove(i)
				simplices[n-1].add(frozenset(d))

	for i in range(1,vertices+1):
		simplices[0].add(frozenset({i}))

	return simplices


def new_vertex_coloring(oldcol): #takes all colorings of n vertices, gives all colorings of n+1 vertices.
	col=[l.copy() for l in oldcol] #Copy the old colorings (they exclude the new vertex)
	col.append([]) #Weight can be one larger than before
	n=len(oldcol[len(oldcol)-1][0]) #the number of old vertices, i.e., the size of the largest coloring 
	# print('oldcol, n are:')
	# print(str(oldcol))
	# print(n)
	for i in range(len(oldcol)):
		for eps in oldcol[i]: #For every old coloring eps with i black vertices
			col[i+1].append(frozenset(set(eps).union(set({n+1})))) #for all old colorings, add a new coloring that also includes 
	return col




if __name__=='__main__':
	main()

