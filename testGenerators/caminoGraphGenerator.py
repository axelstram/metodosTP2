import sys

desde = 10 #tam del grafo mas chico
hasta = 50 #tam grafo mas grande
step = 1 # cuanto salta entre n y n. si step=2 genera con n=2, n=4, n=6...

for n in range(desde,hasta, step):
	file = open('generatedTests/testCaminoGraph_'+str(n)+'_.in','w')
	file.write("# FromNodeId	ToNodeId \n")
	for i in range(1,n):
		file.write(str(i) + " " + str(i+1) + "\n")
