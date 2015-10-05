import sys

desde = 10 #tam del grafo mas chico
hasta = 50 #tam grafo mas grande
step = 1 # cuanto salta entre n y n. si step=2 genera con n=2, n=4, n=6...

for n in range(desde,hasta, step):
	file = open('generatedTests/treeGraph_'+str(n)+'_.in','w')
	file.write("# FromNodeId	ToNodeId \n")
	for j in range(2,n+1):
		file.write(str(j) + " 1 \n")
