import sys
from graphviz import Digraph

comment = 'Arbol n-ario'
grafo = Digraph(comment=comment, format='png')
n = 15
archivo = '2_ario_treeGraph_'+str(n)+'_' # aca le ponen el nombre del archivo a graficar sin el .in

file = open('generatedTests/'+ archivo  + '.in','r')

for line in file:
	if (line[0] == '#'):
		continue
	edge = line.split() ## ej ['1','2']
	grafo.node(edge[0], edge[0]) # el primer parametro es el id del nodo, el segundo el nombre del nodo
	grafo.node(edge[1], edge[1]) # en nuestro caso es id == nombre. si hicieramos algo tipo Equipo 1 == River, se puede hacer que ponga River en el nodo
	grafo.edge(edge[0], edge[1]) # agrego el eje

grafo.render('./graficosGraphs/'+ archivo + '.gv', view=True)