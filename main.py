# Biblioteca feita para rodar em Python 3

import time
from math import sqrt
#from heap import *
import random
from copy import deepcopy

OUTPUT='output.txt' # arquivo onde eu colocarei os resultados

class Grafo:
	def __init__(self,entrada_txt,formato,com_direcao=False,imprimir_propriedades=True):
		assert formato in ['lista','matriz','ambos']
		self.entrada_txt=entrada_txt
		self.formato=formato

		arquivo=open(entrada_txt,'r')
		self.lista_geral=arquivo.read().splitlines()
		arquivo.close()

		self.qtd_vertices=int(self.lista_geral[0])
		
		self.output=open(OUTPUT,'a') # criando/abrindo arquivo no formato 'append byte'

		self.criar_grafos(formato,com_direcao)
		if imprimir_propriedades:
			self.imprimir_propriedades()

	def criar_grafos(self,formato,com_direcao):
		# criando grafo no formato lista de adjacências
		if self.formato in ['lista','ambos']:
			t=time.time()
			self.grafo_lista=[[] for i in range(self.qtd_vertices)]

			for v in self.lista_geral[1:]:
				# exemplo: v="1 2"
				v1,espaco,v2=v.partition(' ')

				try:
					v2,espaco,peso=v2.partition(' ')
					peso=float(peso)
					self.grafo_com_peso=True
				except:
					self.grafo_com_peso=False

				v1=int(v1)-1
				v2=int(v2)-1

				# utilizando técnica de inserção binária para inserir os elementos de forma muito mais rápida
				
				# adicionando na posição de v1
				lista=self.grafo_lista[v1]
				if len(lista):
					indice=self.indice_pra_insercao_binaria(lista,v2,0,len(lista)-1)
				else:
					indice=0

				if self.grafo_com_peso:
					lista.insert(indice,[v2,peso])
				else:
					lista.insert(indice,v2)

				if not com_direcao:
					# adicionando na posição de v2
					lista=self.grafo_lista[v2]
					if len(lista):
						indice=self.indice_pra_insercao_binaria(lista,v1,0,len(lista)-1)
					else:
						indice=0

					if self.grafo_com_peso:
						lista.insert(indice,[v1,peso])
					else:
						lista.insert(indice,v1)

			self.tempo_consumido_grafo_lista=time.time()-t

			self.memoria_consumida_grafo_lista=self.grafo_lista.__sizeof__()
			for v in self.grafo_lista:
				self.memoria_consumida_grafo_lista+=v.__sizeof__()
			self.memoria_consumida_grafo_lista=self.memoria_consumida_grafo_lista/(1024**2)
		else:
			self.tempo_consumido_grafo_lista=0
			self.memoria_consumida_grafo_lista=0

		# criando grafo no formato matriz de adjacências
		if self.formato in ['matriz','ambos']:
			t=time.time()
			self.grafo_matriz=[]
			linha=[0]*self.qtd_vertices
			for i in range(self.qtd_vertices):
				self.grafo_matriz.append(linha[:]) #ao fazer linha[:] em vez de apenas linha, eu obrigo o python a criar matrizes iguais a linha porém sem utilizar o mesmo endereço de memória (o que seria errado)

			for v in self.lista_geral[1:]:
				# exemplo: v="1 2"
				v1,espaco,v2=v.partition(' ')

				try:
					v2,espaco,peso=v2.partition(' ')
					peso=float(peso)
					self.grafo_com_peso=True
				except:
					self.grafo_com_peso=False

				v1=int(v1)-1
				v2=int(v2)-1

				if self.grafo_com_peso:
					valor=peso
				else:
					valor=1

				self.grafo_matriz[v1][v2]=valor
				self.grafo_matriz[v2][v1]=valor

			self.tempo_consumido_grafo_matriz=time.time()-t

			self.memoria_consumida_grafo_matriz=self.grafo_matriz.__sizeof__()
			for v in self.grafo_matriz:
				self.memoria_consumida_grafo_matriz+=v.__sizeof__()
			self.memoria_consumida_grafo_matriz=self.memoria_consumida_grafo_matriz/(1024**2)
		
		else:
			self.tempo_consumido_grafo_matriz=0
			self.memoria_consumida_grafo_matriz=0

	def indice_pra_insercao_binaria(self,lista_ordenada,alvo,inicio,fim):
		u'''É um método para inserção binária de elementos. Criei esse método para que, quando eu for utilizar a visualização dos grafos por listas, para cada vértice vi, seus adjacentes vj estejam ordenados, de modo a otimizar futuras operações com esta lista.'''
		if fim<=inicio:
			if self.grafo_com_peso:
				return inicio+1 if alvo>lista_ordenada[fim][0] else inicio
			else:
				return inicio+1 if alvo>lista_ordenada[fim] else inicio

		meio=int((inicio+fim)/2)

		if self.grafo_com_peso:
			encontrada=lista_ordenada[meio][0]
		else:
			encontrada=lista_ordenada[meio]

		if alvo>encontrada:
			return self.indice_pra_insercao_binaria(lista_ordenada,alvo,meio+1,fim)
		else:
			if alvo<encontrada:
				return self.indice_pra_insercao_binaria(lista_ordenada,alvo,inicio,meio-1)
			return meio+1

	def calcular_numero_arestas(self):
		qtd=0
		if self.formato in ['lista','ambos']:
			for i in self.grafo_lista:
				qtd+=len(i)

		else:
			for i in self.grafo_matriz:
				qtd+=sum(i)
		self.qtd_arestas=int(qtd/2)

	def criar_lista_graus(self):
		self.lista_graus=[0]*self.qtd_vertices
		if self.formato in ['lista','ambos']:
			for i in range(self.qtd_vertices):
				self.lista_graus[i]=len(self.grafo_lista[i])

		else:
			for i in range(self.qtd_vertices):
				self.lista_graus[i]=sum(self.grafo_matriz[i])
		self.grau_min=min(self.lista_graus)
		self.grau_max=max(self.lista_graus)
		self.grau_medio=sum(self.lista_graus)/self.qtd_vertices

		graus=self.lista_graus[:] # cópia
		graus.sort()
		if self.qtd_vertices%2: # qtd ímpar de vértices
			self.mediana_graus=graus[int(self.qtd_vertices/2)]
		else:
			self.mediana_graus=(graus[int(self.qtd_vertices/2)]+graus[int(self.qtd_vertices/2)-1])/2

	def ver_vizinhos(self,vertice):
		# Retorna a lista de vizinhos de um determinado vértice
		if self.formato in ['lista','ambos']:
			vizinhos=self.grafo_lista[vertice]
			vizinhos.sort()
			return vizinhos

		else:
			if self.grafo_com_peso:
				return [[i,v] for i,v in enumerate(self.grafo_matriz[vertice]) if v!=0]
			else:
				return [i for i,v in enumerate(self.grafo_matriz[vertice]) if v==1]

	def gerar_arvore_bfs(self,vertice_inicial,output=False):
		u'''Algoritmo:
		1.Desmarcar todos os vértices
		2.Definir fila "descobertos" inicialmente como vazia
		3.Marcar o "indice_inicial" e inserí-lo na fila "descobertos"
		4.Enquanto "descobertos" não estiver vazia:
			5.Retirar o vértice "v_descoberto" de "descobertos"
			6.Para todo vértice vizinho "v_vizinho" de v_descoberto faça
				7.Se v_vizinho não estiver marcado
					8.marcar v_vizinho
					9.inserir v_vizinho em "descobertos"
		
		#complexidades
		passo 1: O(n)
		passos de 4 a 9: O(m), onde m=qtd de vértices.
		
		resultado: O(n+m)
		'''
		t=time.time()
		vertice_inicial+=-1 # pois a notação do professor começa com o "1"

		# estou utilizando "set" em vez de lista ou outro contêiner porque "set" me diz mais rapidamente se um elemento está contido nele do que a lista
		marcados=set()
		marcados.add(vertice_inicial)

		descobertos=[vertice_inicial]
		caminho=[]
		
		# marquei com "-2" todas os vértices não analisados
		arvore_bfs=[-2]*self.qtd_vertices
		
		# marquei com "-1" o vértice inicial, pois é o que vou usar
		arvore_bfs[vertice_inicial]=-1

		camadas=[-1]*self.qtd_vertices
		camadas[vertice_inicial]=0
		camada_atual=0

		indice_descobertos=0

		while(len(descobertos)>indice_descobertos):
			v_descoberto=descobertos[indice_descobertos]
			for v_vizinho in self.ver_vizinhos(v_descoberto):
				if self.grafo_com_peso:
					v_vizinho=v_vizinho[0]
				if v_vizinho not in marcados:
					marcados.add(v_vizinho)
					descobertos.append(v_vizinho)
					arvore_bfs[v_vizinho]=v_descoberto+1
					camadas[v_vizinho]=camadas[v_descoberto]+1
			caminho.append(descobertos[indice_descobertos]+1)
			indice_descobertos+=1

		dt=time.time()-t
		if output:
				self.output.write('''
\t\tBFS
Tomando o vértice "%s" como ponto de partida, o caminho percorrido pelo algoritmo da BFS foi:
%s

Os pais de cada vértice são: 
%s

As camadas de cada vértice são:
%s
'''%(vertice_inicial+1,\
	str(caminho)[1:-1],\
	str(arvore_bfs)[1:-1],\
	str(camadas)[1:-1]))
		return caminho,arvore_bfs,camadas,dt

	def gerar_arvore_dfs(self,vertice_inicial,output=False):
		u'''Algoritmo:
		1.Desmarcar todos os vérticecs -- O(n)
		2.Definir pilha P com um elemento s
		3.Enquanto P não estiver vazia --------- passos 3 a 8 geram uma ordem igual a O(m), onde m=qtd de vértices.
			4.Remover v_descoberto de P // no topo da pilha
			5.Se v_descoberto não estiver marcado
				6.Marcar v_descoberto como descoberto
				7.Para cada aresta (v_descoberto,v_vizinho) incidente a v_descoberto
					8.Adicionar v_vizinho em P // no topo
		resultado: O(n+m)
		'''
		t=time.time()
		vertice_inicial+=-1 # pois a notação do professor começa com o "1"

		P=[vertice_inicial]
		Pset={vertice_inicial}

		caminho=[]
		caminho_set=set()

		camadas=[-1]*self.qtd_vertices
		camadas[vertice_inicial]=0

		# marquei com "-2" todos os vértices não analisados
		arvore_dfs=[-2]*self.qtd_vertices
		
		# marquei com "-1" o vértice inicial, pois é o que vou usar
		arvore_dfs[vertice_inicial]=-1

		concluidos=set()

		i=0
		while i>=0:
			v=P[i]

			if v in concluidos:
				P.pop()
				Pset.remove(v)
				i+=-1

			else:
				if self.grafo_com_peso:
					vizinhos=[b[0] for b in self.ver_vizinhos(v)]
				else:
					vizinhos=self.ver_vizinhos(v)

				if v+1 not in caminho_set:
					caminho.append(v+1)
					caminho_set.add(v+1)

				for viz in vizinhos:
					if viz not in concluidos:
						if viz not in Pset:
							P.append(viz)
							Pset.add(viz)
							camadas[viz]=camadas[v]+1
							arvore_dfs[viz]=v+1
							i+=1
							break

				if P[-1]==v:
					concluiu=True
					for viz in vizinhos:
						if viz not in concluidos and viz not in Pset:
							concluiu=False
							break
					if concluiu:
						concluidos.add(v)

					P.pop()
					Pset.remove(v)
					i+=-1

		dt=time.time()-t
		if output:
				self.output.write('''
\t\tDFS
Tomando o vértice "%s" como ponto de partida, o caminho percorrido pelo algoritmo da DFS foi:
%s

Os pais de cada vértice são: 
%s

As camadas de cada vértice são:
%s
'''%(vertice_inicial+1,\
	str(caminho)[1:-1],\
	str(arvore_dfs)[1:-1],\
	str(camadas)[1:-1]))
		return caminho,arvore_dfs,camadas,dt

	def componentes_conexas(self):
		t=time.time()
		# cc = componente conexa
		vertices_conexos=set()
		num_cc=0
		lista_cc=[]
		lista_tamanhos_cc=[]
		tamanho_maior_cc=0
		tamanho_menor_cc=self.qtd_vertices

		#Examina cada vértice do grafo
		for v in range(self.qtd_vertices):
			if v not in vertices_conexos:
				caminho,arvore_bfs,camadas,dt=self.gerar_arvore_dfs(v+1,output=False)

				lista_cc.append(caminho)
				num_cc+=1
				if len(caminho)>tamanho_maior_cc:
					tamanho_maior_cc=len(caminho)

				if len(caminho)<tamanho_menor_cc:
					tamanho_menor_cc=len(caminho)

				for i in caminho:
					vertices_conexos.add(i-1)
		dt=time.time()-t
		return num_cc,tamanho_maior_cc,tamanho_menor_cc,lista_cc,dt

	def diametro(self):
		t=time.time()
		diametro = 0

		#Checa se existe algum par de vértices que não esteja conectado; caso exista, o diâmetro é infinito
		num_cc=self.componentes_conexas()[0]
		if num_cc!=1:
			return -1,time.time()-t

		else:
			#Faz uma bfs para cada vértice do grafo
			for i in range(self.qtd_vertices):
				camadas=self.gerar_arvore_bfs(i+1,output=False)[2]

				#Checa o valor da maior camada, e o atribui como maior menor distância do vértice
				maior_distancia=max(camadas)

				#Se a maior menor distâcia deste vértice for maior que o diâmtro, este recebe o novo valor
				if maior_distancia>diametro:
					diametro=maior_distancia
		return diametro,time.time()-t

	def imprimir_propriedades(self):
		self.calcular_numero_arestas()
		self.criar_lista_graus()

		self.output.write('''
\n\n
###########################################
###########################################

\t%s

Memória consumida na construção dos grafos:
	Lista de adjacências: %.5f Mbs
	Matriz de adjacências: %.5f Mbs

Tempo consumido na construção dos grafos:
	Lista de adjacências: %.5f seg
	Matriz de adjacências: %.5f seg

Número de vértices: %s
Número de arestas: %s
Grau mínimo: %s
Grau máximo: %s
Grau médio: %s
Mediana de grau: %s
'''%(self.entrada_txt.upper(),\
	self.memoria_consumida_grafo_lista,\
	self.memoria_consumida_grafo_matriz,\
	self.tempo_consumido_grafo_lista,\
	self.tempo_consumido_grafo_matriz,\
	self.qtd_vertices,\
	self.qtd_arestas,\
	self.grau_min,\
	self.grau_max,\
	self.grau_medio,\
	self.mediana_graus))

	def dijkstra(self, vertice_inicial):

		for i in range(self.qtd_vertices):
			for peso in self.ver_vizinhos(i):
				if peso[1] < 0:
					print("Djikstra não pode ser realizado, pois o grafo possui pesos negativos")
					return

		vertice_inicial = vertice_inicial-1

		# Criação e inicialização do heap que será utilizado
		dist = BHeap()
		dist.buildHeap()

		#Preenche o heap com peso infinito para cada vértice do grafo; elementos do heap tem formato [peso,vertice]
		for vertice in range(self.qtd_vertices):
			dist.add([float('inf'),vertice])

		#Adiciona o indice inicial com peso 0 ao heap
		dist.add([0,vertice_inicial])

		S = [float('inf')]*self.qtd_vertices #Lista de distâncias de cada vértice até o vértice inicial
		S[vertice_inicial] = 0 #Distância do vértice até ele próprio definida como 0
		caminho = [[] for i in range(self.qtd_vertices)] #Lista de caminhos de cada vértice até o vértice inicial
		caminho[vertice_inicial].append(vertice_inicial+1) #Vértice inicial adicionado como primeiro elemento do caminho percorrido

		explorados = set() #Set que contém todos os vértices já explorados

		while len(dist) > 0:

			u = dist.remove() #Remove elemento de menor peso do heap

			#Checa se o elemento já foi explorado; caso não tenha sido, adiciona-se ele ao set explorados
			if u[1] not in explorados:
				explorados.add(u[1])

				#Para cada vizinho não explorado de u
				for v_vizinho in self.ver_vizinhos(u[1]):
					if v_vizinho[0] not in explorados:

						#Verifica se o peso atual do vizinho é maior que a soma do peso para chegar até o vértice u e ir deste até o vizinho;
						#em caso positivo, atualiza-se o peso atual do vizinho, adiciona-se o valor atualizado de sua distância ao heap
						#e seu caminho até o vértice inicial é atualizado com o caminho ao vértice u adicionado do vizinho
						if S[v_vizinho[0]] > u[0] + v_vizinho[1]:
							S[v_vizinho[0]] = u[0] + v_vizinho[1]
							dist.add([S[v_vizinho[0]], v_vizinho[0]])
							caminho[v_vizinho[0]] = caminho[u[1]] + [v_vizinho[0]+1]

		return S, caminho

	def dijkstra_caminho_minimo_distancia(self, vertice_inicial, vertice_desejado):

		valores = self.dijkstra(vertice_inicial)
		return valores[1][vertice_desejado-1], valores[0][vertice_desejado-1]

	def prim(self, vertice_inicial):

		vertice_inicial = vertice_inicial-1

		# Criação e inicialização do heap que será utilizado
		custo = BHeap()
		custo.buildHeap()

		#Preenche o heap com custo infinito para cada vértice do grafo; elementos do heap tem formato [custo,vertice]
		for vertice in range(self.qtd_vertices):
			custo.add([float('inf'),vertice])

		#Adiciona o indice inicial com custo 0 ao heap
		custo.add([0,vertice_inicial])

		caminho = [] #Caminho percorrido pela MST
		arvore_prim = [float('inf')]*self.qtd_vertices #Lista contendo o pai de cada vértice na MST
		peso = 0 #Peso total da MST

		distancias = [float('inf')]*self.qtd_vertices #Lista contendo as menores distâncias necessárias para chegar a cada vértice
		distancias[vertice_inicial] = 0

		arvore_prim[vertice_inicial] = -1
		explorados = set() #Set que contém todos os vértices já explorados

		while len(custo) > 0:
			u = custo.remove() #Remove do heap o elemento de menor custo

			#Caso o custo seja infinito, todos os elementos restando no heap possuem custo infinito; logo, quebra-se o loop
			if u[0] == float('inf'):
				break

			#Caso u não tenha sido explorado, adiciona-se ele aos explorados, ao caminho percorrido pela MST e seu peso ao peso total
			if u[1] not in explorados:
				explorados.add(u[1])
				caminho.append(u[1]+1)
				peso += u[0]

				for v_vizinho in self.ver_vizinhos(u[1]):
					if v_vizinho[0] not in explorados:

						#Verifica se a distância atual para chegar ao vértice é maior que a distância para chegar a este vértice partindo de u;
						#caso ela seja, sua distância é atualizada em distancias e no heap
						if distancias[v_vizinho[0]] > v_vizinho[1]:
							distancias[v_vizinho[0]] = v_vizinho[1]
							custo.add([distancias[v_vizinho[0]], v_vizinho[0]])

							if u[0] < arvore_prim[v_vizinho[0]]:
								arvore_prim[v_vizinho[0]] = u[1]+1

		return peso, caminho, arvore_prim

	def excentricidade(self, vertice):

		excentricidades = self.dijkstra(vertice)[0]
		return max(excentricidades)

	def maiores_graus_prim(self, vertice):

		arvore_prim = self.prim(vertice)[2] #Lista de pais de cada vértice da MST obtida por Prim
		lista_graus = [1]*self.qtd_vertices #Lista contendo o grau de cada vértice da MST; todos começam com grau 1 pois a função assume que todo vértice exceto a raiz da MST possui um pai
		lista_graus[vertice-1] = 0 #Reduz-se o grau da raiz por 1, dado que ela não possui pai na MST
		maiores_graus = [] #Lista contendo os maiores graus e seus respectivos vértices encontrados, no formato [vértice, grau]

		for i in range(len(arvore_prim)):

			#Cada vez que um vértice é considerado como pai na lista de vértices, soma-se 1 à seu grau total
			if arvore_prim[i] != float('inf') and arvore_prim[i] > -1:
				lista_graus[arvore_prim[i]-1] += 1

		#Encontra-se o vértice de maior grau da lista, que é então adicioná-lo a lista maiores_graus; seu valor em lista_graus é atualizado para 0
		#O processo se repete 3 vezes, coletando assim os 3 maiores graus
		for i in range(3):
			indice = lista_graus.index(max(lista_graus))
			maiores_graus.append([indice+1,lista_graus[indice]])
			lista_graus[indice] = 0

		return maiores_graus

	def veririfcar_biparticao(self):
		grupo1={0}
		grupo2=set()
		vertices_e_grupos=[[0,0]]
		while len(vertices_e_grupos):
			v,grupo=vertices_e_grupos[0]
			del vertices_e_grupos[0]

			for viz in self.ver_vizinhos(v):
				if grupo==0:
					if viz in grupo1:
						return False,False # não é bipartido
					elif viz not in grupo2:
						grupo2.add(viz)
						vertices_e_grupos.append([viz,1])
				else:
					if viz in grupo2:
						return False,False # não é bipartido
					elif viz not in grupo1:
						grupo1.add(viz)
						vertices_e_grupos.append([viz,0])
		if len(grupo1)+len(grupo2)==self.qtd_vertices:
			return grupo1,grupo2
		return False,False

	def capturar_emparelhamento_max(self):
		g1,g2=self.veririfcar_biparticao()
		if g1:
			grafo_dict=dict().fromkeys(g1,None)
			for v in g1:
				viz=self.ver_vizinhos(v)
				grafo_dict[v]=viz
			t0=time.time()
			emparelhamento_max=HopcroftKarp(grafo_dict).maximum_matching(keys_only=True)
			dt=time.time()-t0
			lista=[]
			for i,v in emparelhamento_max.items():
				if v!=None:
					lista.append([i,v])
			return len(lista),lista,dt
		return False,False,False

	def bellmanford(self,v,vf):
		g=Bellman(self.qtd_vertices)
		for v in range(self.qtd_vertices):
			for viz in self.ver_vizinhos(v):
				g.addEdge(v,viz[0],viz[1])
		to=time.time()
		dist,neg=g.BellmanFord(v)
		dt=time.time()-to
		return dist[vf],neg,dt



class HopcroftKarp(object):
    def __init__(self, graph):
        self._matching = {}
        self._dfs_paths = []
        self._dfs_parent = {}
        
        self._graph = deepcopy(graph)
        self._left = set(self._graph.keys())
        self._right = set()

        for value in self._graph.values():
            self._right.update(value)
        for vertex in self._left:
            for neighbour in self._graph[vertex]:
                if neighbour not in self._graph:
                    self._graph[neighbour] = set()
                    self._graph[neighbour].add(vertex)
                else:
                    self._graph[neighbour].add(vertex)

    def __bfs(self):
        layers = []
        layer = set()
        for vertex in self._left:  
            if vertex not in self._matching:  # confirma se o vértice está livre
                layer.add(vertex)
        layers.append(layer)
        visited = set()  # para acompanhar o vérices visitados 
        while True:
            # pegamos a camada mais recente no particionamento a cada repetição
            layer = layers[-1]
            new_layer = set()  # nova lista com as camadas subsequentes
            for vertex in layer:
                if vertex in self._left:
                    visited.add(vertex)
                    for neighbour in self._graph[vertex]:
                        if neighbour not in visited and (vertex not in self._matching or neighbour != self._matching[vertex]):
                            new_layer.add(neighbour)
                else:  
                    visited.add(vertex)  
                    for neighbour in self._graph[vertex]:
                        if neighbour not in visited and (vertex in self._matching and neighbour == self._matching[vertex]):
                            new_layer.add(neighbour)
            layers.append(new_layer)  # adicionamos a nova camada ao set de camadas
            # se a nova camada for vazia, temos que parar o while loop da BFS
            if len(new_layer) == 0:
                return layers   # break
            if any(vertex in self._right and vertex not in self._matching for vertex in new_layer):
                return layers  # break
                # break

    def __dfs(self, v, index, layers):
        if index == 0:
            path = [v]
            while self._dfs_parent[v] != v:
                path.append(self._dfs_parent[v])
                v = self._dfs_parent[v]
            self._dfs_paths.append(path)
            return True
        for neighbour in self._graph[v]:  
            if neighbour in layers[index - 1]:
                if neighbour in self._dfs_parent:
                    continue
                if (neighbour in self._left and (v not in self._matching or neighbour != self._matching[v])) or \
                        (neighbour in self._right and (v in self._matching and neighbour == self._matching[v])):
                    self._dfs_parent[neighbour] = v
                    if self.__dfs(neighbour, index-1, layers):
                        return True
        return False

    def maximum_matching(self, keys_only=False):
        while True:
            layers = self.__bfs()
            # Rompemos o loop while se a camada mais recente adicionada às camadas estiver vazia,
            # pois se não houver vértices na camada recente, não será possível encontrar caminhos aumentantes
            if len(layers[-1]) == 0:
                break
            free_vertex = set([vertex for vertex in layers[-1] if vertex not in self._matching])
            del self._dfs_paths[:]
            self._dfs_parent.clear()

            for vertex in free_vertex:  # O(m)
                self._dfs_parent[vertex] = vertex
                self.__dfs(vertex, len(layers)-1, layers)
                
            if len(self._dfs_paths) == 0:
                break
            
            for path in self._dfs_paths:
                for i in range(len(path)):
                    if i % 2 == 0:
                        self._matching[path[i]] = path[i+1]
                        self._matching[path[i+1]] = path[i]
        if keys_only:
            self._matching = {k:v for k,v in self._matching.items() if k in self._left}
        return self._matching

class Bellman: 
	def __init__(self, vertices): 
		self.V = vertices 
		self.graph = [] 
 
	def addEdge(self, u, v, w): 
		self.graph.append([u, v, w]) 
    
	def BellmanFord(self, src): 
		# Inicializa a distância de src para todos os outros vértices como infinito
		dist = [float("Inf")] * self.V 
		dist[src] = 0
		
		# Relaxar as arestas |V| - 1 vezes. 
		for i in range(self.V - 1):  
			for u, v, w in self.graph: 
				if dist[u] != float("Inf") and dist[u] + w < dist[v]: 
						dist[v] = dist[u] + w 

		
                # Checar se tem ciclo negativo
		neg=False
		for u, v, w in self.graph: 
			if dist[u]!=float("Inf") and dist[u] + w < dist[v]:
				neg=True
				break
		
		return dist,neg


# TRABALHO 3
open(OUTPUT,'w') # apenas para limpar o arquivo de resultados (caso ele já exista)
grafos=['grafo_teste_%s.txt'%i for i in range(7,11)]
#grafos = 'grafo_teste_10.txt'
grafos2=['ER_%s.txt'%i for i in [50,100,500,1000,1500]]


for arquivo in grafos:
	texto_output=''
	g=Grafo(arquivo,'lista',imprimir_propriedades=False)
	qtd,emp,dt=g.capturar_emparelhamento_max()
	texto_output+='\n\nGrafo %s'%arquivo
	if qtd:	
		texto_output+='\nQtd emparelhamentos: %s'%qtd
		texto_output+='\nTempo: %s seg'%dt
		#texto_output+='\nEmparelhamento: %s'%emp
	else:
		texto_output+='\nEste grafo não é bipartido!'

	g.output.write(texto_output)
	g.output.close()


for arquivo in grafos2:
	texto_output=''
	g=Grafo(arquivo,'lista',com_direcao=True,imprimir_propriedades=False)
	texto_output+='\n\nGrafo %s'%arquivo
	for v1,v2 in [[1,10],[2,20],[3,30]]:
		dist,neg,dt=g.bellmanford(v1-1,v2-2)
		texto_output+='\nCiclo Negativo: %s, Distância (%s,%s): %s, Tempo: %s'%(['Não','Sim'][int(neg)],v1,v2,dist,dt)

	g.output.write(texto_output)
	g.output.close()


u"""
# TRABALHO 2
grafos=['grafo_%s.txt'%i for i in [1,2,3,4,5]]

texto_output=''
for arquivo in grafos:
	g=Grafo(arquivo,'lista',imprimir_propriedades=False)

	# exercício 1
	for v in [10,20,30,40,50]:
		caminho,distancia=g.dijkstra_caminho_minimo_distancia(1,v)
		texto_output+='\nRelação entre os vértices 1 e %s:\n'%v
		texto_output+='Distância: %s\n'%distancia
		texto_output+='Caminho mínimo:\n%s\n'%caminho

	g.output.write(texto_output)
	# exercício 2
	texto_output+='\n\n'
	for v in [10,20,30,40,50]:
		texto_output+='\nA excentricidade do vértice %s é: %s'%(v,g.excentricidade(v))

	g.output.write(texto_output)
	# exercício 3
	vertices_aleatorios=random.sample(range(g.qtd_vertices),100) # pegando 100 vértices aleatórios
	qtd_repeticoes=1

	t1=time.time()
	for i in range(qtd_repeticoes):
		for v in vertices_aleatorios:
			g.excentricidade(v)
	t2=time.time()
	tempo_medio=(t2-t1)/qtd_repeticoes
	tempo_medio_por_vertice=tempo_medio/len(vertices_aleatorios)
	texto_output+='\n\nTempo médio para o cálculo da excentricidade de cada vértice: %.6f seg\n'%tempo_medio_por_vertice
		

	g.output.write(texto_output)
	# exercício 4
	t1=time.time()
	peso,caminho,arvore_prim=g.prim(1)
	t2=time.time()
	texto_output+='\n\nPrim começando pelo vértice 1:'
	texto_output+='\nPeso: %s'%peso
	texto_output+='\nCaminho: %s'%caminho
	texto_output+='\nÁrvore Prim: %s'%arvore_prim
	texto_output+='\nTempo para execução: %.6f seg'%(t2-t1)
	texto_output+='\n\n\n'

	g.output.write(texto_output)
	g.output.close()


texto_output='\n\n'
colaboradores=[]
with open('rede_colaboracao_vertices.txt','r',encoding="utf8") as f:
	lista=str(f.read()).splitlines()
	for i in lista:
		colaboradores.append(i.split(',')[1])

g=Grafo('rede_colaboracao.txt','lista',imprimir_propriedades=False)

nome_inicial='Edsger W. Dijkstra'
for i,v in enumerate(colaboradores):
	if v=='Edsger W. Dijkstra':
		v_inicial=i+1
		break

for nome in ['Alan M. Turing','J. B. Kruskal','Jon M. Kleinberg','Éva Tardos','Daniel R. Figueiredo']:
	for i,v in enumerate(colaboradores):
		if v==nome:
			w=i+1
			break
	
	caminho,distancia=g.dijkstra_caminho_minimo_distancia(v_inicial,w)
	caminho=[colaboradores[i-1] for i in caminho]
	texto_output+='\n\nSobre a relação entre %s e %s:\n'%(nome_inicial,nome)
	texto_output+='Distância: %s\n'%distancia
	texto_output+='Caminho: %s\n\n'%caminho
v_daniel=w

# capturando os três maiores vértices da MST
tres_maiores=g.maiores_graus_prim(1)
tres_maiores=[[colaboradores[i[0]-1],i[1]] for i in tres_maiores]
texto_output+='\n\nOs três vértices de maiores graus da MST com seus respectivos graus são:\n%s'%tres_maiores


peso,caminho,arvore_prim=g.prim(1)

vizinhos_edsger=[]
vizinhos_daniel=[]
for i,v in enumerate(arvore_prim):
	if v==v_inicial:
		vizinhos_edsger.append(colaboradores[i])

	if v==v_daniel:
		vizinhos_daniel.append(colaboradores[i])

texto_output+='\n\nOs vizinhos de %s na MST são:\n%s'%(nome_inicial,vizinhos_edsger)
texto_output+='\n\nOs vizinhos de Daniel R. Figueiredo na MST são:\n%s'%vizinhos_daniel


g.output.write(texto_output)
g.output.close()


# TRABALHO 1
grafos={
	'teste.txt':{
		'tipo': 'ambos',
		'qtd_vertices_teste': 1000,
		'output_arvores': True,
	},
	'as_graph.txt':{
		'tipo': 'lista',
		'qtd_vertices_teste': 20,
		'output_arvores': False,
	},
	'dblp.txt':{
		'tipo': 'lista',
		'qtd_vertices_teste': 1,
		'output_arvores': False,
	},

}

texto_output=''
for grafo,w in grafos.items():
	print(u'\nTrabalhando com o grafo de %s'%grafo.upper())
	g=Grafo(grafo,w['tipo'])
	g.gerar_arvore_bfs(1,output=w['output_arvores'])
	g.gerar_arvore_dfs(1,output=w['output_arvores'])

	# pegando 1000 vértices
	qtd_vertices_testar=w['qtd_vertices_teste']
	vertices=[i for i in range(1,1+g.qtd_vertices)]
	while len(vertices)<qtd_vertices_testar:
		vertices=vertices*5
	vertices=vertices[:qtd_vertices_testar]
	tempo_medio_bfs=0
	tempo_medio_dfs=0
	for v in vertices:
		a,b,c,t1=g.gerar_arvore_bfs(v,output=False)
		a,b,c,t2=g.gerar_arvore_dfs(v,output=False)
		tempo_medio_bfs+=t1
		tempo_medio_dfs+=t2

	tempo_medio_bfs=tempo_medio_bfs/qtd_vertices_testar
	tempo_medio_dfs=tempo_medio_dfs/qtd_vertices_testar

	texto_output+='\nO tempo médio para rodar a BFS foi de %.7f seg\n'%tempo_medio_bfs
	texto_output+='\nO tempo médio para rodar a DFS foi de %.7f seg\n\n\n'%tempo_medio_dfs

	# pais dos vértices 10, 20 e 30 quando iniciamos o caminho por 1, 2 e 3
	texto='\nOperação\t\t\tFilho\t\t\tVértice Inicial\t\t\tPai'
	for operacao in ['BFS','DFS']:
		for filho in [10,20,30]:
			if g.qtd_vertices>=filho:
				for inicio in [1,2,3]:
					if operacao=='BFS':
						caminho,arvore_bfs,camadas,dt=g.gerar_arvore_bfs(inicio,output=False)
						pai=arvore_bfs[filho-1]
					else:
						caminho,arvore_dfs,camadas,dt=g.gerar_arvore_dfs(inicio,output=False)
						pai=arvore_dfs[filho-1]
					texto+='\n%s\t\t\t\t\t%s\t\t\t\t%s\t\t\t\t\t\t%s'%(operacao,filho,inicio,pai)
	texto_output+=texto

	# calculando distâncias
	texto='\n\n'
	for p1,p2 in [[10,20],[10,30],[20,30]]:
		if g.qtd_vertices>p2:
			caminho,arvore_bfs,camadas,dt=g.gerar_arvore_bfs(p1,output=False)
			distancia=camadas[p2-1]
			texto+='\nDistância (%s,%s) = %s'%(p1,p2,distancia)
	texto_output+=texto

	# capturando as componentes conexas
	num_cc,tamanho_maior_cc,tamanho_menor_cc,lista_cc,dt=g.componentes_conexas()
	texto_output+='''

Número de componentes conexas: %s
Tamanho da maior componente conexa: %s
Tamanho da menor componente conexa: %s

'''%(num_cc,tamanho_maior_cc,tamanho_menor_cc)

	# calculando o diâmetro
	diametro,dt=g.diametro()
	if diametro==-1:
		texto='''
Diâmetro=Infinito
Tempo para efetuar o cálculo: %.7f seg
'''%dt
	else:
		texto='''
Diâmetro=%s
Tempo para efetuar o cálculo: %.7f seg
'''%(diametro,dt)
	texto_output+=texto

	g.output.write(texto)
	g.output.close()
"""
