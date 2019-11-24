# -*- coding: utf-8 -*-

class BHeap:

	def __init(self):
		self.heapList = [0] + []
		self.tamanhoAtual = 0

	def __len__(self):

		return len(self.heapList)-1

	def printHeap(self):

		print(self.heapList)

	def remove(self):

		menor_valor = self.heapList[1]
		self.heapList[1] = self.heapList[self.tamanhoAtual]
		self.tamanhoAtual -= 1
		self.heapList.pop()
		self.descer_item(1)

		return menor_valor

	def add(self, elemento):

		self.heapList.append(elemento)
		self.tamanhoAtual += 1
		self.subir_item(self.tamanhoAtual)

	def descer_item(self, i):

		while (i * 2) <= self.tamanhoAtual:
			mc = self.minChild(i)

			if self.heapList[i] > self.heapList[mc]:
				tmp = self.heapList[i]
				self.heapList[i] = self.heapList[mc]
				self.heapList[mc] = tmp

			i = mc

	def subir_item(self, i):

		while i // 2 > 0:
			if self.heapList[i] < self.heapList[i // 2]:
				tmp = self.heapList[i //2]
				self.heapList[i //2] = self.heapList[i]
				self.heapList[i] = tmp

			i = i // 2

	def minChild(self, i):

		if i * 2 + 1 > self.tamanhoAtual:
			return i * 2

		else:
			if self.heapList[i*2] < self.heapList[i*2+1]:
				return i * 2
			else:
				return i * 2 + 1

	def buildHeap(self):

		self.tamanhoAtual = 0
		self.heapList = [0] + []

if __name__ == '__main__':
	bh = BHeap()
	bh.buildHeap()
	bh.add([0,0])
	bh.add([6,3])
	bh.add([9,0])
	bh.add([5,1])
	bh.add([8,2])
	bh.add([3,5])
	bh.add([1,4])
	bh.add([0,12])

	bh.printHeap()

	print(bh.remove())
	print(bh.remove())
	print(bh.remove())
	print(bh.remove())
	print(bh.remove())
	print(bh.remove())
	print(bh.remove())