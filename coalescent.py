import numpy as np
import click
import time
import matplotlib.pyplot as plt

class Arrayset:
    array: np.ndarray
    length: int
    
    def __init__(self, array: np.ndarray):
        self.array = array
        self.length = len(array)
    
    def __len__(self) -> int:
        return self.length
    
    def __getitem__(self, i: int) -> int:
        return self.array[i]
        
    def remove(self, i: int):
        self.length -= 1
        self.array[i] = self.array[self.length]
    
    def append(self, i: int):
        if self.length >= len(self.array):
            raise IndexError("Arrayset is full")
        self.array[self.length] = i
        self.length += 1
        

class Coalescent:
    generation: int
    starting_times: np.ndarray
    lengths: np.ndarray
    lineages: Arrayset
    parents: np.ndarray
    num_descendants: np.ndarray
    N: int

    def __init__(self, N: int):
        self.N = N
        self.generation = 0
        self.starting_times = np.zeros(self.num_nodes, dtype=np.int64)
        self.lengths = np.zeros(self.num_nodes, dtype=np.int64)
        self.lineages = Arrayset(np.arange(2*N, dtype=np.int64))
        self.parents = -np.ones(self.num_nodes, dtype=np.int64)
        self.num_descendants = np.zeros(self.num_nodes, dtype=np.int64)
        self.num_descendants[:2*N] = 1

    @property
    def num_nodes(self) -> int:
        return 4*self.N - 1

    def remove_lineage(self, i: int, parent: int):
        lineage = self.lineages[i]
        self.parents[lineage] = parent
        self.lengths[lineage] = self.generation - self.starting_times[lineage]
        self.lineages.remove(i)
        self.num_descendants[parent] += self.num_descendants[lineage]
    
    def add_lineage(self, parent: int):
        self.lineages.append(parent)
        self.starting_times[parent] = self.generation

    def get_coalescent_time(self) -> np.ndarray:
        k = len(self.lineages)
        return 4*self.N / (k*(k-1))

    def simulate(self):
        for k in range(2*self.N-1):
            new_node = 2*self.N + k

            coalescent_time = self.get_coalescent_time()
            self.generation += coalescent_time

            i = np.random.randint(len(self.lineages))
            self.remove_lineage(i, new_node)
            j = np.random.randint(len(self.lineages))
            self.remove_lineage(j, new_node)
            self.add_lineage(new_node)

        assert np.sum(self.parents == -1) == 1
        assert self.num_descendants[-1] == 2*self.N

    def compute_sfs(self):
        result = np.zeros(2*self.N, dtype=np.int64)
        for i in range(self.num_nodes-1):
            length, num_descendants = self.lengths[i], self.num_descendants[i]
            result[num_descendants] += length
        return result

@click.command()
@click.option("-N", "--population-size", type=int, default=10)
def main(population_size: int):
    t = time.time()
    replicates = 100
    result = np.zeros(2*population_size, dtype=np.int64)
    for _ in range(replicates):
        coalescent = Coalescent(population_size)
        coalescent.simulate()
        result += coalescent.compute_sfs()
    
    result = result / np.sum(result)
    
    x = np.arange(1, 2*population_size)
    y = 1/x
    y = y / np.sum(y)
    plt.plot(x, np.cumsum(result[1:]), color="blue")
    plt.plot(x, np.cumsum(y), color="red")
    plt.show()
    
    

if __name__ == "__main__":
    main()
    
            
            
        

