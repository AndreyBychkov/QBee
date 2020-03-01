import networkx as nx
import matplotlib.pyplot as plt

if __name__ == '__main__':
    g = nx.Graph()
    g.add_node("A")
    g.add_node("B")
    g.add_node("C")

    g.add_edges_from([("A", "B"), ("A", "C")])
    nx.draw(g, with_labels=True)
    plt.show()