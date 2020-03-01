import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

if __name__ == '__main__':
    log_df = pd.read_csv('log.csv')
    g = nx.from_pandas_edgelist(log_df, 'from', 'name', create_using=nx.DiGraph, edge_attr='replacement')
    pos = nx.spring_layout(g)
    print(g.edges)

    nx.draw(g, pos=pos)

    edge_labels = dict([((fr, to), data['replacement']) for fr, to, data in g.edges(data=True)])
    nx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels)
    plt.show()