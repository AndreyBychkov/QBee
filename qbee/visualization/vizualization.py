import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from typing import Tuple


def visualize(log_file: str,
              figsize: Tuple[float, float] = (20, 20),
              node_size: int = 50,
              edge_width: float = 0.2,
              edge_label_font_size=5,
              ) -> None:
    log_df = pd.read_csv(log_file)
    start_node, end_node, substitutions = log_df.tail(1).values[0]
    substitutions = substitutions.strip('[]').split(', ')
    log_df.drop(log_df.tail(1).index, inplace=True)

    g = nx.from_pandas_edgelist(log_df, 'from', 'name', create_using=nx.DiGraph, edge_attr='substitution')
    pos = nx.spring_layout(g, scale=5, k=2)

    nodes = list(g.nodes)
    node_color = [10] * len(nodes)
    node_color[nodes.index(start_node)] = 0
    node_color[nodes.index(end_node)] = 20

    plt.figure(figsize=figsize)
    nx.draw(g, pos=pos, node_size=node_size, width=edge_width, node_color=node_color)

    shortest_path = nx.shortest_path(g, start_node, end_node)
    shortest_path_edges = list(zip(shortest_path, shortest_path[1:]))
    nx.draw_networkx_edges(g, pos, edgelist=shortest_path_edges, edge_color='r', width=edge_width * 2.5, arrowsize=10 * edge_width * 2.5)

    edge_labels = dict([((fr, to), data['substitution']) for fr, to, data in g.edges(data=True)])
    nx.draw_networkx_edge_labels(g, pos, edge_labels=edge_labels, font_size=edge_label_font_size)

    plt.show()


if __name__ == '__main__':
    visualize('log.csv')
