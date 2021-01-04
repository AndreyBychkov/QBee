import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import holoviews as hv
import tkinter as tk
from bokeh.plotting import show
from pyvis.network import Network

root = tk.Tk()

screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

hv.extension('bokeh')


def make_edges(df: pd.DataFrame):
    df['edge'] = df.apply(lambda x: set(eval(x['to'])).difference(set(eval(x['from']))), axis=1)


def remove_braces(df: pd.DataFrame):
    df['edge'] = df['edge'].apply(lambda x: str(x).replace('{', '').replace('}', ''))
    df['edge'] = df['edge'].apply(lambda x: str(x).replace('\'', '').replace('\'', ''))
    df['from'] = df['from'].apply(lambda x: str(x).replace('[', '').replace(']', ''))
    df['from'] = df['from'].apply(lambda x: str(x).replace('{', '').replace('}', ''))
    df['from'] = df['from'].apply(lambda x: str(x).replace('\'', '').replace('\'', ''))
    df['to'] = df['to'].apply(lambda x: str(x).replace('[', '').replace(']', ''))
    df['to'] = df['to'].apply(lambda x: str(x).replace('{', '').replace('}', ''))
    df['to'] = df['to'].apply(lambda x: str(x).replace('\'', '').replace('\'', ''))


def get_df(log_file: str):
    log_df = pd.read_csv(log_file)
    make_edges(log_df)
    remove_braces(log_df)
    return log_df


def visualize_pyvis(log_file: str):
    df = get_df(log_file)
    g = nx.from_pandas_edgelist(df, "from", "to", edge_attr="edge", create_using=nx.DiGraph)
    for node, attributes in g.nodes.items():
        attributes['title'] = node
    g = nx.relabel_nodes(g, dict(zip(g.nodes.keys(), range(len(g.nodes)))))

    nt = Network()
    nt.add_nodes(list(g.nodes.keys()), title=[v['title'] for v in g.nodes.values()])
    for (f, t), l in zip(g.edges, df['edge']):
        nt.add_edge(f, t, label=l)
    nt.show('quad.html')


def visualize_bokeh(log_file: str):
    df = get_df(log_file)
    g = hv.Graph(df)
    g = g.relabel('Directed Graph').opts(directed=True,
                                         node_size=10,
                                         arrowhead_length=0.01,
                                         width=int(screen_width * 0.8),
                                         height=int(screen_height * 0.8), )
    show(hv.render(g))


if __name__ == '__main__':
    visualize_pyvis('../log/log.csv')
