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


def get_processed_nodes(df: pd.DataFrame):
    return pd.unique(df['from'])


def get_nodes_enumeration_in_process_order(df: pd.DataFrame) -> dict:
    return dict([(v, k) for (k, v) in enumerate(df['from'].unique())])


def get_df(log_file: str):
    log_df = pd.read_csv(log_file)
    make_edges(log_df)
    remove_braces(log_df)
    return log_df


def in_edges_count(G: nx.DiGraph):
    return [G.in_degree(n) for n in G.nodes]


def nodes_with_multiple_in_edges_and_count(G: nx.DiGraph):
    return list(filter(lambda nd: nd[1] > 1, zip(G.nodes, in_edges_count(G))))


def visualize_pyvis(log_file: str,
                    width=int(screen_width * 0.8),
                    height=int(screen_height * 0.8)):
    df = get_df(log_file)
    g = nx.from_pandas_edgelist(df, "from", "to", edge_attr="edge", create_using=nx.DiGraph)
    g = nx.subgraph(g, get_processed_nodes(df))
    for node, attributes in g.nodes.items():
        attributes['title'] = node
    g = nx.relabel_nodes(g, get_nodes_enumeration_in_process_order(df))

    for (node, attributes), n_parents in zip(g.nodes.items(), in_edges_count(g)):
        attributes['n_parents'] = n_parents
        if n_parents > 1:
            attributes['color'] = '#ffff66'  # light yellow
        elif node == 0:
            attributes['color'] = '#32cd32'  # lime green
        else:
            attributes['color'] = '#87cefa'  # light sky bly

    nt = Network(directed=True,
                 height=f"{height}px", width=f"{width}px",
                 heading="Quadratization algorithm visualization")
    nt.add_nodes(list(g.nodes.keys()),
                 title=[v['title'] for v in g.nodes.values()],
                 color=[v['color'] for v in g.nodes.values()])
    for (f, t), l in zip(g.edges, df['edge']):
        nt.add_edge(f, t, label=l, arrowStrikethrough=True)
    nt.show_buttons(filter_=['physics', 'layout'])
    options = """
    var options = {
      "nodes": {
        "borderWidth": 2,
        "borderWidthSelected": 3
      },
      "edges": {
        "color": {
          "inherit": true
        },
        "smooth": false
      },
      "layout": {
        "hierarchical": {
          "enabled": true,
          "levelSeparation": 315,
          "nodeSpacing": 245,
          "treeSpacing": 325
        }
      },
      "interaction": {
        "keyboard": {
          "enabled": true
        },
        "navigationButtons": true,
        "tooltipDelay": 100
      },
      "physics": {
        "hierarchicalRepulsion": {
          "centralGravity": 0,
          "springLength": 180,
          "springConstant": 0.15,
          "nodeDistance": 380
        },
        "minVelocity": 0.75,
        "solver": "hierarchicalRepulsion"
      }
    }
    """
    nt.set_options(r"""
    var options = {
        "configure": {
            "enabled": true
        },
        "edges": {
            "color": {
                "inherit": true
            },
            "smooth": {
                "enabled": false,
                "type": "continuous"
            }
        },
        "layout": {
            "hierarchical": {
              "enabled": true,
              "levelSeparation": 315,
              "nodeSpacing": 245,
              "treeSpacing": 325
            }
        },
        "interaction": {
            "dragNodes": true,
            "hideEdgesOnDrag": false,
            "hideNodesOnDrag": false,
            "keyboard": {
                "enabled": true
            },
            "navigationButtons": true,
            "tooltipDelay": 100
        },
        "physics": {
            "enabled": true,
            "stabilization": {
                "enabled": true,
                "fit": true,
                "iterations": 1000,
                "onlyDynamicEdges": false,
                "updateInterval": 50
            },
            "hierarchicalRepulsion": {
                "centralGravity": 0,
                "springLength": 180,
                "springConstant": 0.15,
                "nodeDistance": 380
            },
            "minVelocity": 0.75,
            "solver": "hierarchicalRepulsion"
        }
    }
    """)
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
