import pickle
import pandas as pd
import numpy as np
import networkx as nx
import tkinter as tk
from pyvis.network import Network

try:
    root = tk.Tk()
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
except Exception as e:
    screen_width = 1920
    screen_height = 1080
    print(f"Can not initialize display. Default screen_width={screen_width}, screen_height={screen_height}")


def make_edges(df: pd.DataFrame):
    df['edge'] = df.apply(lambda x: set(eval(x['to'])).difference(set(eval(x['from']))), axis=1)


def clear_system_str(s: str):
    s = s.replace('{', '').replace('}', '')
    s = s.replace('\'', '').replace('\'', '')
    s = s.replace('[', '').replace(']', '')
    return s


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
    log_df = pd.read_feather(log_file)
    make_edges(log_df)
    remove_braces(log_df)
    return log_df


def in_edges_count(G: nx.DiGraph):
    return [G.in_degree(n) for n in G.nodes]


def nodes_with_multiple_in_edges_and_count(G: nx.DiGraph):
    return list(filter(lambda nd: nd[1] > 1, zip(G.nodes, in_edges_count(G))))


def visualize_pyvis(log_file: str,
                    quad_systems_file: str,
                    output_html="quad.html",
                    width=int(screen_width * 0.8),
                    height=int(screen_height * 0.8)):
    df = get_df(log_file)

    quad_systems = pickle.load(open(quad_systems_file, 'rb'))
    quad_systems = set(map(lambda s: clear_system_str(str(s)), quad_systems))
    print(quad_systems)
    g = nx.from_pandas_edgelist(df, "from", "to", edge_attr="edge", create_using=nx.DiGraph)
    quad_edges = list(filter(lambda e: e[1] in quad_systems, g.edges))
    quad_edges_labels = list(map(lambda e: g.get_edge_data(*e)['edge'], quad_edges))

    g = nx.subgraph(g, get_processed_nodes(df))
    g = nx.DiGraph(g)
    for e, l in zip(quad_edges, quad_edges_labels):
        g.add_edge(*e, edge=l)
    for node, attributes in g.nodes.items():
        attributes['title'] = node
    g = nx.relabel_nodes(g, get_nodes_enumeration_in_process_order(df))
    nodes_in_edges_count = nodes_with_multiple_in_edges_and_count(g)

    print(f"Count of nodes with multiple parents: {len(nodes_in_edges_count)} / {len(g.nodes)} "
          f"= {np.round(len(nodes_in_edges_count) / len(g.nodes) * 100, 1)}%")
    print(nodes_in_edges_count)

    for (node, attributes), n_parents in zip(g.nodes.items(), in_edges_count(g)):
        attributes['n_parents'] = n_parents
        if n_parents > 1:
            attributes['color'] = '#ffff66'  # light yellow
        elif node == 0:
            attributes['color'] = '#32cd32'  # lime green
        else:
            attributes['color'] = '#87cefa'  # light sky bly

        if attributes['title'] in quad_systems:
            attributes['color'] = '#ff0000'  # red

    nt = Network(directed=True,
                 height=f"{height}px", width=f"{width}px",
                 heading="Quadratization algorithm visualization")
    nt.add_nodes(list(g.nodes.keys()),
                 title=[v['title'] for v in g.nodes.values()],
                 color=[v['color'] for v in g.nodes.values()])
    for f, t in g.edges:
        nt.add_edge(f, t, label=g.get_edge_data(f, t)['edge'], arrowStrikethrough=True)
    nt.show_buttons(filter_=['physics', 'layout'])
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
              "treeSpacing": 325,
              "sortMethod": "directed"
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
    nt.show(output_html)


if __name__ == '__main__':
    visualize_pyvis('log/log.feather', 'log/quad_systems.pkl')
