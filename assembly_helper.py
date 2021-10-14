import networkx as nx
import pandas as pd
import copy
from collections import Counter
import numpy as np

def to_adj(T):
    try:
        return pd.DataFrame(nx.adjacency_matrix(T).todense(),index=T.nodes(),columns=T.nodes())
    except:
        print("Cannot convert to adjacency matrix")
    return None

def show(T):
    T = copy.deepcopy(T)
    width_dict = Counter(T.edges())
    edge_width = [ (u, v, {'width': value}) 
                  for ((u, v), value) in width_dict.items()]
    
    G_new = nx.DiGraph()
    G_new.add_edges_from(edge_width)
    pos=nx.kamada_kawai_layout(G_new)
    #pos=nx.spring_layout(G_new)
    nx.draw(G_new, pos)
    edge_labels=dict([((u,v,),d['width'])
                 for u,v,d in G_new.edges(data=True)])
    
    nx.draw(G_new,pos,with_labels=True)
    nx.draw_networkx_edges(G_new, pos=pos)
    nx.draw_networkx_edge_labels(G_new, pos, edge_labels=edge_labels,
                                 label_pos=0.55, font_size=10)
    
def eulerian_cycle(G,start=None):
    # YOUR SOLUTION HERE
    cycle = None
    ## BEGIN SOLUTION
    edges = {}
    cnt = 0
    for u,v in G.edges():
        if u not in edges:
            edges[u] = []
        edges[u].append(v)
        cnt += 1
    if start is None:
        start = list(edges.keys())[0]
    cycle = [start]
    current = start
    ecnt = 0
    while ecnt < cnt:
        while True:
            neighbors = edges[current]
            if len(neighbors) > 0:
                break
            # we need to keep shifting until we get to a point that has unused portions
            cycle = [cycle[-2]] + cycle[0:-2] + [cycle[-2]]
            current = cycle[-1]
            
        current = neighbors.pop()
        cycle.append(current)
        ecnt += 1
    ## END SOLUTION
    return cycle

def de_bruijn(patterns):
    dB = nx.MultiDiGraph()
    # dB.add_edge("AA","AT") # sample edge in case you want to run the code without implementing your solution
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    k = len(patterns[0])
    for kmer in patterns:
        prefix = kmer[:(k-1)]
        suffix = kmer[1:]
        dB.add_edge(prefix,suffix)
    ## END SOLUTION
    return dB

def calc_in_out(G):
    in_deg = {}
    out_deg = {}
    for u,v in G.edges():
        if v not in in_deg:
            in_deg[v] = 0
        if u not in out_deg:
            out_deg[u] = 0
        in_deg[v] += 1
        out_deg[u] += 1
    in_out = pd.Series(in_deg,name="in").to_frame().join(pd.Series(out_deg,name="out").to_frame(),how='outer')
    return in_out.fillna(0).astype(int)

def eulerian_path(G):
    # YOUR SOLUTION HERE
    path = []
    ## BEGIN SOLUTION
    in_out = calc_in_out(G)
    diff = in_out["out"] - in_out["in"]
    end = list(in_out.index[diff < 0])[0]
    diff = in_out["in"] - in_out["out"]
    start = list(in_out.index[diff < 0])[0]
    G2 = copy.deepcopy(G)
    G2.add_edge(end,start)
    cycle = eulerian_cycle(G2)
    while not (cycle[0] == start and cycle[-2] == end):
        cycle = [cycle[-2]] + cycle[0:-2] + [cycle[-2]]
    path = cycle[:-1]
    ## END SOLUTION
    return path

def reconstruct(kmers):
    dB = de_bruijn(kmers)
    path = eulerian_path(dB)
    text = ""
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    prefix = path.pop(0)
    text = [prefix]
    while len(path) > 0:
        suffix = path.pop(0)
        text.append(suffix[-1])
    text = "".join(text)
    ## END SOLUTION
    return text

def read_fasta(file):
    seqs = []
    headers = []
    # implement a read fasta file
    # YOUR SOLUTION HERE
    ## BEGIN SOLUTION
    f = open(file)
    header = None
    sequence = []
    for line in f:
        line = line.strip()
        if line[0] == ">":
            if header is not None:
                headers.append(header)
                seq = "".join(sequence)
                seqs.append(seq.strip())
                sequence = []
            header = line
        else:
            sequence.append(line.strip())
    if header is not None:
        headers.append(header)
        seq = "".join(sequence)
        seqs.append(seq.strip())
        sequence = []
    ## END SOLUTION    
    return headers,seqs

def align_dynamic2(s1,s2,match_score=1,mismatch_score=0,gap_score=0,verbose=False):
    scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    aligned = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    for s2_part in scores.columns:
        scores.loc["-",s2_part] = 0
        if s2_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc["-",s2_part] = ("".join(["-" for i in range(len(s2_part))]),s2_part)
    for s1_part in scores.index:
        scores.loc[s1_part,"-"] = 0
        if s1_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc[s1_part,"-"] = (s1_part,"".join(["-" for i in range(len(s1_part))]))
    if verbose:
        display(aligned)
    
    nrows,ncols = scores.shape
    for i in range(1,nrows):
        for j in range(1,ncols):
            # What are our three options
            opt1_s1 = scores.index[i-1] # remember the rows are representative of s1
            opt1_s2 = scores.columns[j-1] # remember the columns are representative of s2
            score_opt1 = -np.Inf # FIX THIS!
            s1_aligned_opt1 = "" # FIX THIS!
            s2_aligned_opt1 = "" # FIX THIS!
            ## BEGIN SOLUTION
            if scores.index[i][-1]==scores.columns[j][-1]: 
                score = match_score
            else:
                score = mismatch_score
            score_opt1 = scores.loc[opt1_s1,opt1_s2] + score
            s1_aligned_opt1 = aligned.loc[opt1_s1,opt1_s2][0] + scores.index[i][-1]
            s2_aligned_opt1 = aligned.loc[opt1_s1,opt1_s2][1] + scores.columns[j][-1]
            ## END SOLUTION
            
            opt2_s1 = scores.index[i-1]
            opt2_s2 = scores.columns[j]
            score_opt2 = -np.Inf # FIT THIS!
            s1_aligned_opt2 = "" # FIX THIS!
            s2_aligned_opt2 = "" # FIX THIS!
            ## BEGIN SOLUTION
            score_opt2 = scores.loc[opt2_s1,opt2_s2]+gap_score
            s1_aligned_opt2 = aligned.loc[opt2_s1,opt2_s2][0] + scores.index[i][-1]
            s2_aligned_opt2 = aligned.loc[opt2_s1,opt2_s2][1] + "-"
            ## END SOLUTION
            
            opt3_s1 = scores.index[i]
            opt3_s2 = scores.columns[j-1]
            score_opt3 = -np.Inf # FIT THIS!
            s1_aligned_opt3 = "" # FIX THIS!
            s2_aligned_opt3 = "" # FIX THIS!
            ## BEGIN SOLUTION
            score_opt3 = scores.loc[opt3_s1,opt3_s2]+gap_score
            s1_aligned_opt3 = aligned.loc[opt3_s1,opt3_s2][0] + "-"
            s2_aligned_opt3 = aligned.loc[opt3_s1,opt3_s2][1] + scores.columns[j][-1]
            ## END SOLUTION
            
            scores.loc[scores.index[i],scores.columns[j]] = max(score_opt1,score_opt2,score_opt3)
            if max(score_opt1,score_opt2,score_opt3) == score_opt1:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt1,s2_aligned_opt1)
            elif max(score_opt1,score_opt2,score_opt3) == score_opt2:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt2,s2_aligned_opt2)
            else:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt3,s2_aligned_opt3)
    if verbose:
        display(scores)
        display(aligned)
    return scores.loc[s1,s2],aligned.loc[s1,s2][0],aligned.loc[s1,s2][1]

def print_alignment(aligned_s1_,aligned_s2_,num_to_print=100):
    chunks_s1 = [aligned_s1_[i:i+num_to_print] for i in range(0, len(aligned_s1_), num_to_print)]
    chunks_s2 = [aligned_s2_[i:i+num_to_print] for i in range(0, len(aligned_s2_), num_to_print)]

    for aligned_s1,aligned_s2 in zip(chunks_s1,chunks_s2):
        for i in range(len(aligned_s1)):
            print(aligned_s1[i],end="")
        print()
        for i in range(len(aligned_s1)):
            if aligned_s1[i] == aligned_s2[i]:
                print("|",end="")
            else:
                print(" ",end="")
        print()
        for i in range(len(aligned_s2)):
            print(aligned_s2[i],end="")
        print()