import argparse
import subprocess
import os

import sys
import math
from collections import defaultdict
import pandas as pd
logs = sys.stderr

preamble = r'''
\documentclass{standalone}
%\usepackage{fullpage}
\pagestyle{empty}
\usepackage{tikz}
\usepackage{tkz-euclide}
\usepackage{siunitx}
\usetikzlibrary{shapes, shapes.multipart}
\usepackage{xcolor}

\usepackage{verbatim}
\usepackage{lipsum}
\begin{document}
'''

picturepre = r'''
%\hspace{-3cm}
%\resizebox{1.2\textwidth}{!}{
\begin{tikzpicture}
'''

# dataset = "."
MAXLEN = 5650
MINLEN = 0
circular = True

lbs = [
    '(',
    '[',
    '{',
    '<'
]

rbs = [
    ')',
    ']',
    '}',
    '>'
]

nuc2color = {'A': 'black', 'C': 'brown', 'G': 'green', 'U': 'purple'}

counter_clockwise = False # hzhang
rotate = 180 -5

def drawarc_clockwise(a, b, deg, style, length, lengthfix): ## counterclock wise
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    angle_a = angle_a-70 + rotate
    angle_b = angle_b-70 + rotate

    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    if alpha < 170:
        print("\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_b,
                                                                                                    angle_b + 90,
                                                                                                    angle_a + 270,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(palpha/2.),
                                                                                                    a, b, length), file=output
        )
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        print("\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_a,
                                                                                                    angle_a + 90,
                                                                                                    angle_b - 90,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(pbeta/2.),
                                                                                                    a, b, length), file=output
        )
    else:
        print("\\draw[line width = 0.0mm] %s (%d) to [bend left=%.1f] (%d);" % (style, 
                                                                              a, 
                                                                              2*deg,
                                                                              b), file=output)

def drawarc_counterclockwise(a, b, deg, style, length, lengthfix):
    a,b = b,a
    angle_a = 360./(length+lengthfix)*a
    angle_b = 360./(length+lengthfix)*b
    angle_a = 450 - angle_a -20
    angle_b = 450 - angle_b -20

    alpha = (angle_b-angle_a)
    pa = angle_a * math.pi / 180.
    pb = angle_b * math.pi / 180.
    palpha = alpha * math.pi / 180.
    if alpha < 170:
        print("\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_b,
                                                                                                    angle_b + 90,
                                                                                                    angle_a + 270,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(palpha/2.),
                                                                                                    a, b, length), file=output
        )
    elif alpha > 190:
        beta = 360. - alpha
        pbeta = beta * math.pi / 180.
        print("\\draw[line width = 0.0mm] %s ([shift=(%f:10cm)]0,0) arc (%f:%f:%fcm); %% %d %d %d" % (style, angle_a,
                                                                                                    angle_a + 90,
                                                                                                    angle_b - 90,
                                                                                                    # 10*math.sin(palpha/2.) / math.sin(math.pi/2. - palpha/2.),
                                                                                                    10 * math.tan(pbeta/2.),
                                                                                                    a, b, length), file=output
        )
    else:
        print("\\draw[line width = 0.0mm] %s (%d) to [bend left=%.1f] (%d);" % (style, 
                                                                              a, 
                                                                              2*deg,
                                                                              b), file=output)

        
def draw_eclipse(left, right, style):
    radius = (right - left)/2
    print(f"\\draw{style} ({right}, 0) arc (0:180:{radius} and {radius/2});", file=output)
    
    
    
def agree(pres, pref, a, b): ## pres[a] = b
    if pref[a] == b:
        return True
    elif pref.get(a-1,-1) == b or pref.get(a+1,-1) == b:
        return True
    elif pref.get(b-1,-1) == a or pref.get(b+1,-1) == a:
        return True
    else:
        return False

# def pair_agree(pres, pref, a, b):
#     if pres[a] == b or pres.get(a-1,-1) == b or pres.get(a+1,-1) == b or pres.get(b-1,-1) == a or pres.get(b+1,-1) == a:
#         return True
#     else:
#         return False

def pair_agree(a, b, structure):
    return structure[a-1]=="(" and structure[b-1]==")"

def get_pairs(ss): # find the pairs in a secondary structure, return a dictionary
    assert len(ss) > MINLEN
    pairs = []
    stack = []
    for i, s in enumerate(ss):
        if s==".":
            pass
        elif s=="(":
            stack.append(i)
        elif s==")":
            j = stack.pop()
            assert j < i
            pairs.append((j, i))
        else:
            raise ValueError(f'the value of structure at position: {i} is not right: {s}!')
    return pairs


def draw_bpp(seq, ss, tex_file, bpp_file):
    global output
    output=open(tex_file, 'w')
    print(preamble, file=output)
    # bpp_dict = defaultdict(int)
    if bpp_file is not None:
        bpp_dict = []
        for line in open(bpp_file).readlines():
            if line.startswith(">"): continue
            line = line.strip().split()
            if len(line) != 3: continue
            i, j, prob = int(line[0]), int(line[1]), float(line[2])

            if prob > 1e-4: # threshold to avoid "Dimension too large" error
                # bpp_dict[i,j] = prob
                if prob > 1.0:
                    prob = 1.0
                bpp_dict.append((i, j, prob))
        bpp_dict.sort(key=lambda x: x[2])
    
    length = len(seq)
    bases = [''] + list(seq)

    if length > MAXLEN:
        print("%s too long (%d)" % (filename, length), file=logs)
        # continue # too long    

    if length < MINLEN:
        print("%s too short (%d)" % (filename, length), file=logs)
        # continue # too short    

    print(picturepre, file=output)


    lengthfix = int(length/9.0)
    gap = 5
    for i, nuc in enumerate(seq):
        xloc = i
        print(f"\\filldraw [gray, fill={nuc2color[nuc]}] ({xloc}, 0) circle (10pt); ", file=output)
        print(f"\\node[below=0.5cm, {nuc2color[nuc]}, very thick, scale=1.2] at ({xloc},0) {{\Huge \\bf {nuc}}};", file=output)
        if (i+1)%gap == 0 or i==0:
            print(f"\\node[below=0.5cm, black, scale=1.2] at ({xloc}, -1) {{\Huge  {i+1}}};", file=output)

    # add for bpp
    #####################################
    if bpp_file is not None:
        for a, b, prob in bpp_dict:
            color = "red" # not in gold
            if pair_agree(a, b, ss): # in gold
                color = "blue" 
            lw = ", thick"
            prob *= 100
            color += "!" + str(prob)
            style = "[" + color + lw + "]"
            draw_eclipse(a-1, b-1, style)
            
    else:
        pairs = get_pairs(ss)
        for pair in pairs:
            a, b = pair
            color = 'blue'
            lw = ", thick"
            prob = 100
            color += "!" + str(prob)
            style = "[" + color + lw + "]"
            draw_eclipse(a, b, style);

    print("\\end{tikzpicture}", file=output)

    print("\\end{document}", file=output)
    
    output.close()
        
        
def linearpartition(seq, ss, tex_file='bpp.tex', bpp_file='bpp.txt', directory="./", delete=True):
    assert len(seq)==len(ss), f"{len(seq)} {len(ss)}"
    cmd = ['echo', seq]
    echo_process = subprocess.run(cmd, stdout=subprocess.PIPE)
    if bpp_file is not None:
        if os.path.exists(bpp_file):
            os.remove(bpp_file)
        cmd = ['./linearpartition', '-V', '-o', bpp_file]
        lp_process = subprocess.run(cmd, input=echo_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        # print(lp_process)
        if lp_process.returncode != 0:
            print(lp_process)
    draw_bpp(seq, ss, tex_file, bpp_file)
    cmd = ['pdflatex', '-output-directory', directory, tex_file]
    latex_process = subprocess.run(cmd, check=True)
    if latex_process.returncode != 0:
        print(latez_process)
    if delete:
        if os.path.exists(bpp_file):
            os.remove(bpp_file)
        if os.path.exists(tex_file):
            os.remove(tex_file)
        
        
def batch_plot():
    df = pd.read_csv('/nfs/stak/users/zhoutian/acl/repo/RNA-Fold/data/eterna100/nemo_default.csv')
    for i, row in df.iterrows():
        seq = row['rna']
        ss = row['structure']
        linearpartition(seq, ss, tex_file=f'nemo_{i}.tex', bpp_file=f'nemo_{i}.txt', directory='flat/nemo')
        
    df = pd.read_csv('/nfs/stak/users/zhoutian/acl/repo/RNA-Fold/data/eterna100/df_walk_ee_r05.csv')
    for i, row in df.iterrows():
        seq = row['rna']
        ss = row['structure']
        linearpartition(seq, ss, tex_file=f'walk_05_{i}.tex', bpp_file=f'walk_05_{i}.txt', directory='flat/walk_05')
    
if __name__ == '__main__':
    seq = "GAAAAUCGAUGCUCUUGCCGCACGCCCAUUGCUGCUCGCGCACGAUCAUGAUCUGAAAAAGUCUGUCUGGCAGCACAAUGGCGUCGGAAGAGCUCGGCAAAA"
    ss  = ".....((((.((((((.(((.(((((.((((.((((.((.((.((.((.((.((.....))))))))))))))))))))))))))))))))))))))....."
    linearpartition(seq, ss, tex_file='bpp_flat.tex', bpp_file='bpp_flat.txt')
