import argparse
import subprocess
import os

import sys
import math
from collections import defaultdict
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
\begin{tikzpicture}[darkstyle/.style={}, scale=2] %circle,draw,fill=gray!10}]
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
    for i, base in enumerate(bases[1:], 1):
        if circular:
            angle = 360./(length+lengthfix)*i
            if counter_clockwise: # hzhang
                angle = 450-angle -20
            else:
                angle = angle-70 + rotate 

            print("\\node [darkstyle]  (%d) at (%f:10cm) {};" % (i, angle), file=output)
            if length <= 100:
                gap = 5
            elif length <= 200:
                gap = 10
            elif length <= 300:
                gap = 20
            elif length <= 700:
                gap = 50
            elif length <= 2000:
                gap = 100
            elif length <= 3000:
                gap = 200
            elif length <= 5000:
                gap = 300
            else:
                gap = 400

            if i % gap == 0 and i < len(bases)-10:
                print("\\node [scale=2]           (%d,1) at (%f:10.8cm) {\Huge %d};" % (i, angle, i), file=output)
            if i > 1:
                print("\\draw (%d.center) -- (%d.center);" % (i, i-1), file=output)

        else:
            print("\\node [darkstyle]  (%d) at (%d,0) {};" % (i, i), file=output)
            if i % 5 == 0:
                print("\\node []           (%d,1) at (%d,-1) {%d};" % (i, i, i), file=output)

    if circular:

        for j in range(4):
            i += 1
            angle = 360./(length+lengthfix)*(i-4)

            if counter_clockwise: # hzhang
                angle = 450-angle -20
            else:
                angle = angle-70 + rotate 

            print("\\node [darkstyle]  (%d) at (%f:10cm) {};" % (i, angle), file=output)

        angle = 360./(length+lengthfix)
        angle = 450-angle
        if length <= 50:
            angle5 = angle - 30
            angle3 = angle + 30
        elif length <= 100:
            angle5 = angle - 20
            angle3 = angle + 20
        elif length <= 200:
            angle5 = angle - 17
            angle3 = angle + 17
        else:
            angle5 = angle - 13
            angle3 = angle + 13
        print("\\node [scale=2](3prime) at (%f:10cm) {\LARGE \\textbf{%s}};" % (angle3, "5'"), file=output)
        print("\\node [right=9.5cm of 3prime, scale=2] {\LARGE \\textbf{%s}};" % "3'", file=output)


    # add for bpp
    #####################################
    if bpp_file is not None:
        for a, b, prob in bpp_dict:
            color = "red" # not in gold
            if pair_agree(a, b, ss): # in gold
                color = "blue" 
            lw = ",thick"
            prob *= 100
            color += "!" + str(prob)
            style = "[" + color + lw + "]"

            if circular: 
                dist = b - a
                revdist = length+lengthfix - dist

                deg = 90 * (0.5-(dist+.0) / (length+lengthfix+.0))

            else:
                deg = 20 

            if counter_clockwise:
                drawarc_counterclockwise(a, b, deg, style, length, lengthfix)
            else:
                drawarc_clockwise(a, b, deg, style, length, lengthfix)
    else:
        pairs = get_pairs(ss)
        for pair in pairs:
            a, b = pair
            a += 1
            b += 1
            color = 'blue'
            lw = ",thick"
            prob = 100
            color += "!" + str(prob)
            style = "[" + color + lw + "]"
            if circular: 
                dist = b - a
                revdist = length+lengthfix - dist

                deg = 90 * (0.5-(dist+.0) / (length+lengthfix+.0))

            else:
                deg = 20 

            if counter_clockwise:
                drawarc_counterclockwise(a, b, deg, style, length, lengthfix)
            else:
                drawarc_clockwise(a, b, deg, style, length, lengthfix)
                

    print("\\end{tikzpicture}", file=output)

    print("\\end{document}", file=output)
    
    output.close()

    # print("%d out of %d sequences have pseudoknots" % (num_hasknot, index), file=logs) #TODO
    


def linearpartition(seq, ss, tex_file='bpp.tex', bpp_file='bpp.txt', ):
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
    cmd = ['pdflatex', tex_file]
    latex_process = subprocess.run(cmd, check=True)
    if latex_process.returncode != 0:
        print(latez_process)
    

if __name__ == '__main__':
    # seq = "AGCCGUGC"
    # seq = "GAAAAUCGAUGCUCUUGCCGCACGCCCAUUGCUGCUCGCGCACGAUCAUGAUCUGAAAAAGUCUGUCUGGCAGCACAAUGGCGUCGGAAGAGCUCGGCAAAA"
    # ss  = ".....((((.((((((.(((.(((((.((((.((((.((.((.((.((.((.((.....))))))))))))))))))))))))))))))))))))))....."
    import pandas as pd
    # df = pd.read_csv('/nfs/stak/users/zhoutian/acl/repo/RNA-Fold/data/eterna100/rna_inverse_v249_2.csv')
    df = pd.read_csv('/nfs/stak/users/zhoutian/acl/repo/RNA-Fold/data/eterna100/nemo_default.csv')
    seq = df.iloc[68]['rna']
    # ss = df.iloc[68]['structure']
    ss = eval(df.iloc[68]['opts'])[1][1]
    # assert len(seq)==len(ss), f"{len(seq)} {len(ss)}"
    # cmd = ['echo', seq]
    # echo_process = subprocess.run(cmd, stdout=subprocess.PIPE)
    # if os.path.exists(bpp_file):
    #     os.remove(bpp_file)
    # cmd = ['./linearpartition', '-V', '-o', bpp_file]
    # lp_process = subprocess.run(cmd, input=echo_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # print(lp_process.stdout.decode('utf-8'), end="")
    # print(lp_process.stderr.decode('utf-8'), end="")

    # print(lp_process)
    linearpartition(seq, ss, tex_file='nemo68_2.tex', bpp_file=None)
