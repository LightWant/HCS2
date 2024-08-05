import sys
import re

allG0 = 0
total = 0
cG0 = 0
aG0 = 0
bG0 = 0
bcG0 = 0
with open(sys.argv[1], "r") as f:
    lines = f.readlines()
    lines = lines[7:-3]
    # print "lines:", len(lines)
    # print lines[0], lines[-1]
    for line in lines:
        # p = line.index(',')
        # abc = line[p + 2:].split()
        # # d = line[:p]
        # a, b, c = map(int, )

        # # print d, a, b, c
        
        mt = re.findall('\d+', line)
        ns = map(int, mt)
        if len(ns) != 4:
            continue
        d, a, b, c = ns

        total += 1
        if a > 0 and b > 0:
            allG0 += 1
        if c > 0:
            cG0 += 1
        if a > 0:
            aG0 += 1

        if b > 0:
            bG0 += 1
        
        if b > 0 or c > 0:
            bcG0 += 1
    

print aG0, bG0, allG0, cG0, bcG0, total

# ep2_ABC.txt 0 8036253 81936485
# dblp2_ABC.txt 5189 300 0 262 5803
# dblp3_ABC.txt 120216 20072 1856 7115 175976
# as3_ABC.txt 20065 7344 809 2359 30361
# as2_ABC.txt 10234 3172 0 2465 14541