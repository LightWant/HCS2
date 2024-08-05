tmp = []
with open("tmp.txt", "r") as f:
    for line in f.readlines():
        if line.startswith('ut') == False:
            continue
        u, c = map(int, line[3:].split())
        # if c > 0:
        tmp.append((u, c))
tmp2 = []
with open("tmp2.txt", "r") as f:
    for line in f.readlines():
        if line.startswith('ut') == False:
            continue
        u, c = map(int, line[3:].split())
        # if c > 0:
        tmp2.append((u, c))

print len(tmp), len(tmp2)
for i in range(len(tmp)):
    if tmp[i][1] != tmp2[i][1]:
        print tmp[i], tmp2[i]

# import commands

# T = 100
# outF = "data/my/"
# logF = 'logs/diffMy/'

# for i in range(T):
#     f = outF + "{}.txt".format(i)
#     run = "bin/run noUVM -f_txt " + f + " -k 2 -q 3 > " + logF + "{}a.txt".format(i)
#     # print run
#     commands.getstatusoutput(run)
    
#     # binF = outF + "{}.bin".format(i)
#     # run = "./ListPlex/pro1/toBin " + f +" "+ binF
#     # commands.getstatusoutput(run)
#     # run = "./ListPlex/pro1/listPlex " + binF + " 2 3 " + " > " + logF + "{}b.txt".format(i)
#     # commands.getstatusoutput(run)
#     f = outF + "{}_2.txt".format(i)
#     run = "./kplexEnum/src/kplexes "+f+" -k=2 -q=3 -t=1 > " + logF + "{}b.txt".format(i)
#     commands.getstatusoutput(run)

# for i in range(T):
#     c1 = 0
#     c2 = 0

#     with open(logF + "{}a.txt".format(i), "r") as f:
#         lines = f.readlines()

#         for line in lines:
#             if line.startswith('Maximal K-Plex Count:'):
#                 # maxBiCliqueCount: 3
#                 p = line.index(':')

#                 c1 = int(line[p + 1:])
    
#     with open(logF + "{}b.txt".format(i), "r") as f:
#         lines = f.readlines()

#         for line in lines:
#             if line.startswith('Number of plexes    :'):
#                 # maxBiCliqueCount: 3
#                 p = line.index(':')

#                 c2 = int(line[p + 1:])
    
#     # print c1, c2
#     if c1 != c2:
#         print i, c1, c2