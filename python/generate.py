from random import randint
import commands

n = 100
T = 10
outF = "data/my/"

for i in range(T):
    adj = [[] for u in range(n)]
    m = 0

    for u in range(n - 1):
        du = randint(1, n - u - 1)
        m += du

        for j in range(du):
            v = randint(u + 1, n - 1)
            while v in adj[u]:
                v = randint(u + 1, n - 1)
            adj[u].append(v)

        adj[u].sort()

    with open(outF + "{}.txt".format(i), 'w') as f:
        # f.write("{} {} {}\n".format(n1, m))

        for u in range(n):
            for v in adj[u]:
                f.write("{} {}\n".format(u, v))
    
    # with open(outF + "{}D.txt".format(i), 'w') as f:

    #     for u in range(n - 1):
    #         for v in adj[u]:
    #             f.write("{},{}\n".format(u, v))
    #             f.write("{},{}\n".format(v, u))

    with open(outF + "{}_2.txt".format(i), 'w') as f:
        f.write("{} {}\n".format(n, m))
        for u in range(n):
            for v in adj[u]:
                f.write("{} {}\n".format(u, v))

logF = 'logs/diffMy/'
for i in range(T):
    f = outF + "{}.txt".format(i)
    run = "bin/d2kRun noUVM -f_txt " + f + " -k 2 -q 3 > " + logF + "{}a.txt".format(i)
    # print run
    commands.getstatusoutput(run)
    
    # binF = outF + "{}.bin".format(i)
    # run = "./ListPlex/pro1/toBin " + f +" "+ binF
    # commands.getstatusoutput(run)
    # run = "./ListPlex/pro1/listPlex " + binF + " 2 3 " + " > " + logF + "{}b.txt".format(i)
    # commands.getstatusoutput(run)
    f = outF + "{}_2.txt".format(i)
    run = "./kplexEnum/src/kplexes "+f+" -k=2 -q=3 -t=1 > " + logF + "{}b.txt".format(i)
    commands.getstatusoutput(run)

for i in range(T):
    c1 = 0
    c2 = 0

    with open(logF + "{}a.txt".format(i), "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith('Maximal K-Plex Count:'):
                # maxBiCliqueCount: 3
                p = line.index(':')

                c1 = int(line[p + 1:])
    
    with open(logF + "{}b.txt".format(i), "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith('Number of plexes    :'):
                # maxBiCliqueCount: 3
                p = line.index(':')

                c2 = int(line[p + 1:])
    
    # print c1, c2
    if c1 != c2:
        print i, c1, c2