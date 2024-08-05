from random import randint
import commands

n = 10
T = 100
K = 2
Q = 3
Q2 = 3
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
        for u in range(n):
            for v in adj[u]:
                f.write("{} {}\n".format(u, v))

logF = 'logs/diffMy/'
for i in range(T):
    f = outF + "{}.txt".format(i)
    run = "bin/pl noUVM -f_txt " + f + " -k "+str(K)+" -q "+str(Q)+" > " + logF + "{}a.txt".format(i)
    commands.getstatusoutput(run)


for i in range(T):
    c1 = 0
    c2 = 0

    with open(logF + "{}a.txt".format(i), "r") as f:
        lines = f.readlines()

        for line in lines:
            if line.startswith('answer:'):
                # maxBiCliqueCount: 3
                p = line.index(':')

                c1 = int(line[p + 1:])
            if line.startswith('cnt:'):
                p = line.index(':')

                c2 = int(line[p + 1:])
    
    # print c1, c2
    if c1 != c2:
        print i, c1, c2