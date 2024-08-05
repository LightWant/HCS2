from random import randint, random, uniform

def ER(n1, n2,  p, outF):
    #```https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93R%C3%A9nyi_model```
    adj = [[] for u in range(n1)]
    m = 0

    for u in range(n1):
        for v in range(n2):
            if random() <= p:
                adj[u].append(v)
                m += 1

    with open(outF, 'w') as f:
        f.write("{} {} {}\n".format(n1, n2, m))

        for u in range(n1):
            for v in adj[u]:
                f.write("{} {}\n".format(u, v))

# p = 0.1
# while p < 0.5:
#     ER(200, 200, p, "data/ER/1000all/{}.txt".format(p))
#     p += 0.01

T = 30
for i in range(T):
    p = uniform(0.1, 0.4)
    # print "\"" + "data/ER/100all/" + "{}".format(p) + ".txt"+ "\",",
    print "\"{}\"".format(p) + ",",
    ER(100, 100, p, "data/ER/100all/{}.txt".format(p))


# def BA(n1, n2, p, outF):
#     #```https://en.wikipedia.org/wiki/Barab%C3%A1si%E2%80%93Albert_model```
#     adj = [[] for u in range(n1)]
#     degL = [0 for i in range(n1)]
#     degR = [0 for i in range(n2)]
#     m = 0

#     for u in range(100):
#         for v in range(100):
#             if random() <= p:
#                 adj[u].append(v)
#                 m += 1

#     nl = 100
#     nr = 100

#     while nl < n1 or nr < n2:
#         if nl < n1:
#             for v in range(nr):
#                 if random()*m <= degR[v]:
#                     adj[nl].append(v)
#                     m += 1
            
#             nl += 1

#         if nr < n2:
#             for u in range(nl):
#                 if random() * m <= degL[u]:
#                     adj[u].append(nr)
#                     m += 1
#             nr += 1

#     with open(outF, 'w') as f:
#         f.write("{} {} {}\n".format(n1, n2, m))

#         for u in range(n1):
#             for v in adj[u]:
#                 f.write("{} {}\n".format(u, v))

# BA(5000, 5000, 0.5, "data/BA.txt")
