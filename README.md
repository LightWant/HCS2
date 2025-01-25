# The source code and running stripts of "Counting Cohesive Subgraphs with Hereditary Properties".

`src/run are the starters.`

`src/plex are the code of plex counting`

`src/sdc are the code of deftectiver clique counting`

# The full version of the pdf.

`Compared to the submitted version, fullversion has details of listing based algorithms and applications.`


# compile

`clone the resp and cd into HCS2.`

`cmake -S ./`

`make`

We get the bin directory.

# data format: input command

`each edge per line: noUVM -f_txt filepath`

`the number of nodes and edges and then each edge per line: -f_txt filepath`

#  The description of the programs under bin/ (by examples)


pivot-based: count the dclique with s=1 and size in [5, 10] 
`bin/sdcc noUVM -f_txt data/caida.txt -s 1 -q 5 -Q 10`
Output

```
load graph: n 26475 m 106762 maxD 2628
changeToCoreOrder:2ms
coreNumber:22
s:1 q:5 Q:10
sdcCounting2DP::sdcCounting
sdcCounting2DPUB.cpp::runinit Hash
time:192ms
5-sdc:1000434
6-sdc:909487
7-sdc:776923
8-sdc:594667
9-sdc:390974
10-sdc:211550
```

listing-based: count the dclique with s=1 and size 5
`bin/pd noUVM -f_txt data/caida.txt -s 1 -q 5`

pivot-based:count the plex with s=1 and size in [5, 10] 
`bin/pc noUVM -f_txt data/caida.txt -s 1 -q 5 -Q 10`

listing-based:count the plex with s=1 and size 5
`bin/plv2 noUVM -f_txt data/caida.txt -s 1 -q 5`