import sys
from collections import defaultdict

def sim(a,b):
    return float(len(set.intersection(set(a), set(b)))) / float(len(set.union(set(a), set(b))))


if __name__ == "__main__":

    d = defaultdict(list)
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            line = line.strip()
            tokens = line.split("\t")
            d[tokens[0]] = [int(i) for i in tokens[1].split(",")]

    for i in d:
        for j in d:
            if i != j:
                print (i, j, sim(d[i], d[j]))

