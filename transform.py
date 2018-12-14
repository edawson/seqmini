import sys

if __name__ == "__main__":

    k = ""
    prev = None
    curr = []
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            line = line.strip()
            tokens = line.split("\t")
            k = tokens[1]
            if prev is None:
                prev = k

            if k != prev:
                print (k + "\t" + ",".join(curr) )
                curr = []
                prev = k
            elif k == prev:
                curr.append(tokens[0])
