with open("./test-long/parsnp.xmfa") as xmfa:
    line = f.readline()
    ref = open("./test-long/ref.fna", "r")
    test1 = open("./test-long/test1.fna", "r")
    test2 = open("./test-long/test2.fna", "r")
    test3 = open("./test-long/test3.fna", "r")

    while line[0] == '#':
        line = f.readline()
    
    while line:
        if line[0] = '>':
            seq = int(line[1])


    


    ref.close()
    test1.close()
    test2.close()
    test3.close()
    


with open('./test-long/test3.fna') as f:
    contig = 0
    line = f.readline()
    while line:
        # print(start)
        if line[0] == ">":
            contig += 1
        elif start < 80:
            print(line[start:-1], end='')
            line = f.readline()
            print(line)
            break
        else:
            start -= 80

        line= f.readline()

