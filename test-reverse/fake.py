



with open("test.fna", "x") as new:
    with open("/ogref.fna") as ref:
        ref.readline()
        while line:=ref.readline():
            new.write(line)
            if line2:=ref.readline() and line3:=ref.readline():
                