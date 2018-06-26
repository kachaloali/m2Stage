f = open("sp.txt", "r")
for line in f:
    ec = str_list = list(filter(' ', line.strip(' \n\r\t').split('"')[1:])) # fastest
    print(ec)
f.close()