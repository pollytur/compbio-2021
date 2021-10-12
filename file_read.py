
d = {}
with open("dic.txt", "r") as f:
    for line in f.readlines():
        S = line.split()
        d[S[0]] = S[1]

keys = d.keys()

cpt = 0
for k in keys:
    if d[k] == "1":
        cpt += 1
    else:
        break

counter = 0
for k in keys:
    counter += int(d[k])

cpt2=0
for k in keys:
    cpt2 += 1

print(cpt2)
print(cpt/counter)
