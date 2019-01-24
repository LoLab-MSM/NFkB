result = []
with open('_diswitch.txt','r') as inputFile:
    for line in inputFile:
        result.append(float(line))