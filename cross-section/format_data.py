#!/usr/bin/python
f = open('data','r')
energydata = dict()

for line in f:
    points = line.split()
    if points==[]:
        break
    if points[1]=='-1':
        points[1] = 'BG'
    if points[0] in energydata.keys():
        energydata[points[0]][0].append(points[1])
        energydata[points[0]][1].append(points[2])
        energydata[points[0]][2].append(points[3])
    else:
        energydata[points[0]]=([points[1]],[points[2]],[points[3]])

f2 = open('textables','w')
isleft = True
for e in energydata.keys():
    f2.write('\\begin{tabular}{|l|l|l|}\n')
    f2.write('\\hline\n\\multicolumn{2}{|c|}{'+e+'keV}\\\\\n')
    f2.write('\\hline\n')
    f2.write('Thickness(mm) & Count & Time(s)\\\\\n')
    for i in range(len(energydata[e][0])):
        f2.write('$'+energydata[e][0][i]+'$ & $' + energydata[e][1][i] + '$ & $' + energydata[e][2][i] + '$ \\\\\n')
    f2.write('\\end{tabular}\n')
    if isleft:
        f2.write('\\quad\n')
    else:
        f2.write('\n')
    isleft = not isleft
