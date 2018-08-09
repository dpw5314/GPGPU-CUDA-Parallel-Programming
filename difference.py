import os

for n in range (0,11):
    for i in range (1,21):
        for j in range (1,13):
            os.system("diff -w Result_out/InVivo/n_%d_v_%d_j_%d.txt Answer/InVivo/n_%d_v_%d_j_%d.txt" %(n,i,j,n,i,j)) 
