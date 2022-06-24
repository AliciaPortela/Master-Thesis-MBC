import tskit
import json
ts = tskit.load("hgdp_tgp_sgdp_chr12_p.dated.trees")
import sys
orig_stdout = sys.stdout
f = open('out.txt', 'w')
sys.stdout = f
for i in range(3800):
	print(json.loads(ts.individual(i).metadata))
	 
sys.stdout = orig_stdout
f.close()
