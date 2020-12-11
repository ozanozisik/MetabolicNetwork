"""
Created on Tue Apr 21 16:44:01 2020
@author: Maroua, Ozan

This code converts Recon3D_301.mat human metabolism data obtained from
Virtual Metabolic Human database ( https://www.vmh.life/ ) to edge lists.
For genes, NCBI ID - Approved Symbol conversion is done using the file
obtained from HGNC ( https://www.genenames.org/ )
"""

import cobra.test
import cobra
import cobra.core.Reaction
import pandas as pd

df=pd.read_csv('HGNCGeneEntrez.txt', sep='\t', dtype={'Approved symbol':str, 'NCBI Gene ID(supplied by NCBI)':str})
df=df.set_index('NCBI Gene ID(supplied by NCBI)')
entrezToGene=df.to_dict()
entrezToGene=entrezToGene['Approved symbol']

model=cobra.io.load_matlab_model("Recon3D_301.mat")


unipartite=set()
bipartite=set()
tripartite=set()
unipartitePlusGenes=set()
nodeTypes=set()

for r in model.reactions:

	rId=r.id+'_reaction' #To prevent name conflicts
	nodeTypes.add((rId, 'Reaction'))

	for rt in r.reactants:
		rtId=rt.id.split('[')[0]
		nodeTypes.add((rtId, 'Metabolite'))
		for pr in r.products:
			prodId=pr.id.split('[')[0]
			unipartite.add((rtId,prodId))
			unipartitePlusGenes.add((rtId,prodId))
		bipartite.add((rtId,rId))
		tripartite.add((rtId,rId))


	for pr in r.products:
		prodId=pr.id.split('[')[0]
		nodeTypes.add((prodId, 'Metabolite'))
		bipartite.add((rId,prodId))
		tripartite.add((rId,prodId))

	for g in r.genes:
		gID=g.id.split('.')[0]
		if gID in entrezToGene:
			gName=entrezToGene[gID]
			nodeTypes.add((gName,'Gene'))
			tripartite.add((gName,rId))
			for pr in r.products:
				prodId=pr.id.split('[')[0]
				unipartitePlusGenes.add((gName,prodId))


f=open('nodeTypes.tsv','w')
for t in nodeTypes:
	f.write(t[0]+'\t'+t[1]+'\n')
f.close()


#Directed interaction from reactant to product
f=open('unipartite-metab.tsv','w')
for t in unipartite:
	f.write(t[0]+'\t'+t[1]+'\n')
f.close()

#Directed interaction from reactant to reaction to product
f=open('bipartite-metab-reac.tsv','w')
for t in bipartite:
	f.write(t[0]+'\t'+t[1]+'\n')
f.close()


#Directed interaction from reactant and gene to reaction to product
f=open('tripartite-metab-gene-reac.tsv','w')
for t in tripartite:
	f.write(t[0]+'\t'+t[1]+'\n')
f.close()


#Directed interaction from reactant and gene to product
f=open('bipartite-metab-gene.tsv','w')
for t in unipartitePlusGenes:
	f.write(t[0]+'\t'+t[1]+'\n')
f.close()
