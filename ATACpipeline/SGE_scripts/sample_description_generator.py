import sys
import os
def sample_description_generator(maindir):
	outf=open(os.path.join(os.getcwd(), 'sample_description.txt'),'w')
	fields=['flowcell','raw_fastq_name_r1','raw_fastq_name_r2','sample_name','dataset']
	outf.write('\t'.join(fields)+'\n')
	flowcelldirs=os.listdir(maindir)
	for flowcelldir in flowcelldirs:
		filenamelist=[]
		samplenamelist=[]
		if flowcelldir.startswith('.'):
			continue
		flowcellpath=os.path.join(maindir, flowcelldir)
		for filename in os.listdir(flowcellpath):
			if filename.startswith('.'):
				continue
			if filename.endswith('.gz'):
				filenamelist.append(filename)
				samplename=filename.split('.')[0][:-2]
				samplenamelist.append(samplename)
		samplenamelist=list(set(samplenamelist))
		for sample in sorted(samplenamelist):
			samplename=sample.split('_')[-1]
			rowsamplename1=sample+'_1.fastq.gz'
			rowsamplename2=sample+'_2.fastq.gz'
			if rowsamplename1 in filenamelist and rowsamplename2 in filenamelist:
				outf.write('\t'.join([flowcelldir, rowsamplename1,rowsamplename2, samplename, samplename])+'\n')
			else:
				print('{} does not have a pair of reads'.format(sample))
	outf.close()

if __name__=='__main__':
	sample_description_generator(sys.argv[1])
