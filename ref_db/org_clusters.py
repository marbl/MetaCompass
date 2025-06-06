import sys

def only_numerics(seq):
	digits = "".join(filter(str.isdigit, seq))
	return digits
    # return filter(type(seq).isdigit, seq)

def main():

	gene = sys.argv[1]
	# gene = 'ArgS_COG0018'
	# wrkpth = "/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2/RefSeq_V2_db/marker_index/"
	wrkpth = sys.argv[2]

	clusterDict=dict()
	cluster = -1
	with open(wrkpth+gene+'/'+gene+'_clustered.clstr') as f:
		for num, line in enumerate(f):
			seq = ""
			rep = "2"
			val = line.strip().split('\t')
			if (val[0][0] == ">"):
				cluster = only_numerics(val[0])
			else:
				seq = val[1].split(' ')[1][1:-3]
				rep =  val[1].split(' ')[2][0]
				if (rep == "*"):
					clusterDict[cluster]=seq


	final_fr = open(wrkpth+gene+'/'+gene+'_clustered.clusters', 'w')
	final_fr.write('#seq_id\tcluster\trepresentative\n')

	cluster = -1
	with open(wrkpth+gene+'/'+gene+'_clustered.clstr') as f:
		for num, line in enumerate(f):
			seq = ""
			rep = "2"
			val = line.strip().split('\t')
			if (val[0][0] == ">"):
				cluster = only_numerics(val[0])
			else:
				seq = val[1].split(' ')[1][1:-3]
				representative =  clusterDict[cluster]

				# Check if the cluster ID exists in clusterDict before retrieving the representative sequence
				# THIS IS AN ADDITION(2023) -- Need to find root cause of why this is needed
				# if cluster in clusterDict:
				# 	representative = clusterDict[cluster]
				# else:
				# 	representative = "N/A"  # Set a default value if the cluster ID is not found

				# print(str(seq)+'\t'+str(cluster)+'\t'+str(representative)+'\n')
				final_fr.write(str(seq)+'\t'+str(cluster)+'\t'+str(representative)+'\n')

	final_fr.close()

if __name__ == '__main__':
	main()
