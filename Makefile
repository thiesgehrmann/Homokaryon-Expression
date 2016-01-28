all:
	~/scala-2.11.7/bin/scalac fasta.scala
	~/scala-2.11.7/bin/scalac utils.scala
	~/scala-2.11.7/bin/scalac kmertools.scala
	~/scala-2.11.7/bin/scalac unique_probes.scala

run: all
	~/scala-2.11.7/bin/scala uniqueProbes test1.fasta,test2.fasta test1.fasta,test2.fasta 21 out_prefix
	~/scala-2.11.7/bin/scala uniqueProbes mapping_0.fasta,mapping_0.fasta 
	#~/scala-2.11.7/bin/scala -J-Xmx2g uniqueProbes mapping_0_0.fasta,mapping_0_1.fasta mapping_0_0.fasta,mapping_0_1.fasta 21
