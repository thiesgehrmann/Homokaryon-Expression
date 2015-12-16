all:
	~/scala-2.11.7/bin/scalac fasta.scala
	~/scala-2.11.7/bin/scalac utils.scala
	~/scala-2.11.7/bin/scalac kmertools.scala
	~/scala-2.11.7/bin/scalac unique_probes.scala

run: all
	#~/scala-2.11.7/bin/scala unique_probes test1.fasta,test2.fasta mapping_0_0.fasta,mapping_0_1.fasta 21
	~/scala-2.11.7/bin/scala -J-Xmx2g uniqueProbes mapping_0_0.fasta,mapping_0_1.fasta mapping_0_0.fasta,mapping_0_1.fasta 21
