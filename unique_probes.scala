/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs               *
 *                                                                            *
 ******************************************************************************/

import hse.Fasta
import hse.kmerTools
import hse.Utils
import hse.genomeKmers

object uniqueProbes {

  /////////////////////////////////////////////////////////////////////////////

  def main(args: Array[String]) {

    if(args.length < 3) {
      usage("uniqueProbes")
      System.exit(1);
    }

      // Read input
    val inputFastas  = args(0).split(',').toList
    val inputGenomes = args(1).split(',').toList
    val k            = args(2).toInt
    val out_prefix   = args(3)

      // Check input makes sense
    if(inputFastas.length != inputGenomes.length) {
      usage("uniqueProbes")
      System.exit(2);
    } 

      // Read the genes from different alleles
    val fastas  = inputFastas.map(Fasta.read);
      // Read genome FASTA files
    val genomes = inputGenomes.map(Fasta.read);

      // Determine the unique kmers for each group of alleles
    val orthologous_kmers = Utils.suckLessZip(fastas).map( o => kmerTools.orthologUniqueKmers(o, k));
      // Create the HashMaps for the genomes
    val genomeKmerCounts  = new genomeKmers(genomes.flatten, k);

      // Determine the genome-wide unique probes for each allele
    val unique_probes = kmerTools.orthologGenomeUniqueKmers(orthologous_kmers, genomeKmerCounts)

      // Output
    output(unique_probes, inputGenomes.length, out_prefix)

    System.exit(0);

  }

  /////////////////////////////////////////////////////////////////////////////

  def output(probes: List[(List[String], List[Set[String]])], n_genomes: Int, prefix: String): Unit = {

    // Still need to figure out file IO...

    for(genome_i <- 0 until n_genomes) {

     println("--------WUT WUT" + genome_i.toString)

      for(gene_i <- 0 until probes.length){
        val gene_probes = probes(gene_i)._2(genome_i).toList
        for(probe_i <- 0 until gene_probes.length)
          println(probes(gene_i)._1(genome_i) + '\t' + probe_i.toString + '\t' + gene_probes(probe_i))
      }
    }

  }

  /////////////////////////////////////////////////////////////////////////////*/

  def usage(arg0: String): Unit = {
    println("Unique identifier detection")
    println("Usage: " + arg0 + " <sequences> <genomes> <k> <out_prefix>");
    println("  sequences:  Comma seperated list of FASTA files containing gene DNA sequences")
    println("  genomes:    Comma seperated list of FASTA files containing genome DNA sequences")
    println("  k:          The size of kmer to use")
    println("  out_prefix: The prefix of the output files. Will produce *.tsv files for each input file")
  }

  /////////////////////////////////////////////////////////////////////////////*/

}
