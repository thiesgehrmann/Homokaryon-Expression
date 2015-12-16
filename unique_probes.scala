/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs
 * 
 *****************************************************************************/

import hse.Fasta
import hse.kmerTools
import hse.Utils

object uniqueProbes {

  def main(args: Array[String]) {

    if(args.length < 3) {
      usage("uniqueProbes")
      System.exit(1);
    }

    val inputFastas  = args(0).split(',').toList
    val inputGenomes = args(1).split(',').toList
    val k            = args(2).toInt
    
    val fastas  = inputFastas.map(Fasta.read);
    val genomes = inputGenomes.map(Fasta.read);
    
    println(Utils.suckLessZip(fastas).map( o => kmerTools.orthologKmers(o, k)))

    println(Utils.suckLessZip(List(List(1,2,3,4), List(1,2,3,4), List(1,2,3,4), List(1,2,3,4))))

    System.exit(0);

  }

/*****************************************************************************/

  def usage(arg0: String): Unit = {
    println("Unique identifier detection")
    println("Usage: " + arg0 + " <sequences> <genome>");
  }

}
