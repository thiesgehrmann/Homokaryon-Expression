/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs
 * 
 *****************************************************************************/

package hse

import Stream._

///////////////////////////////////////////////////////////////////////////////

object kmerTools {
/*****************************************************************************
 * Tools to handle kmers.                                                    *
 *****************************************************************************/

  def usage(arg0: String): Unit = {
    println("Unique identifier detection")
    println("Usage: " + arg0 + " <sequences> <genomes>");
  }


  /////////////////////////////////////////////////////////////////////////////

  def kmer_at_i(seq: String, i: Int, upst_downst: Int): String = {
    /*************************************************************************
     * Return a substring based at location i, with up and downstream        *
     *   limits                                                              *
     *************************************************************************/
    seq.substring(i-upst_downst, i+upst_downst+1)
  }


  /////////////////////////////////////////////////////////////////////////////

  def allKmersSet(seq: String, k: Int): Set[String] = {
    /*************************************************************************
     * Return all the Kmers in a string seq, given a size of k, as a set     *
     *************************************************************************/

    val upst_downst = ((k.toDouble - 1.0) / 2.0).toInt
    val i_locations = ((upst_downst) until (seq.length() - upst_downst)).toList

    i_locations.map( (i: Int) => kmer_at_i(seq, i, upst_downst)).toSet

  }

  /////////////////////////////////////////////////////////////////////////////

  def allKmersStream(seq: String, k: Int): Stream[String] = {
    /*************************************************************************
     * Return all the Kmers in a string seq, given a size of k, as a stream  *
     *************************************************************************/

    val upst_downst = ((k.toDouble - 1.0) / 2.0).toInt
    val limit       = seq.length - upst_downst - 1

    def allKmersStreamHelper(seq: String, i: Int): Stream[String] = {
      if (i > limit){
       return empty
      }
      return kmer_at_i(seq, i, upst_downst) #:: allKmersStreamHelper(seq, i+1)
    }

    allKmersStreamHelper(seq, upst_downst)

  }

  /////////////////////////////////////////////////////////////////////////////

  def uniqueKmers(kmers: List[Set[String]]): List[Set[String]] = {
    /**************************************************************************
     * Given a list of outputs from allKmersSet, return only those which are  *
     *   unique to one of them.                                               *
     **************************************************************************/
     val kmer_index : List[(Set[String],Int)]   = kmers.zipWithIndex;
     val remove_index = (L: List[(Set[String],Int)], index: Int) =>(L.patch(index, Nil, 1))
     
     def setUnion(L: List[Set[String]]): Set[String] = {
       var U = L(0);
       L.fold(U){ (a,b) => a.union(b) }
     }
     kmer_index.map(elem => elem._1.diff(setUnion(remove_index(kmer_index, elem._2).map(vi => vi._1))))
  }

  /////////////////////////////////////////////////////////////////////////////

  def orthologUniqueKmers(orthologs: List[Fasta.Entry], k: Int) : (List[String], List[Set[String]]) = {
    /*************************************************************************
     * Given a list of orthologs, return sets of unique kmers for each       *
     *   ortholog.                                                           *
     *************************************************************************/
    
    val orthologNames: List[String] = orthologs.map(o => o.description).toList;
    val orthologSeqs: List[String]  = orthologs.map(o => o.sequence).toList;

    (orthologNames, (uniqueKmers(orthologSeqs.map( s => allKmersSet(s,k)))))

  }

  /////////////////////////////////////////////////////////////////////////////


  def genomeUniqueKmers(kmers: Set[String], gKmers: genomeKmers): Set[String] = {
    def genomeUniqueKmerHelper(kmer: String): List[String] = {
      if(gKmers.countOccurences(kmer) > 1){
        Nil
      }else{
        List(kmer)
      }
    }
    kmers.map(genomeUniqueKmerHelper).flatten.toSet
  }

  /////////////////////////////////////////////////////////////////////////////

  def orthologGenomeUniqueKmers(oKmers: List[(List[String], List[Set[String]])], gKmers: genomeKmers): List[(List[String], List[Set[String]])] = {
     oKmers.map(x => (x._1, x._2.map( y => genomeUniqueKmers(y, gKmers))))
  }

  /////////////////////////////////////////////////////////////////////////////

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Just so I can refer to them without all the stuff
import collection.mutable.{HashMap  => mutHashMap}
import collection.immutable.{HashMap => imHashMap}

class genomeKmers(fasta: List[Fasta.Entry], k: Int) {
/*****************************************************************************
 * A class to handle the kmer HashMaps for genomes                           *
 *   Using the allKmersStream function, we prevent the construction of all   *
 *   kmers in memory at once.                                                *
 *****************************************************************************/

  case class hashMapEntry(description: String,  hm: imHashMap[String,Int])

  val tables = fasta.map(kmerHashMap)

  /////////////////////////////////////////////////////////////////////////////

  def kmerHashMap(fasta: Fasta.Entry): hashMapEntry = {
    /**************************************************************************
      *  Create a hashmap for each fasta sequence. Iterate over all kmers in  *
      *   the stream and add them to the hashMap one by one                   *
      *************************************************************************/

    var kmerCounts = new mutHashMap[String,Int]()

    val kmerIterator = kmerTools.allKmersStream(fasta.sequence, k).iterator

    var kmer  = new String("  ")
    var count = 0

    while(kmerIterator.hasNext){
      kmer  = kmerIterator.next()
      count = 0
      if (kmerCounts.contains(kmer)){
        count = kmerCounts(kmer)
      }
      kmerCounts += Tuple2(kmer,count + 1);
    }
    
    new hashMapEntry(fasta.description, new imHashMap[String, Int] ++ kmerCounts)

  }

  /////////////////////////////////////////////////////////////////////////////

  def countOccurences(kmer: String): Int = {
     /**************************************************************************
      *  Look through all the hashmaps in the structure for the number of     *
      *   kmers and sum them all at the end.                                  *
      *************************************************************************/
     def lookupHelper[A,B](hm: imHashMap[A,B], key: A, default: B): B = {
       if (hm.contains(key)){
         hm(key)
       }else{
         default
       }
     }
     
     tables.map( x => lookupHelper(x.hm, kmer, 0)).foldLeft(0)(_ + _)
  }

  /////////////////////////////////////////////////////////////////////////////

}

///////////////////////////////////////////////////////////////////////////////

