/******************************************************************************
 * Determine Unique sequence identifiers based on on two inputs
 * 
 *****************************************************************************/

package hse

import Stream._

/////////////////////////////////////////////////////////////////////////////

case class kmerType(seq: String, index: Int) {
  /*************************************************************************
   * Define a kmertype such that we keep the index together with the kmer  *
   *   sequence. We have to redefine the equality and hashCode functions   *
   *   such that this type works with the sets and hashtables              *
   *************************************************************************/
  override def equals(o: Any) = o match {
    case that: kmerType => that.seq.equalsIgnoreCase(this.seq)
    case _ => false
  }
  override def hashCode = seq.toUpperCase.hashCode
}

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

  def kmer_at_i(seq: String, i: Int, upst_downst: Int): kmerType = {
    /*************************************************************************
     * Return a substring based at location i, with up and downstream        *
     *   limits                                                              *
     *************************************************************************/
    new kmerType(seq.substring(i-upst_downst, i+upst_downst+1), i)
  }


  /////////////////////////////////////////////////////////////////////////////

  def allKmersSet(seq: String, k: Int): Set[kmerType] = {
    /*************************************************************************
     * Return all the Kmers in a string seq, given a size of k, as a set     *
     *************************************************************************/

    val upst_downst = ((k.toDouble - 1.0) / 2.0).toInt
    val i_locations = ((upst_downst) until (seq.length() - upst_downst)).toList

    i_locations.map( (i: Int) => kmer_at_i(seq, i, upst_downst)).toSet

  }

  /////////////////////////////////////////////////////////////////////////////

  def allKmersStream(seq: String, k: Int): Stream[kmerType] = {
    /*************************************************************************
     * Return all the Kmers in a string seq, given a size of k, as a stream  *
     *************************************************************************/

    val upst_downst = ((k.toDouble - 1.0) / 2.0).toInt
    val limit       = seq.length - upst_downst - 1

    def allKmersStreamHelper(seq: String, i: Int): Stream[kmerType] = {
      if (i > limit){
       return empty
      }
      return kmer_at_i(seq, i, upst_downst) #:: allKmersStreamHelper(seq, i+1)
    }

    allKmersStreamHelper(seq, upst_downst)

  }

  /////////////////////////////////////////////////////////////////////////////

  def uniqueKmers(kmers: List[Set[kmerType]]): List[Set[kmerType]] = {
    /**************************************************************************
     * Given a list of outputs from allKmersSet, return only those which are  *
     *   unique to one of them.                                               *
     * Special care if there is only one sample, then there is no need        *
     **************************************************************************/
     if (kmers.length > 1){
       val kmer_index : List[(Set[kmerType],Int)]   = kmers.zipWithIndex;
       val remove_index = (L: List[(Set[kmerType],Int)], index: Int) =>(L.patch(index, Nil, 1))
     
       def setUnion(L: List[Set[kmerType]]): Set[kmerType] = {
         var U = L(0);
         L.fold(U){ (a,b) => a.union(b) }
       }
       kmer_index.map(elem => elem._1.diff(setUnion(remove_index(kmer_index, elem._2).map(vi => vi._1))))
     } else {
       // We don't need to do the above, because it is a lonely species
       kmers
     }
  }

  /////////////////////////////////////////////////////////////////////////////

  def orthologUniqueKmers(orthologs: List[Fasta.Entry], k: Int) : (List[String], List[Set[kmerType]]) = {
    /*************************************************************************
     * Given a list of orthologs, return sets of unique kmers for each       *
     *   ortholog.                                                           *
     *************************************************************************/
    
    val orthologNames: List[String] = orthologs.map(o => o.description).toList;
    val orthologSeqs: List[String]  = orthologs.map(o => o.sequence).toList;

    (orthologNames, (uniqueKmers(orthologSeqs.map( s => allKmersSet(s,k)))))

  }

  /////////////////////////////////////////////////////////////////////////////


  def genomeUniqueKmers(kmers: Set[kmerType], gKmers: genomeKmers): Set[kmerType] = {
    def genomeUniqueKmerHelper(kmer: kmerType): List[kmerType] = {
      if(gKmers.countOccurences(kmer.seq) > 1){
        Nil
      }else{
        List(kmer)
      }
    }
    kmers.map(genomeUniqueKmerHelper).flatten.toSet
  }

  /////////////////////////////////////////////////////////////////////////////

  def orthologGenomeUniqueKmers(oKmers: List[(List[String], List[Set[kmerType]])], gKmers: genomeKmers): List[(List[String], List[Set[kmerType]])] = {
     oKmers.map(x => (x._1, x._2.map( y => genomeUniqueKmers(y, gKmers))))
  }

  /////////////////////////////////////////////////////////////////////////////

  //val oKmers = Set(new kmerType("ABCDE", 2), new kmerType("CDEFG", 4), new kmerType("WXYZA", 20))
  //val kmers = oKmers
  //val kmersList = kmers.toList.sortWith(_.index < _.index)
  //kmersList.foldLeft(List.empty[kmerType]){ (a,b) => {if (a.length == 0 || (a.last.index + k) < b.index){ a :+ b }else{ a }}}.toSet

  def nonRedundantKmers(oKmers: List[(List[String], List[Set[kmerType]])], k: Int): List[(List[String], List[Set[kmerType]])] = {
    def nonRedundantKmersHelper(kmers: Set[kmerType]): Set[kmerType] = {
      val kmersList = kmers.toList.sortWith(_.index < _.index)
      kmersList.foldLeft(List.empty[kmerType]){ (a,b) => {if (a.length == 0 || (a.last.index + k) < b.index){ a :+ b }else{ a }}}.toSet
    }

    oKmers.map(x => (x._1, x._2.map(y => nonRedundantKmersHelper(y))))
  }

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//class genomeKmers(fastas: List[Fasta.Entry], k: Int) {
///*****************************************************************************
// * A class to handle the kmer HashMaps for genomes                           *
// *   Using the allKmersStream function, we prevent the construction of all   *
// *   kmers in memory at once.                                                *
// *****************************************************************************/
//
//  case class hashMapEntry(description: String,  hm: imHashMap[String,Int])
//
//  //val tables = fastas.map(kmerHashMap)
//  val table = kmerHashMapMem(fastas)
//
//  /////////////////////////////////////////////////////////////////////////////
//
//  def kmerHashMap(fasta: Fasta.Entry): hashMapEntry = {
//    /**************************************************************************
//      *  Create a hashmap for each fasta sequence. Iterate over all kmers in  *
//      *   the stream and add them to the hashMap one by one                   *
//      *************************************************************************/
//
//    var kmerCounts = new mutHashMap[String,Int]()
//
//    val kmerIterator = kmerTools.allKmersStream(fasta.sequence, k).iterator
//
//    var kmer  = new String("  ")
//    var count = 0
//
//    while(kmerIterator.hasNext){
//      kmer  = kmerIterator.next()
//      count = 0
//      if (kmerCounts.contains(kmer)){
//        count = kmerCounts(kmer)
//      }
//      kmerCounts += Tuple2(kmer,count + 1);
//    }
//    
//    new hashMapEntry(fasta.description, new imHashMap[String, Int] ++ kmerCounts)
//
//  }
//
//  /////////////////////////////////////////////////////////////////////////////
//
//  def countOccurences(kmer: String): Int = {
//     /**************************************************************************
//      *  Look through all the hashmaps in the structure for the number of     *
//      *   kmers and sum them all at the end.                                  *
//      *************************************************************************/
//     def lookupHelper[A,B](hm: imHashMap[A,B], key: A, default: B): B = {
//       if (hm.contains(key)){
//         hm(key)
//       }else{
//         default
//       }
//     }
//     
//     tables.map( x => lookupHelper(x.hm, kmer, 0)).foldLeft(0)(_ + _)
//  }
//
//  /////////////////////////////////////////////////////////////////////////////
//
//}

// Just so I can refer to them without all the stuff
import collection.mutable.{HashMap  => mutHashMap}
import collection.immutable.{HashMap => imHashMap}


class genomeKmers(fastas: List[Fasta.Entry], k: Int) {
/*****************************************************************************
 * A class to handle the kmer HashMap for genome sequences                   *
 *   Using the allKmersStream function, we prevent the construction of all   *
 *   kmers in memory at once. This is the memory efficient version           *
 *****************************************************************************/

  val table = kmerHashMap(fastas)

  /////////////////////////////////////////////////////////////////////////////

  def kmerHashMap(fastas: List[Fasta.Entry]): imHashMap[String,Int] = {
    /**************************************************************************
      *  Create a hashmap. For each fasta sequence, create a kmerstream.      *
      *   Iterate over all kmers in the stream and add them to the hashMap    *
      *   one by one. This one is more memory efficient                       *
      *************************************************************************/
    
    var kmerCounts = new mutHashMap[String,Int]()
    
    for(fasta <- fastas){
      val kmerIterator = kmerTools.allKmersStream(fasta.sequence, k).iterator
      var kmer  = new kmerType("  ", 0)
      var count = 0
      
      while(kmerIterator.hasNext){
        kmer  = kmerIterator.next()
        count = 0
        if (kmerCounts.contains(kmer.seq)){
          count = kmerCounts(kmer.seq)
        }
        kmerCounts += Tuple2(kmer.seq,count + 1);
      }
    }
    
    new imHashMap[String, Int] ++ kmerCounts

  }

  /////////////////////////////////////////////////////////////////////////////

  def countOccurences(kmer: String): Int = {
    /*************************************************************************
     * Look in the hashMap for the kmer, and return the count. otherwise 0   *
     *************************************************************************/
    if (table.contains(kmer)) {
      table(kmer)
    }else{
      0
    }
  }

}
///////////////////////////////////////////////////////////////////////////////

