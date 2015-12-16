
package hse

//import hse.Fasta

/*****************************************************************************/

object Utils {

/*****************************************************************************/

def suckLessZip[A](L: List[List[A]]): List[List[A]] = {
  /****************************************************************************
   * OK, finally figured out how to do polymorphism in Scala!
   * Now, what about these mutable types coming up here... I don't know what to do about that, yet...
   **************************************************************************/

  def zipHelper[A](previous: List[List[A]], current: List[A]): List[List[A]] = {
    val n_items = current.length
    var next = new collection.mutable.MutableList[List[A]]
    for ( i <- 0 until n_items) {
      next += previous(i) :+ current(i)
    }
    return next.toList
  }
  var previous  = L(0).map( x => List(x))
  for ( i <- 1 until L.length) {
    previous = zipHelper(previous, L(i));
  }
  return previous


}

/*****************************************************************************/

def suckZip(L: List[List[Any]]): List[List[Any]] = {
  /** called suckZip because I mix mutable types into this...
      If I figure out how to do this better, I will.
      Zips an arbitrary number of lists of lists together, rather than zip which does max 3 **/
  
  def zipHelper(previous: List[List[Any]], current: List[Any]): List[List[Any]] = {
    val n_items = current.length
    var next = new collection.mutable.MutableList[List[Any]]
    for ( i <- 0 until n_items) {
      next += previous(i) :+ current(i)
    }
    return next.toList
  }
  var previous  = L(0).map( x => List(x))
  for ( i <- 1 until L.length) {
    previous = zipHelper(previous, L(i));
  }
  return previous
}

/*****************************************************************************/

def suckZipFastas(L: List[List[Fasta.Entry]]): List[List[Fasta.Entry]] = {
  /** Same as suckZip, just for not Fasta.Entry specifically, so that I don't lose the type **/
  
  def zipHelper(previous: List[List[Fasta.Entry]], current: List[Fasta.Entry]): List[List[Fasta.Entry]] = {
    val n_items = current.length
    var next = new collection.mutable.MutableList[List[Fasta.Entry]]
    for ( i <- 0 until n_items) {
      next += previous(i) :+ current(i)
    }
    return next.toList
  }
  var previous  = L(0).map( x => List(x))
  for ( i <- 1 until L.length) {
    previous = zipHelper(previous, L(i));
  }
  return previous
}

/*****************************************************************************/

}
