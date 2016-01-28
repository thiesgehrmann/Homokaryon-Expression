package hse

import scala.util.parsing.combinator._
import java.io._


object Fasta {
  /****************************************************************************
   * https://gist.github.com/paradigmatic/3437345                             *
   ****************************************************************************/

  case class Entry( description: String, sequence: String )

  /////////////////////////////////////////////////////////////////////////////

  def read( fn: String ): List[Entry] = {
    val lines = io.Source.fromFile(fn).getLines.mkString("\n")
    fromString( lines )
  }

  /////////////////////////////////////////////////////////////////////////////

  def fromString( input: String ): List[Entry] =
    Parser.parse(input)

  /////////////////////////////////////////////////////////////////////////////

  private object Parser extends RegexParsers {

    lazy val header = """>.*""".r ^^ { _.tail.trim }  
    lazy val seqLine = """[^>].*""".r ^^ { _.trim }

    lazy val sequence = rep1( seqLine ) ^^ { _.mkString }

    lazy val entry = header ~ sequence ^^ { 
      case h ~ s => Entry(h,s)
    }

    lazy val entries = rep1( entry )

    def parse( input: String ): List[Entry]  = {
      parseAll( entries, input ) match {
        case Success( es , _ ) => es
        case x: NoSuccess =>  throw new Exception(x.toString)
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////////

  def print(fastaList: List[Entry]) = {
    def printSingleFasta(fastaSingle: Entry) = {
      println('>' + fastaSingle.description)
      println(fastaSingle.sequence)
    }
    fastaList.map(printSingleFasta)
  }

  /////////////////////////////////////////////////////////////////////////////

  def write(fastaList: List[Entry], fn: String) = {

    val outfd = new PrintWriter(new FileWriter(fn, false))

    for (e <- fastaList){
      outfd.write(">" + e.description + '\n');
      for ( l <- e.sequence.grouped(80).toList) {
        outfd.write(l + '\n');
      }
    }
    outfd.close()
  }

}
