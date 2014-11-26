#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>

class ReadContainer;

// Typedefs for simpler editing on reconsideration

typedef unsigned char chr_num_t;
typedef int chr_pos_t;
typedef std::vector<std::string> chr_names_t;
typedef std::unordered_map<std::string, chr_num_t> chr_nums_t;



class Genome {
public:
  chr_num_t getChrNum ( const std::string name ) const;
  std::string getChrName ( const chr_num_t num ) const;
  ReadContainer* getReadAt ( const chr_num_t chr, const chr_pos_t pos ) const;
  ReadContainer* getReadAt ( const std::string chr, const chr_pos_t pos ) const; 

  void createChromosome ( const std::string name, const uint32_t length = 0 );
  ReadContainer* createRead ( const chr_num_t chr, const chr_pos_t pos, const short len );
  void registerRead( ReadContainer* rc ); // Don't const this one!
  void read ( std::string& fileName );
  void read ( std::ifstream& input );
  friend std::ifstream& operator >> ( std::ifstream& input, Genome& target );

  void printout ( ReadContainer* element = nullptr );
  
private:
  chr_names_t chrNames; // vector: index# -> chromosome name
  chr_nums_t chrNums;   // map: string -> chromosome index#
  std::vector< unsigned short > chrLens;
  std::vector< std::map<chr_pos_t, ReadContainer*>* > readPos; // [chr#].find(pos) -> read starting at pos on chromosome#
  
  uint calcLength ( const std::string cigar ) const;

  bool parseHeaderLine ( const std::string& line );
  bool parseDataLine ( const std::string& line ); 
  enum SAM_FIELDS {
    FLAG = 1, RNAME = 2, POS = 3, CIGAR = 5, TLEN = 8, QUAL = 10
  };
};


#endif // GENOME_H