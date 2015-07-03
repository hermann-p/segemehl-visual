#ifndef GENOME_H
#define GENOME_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <memory>

class ReadContainer; // resolve dependency circle

// Typedefs for flexibility on changing requirements
typedef unsigned char chr_num_t;
typedef int chr_pos_t;
typedef std::vector<std::string> chr_names_t;
typedef std::unordered_map<std::string, chr_num_t> chr_nums_t;
typedef std::vector< std::map<chr_pos_t, std::shared_ptr<ReadContainer>>* > chr_map_t;
typedef std::vector<std::shared_ptr<ReadContainer>> link_list_t;
typedef std::unordered_set<std::shared_ptr<ReadContainer>> link_set_t;

class Genome {
public:
  Genome ( bool init_m = false, bool init_c = false );
  ~Genome ();
  chr_num_t getChrNum ( const std::string name ) const;
  std::string getChrName ( const chr_num_t num ) const;
  std::shared_ptr<ReadContainer> getReadAt ( const chr_num_t chr, const chr_pos_t pos ) const;
  std::shared_ptr<ReadContainer> getReadAt ( const std::string chr, const chr_pos_t pos ) const; 
  unsigned short getLength ( const chr_num_t n ) const;

  void createChromosome ( const std::string name, const uint32_t length = 0 );
  std::shared_ptr<ReadContainer> createRead ( const chr_num_t chr, const chr_pos_t pos, const short len );
  void registerRead( std::shared_ptr<ReadContainer> rc ); // Don't const this one!
  void read ( std::string& fileName );
  void read ( std::ifstream& input );
  friend std::ifstream& operator >> ( std::ifstream& input, Genome& target );

  link_set_t* multistrand;
  link_set_t* circular;

private:
  chr_names_t chrNames; // vector: index# -> chromosome name
  chr_nums_t chrNums;   // map: string -> chromosome index#
  std::vector< unsigned short > chrLens;
  
  // [chr#].find(pos) -> read starting at pos on chromosome# - on the heap because of size
  std::unique_ptr<chr_map_t> readPos; 
  
  uint calcLength ( const std::string cigar ) const;

  bool parseHeaderLine ( const std::string& line );
  bool parseDataLine ( std::vector<std::string>& splits ); 
  enum SAM_FIELDS {
    FLAG = 1, RNAME = 2, POS = 3, CIGAR = 5, TLEN = 8, QUAL = 10
  };
};


#endif // GENOME_H
