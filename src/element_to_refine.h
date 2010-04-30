// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_ELEMENT_TO_REFINE_H
#define __H2D_ELEMENT_TO_REFINE_H

#include "refinement_type.h"

struct H2D_API ElementToRefine { ///< A record of an element and a selected candidate.
  int id; //ID of element
  int comp; //componet
  int split; //proposed refinement
  int p[H2D_MAX_ELEMENT_SONS]; //Ecoded orders of sons. If V order is zero, V order is equal to U order.
  int q[H2D_MAX_ELEMENT_SONS]; //Encoded H orders of sons. If V order is zero, V order is equal to U order.

public:
  ElementToRefine() : id(-1), comp(-1) {};
  ElementToRefine(int id, int comp) : id(id), comp(comp), split(H2D_REFINEMENT_H) {};
  ElementToRefine(const ElementToRefine &orig) : id(orig.id), comp(orig.comp), split(orig.split) {
    copy_orders(p, orig.p);
    copy_orders(q, orig.q);
  };
  ElementToRefine& operator=(const ElementToRefine& orig);
  int get_num_sons() const { return get_refin_sons(split); }; ///< Returns a number of sons.
  static inline void copy_orders(int* dest, const int* src) { ///< Copies array of orders.
    memcpy(dest, src, sizeof(int) * H2D_MAX_ELEMENT_SONS);
  }
};

extern H2D_API std::ostream& operator<<(std::ostream& stream, const ElementToRefine& elem_ref); ///< Dumps contants of the structure. Used for debugging purposes.

/// Checks whether the input contains a given tag. Throws an exception otherwise.
class HERMES2D_API TagChecker {
  const std::string& tag;
public:
  explicit TagChecker(const std::string& tag) : tag(tag) {};
  const std::string& get_tag() const { return tag; };
};
extern HERMES2D_API std::istream& operator>>(std::istream& stream, const TagChecker& checker); ///< Performs checking

/// \brief Refinement writer and reader.
///
/// If binary, file is always little endian.
class HERMES2D_API ElementToRefineStream {
private:
  static const char* H2DER_START_TAG; ///< Tag which start file.
  static const char* H2DER_BIN_TAG; ///< Tag which defined uncompressed binary contents.
  static const char* H2DER_VECTOR_TAG; ///< Tag which start a vector of refinements.
  static const int H2DER_SIZE_BYTESIZE = 1; ///< A size of all size declaration.

  std::fstream stream; ///< Underlaying stream.
  bool little_endian; ///< True if system is little endian.

  static uint8_t get_byte_size(int value); ///< Returns number of bytes necessary to store a given value.
  void write_bytes(const void* data, int num_bytes); ///< Writes a given number of bytes in a little-edinan form.
  void write_bytes(const int data, int num_bytes); ///< Writes a given number of LSB bytes from intenger in a little-edinan form.
  int read_bytes(int num_bytes); ///< Reads bytes and converts then to a signed integer.
  void write_header(); ///< Writes header.
  void read_header(); ///< Writes header.

public:
  explicit ElementToRefineStream(const char* filename, std::ios_base::openmode mode); ///< Opens the stream.
  void open(const char* filename, std::ios_base::openmode mode); ///< Opens the stream.
  bool is_open() const { return stream.is_open(); }; ///< Returns true if the stream is open.
  bool eof() const { return stream.eof(); }; ///< Returns true if EOF has been reached.
  bool good() const { return stream.good(); }; ///< Returns true if stream is read for operation.
  bool operator!() const { return stream.fail(); }; ///< Returns true if error occured.
  void close() { stream.close(); }; ///< Closes the stream

  friend HERMES2D_API ElementToRefineStream& operator<<(ElementToRefineStream& stream, const std::vector<ElementToRefine>& elem_refs);
  friend HERMES2D_API ElementToRefineStream& operator>>(ElementToRefineStream& stream, std::vector<ElementToRefine>& elem_refs);
};
extern HERMES2D_API ElementToRefineStream& operator<<(ElementToRefineStream& stream, const std::vector<ElementToRefine>& elem_refs); ///< Stores a list of refinements.
extern HERMES2D_API ElementToRefineStream& operator>>(ElementToRefineStream& stream, std::vector<ElementToRefine>& elem_refs); ///< Stores a list of refinements.

#endif
