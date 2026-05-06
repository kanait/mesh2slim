////////////////////////////////////////////////////////////////////
//
// $Id: tokenizer.h 2019/02/09 16:02:02 kanai Exp $
//
// Copyright (c) 2002 Takashi Kanai
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
////////////////////////////////////////////////////////////////////

#ifndef TOKENIZER_H
#define TOKENIZER_H 1

#include <algorithm>
#include <string>
using namespace std;

typedef std::string String;

class tokenizer {
 private:
  String::size_type cur_;
  String str_;
  String del_;
  bool   ret_;
//   void skip();
  void skip() {
    if ( cur_ == str_.length() )
      cur_ = String::npos;
    if ( !ret_ && cur_ != String::npos ) {
      String::size_type tmp = str_.find_first_not_of(del_, cur_);
      if ( tmp != String::npos )
	cur_ = tmp;
    };
}
 public:
  typedef std::pair<String::size_type,String::size_type> range_type;
  tokenizer(const String& str, const String& del, bool ret=false) 
    : cur_(0), str_(str), del_(del), ret_(ret) {};
    
//   bool empty();
  bool empty() {
    skip();
    return cur_ == String::npos;
  };
//   String next();
  String next() {
    range_type range = next_range();
    return str_.substr(range.first, range.second);
  };
  String str() { return str_; };
//   tokenizer::range_type tokenizer::next_range() {
  range_type next_range() {
    skip();
    String::size_type start = cur_;
    //String::size_type tmp = cur_;
    if ( cur_ != String::npos ) cur_ = str_.find_first_of(del_,cur_);
    if ( cur_ == String::npos ) return range_type(start,str_.length()-start);
    if ( ret_ && start == cur_ && del_.find(str_[cur_]) != String::npos ) ++cur_;
    return range_type(start,cur_-start);
  };

//   void set_str(const String& str);
  void set_str(const String& str) {
    str_ = str;
    cur_ = 0;
  };

//   void set_delimiter(const String& del, bool ret=false);
  void set_delimiter(const String& delim, bool ret) {
    del_ = delim;
    ret_ = ret;
  };

//   size_t count() const;
  size_t count() const {
    size_t count = 0;
    String::size_type currpos = cur_;
    while ( currpos != String::npos ) {
      if ( !ret_ ) {
	currpos = str_.find_first_not_of(del_,currpos);
	if ( currpos == String::npos ) { ++count; break; }
      } else if ( currpos == str_.length() ) {
	break;
      }
      String::size_type start = currpos;
      if ( currpos != String::npos ) currpos = str_.find_first_of(del_,currpos);
      if ( currpos == String::npos ) { ++count; break; }
      if ( ret_ && start == currpos && del_.find(str_[currpos]) != String::npos ) ++currpos;
      ++count;
    }
    return count;
  };

};

// sample
// Split "apple,banana,cherry" by ","
//stx::tokenizer<string> tok("apple,banana,cherry", ",");
//while ( !tok.empty() ) {
//  cout << tok.next() << " ";
//}
//cout << endl;
// tokenizer (constructor)
     
// tokenizer::tokenizer(const String& str, const String& del, bool ret)
//   : cur_(0), str_(str), del_(del), ret_(ret) {}
/*  str: input string to tokenize */
/*  del: set of delimiter characters */
/*  ret: true to return delimiters as tokens */
/*  If ret == true, delimiters are also returned as tokens. */
/*  For example: */

// sample
// Split "apple_banana__cherry_" by "_"
/*  bool ret = false; */
/*  stx::tokenizer<string> tok("apple_banana__cherry_", "_", ret); */
/*  while ( !tok.empty() ) { */
/*    cout << tok.next() << " "; */
/*  } */
/*  cout << endl; */
/*  This extracts: "apple", "banana", "cherry". If ret = true, it extracts: */
/*  "apple", "_", "banana", "_", "_", "cherry", "_". */
/*  skip (private) */
     
// void tokenizer::skip() {
//   if ( cur_ == str_.length() )
//     cur_ = String::npos;
//   if ( !ret_ && cur_ != String::npos ) {
//     String::size_type tmp = str_.find_first_not_of(del_, cur_);
//     if ( tmp != String::npos )
//       cur_ = tmp;
//   }
// }

/*  empty */
// bool tokenizer::empty() {
//   skip();
//   return cur_ == String::npos;
// }

/*  Returns true when no more tokens are available. */

/*  next */
// String tokenizer::next() {
//   range_type range = next_range();
//   return str_.substr(range.first, range.second);
// }
/*  Returns the next extracted token. */


/*  next_range */

// tokenizer::range_type tokenizer::next_range() {
//   skip();
//   String::size_type start = cur_;
//   //String::size_type tmp = cur_;
//   if ( cur_ != String::npos ) cur_ = str_.find_first_of(del_,cur_);
//   if ( cur_ == String::npos ) return range_type(start,str_.length()-start);
//   if ( ret_ && start == cur_ && del_.find(str_[cur_]) != String::npos ) ++cur_;
//   return range_type(start,cur_-start);
// }

/*  Extracts the next token and returns its range (position and length). */
/*  Returns std::pair<String::size_type,String::size_type> range. */
/*  range.first is the token start position, range.second is the token length. */

/*  set_str */

// void tokenizer::set_str(const String& str) {
//   str_ = str;
//   cur_ = 0;
// }

/*  Sets the input string to tokenize. */
/*  Resets the current position to 0. */

/*  set_delimiter */

// void tokenizer::set_delimiter(const String& delim, bool ret) {
//   del_ = delim;
//   ret_ = ret;
// }

/*  Sets the delimiter character set. */
/*  Set the second argument to true to treat delimiters as tokens. */

/*  count */

// size_t tokenizer::count() const {
//   size_t count = 0;
//   String::size_type currpos = cur_;
//   while ( currpos != String::npos ) {
//     if ( !ret_ ) {
//       currpos = str_.find_first_not_of(del_,currpos);
//       if ( currpos == String::npos ) { ++count; break; }
//     } else if ( currpos == str_.length() ) {
//       break;
//     }
//     String::size_type start = currpos;
//     if ( currpos != String::npos ) currpos = str_.find_first_of(del_,currpos);
//     if ( currpos == String::npos ) { ++count; break; }
//     if ( ret_ && start == currpos && del_.find(str_[currpos]) != String::npos ) ++currpos;
//     ++count;
//   }
//   return count;
// }

/*  Returns the number of tokens available from the current position. */

/*  string_tokenizer / wstring_tokenizer */

/*  Commonly used tokenizer<std::string> and tokenizer<std::wstring> are typedef'ed */
/*  as string_tokenizer and wstring_tokenizer, respectively. */

#endif // TOKENIZER_H
