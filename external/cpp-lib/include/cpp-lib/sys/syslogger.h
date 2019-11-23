//
// Copyright 2015 KISS Technologies GmbH, Switzerland
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Component: SYSUTIL
//
// TODO: This is UNIX only.  Provide a Windows implementation
//

#ifndef CPP_LIB_SYSLOGGER_H
#define CPP_LIB_SYSLOGGER_H

#include <functional>
#include <ostream>
#include <vector>

#include "boost/lexical_cast.hpp"

#include "cpp-lib/util.h"

namespace cpl {

namespace util {

namespace log {


//
// This assumes that the constants in <syslog.h> are LOG_EMERG, LOG_ALERT etc.,
// starting with 0.  Probably this is the case due to RFC 5424.
//

enum class prio {
  EMERG  , // system is unusable
  ALERT  , // action must be taken immediately
  CRIT   , // critical conditions
  ERR    , // error conditions
  WARNING, // warning conditions
  NOTICE , // normal but significant condition
  INFO   , // informational
  DEBUG	   // debug-level messages
} ;

// The default priority, reset after each flush.
inline prio default_prio() { return prio::INFO; }

// Textual representation of the given value.
const char* to_string(prio);

// Logging priority from string.  Throws if s is not
// any of ERROR, WARNING, etc. (NOT the enum names, 
// but prionames in .cpp).
prio prio_from_string(std::string const& s);

// Flags for setminprio(): Set minimum priority for syslog, echo
// stream or both.
unsigned short constexpr SYSLOG = 1;
unsigned short constexpr ECHO   = 2;
unsigned short constexpr BOTH   = SYSLOG | ECHO;

} // log

} // util


namespace detail_ {

struct prio_setter {
  // The priority to set
  cpl::util::log::prio const p;
  // Set SYSLOG, ECHO or both priorities
  unsigned short const which;
};

//
// Syslog_writer implements the writer part of the reader_writer concept, 
// see util.h
//

struct syslog_writer {

  // Better be explicit
  syslog_writer()                      = delete ;
  syslog_writer(syslog_writer const& ) = delete ;
  syslog_writer(syslog_writer      &&) = delete ;

  syslog_writer& operator=(syslog_writer const& ) = delete ;
  syslog_writer& operator=(syslog_writer      &&) = delete ;

  syslog_writer(std::string const& t, std::ostream* echo);

  // reader_writer implementation of write()
  int write(char const* buf, int n);
  void shutdown_write() {}

  // Minimum level to log on syslog
  cpl::util::log::prio minlevel_syslog;
  // Minimum level to log on echo stream
  cpl::util::log::prio minlevel_echo  ;

  // Currently set log level (e.g. os << cpl::util::CRIT;)
  cpl::util::log::prio currlevel;

  // The tag, NULL-terminated character array
  char const* tag() const { return &tag_[0]; }

  // Sets the echo clock, used for timestamps on the echo stream
  // (syslogger always uses system time).
  // If cl returns a negative value, no timestamp is written in the
  // echo stream.
  void set_echo_clock(std::function<double()> const& cl) {
    echo_clock_ = cl;
  }

  void set_echo_stream(std::ostream* const s) {
    echo_ = s;
  }

private:
  // Prepended to each write call
  std::vector<char> tag_;

  // An ostream to echo the logs, typically cout or cerr.
  std::ostream* echo_;

  // The clock used for echo.
  std::function<double()> echo_clock_;
} ;

} // detail_


namespace util {

namespace log {

//
// Base-from-member idiom doesn't work well with streams which
// use multiple inheritance.
//
// So we first construct an 'empty' ostreambuf and then give
// it our reader_writer.
//
// http://www.boost.org/doc/libs/1_47_0/libs/utility/base_from_member.html
//
// It is needed because the syslog_writer needs to be initialized before
// the 'real' base class cpl::util::ostreambuf<>.
//

typedef cpl::util::ostreambuf< cpl::detail_::syslog_writer > syslogstreambuf ;


//
// This is the main class.  It derives from ostream and can be used
// exchangeably with any standard stream.  Use prio::INFO, ... to
// set the log level (priority).
//
// See syslogger-test.cpp for extensive usage examples.
//
// Note: This class is neither copyable nor moveable.  It could
// potentially made be moveable, but that may be difficult due to
// inheritance complications etc.
//

struct syslogger : cpl::util::file::owning_ostream<syslogstreambuf> { 

  //
  // Converts the given tag to a string using boost::lexical_cast<>.
  // Pass "" for no tag.  Pass &std::cout or similar to echo
  // log messages to the given stream.
  //
  // Uses the given clock for UTC timestamps *on the echo stream*.  
  // This is useful for testing.  Timestamps in the 'real' logs
  // are added by the syslog daemon and cannot be modified.
  //

  template< typename T >
  syslogger(T const& tag, std::ostream* echo = NULL,
      std::function<double()> const& echo_clock = cpl::util::utc)
    : cpl::util::file::owning_ostream<syslogstreambuf>(
         syslogstreambuf(
             std::make_shared<cpl::detail_::syslog_writer>(
                 boost::lexical_cast<std::string>(tag), echo)))
  { buffer().reader_writer().set_echo_clock(echo_clock); }

  syslogger() : syslogger("") {}

  // If an echo stream is set, all messages will be output to the echo
  // stream in addition to *this.
  void set_echo_stream(std::ostream* const s) {
    buffer().reader_writer().set_echo_stream(s);
  }

  // Use the given clock function for timestamps on the echo stream.
  void set_echo_clock(std::function<double()> const& cl) {
    buffer().reader_writer().set_echo_clock(cl);
  }

  // Selects the minimum Log level on *this, the echo stream, or both.
  void set_minlevel(cpl::detail_::prio_setter const& ps) {
    if (ps.which & SYSLOG) 
    { buffer().reader_writer().minlevel_syslog = ps.p; }
    if (ps.which & ECHO  ) 
    { buffer().reader_writer().minlevel_echo   = ps.p; }
  }
} ;

// Sets minimal log priority, e.g.
// sl << setminprio(prio::ERR) to log only important messages.
inline cpl::detail_::prio_setter setminprio(
    prio const p, unsigned short const which = BOTH) {
  return cpl::detail_::prio_setter{p, which};
}

// Sets priority for current message
std::ostream& operator<<(std::ostream&, prio);

// Sets min prio, use setminprio() above.
std::ostream& operator<<(std::ostream&, cpl::detail_::prio_setter const&);

// Logs an error, format:
// ERROR: msg: what
void log_error(std::ostream&, std::string const& msg, std::string const& what);

// A class to redirect a logstream to an echo other stream for testing
// purposes.
// TODO: Store and reset SYSLOG priority
struct testmode_sentry {
  testmode_sentry(
      syslogger   & sl,
      std::ostream& echo,
      const double& echo_time = 0.0) 
  : sl_(sl) {
    sl_.set_echo_stream(&echo);
    sl_.set_echo_clock([echo_time]() { return echo_time; });
    sl_ << setminprio(prio::EMERG, SYSLOG);
  }

  ~testmode_sentry() {
    sl_.set_echo_stream(nullptr);
    sl_ << setminprio(default_prio(), SYSLOG);
  }

  private:
    syslogger& sl_;
};



} // namespace log

} // namespace util

} // namespace cpl

#endif // CPP_LIB_SYSLOGGER_H
