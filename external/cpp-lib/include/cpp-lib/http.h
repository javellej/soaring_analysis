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
// Component: HTTP
//


#ifndef CPP_LIB_HTTP_H
#define CPP_LIB_HTTP_H

#include <iosfwd>
#include <string>

namespace cpl {

namespace http {

/// Deprecated.  Use write_content_type() instead.
std::string json_header();

/// The constant "\r\n";
extern const char* const endl;

/// @return Default timeout for HTTP connections [s]
inline double default_timeout() { return 60; }

/// @return Default server identification
std::string default_server_identification();

/// Writes 'Content-type: <content_type>; charset = <charset>' and an empty line
/// to os.
void write_content_type(
    std::ostream& os, 
    const std::string& content_type,
    const std::string& charset = "utf-8");

/// Equivalent to write_content_type(os, "application/json", charset);
void write_content_type_json(
    std::ostream& os,
    const std::string& charset = "utf-8");

/// Equivalent to write_content_type(os, "text/plain", charset);
void write_content_type_text(
    std::ostream& os,
    const std::string& charset = "utf-8");

/// Equivalent to write_content_type(os, "text/csv", charset);
void write_content_type_csv(
    std::ostream& os,
    const std::string& charset = "utf-8");

///
/// @return Content type associated with the extension of name:
/// .html -> text/html
/// .txt  -> text/plain
/// others -> application/octet-stream
/// @todo This is massively incomplete
///

std::string content_type_from_file_name(const std::string& name);

/// Writes an HTTP 'Date:' header
/// @param time Time since UNIX epoch [s] (UTC.)
/// If < 0, uses cpl::util::utc().
void write_date(std::ostream& os, double now = -1);

/// Writes an HTTP 'Connection:' header, e.g. 'Connection: close'
void write_connection(std::ostream& os, const std::string& connection);

/// Writes an HTTP 'Server:' header
void write_server(
    std::ostream& os, 
    const std::string& server = default_server_identification());

///
/// Writes an HTTP header for response code 200 (i.e. OK).
/// After that header, the payload data must be transmitted.
/// @param time Time since UNIX epoch [s] (UTC) for the 'Date:' header,
/// semantics see write_date().
/// @param content_type The 'Content-Type:' header
/// Remarks:
/// * Doesn't send Content-Length and sends 'Connection: close'
///

void write_http_header_200(
    std::ostream& os,
    const std::string& content_type,
    double now = -1,
    const std::string& server = default_server_identification());

///
/// Writes an HTTP header for response code 404 (i.e. not found).
/// No payload is expected after the header.
///
/// @param reason A string to add after the 404 code
/// @todo Add an HTML payload with error information
///

void write_http_header_404(
    std::ostream& os,
    const std::string& reason,
    double now = -1,
    const std::string& server = default_server_identification());

/// Gets the specified URL and pipes the result to os.  Logs to log.
void wget( std::ostream& log , std::ostream& os , std::string url , 
           double timeout = default_timeout() ) ;


/// Data of an HTTP GET request
/// https://www.w3.org/Protocols/rfc2616/rfc2616-sec5.html
struct get_request {
  /// The HTTP version (1.0, 1.1 etc.)
  std::string version;

  /// The host field of the URI
  std::string host;

  /// The absolute path of the URI
  std::string abs_path;

  /// Contents of the 'User-Agent:' field
  std::string user_agent;

  /// Contents of the 'Accept:' field
  std::string accept;
};

/// Parses a GET request from the given istream.  Parses up to the first empty
/// line.
/// @param first_line Contains the actual request, e.g. "GET /foobar HTTP/1.1"
/// @param is Contains the rest of the request, i.e. the headers, e.g.
/// Host:, Accept: etc.
/// @param log If not nullptr, logs some conditions to it
get_request parse_get_request(
    const std::string& first_line,
    std::istream& is,
    std::ostream* log = nullptr);

} // namespace http

} // namespace cpl


#endif // CPP_LIB_HTTP_H
